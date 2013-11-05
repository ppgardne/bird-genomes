#!/usr/bin/perl 

#Given multiple gff files, merge them together and combine entries with the same coordinates. 

use warnings;
use strict;
use Getopt::Long;
use IO::File;

my (@gffs, $clanInfo, $verbose, $help);
my $overlapThreshold=10.0;
my $overlapDir=".";
&GetOptions( 
    "ot|overlapthresh=s"  => \$overlapThreshold,
    "g|gff=s@"            => \@gffs,
    "cl|claninfo=s"       => \$clanInfo,
    "od|overlapdir=s"     => \$overlapDir,
    "v|verbose"           => \$verbose,
    "h|help"              => \$help
    );

if( $help ) {
    &help();
    exit(1);
}
elsif (@gffs == 0){
    print "FATAL: no gff files given\n";
    &help();
    exit(1);
}

my $gffs = join(' ', @gffs);

my %clans;
if(defined($clanInfo) && (-s $clanInfo)){
    open(CL, "< $clanInfo") or die "FATAL: failed to open [$clanInfo]\n[$!]"; 
    while(my $cl=<CL>){
	my @cl=split(/\t/, $cl);
	#rfamacc        clanacc
	$clans{$cl[2]}=$cl[0];
    }
    close(CL);
}
else{
    print "WARNING: no Clan information! (which annotations should be merged?)\n"
}

#rethreshold some problematic families
my %rethreshold=(
    'RF00254'=>70, #mir-16
    'RF00606'=>50 #SNORD93
    );

my $fh = IO::File->new();
# sort features so they can be merged in coordinate
# order
$fh->open( "sortGffs.pl $gffs | uniq | " );
print "#running [sortGffs.pl $gffs |]\n" if(defined $verbose);
#$gffs = join('_', @gffs);
open(OV, "> $overlapDir/$$\_unclassified_overlaps_file.overlaps");
print OV "#GFFS: [$gffs]\n"; 
my ($cur, $prev, $ext, $rfid, $prfid, %uniq);
my $cnt=0;
while($cur = <$fh>) {
    chomp($cur); 
    next if ($cur =~ /^#/);
    if(defined($prev)){
	my @cur  = split(/\t/,  $cur);
	my @prev = split(/\t/, $prev);
	next if (not defined $cur[3] || not defined $cur[4] || not defined $cur[6]); #not a GFF
	$rfid="undef";	
	if($cur[8]=~/rfam-acc=(RF\d+)/){
	    $rfid=$1;
	    if( defined($rethreshold{$rfid})  && ($cur[5] < $rethreshold{$rfid})  ){
		printf STDERR "skipping [$gffs[0]] {[$cur[5]] > [$rethreshold{$rfid}]} [$cur]\n";
		next;
	    }
	}
	
	if ($cur[0] ne $prev[0]){#not the same sequence
	    print $prev . "\n";
	    ($prev,$prfid) = ($cur,$rfid); 
	    next;
	}
	
	$ext  = 0;
	$ext     = overlapExtent($cur[3], $cur[4], $prev[3], $prev[4]  );
	if ($ext>=$overlapThreshold){# FOUND AN OVERLAP
	    
	    my $sameStrand=1; 
	    $sameStrand = 0 if($cur[6] ne $prev[6]);
	    
	    if($sameStrand && ($cur[1] ne "Rfam" && $prev[1] ne "Rfam") && $cur[1] eq $prev[1]){#COMPETE IF PRODUCED BY THE SAME NON-RFAM METHOD 
		my $curp  = compete(\@cur, \@prev);
		$cur = join("\t", @{$curp});
		print "#COMPETE NON-RFAM [$cnt] [$cur[1]]vs[$prev[1]]\n" if(defined $verbose);
	    }
	    elsif($sameStrand  && ($cur[1] eq "Rfam") && ($prev[1] eq "tRNAscan-SE")){#tRNAscan-SE OUTPERFORMS RFAM (1)
		$cur = $prev;
		print "#COMPETE tRNAscan-SE >> RFAM (1) [$cnt]\n" if(defined $verbose);
	    }
	    elsif($sameStrand  && ($cur[1] eq "tRNAscan-SE") && ($prev[1] eq "Rfam")){#tRNAscan-SE OUTPERFORMS RFAM (2)
		$prev = $cur; 
		print "#COMPETE tRNAscan-SE >> RFAM (2) [$cnt]\n" if(defined $verbose);
	    }
	    elsif($sameStrand  && ($cur[1] eq "Rfam") && ($prev[2] eq "microRNA")){#miRBase OUTPERFORMS RFAM (1)
		$cur = $prev;
		print "#COMPETE miRBase >> RFAM (1) [$cnt] [$prev[8]]\n" if(defined $verbose);
	    }
	    elsif($sameStrand  && ($cur[2] eq "microRNA") && ($prev[1] eq "Rfam")){#miRBase OUTPERFORMS RFAM (2)
		$prev = $cur; 
		print "#COMPETE miRBase >> RFAM (2) [$cnt] [$cur[8]]\n" if(defined $verbose);
	    }
	    elsif($sameStrand  && defined($prfid) && defined($clans{$rfid}) && defined($clans{$prfid}) && ($clans{$prfid} eq $clans{$rfid})){#COMPETE IF IN THE SAME CLAN AND ON THE SAME STRAND
		my $curp  = compete(\@cur, \@prev);
		$cur = join("\t", @{$curp}); 
		print "#COMPETE A CLAN [$cnt] clan[$clans{$prfid}]:[$rfid]&[$prfid]\n" if(defined $verbose); 
	    }
	    elsif(!$sameStrand && defined($prfid) && defined($rfid) && ($prfid eq $rfid)){#compete if on the opposite strand and are the same family
		my $curp  = compete(\@cur, \@prev);
		$cur = join("\t", @{$curp});
		print "#COMPETE OPPOSITE STRAND OVERLAP FROM THE SAME FAMILY (miRNA)  [$cnt] [$rfid]vs[$prfid]\n" if(defined $verbose); 
	    }
	    else{#An unclassified overlap, print lowest scoring one to an overlaps file
		
		my $curp  = anticompete(\@cur, \@prev);
		my $overlow = join("\t", @{$curp}); 
		print "#UNCLASSIFIED OVERLAP [$prfid] != [$rfid] [$cnt]\ncur[$cur]\nprv[$prev]\n" if(defined $verbose);
		$curp  = compete(\@cur, \@prev);
		$cur = join("\t", @{$curp}); 
		print OV "#unclassified overlap:[$rfid][$cur[3]/$cur[3]-$cur[4]] vs [$prfid][$prev[3]/$prev[3]-$prev[4]]\n$overlow\n$cur\n" if(defined($prfid)); 
	    }
	    
	}
	else {
	    print $prev . "\n";
	}
    }
    $prev = $cur; 
    $prfid=$rfid;
    $cnt++; 
}
$fh->close;
close(OV);
 
print $prev . "\n";

exit(0);
######################################################################
#compete: return the winner (highest scoring) 
sub compete {
    my ($gff1, $gff2)=@_;
    
    if(isNumeric($gff1->[5]) && isNumeric($gff2->[5])){
	if ($gff1->[5] > $gff2->[5]){     
	    return $gff1;
	}
	else {
	    return $gff2;	
	}
    }
    elsif($gff1->[2] eq "microRNA" or $gff1->[1] eq "tRNAscan-SE"){#more specific methods trump Rfam
	return $gff1;
    }
    elsif($gff2->[2] eq "microRNA" or $gff2->[1] eq "tRNAscan-SE"){#more specific methods trump Rfam
	return $gff2;
    }
    elsif(isNumeric($gff1->[5])){
	return $gff1;
    }
    elsif(isNumeric($gff2->[5])){
	return $gff2;
    }
    else{
	return $gff1;
    }
}

######################################################################
#anticompete: return the loser (lowest scoring) 
sub anticompete {
    my ($gff1, $gff2)=@_;
    
    if(isNumeric($gff1->[5]) && isNumeric($gff2->[5])){
	if ($gff1->[5] < $gff2->[5]){     
	    return $gff1;
	}
	else {
	    return $gff2;	
	}
    }
    elsif(isNumeric($gff1->[5])){
	return $gff1;
    }
    elsif(isNumeric($gff2->[5])){
	return $gff2;
    }
    else{
	return $gff1;
    }

}
######################################################################
sub isNumeric {
    my $num = shift;
    if ($num=~/^-?\d+\.?\d*$/) { 
	return 1; 
    }
    else {
	return 0;
    }
}

######################################################################
#Returns the extent of overlap between two regions A=($x1, $y1) and B=($x2, $y2)
#using the following metric:
#
# D1 = 2*|A n B|/(|A|+|B|)
#
#An alternative, more permissive, metric is:
#
# D2 = |A n B|/min(|A|,|B|)
#
sub overlapExtent {
    my($x1, $y1, $x2, $y2) = @_;
    
    ($x1, $y1) = reorder($x1, $y1);
    ($x2, $y2) = reorder($x2, $y2);
    # 1.
    # x1                   y1
    # |<---------A--------->|
    #    |<------B------>|
    #    x2             y2
    #    XXXXXXXXXXXXXXXXX
    # 2.  x1                     y1
    #     |<---------A----------->|
    # |<-------------B------>|
    # x2                    y2
    #     XXXXXXXXXXXXXXXXXXXX
    # 3. x1             y1
    #    |<------A------>|
    # |<---------B--------->|
    # x2                   y2
    #    XXXXXXXXXXXXXXXXX
    # 4. x1                    y1
    #    |<-------------A------>|
    #        |<---------B----------->|
    #        x2                     y2
    #        XXXXXXXXXXXXXXXXXXXX
    my $D=0;
    my $int=0;
    my $L1=$y1-$x1+1;
    my $L2=$y2-$x2+1;
    my $minL = min($L1,$L2);
    if ( ($x1<=$x2 && $x2<=$y1) && ($x1<=$y2 && $y2<=$y1) ){    #1.
	$D = $L2;
    }
    elsif ( ($x2<=$x1) && ($x1<=$y2 && $y2<=$y1) ){             #2.
	$D = $y2-$x1+1;
    }
    elsif ( ($x2<=$x1 && $x1<=$y2) && ($x2<=$y1 && $y1<=$y2) ){ #3.
	$D = $L1;
    }
    elsif ( ($x1<=$x2 && $x2<=$y1) && ($y1<=$y2) ){             #4.
	$D = $y1-$x2+1;
    }
#D1:
    return 100*(2*$D)/($L1+$L2);
#D2:
#    return $D/$minL;
}

######################################################################
#reorder: given 2 integers, return the smallest first & the largest last:
sub reorder {
    my ($x,$y)=@_;
    
    if ($y<$x){
	my $tmp = $x;
	$x = $y;
	$y = $tmp;
    }
    return ($x,$y);
}

######################################################################
#Max and Min
#max
sub max {
  return $_[0] if @_ == 1;
  $_[0] > $_[1] ? $_[0] : $_[1]
}

#min
sub min {
  return $_[0] if @_ == 1;
  $_[0] < $_[1] ? $_[0] : $_[1]
}

######################################################################
sub help {
    print STDERR <<EOF;

compete_clans.pl: Given multiple gff files, resolve overlapping entries by:
                     a) selecting the entry with the highest bit-score if the families are in the same clan
                     b) selecting the entry with the *best* methods (e.g. tRNAscan >> Rfam)
                     c) opposite strand overlaps that are in the same family

Usage:   compete_clans.pl -g <gff1> -g <gff2> -g <gff3> ... 
Options:       -h|--help                     Show this help.

               -ot|--overlapthresh <val>     Merge features overlapping by <val>\% or more. default=[$overlapThreshold]
                                             NOTE: THIS IS A PERCENTAGE!
	       -g|--gff <file>               Give GFF file names as input. For multiple gff files use additional -g\'s.

EOF
}

