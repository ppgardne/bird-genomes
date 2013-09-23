#!/usr/bin/perl 

#1. Filter contaminants
#2. Select families conserved in >10% of birds
##   --print to GFF
#3. Create Heatmaps
#4. Analyse repeated/pseudogene families
#     --"calibrate" with tRNAscan predictions
#     

#split the tRNAs out. 

use warnings;
use strict;
use Getopt::Long;
#use Statistics::Descriptive;


my ($gffsDir, $outDir, $rDir) = ("data/merged-annotations","data/conserved-merged-annotations", "data/R");

my ($verbose, $help);
&GetOptions( 
    "gd|gffDir=s"         => \$gffsDir,
    "od|outDir"           => \$outDir,
    "rd|rDir"             => \$rDir,
    "v|verbose"           => \$verbose,
    "h|help"              => \$help
    );

my @gffs = glob("$gffsDir/*gff"); 

if( $help ) {
    &help();
    exit(1);
}
elsif ((not defined($gffsDir)) or (@gffs)==0){
    print "FATAL: no gff files given\n";
    &help();
    exit(1);
}

my %whitelist; 
if (-s "$rDir/allRNA.dat"){
    open(W, "< $rDir/allRNA.dat");
    while(my $w=<W>){
	chomp($w);
	my @w=split(/\t/, $w);
	my $w=pop(@w);
	$whitelist{$w}=1;
    }
    close(W);
}
$whitelist{'SNORD34'}=1;
#replace with a family2type hash:
my %lncRNAlist;
if (-s "data/rfam11_lncRNAs.txt"){
    open(W, "< data/rfam11_lncRNAs.txt");    
    while(my $w=<W>){
	chomp($w);
	$lncRNAlist{$w}=1;
    }
    close(W);
}

my %blacklist=(
    mraW=>1,    
    );

print "Reading GFFs\n" if (defined($verbose));
my (%HoHspeciesRfamCount, %familiesSpeciesCounts, %familiesTotalCounts, %statistics, @bitscores, %tRNA, %tRNAfamilies); 
foreach my $f (@gffs){

    print "$f\n" if(defined($verbose));
    my $species;
    if ($f=~/$gffsDir\/(\S+?_\S+?)\.gff/){
	$species=$1; 
    }
    else{
	print "WARNING: no species name parsed from [$f], skipping!\n";
	next;
    }
    
    open(GFF, "< $f");
    open(OUT, "> $outDir/$species\.gff");
    while(my $g=<GFF>){
	chomp($g);
	my @g=split(/\t/, $g);
	my $family; 
	if($g[8]=~/rfam-id=(\S+);evalue/ or $g[8]=~/Alias=(\S+);Note/){#Rfam
	    $family=$1;
	    #push(@bitscores, $g[5]) if (defined($whitelist{$family}) && $g[5]=~/^\d+\.\d+$/); 
	}
	elsif($g[8]=~/ID=(\S+)\_\d+/){#miRBase
	    $family=$1;
	}
	elsif($g[8]=~/type=Pseudo;/){#tRNAscan
	    $family="tRNA-pseudogene";
	    $tRNA{$species}{"Pseudo"}=0 if(not defined($tRNA{$species}{"Pseudo"}));
	    $tRNA{$species}{"Pseudo"}++;
	    $tRNAfamilies{"Pseudo"}=0 if(not defined($tRNAfamilies{"Pseudo"}));
	    $tRNAfamilies{"Pseudo"}++;
	    #print "$g[8]\n";
	}
	elsif($g[8]=~/type=(\S+?);/){#tRNAscan
	    $family="tRNA";
	    $tRNA{$species}{$1}=0 if(not defined($tRNA{$species}{$1}));
	    $tRNA{$species}{$1}++;
	    $tRNAfamilies{$1}=0 if(not defined($tRNAfamilies{$1}));
	    $tRNAfamilies{$1}++;
	    #print "$g[8]\n";
	}
	else{#WTF?
	    $family=$g[2];
	    #print "$g[2]\n";
	}
	
	print "$species\t$family\n" if defined($verbose);
	if ( not defined($HoHspeciesRfamCount{$species}{$family}) ){
	    $HoHspeciesRfamCount{$species}{$family}=0;
	    $familiesSpeciesCounts{$family}=0 if (not defined($familiesSpeciesCounts{$family}));
	    $familiesSpeciesCounts{$family}++;
	}
	$familiesTotalCounts{$family}=0 if (not defined($familiesTotalCounts{$family}));
	$familiesTotalCounts{$family}++;
	$HoHspeciesRfamCount{$species}{$family}++;
	print OUT "$g\n" if defined($whitelist{$family});
    }
    close(GFF);
    close(OUT);
}

#slurp species in phylogenetic order:
my (@speciesPhyloOrder,@speciesCommonPhyloOrder); 
open(SPC, "< data/species_list-phyloOrder.txt");
while(my $s=<SPC>){
    if($s=~/^(\S+)\t(\S+)/){
	push(@speciesPhyloOrder,$1); 
	push(@speciesCommonPhyloOrder,$2); 
    }
}
close(SPC); 

print "Printing R-dat headers\n" if (defined($verbose));
foreach my $ut (("$rDir/snoRNA.dat", "$rDir/miRNA.dat", "$rDir/RNA.dat", "$rDir/allRNA.dat", "$rDir/tRNA.dat", "$rDir/lncRNA.dat")){
    open(UT, "> $ut");
    my $cnt=0;
    foreach my $species (@speciesCommonPhyloOrder){#@speciesPhyloOrder ) {
	print UT "$species\t";
	#print "$species\t$statistics{$speciesPhyloOrder[$cnt]}" if ($ut=~/all/);
	$cnt++;
    }
    print UT "family\n";
    close(UT);
}

my $cnt=0;

print "Sorting family IDs\n" if (defined($verbose));
my @family = sort {
    if($a=~/^mir\-\d+$/ && $b=~/^mir\-\d+$/){
	(my $num_a = $a) =~ s/^mir\-(\d+)$/$1/;
	(my $num_b = $b) =~ s/^mir\-(\d+)$/$1/;
	return $num_a <=> $num_b;
    }
    if($a=~/^mir\-\d+\_\d+$/ && $b=~/^mir\-\d+\_\d+$/){
	(my $num_a = $a) =~ s/^mir\-(\d+)\_\d+$/$1/;
	(my $num_b = $b) =~ s/^mir\-(\d+)\_\d+$/$1/;
	return $num_a <=> $num_b;
    }
    elsif($a=~/^MIR\d+$/ && $b=~/^MIR\d+$/){
	(my $num_a = $a) =~ s/^MIR(\d+)$/$1/;
	(my $num_b = $b) =~ s/^MIR(\d+)$/$1/;
	return $num_a <=> $num_b;
    }
    elsif($a=~/^MIR\d+\_\d+$/ && $b=~/^MIR\d+\_\d+$/){
	(my $num_a = $a) =~ s/^MIR(\d+)\_\d+/$1/;
	(my $num_b = $b) =~ s/^MIR(\d+)\_\d+/$1/;
	return $num_a <=> $num_b;
    }
    elsif($a=~/^SNORA\d+$/ && $b=~/^SNORA\d+$/){
	(my $num_a = $a) =~ s/^SNORA(\d+)$/$1/;
	(my $num_b = $b) =~ s/^SNORA(\d+)$/$1/;
	return $num_a <=> $num_b;
    }
    elsif($a=~/^SNORD\d+$/ && $b=~/^SNORD\d+$/){
	(my $num_a = $a) =~ s/^SNORD(\d+)$/$1/;
	(my $num_b = $b) =~ s/^SNORD(\d+)$/$1/;
	return $num_a <=> $num_b;
    }
    elsif($a=~/^SCARNA\d+$/ && $b=~/^SCARNA\d+$/){
    	(my $num_a = $a) =~ s/^SCARNA(\d+)$/$1/;
    	(my $num_b = $b) =~ s/^SCARNA(\d+)$/$1/;
    	return $num_a <=> $num_b;
    }
    elsif($a=~/^snoZ\d+$/ && $b=~/^snoZ\d+$/){
    	(my $num_a = $a) =~ s/^snoZ(\d+)$/$1/;
    	(my $num_b = $b) =~ s/^snoZ(\d+)$/$1/;
    	return $num_a <=> $num_b;
    }
    elsif($a=~/^snoU\d+$/ && $b=~/^snoU\d+$/){
    	(my $num_a = $a) =~ s/^snoU(\d+)$/$1/;
    	(my $num_b = $b) =~ s/^snoU(\d+)$/$1/;
    	return $num_a <=> $num_b;
    }
    elsif($a=~/^U\d+/ && $b=~/^U\d+/){
    	(my $num_a = $a) =~ s/^U(\d+)\w*/$1/;
    	(my $num_b = $b) =~ s/^U(\d+)\w*/$1/;
    	return $num_a <=> $num_b;
    }
    else { 
	return lc($a) cmp lc($b); 
    } 
} (keys %familiesTotalCounts);  #{ $HoHspeciesRfamCount{$speciesPhyloOrder[0]} } );

my @tRNAfamily = sort {
    lc($a) cmp lc($b); 
}(keys %tRNAfamilies);

print "Printing R-dat files\n" if (defined($verbose));
my %exceptions = (SNORD34 => 1);
foreach my $family ( @family ) {
    
    #print "$family\n";# if ($family=~/Pseudo/i);
    
    #Family must be found in >10% of all species
    #printf "$family: $familiesSpeciesCounts{$family}/%d %0.2f\n", scalar(@speciesPhyloOrder), $familiesSpeciesCounts{$family}/scalar(@family) if($family=~/RNase_MRP/);
    next if( ($familiesSpeciesCounts{$family}/scalar(@speciesPhyloOrder)) < 0.1);
    next if(defined($blacklist{$family}));
    open(AL, ">> $rDir/allRNA.dat");
    if ($family=~/^sno/i or $family=~/^SCA/ or $family=~/^ACA/){
	next if( (($familiesSpeciesCounts{$family}/scalar(@speciesPhyloOrder)) < 0.10 or $familiesTotalCounts{$family}<5) and not defined($exceptions{$family}) );
	open(UT, ">> $rDir/snoRNA.dat");
    }
    elsif ($family=~/^mir/i or $family=~/^let-7/ or $family=~/^lin-4/ ){
	next if( ($familiesSpeciesCounts{$family}/scalar(@speciesPhyloOrder)) < 0.10 or $familiesTotalCounts{$family}<5);
	open(UT, ">> $rDir/miRNA.dat");
    }
    elsif(defined($lncRNAlist{$family})){
	open(UT, ">> $rDir/lncRNA.dat");
    }
    else{
	open(UT, ">> $rDir/RNA.dat");
    }
    
    foreach my $species ( @speciesPhyloOrder ) {
	if(defined($HoHspeciesRfamCount{$species}{$family})){
	    print UT "$HoHspeciesRfamCount{$species}{$family}\t";
	    print AL "$HoHspeciesRfamCount{$species}{$family}\t";
	}
	else {
	    print UT "0\t";
	    print AL "0\t";
	}
    }
    print UT "$family\n";
    print AL "$family\n";
    close(UT);
    close(AL);
    $cnt++;
}

print "Printing tRNA R-dat files\n" if (defined($verbose));
open(TR, ">> $rDir/tRNA.dat");
foreach my $family ( @tRNAfamily ) {
    foreach my $species ( @speciesPhyloOrder ) {

	if(defined($tRNA{$species}{$family})){
	    print TR "$tRNA{$species}{$family}\t";
	}
	else {
	    print TR "0\t";
	}
    }
    print TR "$family\n";
}
close(TR);

print "System calls: creating [$rDir/snoRNA-human-yeast-correspondences.dat] & running R script\n" if (defined($verbose));

#egrep '^snoR38;|SNORA13;|SNORA16;|SNORA2;|SNORA21;|SNORA26;|SNORA27;|SNORA28;|SNORA3;|SNORA36;|SNORA4;|SNORA44;|SNORA48;|SNORA5;|SNORA50;|SNORA52;|SNORA56;|SNORA58;|SNORA62;|SNORA64;|SNORA65;|SNORA66;|SNORA69;|SNORA7;|SNORA74;|SNORA76;|SNORA8;|SNORA9;|SNORD12;|SNORD14;|SNORD15;|SNORD16;|SNORD17;|SNORD18;|SNORD2;|SNORD24;|SNORD27;|SNORD29;|SNORD31;|SNORD33;|SNORD34;|SNORD35;|SNORD36;|SNORD38;|SNORD41;|SNORD43;|SNORD46;|SNORD51;|SNORD52;|SNORD57;|SNORD59;|SNORD60;|SNORD62;|SNORD65;|SNORD74;|SNORD77;|SNORD88;|SNORND104;|snosnR60_Z15' clans_competed/*gff | perl -lane 'if(/\/(\S+?)\-.*gff:(\S+).*\-id=(\S+);eval/ or /\/(\S+?)\-.*gff:(\S+).*Alias=(\S+);Not/){print "$1\t$2\t$3\t$F[6]"}'
#GAS5 snoRNAs: \|SNORD81\$\|SNORD47\$\|SNORD80\$\|SNORD79\$\|SNORD78\$\|SNORD44\$\|SNORD77\$\|SNORD76\$\|SNORD75\$\|SNORD74\$
#DAMN THIS IS NASTY!
system("egrep \47^human\|snoR38\$\|SNORA13\$\|SNORA16\$\|SNORA2\$\|SNORA21\$\|SNORA26\$\|SNORA27\$\|SNORA28\$\|SNORA3\$\|SNORA36\$\|SNORA4\$\|SNORA44\$\|SNORA48\$\|SNORA5\$\|SNORA50\$\|SNORA52\$\|SNORA56\$\|SNORA58\$\|SNORA62\$\|SNORA64\$\|SNORA65\$\|SNORA66\$\|SNORA69\$\|SNORA7\$\|SNORA74\$\|SNORA76\$\|SNORA8\$\|SNORA9\$\|SNORD12\$\|SNORD14\$\|SNORD15\$\|SNORD16\$\|SNORD17\$\|SNORD18\$\|SNORD2\$\|SNORD24\$\|SNORD27\$\|SNORD29\$\|SNORD31\$\|SNORD33\$\|SNORD34\$\|SNORD35\$\|SNORD36\$\|SNORD38\$\|SNORD41\$\|SNORD43\$\|SNORD46\$\|SNORD51\$\|SNORD52\$\|SNORD57\$\|SNORD59\$\|SNORD60\$\|SNORD62\$\|SNORD65\$\|SNORD74\$\|SNORD77\$\|SNORD88\$\|SNORND104\$\|snosnR60_Z15\$\47 $rDir/snoRNA.dat > $rDir/snoRNA-human-yeast-correspondences.dat");
#$rDir/allRNA.dat
#egrep -i '^human|HOXA11-AS1|PCA3|RMST|Six3os1|SOX2OT|ST7-OT3|HOTAIRM1|HOTTIP|RMST|DLEU2|NBR2|SNORD93|miR-15$|miR-16$' ../data/R/allRNA.dat

system("R CMD BATCH --no-save scripts/heatmaps.R");

print "Finished!\n" if (defined($verbose));

exit(0);

######################################################################
sub sortFamilies {
    #my ($a,$b)=@_;
    if($a=~/^mir\-\d+/ && $b=~/^mir\-\d+/){
	(my $num_a = $a) =~ s/^mir\-(\d+)/$1/;
	(my $num_b = $b) =~ s/^mir\-(\d+)/$1/;
	return $num_a <=> $num_b;
    }
    elsif($a=~/^SNORA\d+/ && $b=~/^SNORA\d+/){
	(my $num_a = $a) =~ s/^SNORA(\d+)/$1/;
	(my $num_b = $b) =~ s/^SNORA(\d+)/$1/;
	return $num_a <=> $num_b;
    }
    elsif($a=~/^SNORD\d+/ && $b=~/^SNORD\d+/){
	(my $num_a = $a) =~ s/^SNORD(\d+)/$1/;
	(my $num_b = $b) =~ s/^SNORD(\d+)/$1/;
	return $num_a <=> $num_b;
    }
    elsif(/^SCARNA\d+/){
    	(my $num_a = $a) =~ s/^SCARNA(\d+)/$1/;
    	(my $num_b = $b) =~ s/^SCARNA(\d+)/$1/;
    	return $num_a <=> $num_b;
    }
    elsif(/^SNORD\d+/){
    	(my $num_a = $a) =~ s/^SNORD(\d+)/$1/;
    	(my $num_b = $b) =~ s/^SNORD(\d+)/$1/;
    	return $num_a <=> $num_b;
    }
    elsif(/^snoZ\d+/){
    	(my $num_a = $a) =~ s/^snoZ(\d+)/$1/;
    	(my $num_b = $b) =~ s/^snoZ(\d+)/$1/;
    	return $num_a <=> $num_b;
    }
    elsif(/^snoU\d+/){
    	(my $num_a = $a) =~ s/^snoU(\d+)/$1/;
    	(my $num_b = $b) =~ s/^snoU(\d+)/$1/;
    	return $num_a <=> $num_b;
    }
    elsif(/^U\d+/){
    	(my $num_a = $a) =~ s/^U(\d+)/$1/;
    	(my $num_b = $b) =~ s/^U(\d+)/$1/;
    	return $num_a <=> $num_b;
    }
    else {
    	return lc $a cmp lc $b;
    }
}


######################################################################
sub help {
    print STDERR <<EOF;

gffs2heatmaps.pl: 1. Select families conserved in >10% of birds
                   --print to GFF & to R "dat" files
                  2. Create Heatmaps with R-scripts

Usage:   gffs2heatmaps.pl 
Options:       -h|--help                     Show this help.
               -v|--verbose                  Print lots of stuff.

               -gd|--gffDir <dir>            Directory contain gffs of ncRNA annotation [default:data/merged-annotations]
               -od|--outDir <dir>            Output directory, print gffs of conserved ncRNAs there [default:data/conserved-merged-annotations]
	       -rd|--rDir   <dir>            Output directory, print R-data files there [default:data/R]
	       
EOF
}
