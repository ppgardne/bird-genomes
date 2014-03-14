#!/usr/bin/perl 

use warnings;
use strict;

my %rfam2type;
if (-s "data/rfam2type.txt"){
    open(W, "< data/rfam2type.txt");    
    while(my $w=<W>){
	chomp($w);
	my @w = split(/\t/,$w);
	$rfam2type{$w[1]}=$w[2];
	if($w=~/Cis-reg/){
	    $rfam2type{$w[1]}='Cis-regulatory element';
	}
	elsif($w=~/HACA-box/){
	    $rfam2type{$w[1]}='H/ACA box snoRNA';
	}
	elsif($w=~/CD-box/){
	    $rfam2type{$w[1]}='C/D box snoRNA';
	}
	elsif($w=~/scaRNA/){
	    $rfam2type{$w[1]}='Small cajal body RNA';
	}
	elsif($w=~/rRNA/){
	    $rfam2type{$w[1]}='Ribosomal RNA';
	}
	elsif($w=~/miRNA/ or $w=~/mir-\d+/){
	    $rfam2type{$w[1]}='microRNA';
	}
	elsif($w=~/lncRNA/){
	    $rfam2type{$w[1]}='Long non-coding RNA';
	}
	elsif($w=~/tRNA/){
	    $rfam2type{$w[1]}='Transfer RNA';
	}
	elsif($w[1]=~/U1/ or $w[1]=~/U2/ or $w[1]=~/U4$/ or $w[1]=~/U5/ or $w[1]=~/U6$/){
	    $rfam2type{$w[1]}='Major spliceosomal RNA';
	}
	elsif($w[1]=~/U11/ or $w[1]=~/U12/ or $w[1]=~/U4atac$/ or $w[1]=~/U6atac$/){
	    $rfam2type{$w[1]}='Minor spliceosomal RNA';
	}
	elsif($w=~/RNase/){
	    $rfam2type{$w[1]}='RNase P/MRP RNA';
	}
	elsif($w=~/SRP/){
	    $rfam2type{$w[1]}='SRP RNA';
	}
	elsif($w=~/7SK/){
	    $rfam2type{$w[1]}='7SK RNA';
	}
	elsif($w=~/Y_RNA/){
	    $rfam2type{$w[1]}='Y RNA';
	}
	elsif($w=~/Telomerase/){
	    $rfam2type{$w[1]}='Telomerase RNA';
	}
	elsif($w=~/Vault/){
	    $rfam2type{$w[1]}='Vault RNA';
	}
#	else {
#	    print "UNCLASSIFIED: [$w]\n";
#	}
	
    }
    close(W);
    $rfam2type{'tRNA-pseudogene'}='Transfer RNA pseudogene';
    $rfam2type{'CD-snoRNA'}      ='C/D box snoRNA';
    $rfam2type{'HACA-snoRNA'}    ='H/ACA box snoRNA';
}

my (@species, @families, %countsTot, %countsSpecies, %countsArray, %countsFams, %expressedChicken); 
if (-s "data/R/allRNA.dat"){
    open(W, "< data/R/allRNA.dat");
    while(my $w=<W>){
	chomp($w);
	my @w=split(/\t/, $w);
	my $w=pop(@w);
	
	if($w=~/family/){
	    @species=@w;
	}
	else{
	    
	    if(not defined($rfam2type{$w}) and $w=~/^mir-\d+$/){
		$rfam2type{$w}='microRNA';
	    }
	    
	    if(defined($rfam2type{$w})){
		
		$countsTot{"human"}{$rfam2type{$w}}=0 if(not defined($countsTot{"human"}{$rfam2type{$w}}));
		$countsTot{"human"}{$rfam2type{$w}}+=$w[0];
		
		$countsTot{"chicken"}{$rfam2type{$w}}=0 if(not defined($countsTot{"chicken"}{$rfam2type{$w}}));
		$countsTot{"chicken"}{$rfam2type{$w}}+=$w[7];
		
		$expressedChicken{$rfam2type{$w}}+=expressedInChicken($w); 
		
		#print "$w:\t$rfam2type{$w}\n";
		$countsTot{"all-birds"}{$rfam2type{$w}}=0 if(not defined($countsTot{"all-birds"}{$rfam2type{$w}}));
		
		for (my $i=3; $i<scalar(@w); $i++){
		    $countsTot{"all-birds"}{$rfam2type{$w}}+=$w[$i];
		    #print "[$species[$i]]\t";
		    $countsSpecies{$species[$i]}{$rfam2type{$w}}+=$w[$i];
		}
		#print "\n";
	    }
#	    else {
#		print "What the hell is this? [$w] not in rfam2type\n";
#	    }
	}
    }
    close(W);
}

my @printNames = split(/,\n/,
'Long non-coding RNA,
microRNA,
C/D box snoRNA,
H/ACA box snoRNA,
Small cajal body RNA,
Major spliceosomal RNA,
Minor spliceosomal RNA,
Cis-regulatory element,
7SK RNA,
Telomerase RNA,
Vault RNA,
Y RNA,
Transfer RNA,
Transfer RNA pseudogene,
SRP RNA,
Ribosomal RNA,
RNase P/MRP RNA');

print '
\begin{tabular}{|r|r|r|r|l|}
\hline 
\multicolumn{5}{|l|}{{\bf ncRNA genes in human, chicken and all bird genomes}}\\\\
\hline 
Number in human & median(48 birds) & Number in chicken & Number in chicken  & RNA type\\\\
                &            &                   & confirmed with RNA-seq & \\\\
\hline
';

#foreach my $k (keys %{$countsTot{"all-birds"}}){
foreach my $k (@printNames){
    my @allBirds; 
    for (my $i=3; $i<scalar(@species); $i++){
	push(@allBirds, $countsSpecies{$species[$i]}{$k});
    }
#    printf "%d&%0.1f&%d&%d (%0.0f)&$k\\\\ \n", $countsTot{"human"}{$k}, $countsTot{"all-birds"}{$k}/48, $countsTot{"chicken"}{$k}, $expressedChicken{$k}, 100*$expressedChicken{$k}/$countsTot{"chicken"}{$k};    
    printf "%d&%0.1f&%d&%d (%0.1f\\%%)&$k\\\\ \n", $countsTot{"human"}{$k}, median(@allBirds), $countsTot{"chicken"}{$k}, $expressedChicken{$k}, 100*$expressedChicken{$k}/$countsTot{"chicken"}{$k};    
}

print '\hline
  \end{tabular}
';

######################################################################
sub expressedInChicken {
    
    my $rnaId=shift; 
    
    my $gff = "";
    if ($rnaId=~/^mir/ or $rnaId=~/^let/){
	$gff = `grep \'ID=$rnaId\_[0-9]\\+\$\' data/conserved-merged-annotations/Gallus_gallus.gff`;
    }
    elsif($rnaId eq "tRNA"){
	$gff = `grep "tRNAscan-SE" data/conserved-merged-annotations/Gallus_gallus.gff | grep -v Pseudo`; 
    }
    elsif($rnaId eq "tRNA-pseudogene"){
	$gff = `grep "tRNAscan-SE" data/conserved-merged-annotations/Gallus_gallus.gff | grep    Pseudo`; 
    }
    else {	
	$gff = `grep "Alias=$rnaId;Note" data/conserved-merged-annotations/Gallus_gallus.gff`; 
    }
    
    return 0 if (length($gff)==0);
    #print "gff:[$gff]\n";
    my @gff = split(/\n/, $gff);
    my $numExp=0;
    foreach my $g (@gff){
	
	#print "g:  [$g]\n";

	my @g = split(/\t/, $g);
	
	next if (scalar(@g) != 9); 
	my $ex=`grep "$g[0]_$g[3]" data/RNA-seq/*dat`; 
	
	$ex=~s/\,/\./g; #De-Germanify the numbers
	my @ex = split(/[\n;]/, $ex);
	
	my $sumEx=0;
	foreach my $ex (@ex){
	    
	    next if $ex=~/^data/;
	    $sumEx+=$ex; 
	}
	
	#print "$rnaId\t$g[0]/$g[3]-$g[4]\t$sumEx\n";
	$numExp++ if ($sumEx>10); 
	
    }    

    return $numExp; 
}

######################################################################
#http://www.perlmonks.org/?node_id=474564
sub median
{
    my @vals = sort {$a <=> $b} @_;
    my $len = @vals;
    if($len%2) #odd
    {
        return $vals[int($len/2)];
    }
    else #even
    {
        return ($vals[int($len/2)-1] + $vals[int($len/2)])/2;
    }
}
