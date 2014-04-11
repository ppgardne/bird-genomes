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

my (@species, @families, %countsTot, %countsSpecies, %countsConditions, %countsFams, %expressedChickenPaul, %expressedChickenMario, %expressedChickenPaulRand, %expressedChickenMarioRand); 
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
		$countsTot{"chicken"}{$rfam2type{$w}}+=$w[5];
		
		if(not defined($expressedChickenPaul{$rfam2type{$w}})){
		    $expressedChickenPaul{$rfam2type{$w}}=0;
		    $expressedChickenMario{$rfam2type{$w}}=0;
		    $expressedChickenPaulRand{$rfam2type{$w}}=0;
		    $expressedChickenMarioRand{$rfam2type{$w}}=0;
		}
		my ($paul, $mario, $paulRand, $marioRand)=expressedInChicken($w); 
		$expressedChickenPaul{$rfam2type{$w}}     +=$paul;
		$expressedChickenMario{$rfam2type{$w}}+=$mario if (isNumeric($mario));
		$expressedChickenPaulRand{$rfam2type{$w}}     +=$paulRand;
		$expressedChickenMarioRand{$rfam2type{$w}}+=$marioRand if (isNumeric($marioRand));
		#print "mario:[$mario] paul:[$paul]\n"; 
		
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
                &                  &                   & Chicken ncRNAs & \\\\
                &                  &                   & confirmed with RNA-seq & \\\\
\hline
Number in human & median(48 birds) & Number in chicken & max(RNA$_i$)$>13.0$ & RNA type\\\\
\hline
';

#foreach my $k (keys %{$countsTot{"all-birds"}}){
my ($totalHuman, $totalMedianBird, $totalChicken, $totalExpressedChickenPaul, $totalExpressedChickenMario, $totalExpressedChickenPaulRand, $totalExpressedChickenMarioRand )=(0,0,0,0,0,0,0); 
foreach my $k (@printNames){
    my @allBirds; 
    for (my $i=3; $i<scalar(@species); $i++){
	push(@allBirds, $countsSpecies{$species[$i]}{$k});
    }
#    printf "%d&%0.1f&%d&%d (%0.0f)&$k\\\\ \n", $countsTot{"human"}{$k}, $countsTot{"all-birds"}{$k}/48, $countsTot{"chicken"}{$k}, $expressedChicken{$k}, 100*$expressedChicken{$k}/$countsTot{"chicken"}{$k};
    my $medBirds=median(@allBirds);
    $countsTot{"chicken"}{$k} = 1 if (not defined($countsTot{"chicken"}{$k}) or $countsTot{"chicken"}{$k} == 0);
    printf "%d&%0.1f&%d&%d (%0.1f\\%%) &$k\\\\ \n", $countsTot{"human"}{$k}, $medBirds, $countsTot{"chicken"}{$k}, 
    $expressedChickenMario{$k}, 100*$expressedChickenMario{$k}/$countsTot{"chicken"}{$k};

    $totalHuman+=$countsTot{"human"}{$k};
    $totalMedianBird+=$medBirds; 
    $totalChicken+=$countsTot{"chicken"}{$k}; 
    $totalExpressedChickenPaul+=$expressedChickenPaul{$k}; 
    $totalExpressedChickenMario+=$expressedChickenMario{$k}; 
    $totalExpressedChickenPaulRand+=$expressedChickenPaulRand{$k}; 
    $totalExpressedChickenMarioRand+=$expressedChickenMarioRand{$k}; 
}
print '\hline' . "\n";
printf "%d&%0.1f&%d&%d (%0.1f\\%%) &Total\\\\ \n", $totalHuman, $totalMedianBird, $totalChicken, 
    $totalExpressedChickenMario, 100*$totalExpressedChickenMario/$totalChicken;
print '\hline
\end{tabular}
';

printf "\n\n\n\n totalExpressed[$totalExpressedChickenMario] totalNegControlsExpressed[$totalExpressedChickenMarioRand] totalGenes[$totalChicken] FPR[%0.1f]\n", 100*$totalExpressedChickenMarioRand/$totalChicken;

# print '
# \begin{tabular}{|r|r|r|r|r|l|}
# \hline 
# \multicolumn{5}{|l|}{{\bf ncRNA genes in human, chicken and all bird genomes}}\\\\
# \hline 
#                 &                  &                   & \multicolumn{2}{|l|}{Chicken ncRNAs confirmed with RNA-seq} & \\\\
# \hline
# Number in human & median(48 birds) & Number in chicken & sum(RNA$_i$)$>10.0$      & max(RNA$_i$)$>5.0$ & RNA type\\\\
# \hline
# ';

# #foreach my $k (keys %{$countsTot{"all-birds"}}){
# my ($totalHuman, $totalMedianBird, $totalChicken, $totalExpressedChickenPaul, $totalExpressedChickenMario, $totalExpressedChickenPaulRand, $totalExpressedChickenMarioRand )=(0,0,0,0,0,0,0); 
# foreach my $k (@printNames){
#     my @allBirds; 
#     for (my $i=3; $i<scalar(@species); $i++){
# 	push(@allBirds, $countsSpecies{$species[$i]}{$k});
#     }
# #    printf "%d&%0.1f&%d&%d (%0.0f)&$k\\\\ \n", $countsTot{"human"}{$k}, $countsTot{"all-birds"}{$k}/48, $countsTot{"chicken"}{$k}, $expressedChicken{$k}, 100*$expressedChicken{$k}/$countsTot{"chicken"}{$k};
#     my $medBirds=median(@allBirds);
#     printf "%d&%0.1f&%d&%d (%0.1f\\%%) [FPR:%2.1f\\%%]&%d (%0.1f\\%%) [FPR:%2.1f\\%%]&$k\\\\ \n", $countsTot{"human"}{$k}, $medBirds, $countsTot{"chicken"}{$k}, 
#     $expressedChickenPaul{$k},  100*$expressedChickenPaul{$k}/$countsTot{"chicken"}{$k},  100*$expressedChickenPaulRand{$k}/$countsTot{"chicken"}{$k}, 
#     $expressedChickenMario{$k}, 100*$expressedChickenMario{$k}/$countsTot{"chicken"}{$k}, 100*$expressedChickenMarioRand{$k}/$countsTot{"chicken"}{$k};    

#     $totalHuman+=$countsTot{"human"}{$k};
#     $totalMedianBird+=$medBirds; 
#     $totalChicken+=$countsTot{"chicken"}{$k}; 
#     $totalExpressedChickenPaul+=$expressedChickenPaul{$k}; 
#     $totalExpressedChickenMario+=$expressedChickenMario{$k}; 
#     $totalExpressedChickenPaulRand+=$expressedChickenPaulRand{$k}; 
#     $totalExpressedChickenMarioRand+=$expressedChickenMarioRand{$k}; 
# }
# print '\hline' . "\n";
# printf "%d&%0.1f&%d&%d (%0.1f\\%%) [FPR:%2.1f\\%%]&%d (%0.1f\\%%) [FPR:%2.1f\\%%]&Total\\\\ \n", $totalHuman, $totalMedianBird, $totalChicken, 
#     $totalExpressedChickenPaul,  100*$totalExpressedChickenPaul/$totalChicken,  100*$totalExpressedChickenPaulRand/$totalChicken,
#     $totalExpressedChickenMario, 100*$totalExpressedChickenMario/$totalChicken, 100*$totalExpressedChickenMarioRand/$totalChicken;
# print '\hline
# \end{tabular}
# ';

######################################################################
sub expressedInChicken {
    
    my $rnaId=shift; 
    
    my ($gff,$gffRand) = ("","");
    if ($rnaId=~/^mir/ or $rnaId=~/^let/){
	$gff     = `grep \'ID=$rnaId\_[0-9]\\+\$\' data/conserved-merged-annotations/Gallus_gallus.gff`;
	$gffRand = `grep \'ID=$rnaId\_[0-9]\\+\$\' data/RNA-seq/Gallus_gallus.randomized.gff`;
    }
    elsif($rnaId eq "tRNA"){
	$gff     = `grep "tRNAscan-SE" data/conserved-merged-annotations/Gallus_gallus.gff | grep -v Pseudo`; 
	$gffRand = `grep "tRNAscan-SE" data/RNA-seq/Gallus_gallus.randomized.gff           | grep -v Pseudo`; 
    }
    elsif($rnaId eq "tRNA-pseudogene"){
	$gff     = `grep "tRNAscan-SE" data/conserved-merged-annotations/Gallus_gallus.gff | grep    Pseudo`; 
	$gffRand = `grep "tRNAscan-SE" data/RNA-seq/Gallus_gallus.randomized.gff           | grep    Pseudo`; 
    }
    else {	
	$gff     = `grep "Alias=$rnaId;Note" data/conserved-merged-annotations/Gallus_gallus.gff`; 
	$gffRand = `grep "Alias=$rnaId;Note" data/RNA-seq/Gallus_gallus.randomized.gff`; 
    }
    
    return (0,0,0,0) if (length($gff)==0);
    #print "gff:[$gff]\n";
    
    my ($numExp,    $numExpMario)    =numExpressed($gff,     0);
    my ($numExpRand,$numExpMarioRand)=numExpressed($gffRand, 1);

    return ($numExp,$numExpMario,$numExpRand,$numExpMarioRand); 
}

######################################################################
#numExpressed: input is a string of gff formatted coordinates, return the number of "expressed" regions
sub numExpressed {
    
    my ($gff, $randomDBs) = @_; 
    my @gff = split(/\n/, $gff);
    my ($numExp,$numExpMario) =(0,0);
    my @expArray=(); 
    my @databases = qw(data/RNA-seq/mccarthy_expression_tissue.dat            data/RNA-seq/ulitsky_expression_tissue.dat); 
    my @thresholds = (5, 12); 
    @databases = qw(data/RNA-seq/mccarthy_expression_tissue.randomized.dat data/RNA-seq/ulitsky_expression_tissue.randomized.dat) if ($randomDBs);
    
    
    
    foreach my $g (@gff){
	my $seen=0;
	for (my $i=0; $i<2; $i++){
	    my $databases = $databases[$i];
	    
	    #print "g:  [$g]\n";
	    
	    my @g = split(/\t/, $g);
	    
	    next if (scalar(@g) != 9); 
	    my $ex=`grep "$g[0]_$g[3]" $databases`;
	    
	    $ex=~s/\,/\./g; #De-Germanify the numbers
	    my @ex = split(/[\n;]/, $ex);
	    
	    @expArray=(); 
	    my $sumEx=0;
	    foreach my $ex (@ex){
		
		#next if $ex=~/^data/;
		next if (not isNumeric($ex)); 
		
		$sumEx+=$ex;
		push(@expArray, $ex);
	    }
	
	
	    my $maxExp = 0;
	    $maxExp = maxA(@expArray) if (scalar(@expArray)); 
	    
	    #print "$rnaId\t$g[0]/$g[3]-$g[4]\t$sumEx\n";
	    $numExp++      if ($sumEx>10); 
	    if ($maxExp>$thresholds[$i] and not $seen){
		$numExpMario++;
		$seen++;
	    }
	}
	#print "$rnaId: mario [cnt [$numExpMario] ($maxExp>5)] paul [cnt [$numExp] ($sumEx>10)]   [$g[0]_$g[3]]\n";
	
	
    }
    
    return ($numExp,$numExpMario); 
    
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
#Max and Min for arrays:
#max
sub maxA {
    my $max = $_[0];
    foreach my $a (@_){
	$max = max($max, $a) if isNumeric($a);
    }
    return $max;
}

#min
sub minA {
    my $min = $_[0];
    foreach my $a (@_){
	$min = min($min, $a) if isNumeric($a);
    }
    return $min;
}
######################################################################
sub isNumeric {
    my $num = shift;
    
    return 0 if (not defined($num));
    if ($num=~/^-?\d+\.?\d*$/) { 
	return 1; 
    }
    else {
	return 0;
    }
}
