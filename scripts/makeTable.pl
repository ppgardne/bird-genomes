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
	elsif($w=~/miRNA/){
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
	
    }
    close(W);
    $rfam2type{'tRNA−pseudogene'}='Transfer RNA pseudogene';
}

my (@species, @families, %countsTot, %countsFams); 
if (-s "data/R/allRNA.dat"){
    open(W, "< data/R/allRNA.dat");
    while(my $w=<W>){
	chomp($w);
	my @w=split(/\t/, $w);
	my $w=pop(@w);
	
	if($w=~/human/){
	    @species=@w;
	}
	else{
	    
	    if(defined($rfam2type{$w})){
		
		$countsTot{"human"}{$rfam2type{$w}}=0 if(not defined($countsTot{"human"}{$rfam2type{$w}}));
		$countsTot{"human"}{$rfam2type{$w}}+=$w[0];
		
		$countsTot{"chicken"}{$rfam2type{$w}}=0 if(not defined($countsTot{"chicken"}{$rfam2type{$w}}));
		$countsTot{"chicken"}{$rfam2type{$w}}+=$w[7];
		
		$countsTot{"all-birds"}{$rfam2type{$w}}=0 if(not defined($countsTot{"all-birds"}{$rfam2type{$w}}));

		for (my $i=3; $i<scalar(@w); $i++){
		    $countsTot{"all-birds"}{$rfam2type{$w}}+=$w[$i];
		}
		
		
	    }
	    
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
SRP RNA,
Ribosomal RNA,
RNase P/MRP RNA');

print "human\tave(bird)\tchicken\tRNA-seq\tRNA type\n";

#foreach my $k (keys %{$countsTot{"all-birds"}}){
foreach my $k (@printNames){
     
    printf "%d\t%0.1f\t%d\t\t??\t$k\n", $countsTot{"human"}{$k}, $countsTot{"all-birds"}{$k}/48, $countsTot{"chicken"}{$k};
    
}



