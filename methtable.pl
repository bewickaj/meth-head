#!/usr/bin/perl

use strict;
use warnings;
use Math::CDF;
use Statistics::Multtest qw(:all);
use Getopt::Long;

my $intersect;
my $out;
my $fdr;
my $help;
my %meth;
my %t;
my $summCG;
my $sumtCG;
my $pCG;
my $summCHG;
my $sumtCHG;
my $pCHG;
my $summCHH;
my $sumtCHH;
my $pCHH;
my %hoa;

if (@ARGV > 0) {
	GetOptions (
		'intersect_infile:s' => \$intersect,
		'outfile:s' => \$out,
		'fdr:s' => \$fdr,
		'help' => \$help
	);
}
if ($help) {
	print "Usage: --intersect_infile \"name of outfile from bedtools intersect\" --outfile \"file name of outfile to write to\" --fdr \"yes/no to perform FDR adjustment on p-values\"\n";
	exit;
}

open INTERSECT, "<$intersect" or die;
open OUT, ">$out" or die;

warn "Storing intersect outfile from bedtools, and calculating probability threshold across all CG, CHG, and CHH sites for binomial test...\n";

while (<INTERSECT>) {
	chomp;
	my @array = split "\t";
	my $con = $array[3];
	my $meth = $array[4];
	my $t = $array[5];
	my $call = $array[6];
	my $gm = $array[10];
	if ($t >= 3) {
		$t{$gm}{$con}++;
		$meth{$gm}{$con} += $call;
		if ($con eq "CG") {
			$summCG += $call; 
			$sumtCG++;
			$pCG = $summCG/$sumtCG;
		}
		elsif ($con eq "CHG") {
			$summCHG += $call; 
			$sumtCHG++;
			$pCHG = $summCHG/$sumtCHG;
		}
		elsif ($con eq "CHH") {
			$summCHH += $call; 
			$sumtCHH++;
			$pCHH = $summCHH/$sumtCHH;
		}
	}	
	else {
		next;
	}
}

warn "Checking to see if every methylation context is present in each gene model. If absent, a 0 will be inserted...\n";

my @con = ("CG", "CHG", "CHH");

foreach my $gm (keys %t) {
	foreach (@con) {
		if (exists $t{$gm}{$_}) {
			next;
		}
		else {
			$t{$gm}{$_} = "0";
		}
	}
}

foreach my $gm (keys %t) {
	foreach (@con) {
		if (exists $meth{$gm}{$_}) {
			next;
		}
		else {
			$meth{$gm}{$_} = "0";
		}
	}
}
	

warn "Estimating probaility of per methylation context state per gene model, and writing to outfile...\n";

if ($fdr eq "yes" || "no") {

	print OUT "gene_model\ttotal_CG\tmeth_CG\tpvalue_CG\ttotal_CHG\tmeth_CHG\tpvalue_CHG\ttotal_CHH\tmeth_CHH\tpvalue_CHH\n";

	foreach my $gm (keys %t) {
		print OUT "$gm\t";
		foreach my $con (sort keys %{$t{$gm}}) {
	    	if (exists $meth{$gm}{$con}) {
				if (($meth{$gm}{$con} == 0 or $t{$gm}{$con} == 0) && $con eq "CG") {
					print OUT "$t{$gm}{$con}\t", "$meth{$gm}{$con}\t", "1", "\t";
				}
				elsif ($meth{$gm}{$con} > 0 && $con eq "CG") {
					my $prob = 1-&Math::CDF::pbinom($meth{$gm}{$con}-1, $t{$gm}{$con}, $pCG);
					print OUT "$t{$gm}{$con}\t", "$meth{$gm}{$con}\t", $prob, "\t";
				}
				if (($meth{$gm}{$con} == 0 or $t{$gm}{$con} == 0) && $con eq "CHG") {
					print OUT "$t{$gm}{$con}\t", "$meth{$gm}{$con}\t", "1", "\t";
				}
				elsif ($meth{$gm}{$con} > 0 && $con eq "CHG") {
					my $prob = 1-&Math::CDF::pbinom($meth{$gm}{$con}-1, $t{$gm}{$con}, $pCHG);
					print OUT "$t{$gm}{$con}\t", "$meth{$gm}{$con}\t", $prob, "\t";
				}
				if (($meth{$gm}{$con} == 0 or $t{$gm}{$con} == 0) && $con eq "CHH") {
					print OUT "$t{$gm}{$con}\t", "$meth{$gm}{$con}\t", "1";
				}
				elsif ($meth{$gm}{$con} > 0 && $con eq "CHH") {
					my $prob = 1-&Math::CDF::pbinom($meth{$gm}{$con}-1, $t{$gm}{$con}, $pCHH);
					print OUT "$t{$gm}{$con}\t", "$meth{$gm}{$con}\t", $prob;
				}
			}
		}
		print OUT "\n";
	}
}

close OUT;

if ($fdr eq "yes") {

	warn "Printing additional file with FDR adjusted p-values...\n";
	
	open OUT, "<$out" or die;

	my ($file, undef) = split /\.[\w]+$/, $out;

	open FDR_OUT, ">$file\_fdr.txt" or die;

	my $CG_href;
	my $CHG_href;
	my $CHH_href;

	while (<OUT>) {
		chomp;
		next if $. == 1;
		my @array = split "\t";
		$CG_href->{$_} = $array[3];
		$CHG_href->{$_} = $array[6];
		$CHH_href->{$_} = $array[9];
	}

	my $CG_res = BH($CG_href);
	my $CHG_res = BH($CHG_href);
	my $CHH_res = BH($CHH_href);

	print FDR_OUT "gene_model\ttotal_CG\tmeth_CG\tpvalue_CG\tfdr_CG\ttotal_CHG\tmeth_CHG\tpvalue_CHG\tfdr_CHG\ttotal_CHH\tmeth_CHH\tpvalue_CHH\tfdr_CHH\n";

	foreach (sort keys %$CG_res) {
		if (exists $CHG_res->{$_} && exists $CHH_res->{$_}) {
			my @array = split "\t";
			print FDR_OUT "$array[0]\t$array[1]\t$array[2]\t$array[3]\t$CG_res->{$_}\t$array[4]\t$array[5]\t$array[6]\t$CHG_res->{$_}\t$array[7]\t$array[8]\t$array[9]\t$CHH_res->{$_}\n";
		}
	}
unlink $out;
}

warn "Done!\n";
