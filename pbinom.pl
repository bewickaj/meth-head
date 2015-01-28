#!/usr/bin/perl

use strict;
use warnings;
use Math::CDF;
use Getopt::Long;

my $bis;
my $out;
my $p;
my $c;
my $help;

if (@ARGV > 0) {
	GetOptions (
		'bismark_CX_report:s' => \$bis,
		'outfile:s' => \$out,
		'probability:s' => \$p,
		'confidence:s' => \$c,
		'help' => \$help
	);
}
if ($help) {
	print "Usage: --bismark_CX_report \"name of bismark CX report\" --outfile \"file name of vcf output file\" --probability \"integer indicating the hypothesized probability of success i.e., the non-conversion rate of the chloroplast genome, 0 to 1\" --confidence \"integer indicating confidence level i.e., 0.05 for 95%\"\n";
	exit;
}
if (($p < 0) or (1 < $p)) {
	die "Probability of success must be between 0 and 1"; 
}
if (($c < 0) or (1 < $c)) {
	die "Confidence level must be between 0 and 1"; 
}

open BIS, "<$bis" or die;
open OUT, ">$out" or die;

print OUT "##fileformat=VCFv4.1\n";
print OUT "#chr\tpos\tstrand\tmc_class\tmc_count\ttotal\tmethylated\n";

while (<BIS>) {
	chomp;
	my ($chr, $pos, $strand, $meth, $noMeth, $con, $tri) = split "\t";	
	if ($meth + $noMeth >= 3) {
		my $x = $meth;
		my $n = $meth + $noMeth;
		my $prob;
		if ($x == 0) {
			$prob = 1;
		}
		else {
			$prob = 1-&Math::CDF::pbinom($x-1, $n, $p);
		}
		if ($prob < $c) {
			print OUT "$chr\t$pos\t$strand\t$con\t$meth\t", $meth+$noMeth, "\t1\n";
		}
		else {
			print OUT "$chr\t$pos\t$strand\t$con\t$meth\t", $meth+$noMeth, "\t0\n";
		}
	}
	else {
		print OUT "$chr\t$pos\t$strand\t$con\t$meth\t", $meth+$noMeth, "\t0\n";
	}
}