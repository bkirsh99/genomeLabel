#!/usr/bin/env perl
#
# Use script to get summary statistics for labelled genomic elements.

use strict;
use warnings;
use File::Which;
use File::Basename;

my $usage = "Usage: $0 <chr:start-end>\n";
my $region = shift or die $usage;
my $chr = '';
my $start = 0;
my $end = 0;
 
#get rid of commas from the input region
$region =~ s/,//g;
 
#check region format
if ($region =~ /^(\w+):(\d+)-(\d+)$/){
   $chr = $1;
   $start = $2;
   $end = $3;
} else {
   die "Could not recognise $region\n";
}

my $labelFileDir = join("","out/",$region,"/labelFiles");
my $regionDir = join("","out/",$region);
my $statsFile = join("",$regionDir,"/summary.txt");

if(-e $statsFile){
	warn "\nSummary statistics have already been generated for $region.\nPlease verify '$regionDir.'\n\n";
	die;
}

if (!-d $labelFileDir) {
	warn "\nUnable to locate program-generated label files for the input region.\n\nPlease run makeLabels.pl or ensure that previous outputs were not moved from the default regional directory.\n";
	die;
}

my @files = <$labelFileDir/*.bed.gz>;
my %hash;
my (@file, @avg, @cvg);
push @{ $hash{ $file[$_] } }, [ $avg[$_], $cvg[$_] ] for (0 .. $#file);

my $totalRegion = 0;
my $totalExon = 0;
my $totalRepeat = 0; 


foreach my $file (@files) {
  if(-e $file){
	my ($avg, $cvg) = stats($file);
	@{$hash{$file}}[0]=$avg;
	@{$hash{$file}}[1]=$cvg;
	if($file=~/intron/ || $file=~/exon/ || $file=~/intergenic/){
		$totalRegion+=$cvg;
		if($file=~/exon/){
			$totalExon+=$cvg;
		}
	}
	if($file=~/LINE/ || $file=~/SINE/ || $file=~/Alu/){
		$totalRepeat+=$cvg;
	}
  }
  else{
	warn "\nUnable to locate $file.\n\n.";
	die;
  }
}

open(IN, '>', $statsFile) or die "Could not open file '$statsFile' $!";

print IN "\nCoverage summary per genomic element in $region (percentage)\n\n";
printf IN "Exon: %.2f\n", $totalExon*100/$totalRegion;
printf IN "- Coding exon: %.2f\n", @{$hash{join("",$labelFileDir,"/coding.bed.gz")}}[1]*100/$totalExon;
printf IN "- Noncoding exon: %.2f\n", @{$hash{join("",$labelFileDir,"/noncoding.bed.gz")}}[1]*100/$totalExon;
printf IN "Intron: %.2f\n", @{$hash{join("",$labelFileDir,"/intron.bed.gz")}}[1]*100/$totalRegion;
printf IN "Intergenic: %.2f\n", @{$hash{join("",$labelFileDir,"/intergenic.bed.gz")}}[1]*100/$totalRegion;
printf IN "Repeat: %.2f\n", $totalRepeat*100/$totalRegion;
printf IN "- LINE: %.2f\n", @{$hash{join("",$labelFileDir,"/LINE.bed.gz")}}[1]*100/$totalRepeat;
printf IN "- SINE: %.2f\n", @{$hash{join("",$labelFileDir,"/SINE.bed.gz")}}[1]*100/$totalRepeat;
printf IN "- Alu: %.2f\n", @{$hash{join("",$labelFileDir,"/Alu.bed.gz")}}[1]*100/$totalRepeat;
printf IN "Promoter: %.2f\n", @{$hash{join("",$labelFileDir,"/promoter.bed.gz")}}[1]*100/$totalRegion;

print IN "\nAverage length summary per genomic element in $region (bp)\n\n";
print IN "Exon: @{$hash{join('',$labelFileDir,'/exon.bed.gz')}}[0]\n";
print IN "- Coding exon: @{$hash{join('',$labelFileDir,'/coding.bed.gz')}}[0]\n";
print IN "- Noncoding exon: @{$hash{join('',$labelFileDir,'/noncoding.bed.gz')}}[0]\n";
print IN "Intron: @{$hash{join('',$labelFileDir,'/intron.bed.gz')}}[0]\n";
print IN "Intergenic: @{$hash{join('',$labelFileDir,'/intergenic.bed.gz')}}[0]\n";
print IN "- LINE: @{$hash{join('',$labelFileDir,'/LINE.bed.gz')}}[0]\n";
print IN "- SINE: @{$hash{join('',$labelFileDir,'/SINE.bed.gz')}}[0]\n";
print IN "- Alu: @{$hash{join('',$labelFileDir,'/Alu.bed.gz')}}[0]\n";
print IN "Promoter: @{$hash{join('',$labelFileDir,'/promoter.bed.gz')}}[0]\n\n";

close(IN);

sub stats {

   my ($geneFile) = @_;
   my $coverage = 0;
   my $total = 0;
   my $average = 0;

   open(IN, '-|' ,"gunzip -c $geneFile") || die "Could not open $geneFile: $!\n";
   while(<IN>){
      chomp;
      ++$total;
      my ($chrom, $chromStart, $chromEnd) = split(/\t/);
      my $c = $chromEnd - $chromStart;
      $coverage += $c;
   }
   close(IN);

   $average = sprintf("%.2f", $coverage / $total);

   return($average, $coverage);
}

exit(0);
