#!/usr/bin/env perl
#
# Use this script to get summary statistics for labelled genomic elements and tracks.

use strict;
use warnings;
use File::Which;
use File::Basename;

#check for bedtools installation
my $bedtools = which('bedtools');

if (!$bedtools){
   print STDERR "bedtools was not found in your path:\n\n$ENV{PATH}\n\nPlease install bedtools and add it to your path:\n\n";
   print STDERR "git clone https://github.com/arq5x/bedtools2.git\ncd bedtools2\nmake clean\nmake all\n\n";
   exit(1);
}

#take in input region as command line argument
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

#define required directories
my $regionDir = join("","out/",$region);
my $labelFileDir = join("",$regionDir,"/labelFiles");
my $trackFileDir = join("",$regionDir,"/trackFiles");
my $intersectFileDir = join("",$regionDir,"/intersectFiles");
my $statsFile = join("",$regionDir,"/summary.txt");

if(-e $statsFile){
	warn "\nSummary statistics have already been generated for $region.\nPlease verify '$regionDir.'\n\n";
	die;
}
else {
	mkdir $intersectFileDir or die "Failed to create path";
}

if (!-d $labelFileDir) {
	warn "\nUnable to locate program-generated label files for the input region.\n\nPlease run makeLabels.pl or ensure that previous outputs were not moved from the default regional directory.\n";
	die;
}

if (!-d $trackFileDir) {
	warn "\nUnable to locate program-generated track files for the input region.\n\nPlease run makeTracks.pl or ensure that previous outputs were not moved from the default regional directory.\n";
	die;
}

my @files = <$labelFileDir/*.bed.gz>;
my @repeatFiles;
my %hash;
my %HoH;
my (@file, @avg, @cvg, @reg, @pct);
push @{ $hash{ $file[$_] } }, [ $avg[$_], $cvg[$_] ] for (0 .. $#file);
push @{ $HoH{ $file[$_] } }, [ $reg[$_], $pct[$_] ] for (0 .. $#file);

my $totalRegion = 0;
my $totalExon = 0;
my $totalRepeat = 0; 

#extract and store statistics
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
		my @string = split '/', $file;
		push @repeatFiles, $file;
		my $intersect_file = join("",$intersectFileDir,"/",$string[3]);
 		system "gunzip -c $file | bedtools intersect -wa -wb -a $trackFileDir/genomicCodingTrack.bed -b stdin | gzip > $intersect_file";
		my %hash2 = breakdown($intersect_file);
                foreach my $key (sort keys %hash2){ 
    		    $HoH{$file}{$key} = $hash2{$key}; 
                }
		$totalRepeat+=$cvg;
	}
  }
  else{
	warn "\nUnable to locate $file.\n\n.";
	die;
  }
}

#write 'summary.txt' file
open(IN, '>', $statsFile) or die "Could not open file '$statsFile' $!";

print IN "\nCoverage summary per genomic element in $region (percentage)\n\n";
printf IN "Exon: %.2f\n", $totalExon*100/$totalRegion;
if($totalExon == 0){
   $totalExon = 1;
}
printf IN "- Coding exon: %.2f\n", @{$hash{join("",$labelFileDir,"/coding.bed.gz")}}[1]*100/$totalExon;
printf IN "- Noncoding exon: %.2f\n", @{$hash{join("",$labelFileDir,"/noncoding.bed.gz")}}[1]*100/$totalExon;
printf IN "Intron: %.2f\n", @{$hash{join("",$labelFileDir,"/intron.bed.gz")}}[1]*100/$totalRegion;
printf IN "Intergenic: %.2f\n", @{$hash{join("",$labelFileDir,"/intergenic.bed.gz")}}[1]*100/$totalRegion;
printf IN "Repeat: %.2f\n", $totalRepeat*100/$totalRegion;
printf IN "- LINE: %.2f\n", @{$hash{join("",$labelFileDir,"/LINE.bed.gz")}}[1]*100/$totalRepeat;
printf IN "- SINE: %.2f\n", @{$hash{join("",$labelFileDir,"/SINE.bed.gz")}}[1]*100/$totalRepeat;
printf IN "- Alu: %.2f\n", @{$hash{join("",$labelFileDir,"/Alu.bed.gz")}}[1]*100/$totalRepeat;
printf IN "Promoter: %.2f\n", @{$hash{join("",$labelFileDir,"/promoter.bed.gz")}}[1]*100/$totalRegion;
printf IN "Enhancer: %.2f\n", @{$hash{join("",$labelFileDir,"/enhancer.bed.gz")}}[1]*100/$totalRegion;

print IN "\nAverage length summary per genomic element in $region (bp)\n\n";
print IN "Exon: @{$hash{join('',$labelFileDir,'/exon.bed.gz')}}[0]\n";
print IN "- Coding exon: @{$hash{join('',$labelFileDir,'/coding.bed.gz')}}[0]\n";
print IN "- Noncoding exon: @{$hash{join('',$labelFileDir,'/noncoding.bed.gz')}}[0]\n";
print IN "Intron: @{$hash{join('',$labelFileDir,'/intron.bed.gz')}}[0]\n";
print IN "Intergenic: @{$hash{join('',$labelFileDir,'/intergenic.bed.gz')}}[0]\n";
printf IN "Repeat: %.2f\n", (@{$hash{join('',$labelFileDir,'/LINE.bed.gz')}}[0]+@{$hash{join('',$labelFileDir,'/SINE.bed.gz')}}[0]+@{$hash{join('',$labelFileDir,'/Alu.bed.gz')}}[0]);
print IN "- LINE: @{$hash{join('',$labelFileDir,'/LINE.bed.gz')}}[0]\n";
print IN "- SINE: @{$hash{join('',$labelFileDir,'/SINE.bed.gz')}}[0]\n";
print IN "- Alu: @{$hash{join('',$labelFileDir,'/Alu.bed.gz')}}[0]\n";
print IN "Promoter: @{$hash{join('',$labelFileDir,'/promoter.bed.gz')}}[0]\n";
print IN "Enhancer: @{$hash{join('',$labelFileDir,'/enhancer.bed.gz')}}[0]\n";

print IN "\nRegion overlap per repetitive element in $region (percentage)\n\n";

foreach my $repeat (@repeatFiles) {
    my @str1 = split '/', $repeat;
    my @str2 = split '.bed.gz', $str1[3];
    print IN "$str2[0]:\n";
    foreach my $key (sort keys %{$HoH{$repeat}}) {
        print IN "$key - $HoH{$repeat}{$key}\n";
    }
    print IN "\n";
}

close(IN);

#subroutine to calculate percent composition and average length of genomic elements
sub stats {

   my ($infile) = @_;
   my $coverage = 0;
   my $total = 0;
   my $average = 0;

   open(IN, '-|' ,"gunzip -c $infile") || die "Could not open $infile: $!\n";
   while(<IN>){
      chomp;
      ++$total;
      my ($chrom, $chromStart, $chromEnd) = split(/\t/);
      my $c = $chromEnd - $chromStart;
      $coverage += $c;
   }
   close(IN);
   if($total !=0){
   	$average = sprintf ("%.2f", $coverage / $total);
   	return($average, $coverage);
   }
   else{
        return(0,0);
   }
}

#subroutine to calculate percent overlap of certain elements with different genomic regions 
sub breakdown {

   my ($infile) = @_;
   my $total = 0;
   my %hash;
   my (@region, @overlap);
   push @{ $hash{ $region[$_] } }, [ $overlap[$_] ] for (0 .. $#region);

   open(IN, '-|' ,"gunzip -c $infile") || die "Could not open $infile: $!\n";
   while(<IN>){
      chomp;
      ++$total;
      my @spl = split(/\t/);
      push @region, $spl[3];
    }
    $hash{$_}++ for @region;
    close(IN);

    foreach my $key (sort keys %hash ){
        if($total !=0){
	    $hash{$key} = sprintf ("%.2f", $hash{$key} / $total * 100);
        }
        else{
            $hash{$key} = 0;
        }
    }
    return(%hash);
}

#check file generation
if(-e $statsFile){
	print "\nSummary statistics successfully extracted and saved as '.txt' file to '$regionDir'.\n";
	print "\nWould you like to print the results (Y/N)?\n";
	my $printAns = <STDIN>;
	if($printAns=~m/^[Y]$/i){
		system "cat $statsFile";
	}
}

exit(0);
