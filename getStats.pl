#!/usr/bin/env perl
#
# Use this script to get summary statistics for labelled genomic elements and tracks.

package UndefMod;

use strict;
use warnings;
use 5.010;

no warnings 'uninitialized';

use File::Which;
use File::Basename;
use experimental 'smartmatch';
use Tie::IxHash;

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

my $exon_file = join("",$labelFileDir,"/exon.bed.gz");
my $intergenic_file = join("",$labelFileDir,"/intergenic.bed.gz");
my $intron_file = join("",$labelFileDir,"/intron.bed.gz");
my $coding_file = join("",$labelFileDir,"/coding.bed.gz");
my $noncoding_file = join("",$labelFileDir,"/noncoding.bed.gz");
my $LINE_file = join("",$labelFileDir,"/LINE.bed.gz");
my $SINE_file = join("",$labelFileDir,"/SINE.bed.gz");
my $Alu_file = join("",$labelFileDir,"/Alu.bed.gz");
my $promoter_file = join("",$labelFileDir,"/promoter.bed.gz");
my $enhancer_file = join("",$labelFileDir,"/enhancer.bed.gz");

my $t1 = tie(my %files, 'Tie::IxHash', $exon_file => 'Exonic', $coding_file => 'Coding exonic', $noncoding_file => 'Noncoding exonic', $intron_file => 'Intronic', $intergenic_file => 'Intergenic', $LINE_file => 'LINE repeats', $SINE_file => 'SINE repeats', $Alu_file => 'Alu repeats', $promoter_file => 'Promoter', $enhancer_file => 'Enhancer (active, transcribed in-vivo)');
my $t2 = tie(my %hash, 'Tie::IxHash');
my $t3 = tie(my %HoH, 'Tie::IxHash');
my (@file, @avg, @cvg, @reg, @pct);
my ($file, $avg, $cvg, $reg, $pct);
my @repeatFiles;
my @codingFiles;
push @{ $hash{ $file[$_] } }, [ $avg[$_], $cvg[$_] ] for (0 .. $#file);
push @{ $HoH{ $file[$_] } }, [ $reg[$_], $pct[$_] ] for (0 .. $#file);

my $totalRegion = 0;
my $totalExon = 0;
my $totalRepeat = 0;

#extract and store statistics
foreach $file (keys %files) {
  if(-e $file){
        ($avg, $cvg) = stats($file);
        $hash{$file}[0]=$files{$file};
        $hash{$file}[1]=$avg;
        $hash{$file}[2]=$cvg;
        if($file=~/intron/ || $file=~/exon/ || $file=~/intergenic/){
                $totalRegion+=$cvg;
                if($file=~/exon/){
                        $totalExon+=$cvg;
                }
        }
        if($file=~/LINE/ || $file=~/SINE/ || $file=~/Alu/){
                push @repeatFiles, $file;
                my @string = split '/', $file;
                my $intersect_file = join("",$intersectFileDir,"/",$string[3]);
                system "gunzip -c $file | bedtools intersect -wa -wb -a $trackFileDir/genomicCodingTrack.bed -b stdin | gzip > $intersect_file";
                my %hash2 = breakdown($intersect_file);
                foreach my $key (keys %hash2){
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
foreach my $key (keys %hash) {
    if(@{$hash{$key}}[1]>0){
        if($key=~/coding/){
                printf IN "\t- $files{$key}: %.2f%% of exons (or %.2f%% of region)\n", @{$hash{$key}}[2]*100/$totalExon, @{$hash{$key}}[2]*100/$totalRegion;
        }
        elsif($key=~/LINE/ || $key=~/SINE/ || $key=~/Alu/){
                printf IN "\t- $files{$key}: %.2f%% of repeats (or %.2f%% of region)\n", @{$hash{$key}}[2]*100/$totalRepeat, @{$hash{$key}}[2]*100/$totalRegion;
        }
        else{
                printf IN "$files{$key}: %.2f%% of region\n", @{$hash{$key}}[2]*100/$totalRegion;
                if($key=~/intergenic/){
                        printf IN "Repeat: %.2f%% of region\n", $totalRepeat*100/$totalRegion;
                }
        }
    }
}

print IN "\nAverage length summary per genomic element in $region (bp)\n\n";
foreach my $key (keys %hash) {
    if(@{$hash{$key}}[2]>0){
        printf IN "$files{$key}: %.0f bps\n", @{$hash{$key}}[1];
        if($key=~/intergenic/){                                                                                                                                              printf IN "Repeat: %.2f%%\n of region", $totalRepeat*100/$totalRegion;
                printf IN "Repeat: %.0f\n", (@{$hash{join('',$labelFileDir,'/LINE.bed.gz')}}[1]+@{$hash{join('',$labelFileDir,'/SINE.bed.gz')}}[1]+@{$hash{join('',$labelFileDir,'/Alu.bed.gz')}}[1]);
        }
    }
}

print IN "\nRegion overlap per repetitive element in $region (percentage)\n\n";

foreach my $repeatName (@repeatFiles) {
    my @str1 = split '/', $repeatName;
    my @str2 = split '.bed.gz', $str1[3];
    print IN "$str2[0]:\n";
    foreach my $key (keys %{$HoH{$repeatName}}) {
        print IN "\t- $key - $HoH{$repeatName}{$key}%%\n";
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

    foreach my $key (keys %hash ){
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
