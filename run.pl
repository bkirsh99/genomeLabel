#!/usr/bin/env perl
#
# Use this script to pull and filter data for the genomeLabel tool.

use strict;
use warnings;
use File::Which;
use File::Basename;

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

#create basic directories 
my $outDir = "out";
my $dataDir = "data";
my $regionDir = join("",$outDir,"/",$region);
my $regionFileDir = join("",$regionDir,"/dataFiles");

print "\nChecking for previously downloaded raw data [...]\n\n";
if (!-d $outDir){
	mkdir $outDir or die "Failed to create path";
}

if (!-d $dataDir){
	mkdir $dataDir or die "Failed to create path";
	print "\nNo data found.\nFetching raw data [...]";
}


#paths to store downloaded raw data 
my $geneFile = join("",$dataDir,"/hg38.ncbiRefSeq.gtf.gz");
my $repeatFile = join("",$dataDir,"/rmsk.txt.gz");
my $promoterFile = join("",$dataDir,"/hg38_fair+new_CAGE_peaks_phase1and2.bed.gz");
my $enhancerFile = join ("",$dataDir,"/F5.hg38.enhancers.bed.gz");
my $genome = join("",$dataDir,"/hg38.chrom.sizes");

if (!-e $geneFile){
	system "wget -P $dataDir https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz";
}

if (!-e $repeatFile){
	system "wget -P $dataDir https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz";
}

if (!-e $promoterFile){
	system "wget -P $dataDir https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest/extra/CAGE_peaks/hg38_fair+new_CAGE_peaks_phase1and2.bed.gz";
}

if (!-e $enhancerFile){
	system "wget -P $dataDir https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest/extra/enhancer/F5.hg38.enhancers.bed.gz";
}

if (!-e $genome){
	system "wget -P $dataDir https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes";
}

print "\nChecking for previously filtered regional data [...]\n\n";
if (!-d $regionDir) {
	mkdir $regionDir or die "Failed to create path";
	print "\nNo data found.\nFiltering raw data [...]";
}

if (!-d $regionFileDir) {
	mkdir $regionFileDir or die "Failed to create path";
}

#paths to store filtered data for input region
my $regionGeneFile = join("",$regionFileDir,"/hg38.ncbiRefSeq.gtf.gz");
my $regionRepeatFile = join("",$regionFileDir,"/rmsk.txt.gz");
my $regionPromoterFile = join("",$regionFileDir,"/hg38_fair+new_CAGE_peaks_phase1and2.bed.gz");
my $regionEnhancerFile = join("",$regionFileDir,"/F5.hg38.enhancers.bed.gz");
my $regionGenome = join("",$regionFileDir,"/hg38.chrom.sizes");

if (!-e $regionGeneFile){
        my $command = "gunzip -c $geneFile | awk -F \"\\t\" 'BEGIN{OFS=\"\\t\";} ((\$9~/NM/ || \$9~/NR/) && \$1==\"$chr\" && ((\$5 >= $start && \$5 <= $end) || (\$4 >= $start && \$4 <= $end)))' | gzip > tmp && mv tmp $regionGeneFile";
        system($command);
}

if (!-e $regionRepeatFile){
	my $command = "gunzip -c $repeatFile | awk -F \"\\t\" 'BEGIN{OFS=\"\\t\";} (\$6==\"$chr\" && ((\$8 >= $start && \$8 <= $end) || (\$7 >= $start && \$7 <= $end)))' | gzip > tmp && mv tmp $regionRepeatFile";
	system($command);
}

if (!-e $regionPromoterFile){
	my $command = "gunzip -c $promoterFile | awk -F \"\\t\" 'BEGIN{OFS=\"\\t\";} (\$1==\"$chr\" && ((\$3 >= $start && \$3 <= $end) || (\$2 >= $start && \$2 <= $end)))' | gzip > tmp && mv tmp $regionPromoterFile";
	system($command);
}

if (!-e $regionEnhancerFile){
	my $command = "gunzip -c $enhancerFile | awk -F \"\\t\" 'BEGIN{OFS=\"\\t\";} (\$1==\"$chr\" && ((\$3 >= $start && \$3 <= $end) || (\$2 >= $start && \$2 <= $end)))' | gzip > tmp && mv tmp $regionEnhancerFile";
	system($command);
}

if (!-e $regionGenome){
	my $command = "cat $genome | awk -F \"\\t\" 'BEGIN{OFS=\"\\t\";} (\$1==\"$chr\")' > tmp && mv tmp $regionGenome";
	system($command);
}

#paths to store program-generated annotation files
my $labelFileDir = join("",$regionDir,"/labelFiles");
my $trackFileDir = join("",$regionDir,"/trackFiles");
my $statsFile = join("",$regionDir,"/summary.txt");

#additional functionalities of 'run.pl'
print "\nWould you like to generate labels for the features in $region (Y/N)?\n";
my $labelAns = <STDIN>;
if ($labelAns =~ m/^[Y]$/i){
	system "./makeLabels.pl $region";
}

print "\nWould you like to generate track files for data visualization (Y/N)?\n";
my $trackAns = <STDIN>;
if ($trackAns =~ m/^[Y]$/i){
	system "./makeTracks.pl $region";
}
print "\nWould you like to extract summary statistics for data analysis (Y/N)?\n";
my $summaryAns = <STDIN>;
if ($summaryAns =~ m/^[Y]$/i){
	system "./getStats.pl $region"; 
}

exit(0);

