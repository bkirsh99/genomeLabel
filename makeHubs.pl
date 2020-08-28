#!/usr/bin/env perl
#
# Use script to make a track hub from the BED files created using "makeTracks.pl".
#
# To run the TrackHub Generator enter the following command prompt in a unix shell:
# ./makeTrackHub.pl chr:start-end UCSC EMAIL
#
# Positional arguments:
# chr:start-end		Region for which track hub will be created (needed to retrieve track files)
# UCSC			Path to the folder containing the bedToBigBed binary and the fetchChromSizes script downloaded from UCSC
# EMAIL			Email address to be used as contact for people using the track hub
#
#
# Dependencies:
# The TrackHub Generator relies on two tools provided through UCSC Genomes which can be downloaded from http://hgdownload.soe.ucsc.edu/admin/exe/
# - bedToBigBed			binary tool to create binary BigBed files from input bed files
# - fetchChromSizes		script to download chromosome sizes for indexing BigBed files during generation through bedToBigBed
#
#
# For more information please refere to www.sanger.ac.uk/science/tools/trackhub-generator

use strict;
use warnings;
use File::Which;
use File::Basename;
use File::Find::Rule;

my $usage = "Usage: $0 <chr:start-end> <toolfiles> <email>\n";
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

my $hubDir = "hubs";
my $hubname = join("",$hubDir,"/myHub_",$region);
my $genome = "hg38";
my $bedfiles = join("","out/",$region,"/trackFiles");
my $toolfiles = shift or die $usage;
my $email = shift or die $usage;
my $shortlabel = "Genomic Label Hub ($region)";
my $longlabel = "Genomic Label Hub for Identifying Cis-Regulatory Elements ($region)";
my $genomepath = join("",$hubname,"/",$genome);

if (!-d $hubDir) {
	mkdir $hubDir or die "Failed to create path for hubs";
}
if(!-d $hubname) {
	mkdir $hubname or die "Failed to create path for hub in $region";
	mkdir $genomepath or die "Failed to create path for genome in $region";
}


if(!-d $bedfiles){
	warn "Unable to locate program-generated files for the input region.\n\nPlease run makeLabels.pl and makeTracks.pl or ensure that previous outputs were not moved from the default regional directory.\n";
	die;
}

my $hubfile = join("",$hubname,"/hub.txt");
open(OUT, '>', $hubfile) or die "Could not open file $hubfile";
print OUT "hub ${hubname}\nshortLabel ${shortlabel}\nlongLabel ${longlabel}\ngenomesFile genomes.txt\nemail ${email}";
close OUT;

my $genomesfile = join("",$hubname,"/genomes.txt");
open(OUT, '>', $genomesfile) or die "Could not open file $genomesfile";
print OUT "genome ${genome}\ntrackDb ${genome}/trackDb.txt\n\n";
close OUT;

my $chromsizes =join("","data/",$genome,".chrom.sizes");

if(!-f $chromsizes) {
	warn "Unable to locate program-generated data files. \n\nPlease ensure that they were not moved from the default data directory.\n";
	die;
}

my $bedtobigbedfile = join("",$toolfiles,"/bedToBigBed");
if(-f $bedtobigbedfile) {
	my $trackDbstring = "";
	my @bedfiles = glob "$bedfiles/*.bed";
	for my $file (@bedfiles) {
		my $newtempbed = join("",$file,".2");
		my $result = `sort -k1,1 -k2,2n "${file}" > "${newtempbed}"`;
		my $filename = basename($file);
		$filename =~ s/\.bed//g;
		my $shortfilename = substr($filename,0,17);
		$filename = substr($filename,0,80);
		my $bbfile = join("",$genomepath,"/",$filename,".bb");
		print "Converting $file to bigBed format\n";
		$result = `${bedtobigbedfile} "${newtempbed}" ${chromsizes} "${bbfile}"`;
		$result = `rm "${newtempbed}"`;
		my $trackname = $filename;
		$trackname =~ tr/ //ds;
		$trackDbstring = $trackDbstring . "track ".$trackname."\nbigDataUrl ".$filename.".bb\nshortLabel ".$shortfilename."\nlongLabel ".$filename."\ntype bigBed 12\nvisibility dense\nitemRgb on\n\n";
	}

	my $trackDbfile = join("",$genomepath,"/","trackDb.txt");
	open(OUT, '>', $trackDbfile) or die "Could not open file $trackDbfile";
	print OUT $trackDbstring;
	close OUT;
} else {
	die "Tool \"bedToBigBed\" not found: $bedtobigbedfile";
}

exit(0);
