#!/usr/bin/env perl
#
# Use script to create tracks for genomic elements using the labels created with "makeLabels.pl".
#

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

if (!-d $labelFileDir) {
	warn "Label files have not been created for that region.\nPlease run makeLabels.pl first and then try again.\n";
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

my $trackDir = join("","out/",$region,"/trackFiles");

if (!-d $trackDir) {
	mkdir $trackDir or die "Failed to create path";	
}

my $genomic_track = join("",$trackDir,"/genomicTrack.bed");
my $coding_track = join("",$trackDir,"/codingTrack.bed");
my $combined_track = join("",$trackDir,"/genomicCodingTrack.bed");
my $promoter_track = join("",$trackDir,"/promoterTrack.bed");
my $repeat_track = join("",$trackDir,"/repeatTrack.bed");

if (!-e $genomic_track){
   warn "Creating genomic track (i.e. exon, intron, intergenic)\n";
   my $command = "zcat $exon_file $intron_file $intergenic_file > $genomic_track";
   system($command);
} else {
   warn "$genomic_track already exists; skipping step\n";
}

if (!-e $coding_track){
   warn "Creating coding track (i.e. coding exon, non-coding exon)\n";
   my $command = "zcat $coding_file $noncoding_file > $coding_track";
   system($command);
} else {
   warn "$coding_track already exists; skipping step\n";
}

if (!-e $combined_track){
   warn "Creating genomic coding track (i.e. coding exon, non-coding exon, intron, intergenic)\n";
   my $command = "zcat $coding_file $noncoding_file $intron_file $intergenic_file > $combined_track";
   system($command);
} else {
   warn "$combined_track already exists; skipping step\n";
}

if (!-e $promoter_track){
   warn "Creating promoter track (i.e. immediate promoters)\n";
   my $command = "zcat $promoter_file > $promoter_track";
   system($command);
} else {
   warn "$promoter_track already exists; skipping step\n";
}

if (!-e $repeat_track){
   warn "Creating repetitive element track (i.e. LINEs, SINEs, Alu)\n";
   my $command = "zcat $LINE_file $SINE_file $Alu_file > $repeat_track";
   system($command);
} else {
   warn "$repeat_track already exists; skipping step\n";
}

exit(0);
