#!/usr/bin/env perl
#
# Use this script to create tracks for genomic elements using the labels created with 'makeLabels.pl'.
# To run, enter the following command prompt in a unix shell:
# ./makeTracks.pl <chr:start-end> <utilities>

use strict;
use warnings;
use File::Which;
use File::Basename;

#take in input region as command line argument
my $usage = "Usage: $0 <chr:start-end> <utilities>\n";
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

#define and create required variables
my $toolfiles = shift or die $usage;

#define required directories
my $labelFileDir = join("","out/",$region,"/labelFiles");
my $trackFileDir = join("","out/",$region,"/trackFiles");
my $hubDir = "hubs";
my $hubname = join("","myHub_",$region);
my $hubRegionDir = join("",$hubDir,"/",$hubname);

if (-d $trackFileDir) {
        print "\nTracks have already been generated for $region.\nPlease verify '$trackFileDir.'\n";
        if (!-d $hubRegionDir){
                print "\nWould you like to create a UCSC track hub to using these results (Y/N)?\n";
                my $hubAns = <STDIN>;
                if($hubAns=~m/^[Y]$/i){
                        print "\nPlease provide an email address.\n";
                        my $email = <STDIN>;
                        chomp($email);
                        system "./makeHubs.pl $region $toolfiles $email";
                }
        die;
        }
        else{
                print "\nA hub has also already been generated for $region.\nPlease verify '$hubRegionDir'.\n\n";
                die;
        }
}
else {
        if (-d $labelFileDir) {
                mkdir $trackFileDir or die "Failed to create path";
        }
        else {
                warn "\nUnable to locate program-generated label files for the input region.\n\nPlease run 'makeLabels.pl' or ensure that previous outputs were not moved from the default directory.\n";
                die;
        }
}

#define required and generated files
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
my $functional_element_file = join("",$labelFileDir,"/functional_element.bed.gz");
my $crm_file = join("",$labelFileDir,"/crm.bed.gz");

my $genomic_track = join("",$trackFileDir,"/genomicTrack.bed");
my $coding_track = join("",$trackFileDir,"/codingTrack.bed");
my $combined_track = join("",$trackFileDir,"/genomicCodingTrack.bed");
my $promoter_track = join("",$trackFileDir,"/promoterTrack.bed");
my $enhancer_track = join("",$trackFileDir,"/enhancerTrack.bed");
my $repeat_track = join("",$trackFileDir,"/repeatTrack.bed");
my $functional_element_track = join("",$trackFileDir,"/functional_elementTrack.bed");
my $crm_track = join("",$trackFileDir,"/crmTrack.bed");

if (!-e $genomic_track){
   warn "\nCreating genomic track (i.e. exon, intron, intergenic)\n";
   my $command = "zcat $exon_file $intron_file $intergenic_file > $genomic_track";
   system($command);
} else {
   warn "\n$genomic_track already exists; skipping step\n";
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

if (!-e $enhancer_track){
   warn "Creating enhancer track (i.e. active, transcribed in-vivo)\n";
   my $command = "zcat $enhancer_file > $enhancer_track";
   system($command);
} else {
   warn "$enhancer_track already exists; skipping step\n";
}

if (!-e $repeat_track){
   warn "Creating repetitive element track (i.e. LINEs, SINEs, Alu)\n";
   my $command = "zcat $LINE_file $SINE_file $Alu_file > $repeat_track";
   system($command);
} else {
   warn "$repeat_track already exists; skipping step\n";
}

if (!-e $functional_element_track){
   warn "Creating functional element track (i.e. cell type-agnostic promoter, enhancer, transcribed, bivalent, quiescent regions etc.)\n";
   my $command = "zcat $functional_element_file > $functional_element_track";
   system($command);
} else {
   warn "$functional_element_track already exists; skipping step\n";
}

if (!-e $crm_track){
   warn "Creating cis-regulatory module track (i.e. cell type and transcriptor regulator-agnostic binding sites)\n";
   my $command = "zcat $crm_file > $crm_track";
   system($command);
} else {
   warn "$crm_track already exists; skipping step\n";
}

#check file generation
if(-d $trackFileDir){
        print "\nTracks successfully generated and saved as '.bed' files to '$trackFileDir'.\n";
        print "\nWould you like to create a UCSC track hub to using these results (Y/N)?\n";
        my $hubAns = <STDIN>;
        if($hubAns=~m/^[Y]$/i){
                print "\nPlease provide an email address.\n";
                my $email = <STDIN>;
                chomp($email);
                system "./makeHubs.pl $region $toolfiles $email";
        }
}


exit(0);
