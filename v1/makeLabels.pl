#!/usr/bin/env perl
#
# Use this script to create labels for genomic elements.
# To run, enter the following command prompt in a unix shell:
# ./makeLabels.pl <chr:start-end>
#
# Positional arguments:
# chr:start-end         Region for which labels will be created

use strict;
use warnings;
use File::Which;
use File::Basename;
use List::Priority;
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
my $outDir = "out";
my $dataDir = "data";
my $regionDir = join("",$outDir,"/",$region);
my $regionFileDir = join("",$regionDir,"/dataFiles");
my $labelFileDir = join("",$regionDir,"/labelFiles");

if (-d $labelFileDir) {
        print "\nLabels have already been generated for $region.\nPlease verify '$labelFileDir'.\n";
        die;
}
else {
        if (-d $outDir && -d $dataDir && -d $regionDir && -d $regionFileDir) {
                mkdir $labelFileDir or die "Failed to create path";
        }
        else {
                warn "\nUnable to locate required data files for the input region.\nPlease run 'run.pl' or ensure that previous outputs were not moved from the default directory.\n";
                die;
        }
}

#define required and generated files
my $regionGeneFile = join("",$regionFileDir,"/hg38.ncbiRefSeq.gtf.gz");
my $regionRepeatFile = join("",$regionFileDir,"/rmsk.txt.gz");
my $regionPromoterFile = join("",$regionFileDir,"/hg38_fair+new_CAGE_peaks_phase1and2.bed.gz");
my $regionEnhancerFile = join("",$regionFileDir,"/F5.hg38.enhancers.bed.gz");
my $regionSegwayFile = join("",$regionFileDir,"/segway_encyclopedia.bed.gz");
my $regionReMapFile = join("",$regionFileDir,"/remap2020_crm_macs2_hg38_v1_0.bed.gz");
my $regionGenome = join("",$regionFileDir,"/hg38.chrom.sizes");

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

if (!-e $exon_file){
   warn "Creating exonic regions\n";
   my $command = "gunzip -c $regionGeneFile | awk -F \"\\t\" 'BEGIN{OFS=\"\\t\";}; (\$3==\"exon\") {print \$1,\$4-1,\$5}' | bedtools sort | bedtools merge -i stdin | awk 'BEGIN{OFS=\"\\t\";}; {if(\$2<$start) \$2=$start}; {if(\$3>$end) \$3=$end}; \$3>=$start && \$2<=$end' | awk 'BEGIN{OFS=\"\\t\";} {print \$1,\$2,\$3,\"exon\",0,\".\",\$2,\$3,\"0,255,0\"}' | uniq | gzip > $exon_file";
   system($command);
} else {
   warn "$exon_file already exists; skipping exon step\n";
}

if (!-e $intron_file){
   warn "Creating intronic regions\n";
   my $command = "gunzip -c $regionGeneFile | awk -F \"\\t\" 'BEGIN{OFS=\"\\t\";}; (\$3==\"transcript\") {print \$1,\$4-1,\$5}' | bedtools sort | bedtools merge -i stdin | bedtools subtract -a stdin -b $exon_file | awk 'BEGIN{OFS=\"\\t\";}; {if(\$2<$start) \$2=$start}; {if(\$3>$end) \$3=$end}; \$3>=$start && \$2<=$end' | awk 'BEGIN{OFS=\"\\t\";} {print \$1,\$2,\$3,\"intron\",0,\".\",\$2,\$3,\"255,0,0\"}' | gzip > $intron_file";
   system($command);
} else {
   warn "$intron_file already exists; skipping intron step\n";
}

if (!-e $intergenic_file){
   warn "Creating intergenic regions\n";
   my $command = "gunzip -c $regionGeneFile | awk -F \"\\t\" 'BEGIN{OFS=\"\\t\";}; (\$3==\"transcript\") {print \$1,\$4-1,\$5}' | bedtools sort | bedtools complement -i stdin -g $regionGenome | awk 'BEGIN{OFS=\"\\t\"}; {if(\$2<$start) \$2=$start}; {if(\$3>$end) \$3=$end}; \$3>=$start && \$2<=$end' | awk 'BEGIN{OFS=\"\\t\";} {print \$1,\$2,\$3,\"intergenic\",0,\".\",\$2,\$3,\"0,0,255\"}' | gzip > $intergenic_file";
   system($command);
} else {
   warn "$intergenic_file already exists; skipping intergenic step\n";
}

if (!-e $coding_file){
   warn "Creating coding exonic regions\n";
   my $command = "gunzip -c $regionGeneFile | awk -F \"\\t\" 'BEGIN{OFS=\"\\t\";}; (\$3==\"CDS\" || \$3==\"start_codon\" || \$3==\"stop_codon\") {print \$1,\$4-1,\$5}' |  bedtools sort | bedtools merge -i stdin | awk 'BEGIN{OFS=\"\\t\";}; {if(\$2<$start) \$2=$start}; {if(\$3>$end) \$3=$end}; \$3>=$start && \$2<=$end' | awk 'BEGIN{OFS=\"\\t\";} {print \$1,\$2,\$3,\"coding\",0,\".\",\$2,\$3,\"0,102,0\"}' | uniq | gzip > $coding_file";
   system($command);
} else {
   warn "$coding_file already exists; skipping coding step\n";
}

if (!-e $noncoding_file){
   warn "Creating non-coding exonic regions\n";
   my $command = "gunzip -c $regionGeneFile | awk -F \"\\t\" 'BEGIN{OFS=\"\\t\";}; ((\$3~/UTR/) || (\$9~/NR/ && \$3==\"exon\")) {print \$1,\$4-1,\$5}' | bedtools sort | bedtools subtract -a stdin -b $coding_file | sort -k1,1 -k2,2n | bedtools merge -i stdin | awk 'BEGIN{OFS=\"\\t\";}; {if(\$2<$start) \$2=$start}; {if(\$3>$end) \$3=$end}; \$3>=$start && \$2<=$end' | awk 'BEGIN{OFS=\"\\t\";} {print \$1,\$2,\$3,\"noncoding\",0,\".\",\$2,\$3,\"153,255,153\"}' | uniq | gzip > $noncoding_file";
   system($command);
} else {
   warn "$noncoding_file already exists; skipping non-coding step\n";
}

if (!-e $LINE_file){
   warn "Extracting LINE repeats\n";
   my $command = "gunzip -c $regionRepeatFile | awk -F \"\\t\" 'BEGIN{OFS=\"\\t\";} (\$12==\"LINE\") {print \$6,\$7,\$8}' | bedtools sort | bedtools merge -i stdin | awk 'BEGIN{OFS=\"\\t\";}; {if(\$2<$start) \$2=$start}; {if(\$3>$end) \$3=$end}; \$3>=$start && \$2<=$end' | awk 'BEGIN{OFS=\"\\t\";} {print \$1,\$2,\$3,\"LINE\",0,\".\",\$2,\$3,\"255,255,0\"}' | uniq | gzip > $LINE_file";
   system($command);
} else {
   warn "$LINE_file already exists; skipping LINE step\n";
}

if (!-e $SINE_file){
   warn "Extracting non-Alu SINE repeats\n";
   my $command = "gunzip -c $regionRepeatFile | awk -F \"\\t\" 'BEGIN{OFS=\"\\t\";} (\$12==\"SINE\" && \$13!=\"Alu\") {print \$6,\$7,\$8}' | bedtools sort | bedtools merge -i stdin | awk 'BEGIN{OFS=\"\\t\";}; {if(\$2<$start) \$2=$start}; {if(\$3>$end) \$3=$end}; \$3>=$start && \$2<=$end' | awk 'BEGIN{OFS=\"\\t\";} {print \$1,\$2,\$3,\"SINE\",0,\".\",\$2,\$3,\"255,204,153\"}' | uniq | gzip > $SINE_file";
   system($command);
} else {
   warn "$SINE_file already exists; skipping SINE step\n";
}

if (!-e $Alu_file){
   warn "Extracting Alu repeats\n";
   my $command = "gunzip -c $regionRepeatFile | awk -F \"\\t\" 'BEGIN{OFS=\"\\t\";} (\$13==\"Alu\") {print \$6,\$7,\$8}' | bedtools sort | bedtools merge -i stdin | awk 'BEGIN{OFS=\"\\t\";}; {if(\$2<$start) \$2=$start}; {if(\$3>$end) \$3=$end}; \$3>=$start && \$2<=$end' | awk 'BEGIN{OFS=\"\\t\";} {print \$1,\$2,\$3,\"Alu\",0,\".\",\$2,\$3,\"204,102,0\"}' | uniq | gzip > $Alu_file";
   system($command);
} else {
   warn "$Alu_file already exists; skipping Alu step\n";
}

if (!-e $promoter_file){
   warn "Extracting immediate promoters\n";
   my $command = "gunzip -c $regionPromoterFile | bedtools sort | bedtools merge -i stdin | awk 'BEGIN{OFS=\"\\t\";}; {if(\$2<$start) \$2=$start}; {if(\$3>$end) \$3=$end}; \$3>=$start && \$2<=$end' | awk 'BEGIN{OFS=\"\\t\";} {print \$1,\$2,\$3,\"promoter\",0,\".\",\$2,\$3,\"127,0,255\"}' | uniq | gzip > $promoter_file";
   system($command);
} else {
   warn "$promoter_file already exists; skipping promoter step\n";
}

if (!-e $enhancer_file){
   warn "Extracting active, in-vivo transcribed enhancers\n";
   my $command = "gunzip -c $regionEnhancerFile | bedtools sort | bedtools merge -i stdin | awk 'BEGIN{OFS=\"\\t\";}; {if(\$2<$start) \$2=$start}; {if(\$3>$end) \$3=$end}; \$3>=$start && \$2<=$end' | awk 'BEGIN{OFS=\"\\t\";} {print \$1,\$2,\$3,\"enhancer\",0,\".\",\$2,\$3,\"127,0,127\"}' | uniq | gzip > $enhancer_file";
   system($command);
} else {
   warn "$enhancer_file already exists; skipping enhancer step\n";
}

if (!-e $functional_element_file){
   warn "Extracting putative functional elements\n";
   my $command = "gunzip -c $regionSegwayFile | bedtools sort | bedtools merge -i stdin | awk 'BEGIN{OFS=\"\\t\";}; {if(\$2<$start) \$2=$start}; {if(\$3>$end) \$3=$end}; \$3>=$start && \$2<=$end' | awk 'BEGIN{OFS=\"\\t\";} {print \$1,\$2,\$3,\"functional_element\",0,\".\",\$2,\$3,\"102,178,255\"}' | uniq | gzip > $functional_element_file";
   system($command);
} else {
   warn "$functional_element_file already exists; skipping functional element step\n";
}

if (!-e $crm_file){
   warn "Extracting putative cis-regulatory modules\n";
   my $command = "gunzip -c $regionReMapFile | bedtools sort | bedtools merge -i stdin | awk 'BEGIN{OFS=\"\\t\";}; {if(\$2<$start) \$2=$start}; {if(\$3>$end) \$3=$end}; \$3>=$start && \$2<=$end' | awk 'BEGIN{OFS=\"\\t\";} {print \$1,\$2,\$3,\"crm\",0,\".\",\$2,\$3,\"255,0,127\"}' | uniq | gzip > $crm_file";
   system($command);
} else {
   warn "$crm_file already exists; skipping cis-regulatory module step\n";
}

#check file generation
if(-d $labelFileDir){
        print "\nLabels successfully generated and saved as '.bed.gz' files to '$labelFileDir.'\n";
}

exit();
