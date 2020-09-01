#!/usr/bin/env perl
#
# Use script to create labels for genomic elements.
#

use strict;
use warnings;
use File::Which;
use File::Basename;

my $bedtools = which('bedtools');

if (!$bedtools){
   print STDERR "bedtools was not found in your path:\n\n$ENV{PATH}\n\nPlease install bedtools and add it to your path:\n\n";
   print STDERR "git clone https://github.com/arq5x/bedtools2.git\ncd bedtools2\nmake clean\nmake all\n\n";
   exit(1);
}


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

my $outDir = "out";
my $dataDir = "data";
my $regionDir = join("",$outDir,"/",$region);
my $regionFileDir = join("",$regionDir,"/dataFiles");
my $labelFileDir = join("",$regionDir,"/labelFiles");

if (-d $labelFileDir) {
	print "\nLabels for $region have already been generated.\nPlease verify the directory '$labelFileDir'.\n\n"; 
	die;
}
else {
	if (!-d $outDir) {
		mkdir $outDir or die "Failed to create path";
	}

	if (!-d $dataDir) {
		mkdir $dataDir or die "Failed to create path";
	}


	if (!-d $regionDir) {
		mkdir $regionDir or die "Failed to create path";
	}


	if (!-d $regionFileDir) {
		mkdir $regionFileDir or die "Failed to create path";
	}

	mkdir $labelFileDir or die "Failed to create path";
}

my $geneFile = join("",$dataDir,"/hg38.ncbiRefSeq.gtf.gz");
my $repeatFile = join("",$dataDir,"/rmsk.txt.gz");
my $promoterFile = join("",$dataDir,"/hg38_fair+new_CAGE_peaks_phase1and2.bed.gz");
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

if (!-e $genome){
	system "wget -P $dataDir https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes";
}

my $regionGeneFile = join("",$regionFileDir,"/hg38.ncbiRefSeq.gtf.gz");
my $regionRepeatFile = join("",$regionFileDir,"/rmsk.txt.gz");
my $regionPromoterFile = join("",$regionFileDir,"/hg38_fair+new_CAGE_peaks_phase1and2.bed.gz");
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

if (!-e $regionGenome){
	my $command = "cat $genome | awk -F \"\\t\" 'BEGIN{OFS=\"\\t\";} (\$1==\"$chr\")' > tmp && mv tmp $regionGenome";
	system($command);
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

if (!-e $exon_file){
   warn "Creating exonic regions\n";
   my $command = "gunzip -c $regionGeneFile | awk -F \"\\t\" 'BEGIN{OFS=\"\\t\";} (\$3==\"exon\") {print \$1,\$4,\$5}' | bedtools sort | bedtools merge -i stdin | awk 'BEGIN{OFS=\"\\t\";}; {if(\$2<$start) \$2=$start}; {if(\$3>$end) \$3=$end}; \$3>=$start && \$2<=$end' | awk 'BEGIN{OFS=\"\\t\";} {print \$1,\$2,\$3,\"exon\",0,\".\",\$2,\$3,\"0,255,0\"}' | uniq | gzip > $exon_file";
   system($command);
} else {
   warn "$exon_file already exists; skipping exon step\n";
}

if (!-e $intron_file){
   warn "Creating intronic regions\n";
   my $command = "gunzip -c $regionGeneFile | awk -F \"\\t\" 'BEGIN{OFS=\"\\t\";} (\$3==\"transcript\") {print \$1,\$4,\$5}' | bedtools sort | bedtools merge -i stdin | bedtools subtract -a stdin -b $exon_file | awk 'BEGIN{OFS=\"\\t\";}; {if(\$2<$start) \$2=$start}; {if(\$3>$end) \$3=$end}; \$3>=$start && \$2<=$end' | awk 'BEGIN{OFS=\"\\t\";} {print \$1,\$2,\$3,\"intron\",0,\".\",\$2,\$3,\"255,0,0\"}' | gzip > $intron_file";
   system($command);
} else {
   warn "$intron_file already exists; skipping intron step\n";
}

if (!-e $intergenic_file){
   warn "Creating intergenic regions\n";
   my $command = "gunzip -c $regionGeneFile | awk -F \"\\t\" 'BEGIN{OFS=\"\\t\";} (\$3==\"transcript\") {print \$1,\$4,\$5}' | bedtools sort | bedtools complement -i stdin -g $regionGenome | awk 'BEGIN{OFS=\"\\t\"}; {if(\$2<$start) \$2=$start}; {if(\$3>$end) \$3=$end}; \$3>=$start && \$2<=$end' | awk 'BEGIN{OFS=\"\\t\";} {print \$1,\$2,\$3,\"intergenic\",0,\".\",\$2,\$3,\"0,0,255\"}' | gzip > $intergenic_file";
   system($command);
} else {
   warn "$intergenic_file already exists; skipping intergenic step\n";
}

if (!-e $coding_file){
   warn "Creating coding exonic regions\n";
   my $command = "gunzip -c $regionGeneFile | awk -F \"\\t\" 'BEGIN{OFS=\"\\t\";} (\$3==\"CDS\" || \$3==\"start_codon\" || \$3==\"stop_codon\")' |  bedtools sort | bedtools merge -i stdin | awk 'BEGIN{OFS=\"\\t\";}; {if(\$2<$start) \$2=$start}; {if(\$3>$end) \$3=$end}; \$3>=$start && \$2<=$end' | awk 'BEGIN{OFS=\"\\t\";} {print \$1,\$2,\$3,\"coding\",0,\".\",\$2,\$3,\"0,102,0\"}' | uniq | gzip > $coding_file";
   system($command);
} else {
   warn "$coding_file already exists; skipping coding step\n";
}

if (!-e $noncoding_file){
   warn "Creating non-coding exonic regions\n";
   my $command = "gunzip -c $regionGeneFile | awk -F \"\\t\" 'BEGIN{OFS=\"\\t\";} ((\$3~/UTR/) || (\$9~/NR/ && \$3==\"exon\")) {print \$1,\$4,\$5}' | bedtools sort | bedtools subtract -a stdin -b $coding_file | sort -k1,1 -k2,2n | bedtools merge -i stdin | awk 'BEGIN{OFS=\"\\t\";}; {if(\$2<$start) \$2=$start}; {if(\$3>$end) \$3=$end}; \$3>=$start && \$2<=$end' | awk 'BEGIN{OFS=\"\\t\";} {print \$1,\$2,\$3,\"noncoding\",0,\".\",\$2,\$3,\"153,255,153\"}' | uniq | gzip > $noncoding_file";
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

if(-d $labelFileDir){
	print "\nLabels for $region have been succesfully created and saved to '$labelFileDir'.\n\n";
}

exit();
