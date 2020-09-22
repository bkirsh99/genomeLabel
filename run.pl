#!/usr/bin/perl

use JSON;
use strict;
use warnings;
use Data::Dumper;
use Bio::ToolBox::parser::gff;
use 5.010;
use 5.005;
use Getopt::Long;
use Pod::Usage;
use Bio::Tools::Run::BEDTools;


use vars qw($region $chr $start $end);
use constant LABELS => {
        EXON => ['exon.bed', 'Exonic Region', '0,255,0', 'Green'],
        INTRON => ['intron.bed', 'Intronic Region', '255,0,0', 'Red'],
        INTERGENIC => ['intergenic.bed', 'Intergenic Region', '0,0,255', 'Blue'],
        CODING => ['coding.bed', 'Coding Exonic Region', '0,102,0', 'Dark Green'],
        NONCODING => ['noncoding.bed', 'Non-Coding Exonic Region', '153,255,153', 'Light Green'],
        LINE => ['LINE Element', '255,255,0', ''],
        REPEAT => ['repeat.bed'],
        SINE => ['SINE Element', '255,204,153', ''],
        ALU => ['Alu Element', '204,102,0', ''],
        PROMOTER => ['promoter.bed', 'Promoter Element', '127,0,255', ''],
        FUNC_ELEMENT => ['f_element.bed', 'Functional Element', '102,178,255', ''],
        CIS_REG_MODULE => ['cr_module.bed', 'Cis-Regulatory Module', '255,0,127', ''],
        TRANSCRIPT => ['transcript.bed'],
};

use lib '.';
use fetchFeature qw($region $chr $start $end);
use bedFeature;

GetOptions(
  'biotypes:s' => \ ( my $biotype ),
  'regulator:s' => \ ( my $regulator ),
  'repClass:s' => \ ( my $repClass ),
  'repFamily:s' => \ (my $repFamily ),
  'repName:s' => \ (my $repName ),
  'makeTracks' => \ ( my $makeTracks = 0 ),
  'getStats' => \ ( my $getStats = 0 ),
  'help' => \ ( my $help = 0 ),
) or die "Invalid options passed to $0\n";

pod2usage(2) if $help;

die "ERROR: $0 requires user input <chr:start-stop>\n" unless @ARGV;
($region, $chr, $start, $end) = getRegion(@ARGV);

$regulator = [split(/,/, $regulator)] if length $regulator;
$repClass = [split(/,/, $repClass)] if length $repClass;
$repFamily = [split(/,/, $repFamily)] if length $repFamily;
$repName = [split(/,/, $repName)] if length $repName;
my $repeat = { repClass => $repClass, repFamily => $repFamily, repName => $repName };

my $options = { region => $region, chr => $chr, start => $start,
        end => $end, biotype => $biotype, regulator => $regulator,
        repeat => $repeat };

my (@exon, @transcript, @intron, @intergenic, @coding, @noncoding, @promoter, @repeat, @f_element, @cr_module); #arrays of parsed bedFeatures
my %file = ( LABELS->{EXON}[0] => \@exon, LABELS->{TRANSCRIPT}[0] => \@transcript,
                LABELS->{INTRON}[0] => \@intron, LABELS->{INTERGENIC}[0] => \@intergenic,
                LABELS->{CODING}[0] => \@coding, LABELS->{NONCODING}[0] => \@noncoding,
                LABELS->{PROMOTER}[0] => \@promoter, LABELS->{REPEAT}[0] => \@repeat,
                LABELS->{FUNC_ELEMENT}[0] => \@f_element, LABELS->{CIS_REG_MODULE}[0] => \@cr_module, );


#Flow:
#1. Create two hashes of arrays, one with original entries and one with bedFeature objects (keys are the source)
#The bedFeature objects are already filtered according to user options, and only the desired ones are returned
#my ($hashOriginal, $hashbedFeatures) = collect_All();

#2. Create bed files from the hash of bedFeatures
#make_Labels(%$hashbedFeatures);

sub getRegion{
        my $region = shift;
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
                die "ERROR: Invalid input region <chr:start-stop>\n";
        }
        return ($region, $chr, $start, $end);
}

sub collect_All{
        my %HoOriginal;
        my %HobedFeatures;

        my $chromSizes = fetchFeature->new({
        name => 'chromSizes',
        url => 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes',
        fixed => 1,
        });

        my $liftOver = fetchFeature->new({
        name => 'liftOver',
        url => 'http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz',
        fixed => 1,
        });

        my $ncbiRefSeq = fetchFeature->new({
        name => 'ncbiRefSeq',
        url => 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz',
        });

        my $RepeatMasker = fetchFeature->new({
        name => 'RepeatMasker',
        trackName => 'rmsk',
        });

        my $ReMap = fetchFeature->new({
        name => 'ReMap',
        trackName => 'ReMap',
        hub => 'http://remap.univ-amu.fr/storage/public/hubReMap2020UCSC/hub.txt',
        });

        my $FANTOM5 = fetchFeature->new({
        name => 'FANTOM5',
        trackName => 'robustPeaks',
        hub => 'https://fantom.gsc.riken.jp/5/datahub/hub.txt',
        });

        my $SegWay = fetchFeature->new({
        name => 'SegWay',
        trackName => join('',$biotype, '_Segway_segmentation'),
        hub => 'http://ftp.ebi.ac.uk/pub/databases/ensembl/encode/integration_data_jan2011/hub.txt',
        genome => 'hg19',
        });

        push my @fetchFeatures, ($ncbiRefSeq, $RepeatMasker, $ReMap, $FANTOM5, $SegWay);

        foreach my $feature (@fetchFeatures){
                my @AoOriginal;
                my @AobedFeatures;
                foreach ( @{ $feature->{array} } ) {
                        push @AoOriginal, $_;
                        push @AobedFeatures, bedFeature->new($_, $options);
                }
                $HoOriginal{ $feature->{name} } = [ @AoOriginal ];
                $HobedFeatures{ $feature->{name} } = [ @AobedFeatures ];
        }
        return ( \%HoOriginal, \%HobedFeatures );

}


sub to_Bed{
        my $file = shift;
        my $array = shift;

        open(FH, '>', $file) or die "Could not open file '$file'\n";
        foreach (@{ $array }){
                print FH $_->bed_string;
        }
        close(FH);
}

sub from_Bed{
        my $file = shift;
        my $array;

        open(FH, '<', $file) or die "Could not open file '$file'\n";
        while(<FH>){
                my $feature = bedFeature->new();
                $feature->from_bed_string($_);
                push @{ $array }, $feature;
        }
        close(FH);
        return $array;
}


sub make_Labels{
        my %hash = @_; #hash of arrays of bedFeatures

        foreach my $key (keys %hash){
                my @array = @{ $hash{$key} };
                foreach (@array){
                        if ($key eq 'ncbiRefSeq'){
                                if ($_->{original}->{TYPE} eq 'exon'){
                                        push @exon, $_;
                                } elsif ($_->{original}->{TYPE} eq 'transcript'){
                                        push @transcript, $_;
                                } elsif ($_->{original}->{TYPE} =~ /CDS|codon/i){
                                        push @coding, $_;
                                } elsif ($_->{original}->{TYPE} eq 'UTR' || ( $_->{original}->{TYPE} eq 'exon' && $_->{original}->{ATTR}->{Parent} =~ /NR/i) ){
                                        push @noncoding, $_;
                                }
                        } elsif ($key eq 'FANTOM5'){
                                push @promoter, $_;
                        } elsif ($key eq 'RepeatMasker'){
                                push @repeat, $_;
                        } elsif ($key eq 'SegWay'){
                                push @f_element, $_;
                        } elsif ($key eq 'ReMap'){
                                push @cr_module, $_;
                        }
                }
        }

        while (my ($key, $value) = each (%file)){
                to_Bed($key, $value);
        }

        #apply_bedTools(keys %file);
}




sub next_bedFeature{
        my @AobedFeatures = @_;

        foreach(@AobedFeatures){
                return $_;
        }
}



sub apply_bedTools{
        my $files = shift;

        my $bedtools_sort = Bio::Tools::Run::BEDTools->new( -command => 'sort' );
        my $bedtools_merge = Bio::Tools::Run::BEDTools->new( -command => 'merge' );
        my $bedtools_subtract = Bio::Tools::Run::BEDTools->new( -command => 'subtract' );
        my $bedtools_complement = Bio::Tools::Run::BEDTools->new( -command => 'complement' );
        my $bedtools_coverage = Bio::Tools::Run::BEDTools->new( -command => 'coverage' );

        foreach (@{ $files}){
                $bedtools_sort->run( -bgv => $_, -out => $_);
                if ($_ =~ /exon/ || /transcript/ || /coding/ ){
                        $bedtools_merge->run( -bgv => $_, -out => $_);
                }
        }

$bedtools_subtract->run( -bgv1 => 'noncoding.bed', -bgv2 => 'coding.bed', -out => 'noncoding.bed' );

$bedtools_subtract->run( -bgv1 => 'transcripts.bed', -bgv2 => 'exons.bed', -out => 'intron.bed' );
push @{ $files}, LABELS->{INTRON}[0];
$bedtools_complement->run( -bgv => 'transcripts.bed', -genome => 'hg38.chrom.sizes', -out => 'intergenic.bed');
push @{ $files}, LABELS->{INTERGENIC}[0];

}

__END__

=head1 NAME

genomeLabel - Label genomic features in a region.

=head1 SYNOPSIS

        ./genomeLabel.pl [--biotype] [--regulator] [--repeat] [--makeTracks] [--makeHubs] [--getStats] chr:start-stop

        Options:
        --biotype       One of the six human cell types used by the ENCODE Consortium: GM12878, H1-hESC, K562, HeLa-S3, HepG2, HUVEC.
                        A description of these cell types is available at L<http://genome.ucsc.edu/cgi-bin/hgEncodeVocab?term=GM12878,K562,H1-hESC,HeLa-S3,HepG2,HUVEC>

        --regulator     A list of one or more of the 960 DNA-binding proteins, including transcription factors (TFs), transcription co-activators (TCFs) and chromatin-remodeling factors (CRFs) used by the ReMap Atlas.
                        Access to ReMap documentation is available at L<http://remap.univ-amu.fr>

        --repeat        A list of one or more of the families, classes, or names of repeats used by RepeatMasker.
                        Access to RepeatMasker documentation is available at L<http://www.repeatmasker.org/webrepeatmaskerhelp.html>

        --makeTracks    Generate Genome Browser tracks.

        --getStats      Generate summary statistics.

        chr:start-stop  Region used in feature annotation.

=head1 OPTIONS

=over 4

=item B<--biotype>

Select cell line or default to cell type-agnostic annotation.

=item B<--regulator>

Select transcriptional regulator or default to regulator-agnostic annotation.

=item B<--repeat>

Select repeat or default to repeat type-agnostic annotation.

=item B<--makeTracks>

Generate tracks.

=item B<--getStats>

Generate statistics.

=itemB<chr:start-stop>

Input chosen region.

=back

=head1 DESCRIPTION

B<genomeLabel> will annotate the genomic features within the given region and work with the contents thereof.

=cut
