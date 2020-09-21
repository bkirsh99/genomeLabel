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
use experimental;
use once;
use Scalar::Util 'blessed';
#use Bio::Tools::Run::BEDTools;

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
#my @regulator = split(/,/,join(',',$regulator)) if length $regulator;
$repClass = [split(/,/, $repClass)] if length $repClass;
$repFamily = [split(/,/, $repFamily)] if length $repFamily;
$repName = [split(/,/, $repName)] if length $repName;
#print Dumper @regulator;
#print Dumper @repClass;
my $repeat = { repClass => $repClass, repFamily => $repFamily, repName => $repName }; 
#print Dumper $repeat;

my $options = { region => $region, chr => $chr, start => $start, 
	end => $end, biotype => $biotype, regulator => $regulator, 
	repeat => $repeat }; 
#print Dumper $options;

#my $options = { region => $region, chr => $chr, start => $start,
 #     end => $end, biotype => $biotype, regulator => @regulator,
  #    repClass => @repClass, repFamily => @repFamily, repName => @repName };

#collect_All();
my ($allOriginal, $allbedFeatures) = collect_All();
#print Dumper %$allbedFeatures;
#peek(%$allOriginal);
#peek(%$allbedFeatures);

my @ncbi = collect_by_Source('ncbiRefSeq', %$allbedFeatures);
#print Dumper @ncbi;
parse_ncbiRefSeq(@ncbi);
#make_Bed(@ncbi);

#my $test = bedFeature->new(%options);
#print Dumper %$options;
#print Dumper @$options;
#my $sampleoriginal = sampleoriginal(%$allOriginal);
#print Dumper $sampleoriginal;
#my $name = "remap";
#recieve($options, $name);

sub recieve{
	my $options = shift;
	if (@_){
		print Dumper @_;
	}
}
#my $test = bedFeature->new($sampleoriginal, $options, $name);

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
=begin comment

	my $chromSizes = fetchFeature->new({
	name => 'chromSizes',
	url => 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes',
	fixed => 1,
	});
print "done\n";
	my $liftOver = fetchFeature->new({
	name => 'liftOver',
	url => 'http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz',
	fixed => 1,
	});
print "done\n";
=end comment
=cut
	my $ncbiRefSeq = fetchFeature->new({
        name => 'ncbiRefSeq',
        url => 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz',
	});
print "done\n";

=begin comment
	my $RepeatMasker = fetchFeature->new({
        name => 'RepeatMasker',
        trackName => 'rmsk',
	});
print "done\n";
	my $ReMap = fetchFeature->new({
        name => 'ReMap',
        trackName => 'ReMap',
        hub => 'http://remap.univ-amu.fr/storage/public/hubReMap2020UCSC/hub.txt',
	});
print "done\n";
	my $FANTOM5 = fetchFeature->new({
        name => 'FANTOM5',
        trackName => 'robustPeaks',
        hub => 'https://fantom.gsc.riken.jp/5/datahub/hub.txt',
	});
print "done\n";
	my $SegWay = fetchFeature->new({
        name => 'SegWay',
        trackName => join('',$biotype, '_Segway_segmentation'),
        hub => 'http://ftp.ebi.ac.uk/pub/databases/ensembl/encode/integration_data_jan2011/hub.txt',
        genome => 'hg19',
	});	 
	print "done\n";
=end comment
=cut

#=begin comment	
#	push my @fetchFeatures, ($ncbiRefSeq, $RepeatMasker, $ReMap, $FANTOM5, $SegWay);
	push my @fetchFeatures, ($ncbiRefSeq);
	foreach my $feature (@fetchFeatures){
		$feature->make_AoH;
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
#=end comment
#=cut

}


sub collect_by_Source{
        my ($source, %hash) = @_;
        foreach my $key (keys %hash){
                if ($key eq $source){
			return @{ $hash{$key} };
        	}
	}
}

sub make_Bed{
	my $file = shift;
	my $array = shift;
	#my $file = shift;

	open(FH, '>', $file) or die "Could not open file '$file'\n";
	foreach (@{ $array }){
		print FH $_->bed_string;
	}
	close(FH);
}

sub print_allSources{
	my %HobedFeatures = @_;
	foreach my $key (keys %HobedFeatures){
		print "$key\n";
	}
}

=begin comment
sub make_Labels{
        my $array1 = shift; #original hash array
        my $options = shift;
	my @array2; #bedFeature hash array

        foreach (@AoH) {
                push (@AobedFeatures, bedFeature->new($_));
        }
        
	foreach (@AobedFeatures){
                print Dumper @AobedFeatures unless defined $_->get_original_entry->{SRC};
		if ($_->get_original_entry->{SRC} eq 'ncbiRefSeq'){
                        print "refseq\n";
                }
        }
}
=end comment
=cut

sub peek{
	my %hash = @_;
        foreach my $key (keys %hash){
		if ( blessed %hash && %hash->isa("bedFeature") ) { #not working
			print "This is what the $key bedFeature object looks like:\n";
		} else {
			print "This is what the orginal $key entry looks like:\n";
		}
		print Dumper $hash{$key}[0]; 
        }
}

sub sampleoriginal{
        my %hash = @_;
        #return $hash{ncbiRefSeq}[0];
#	return $hash{ReMap}[0];
#	return $hash{RepeatMasker}[0];
	#return $hash{FANTOM5}[0];
	#return $hash{SegWay}[0];
}

sub next_bedFeature{
	my @AobedFeatures = @_; 
}

sub parse_ncbiRefSeq{
        my @array = @_; #array of bedFeatures
	my (@exon, @transcript, @coding, @noncoding); #array of parsed bedFeatures
	my %hash = ( LABELS->{EXON}[0] => \@exon, LABELS->{TRANSCRIPT}[0] => \@transcript,
		LABELS->{CODING}[0] => \@coding, LABELS->{NONCODING}[0] => \@noncoding );	

	foreach ( @array ) {
		if ($_->{original}->{TYPE} eq 'exon'){
			push @exon, $_;
			#$file = LABELS->{EXON}[0];
		}
                elsif ($_->{original}->{TYPE} eq 'transcript'){
			push @transcript, $_;
			#$file = LABELS->{TRANSCRIPT}[0];
#make transcript bed
#call bedtools here to sort exon, transcript, and make intron b$#
                }
		elsif ($_->{original}->{TYPE} =~ /CDS|codon/i){
			push @coding, $_;
			#$file
#make coding bed
                }
               elsif ($_->{original}->{TYPE} eq 'UTR' || ( $_->{original}->{TYPE} eq 'exon' && $_->{original}->{ATTR}->{Parent} =~ /NR/i) ){
			push @noncoding, $_;
#make noncoding bed
                }
	}

	while (my ($key, $value) = each (%hash)){
		make_Bed($key, $value);
#		print Dumper $key;
#		print Dumper $value;
	}

}

=begin comment
sub apply_bedTools{
        my $command = shift;
        my $file1 = shift;
        if (@_) {
		my $file2 = $_[0];
        }

#filter self->original->type = exon, makebed (perm) -> bedtools sort, merge
#filter self->original->type =5tramscvript, make bed (temp_) -> bedtools sort, merge, subtract -a transcript -b exon -> makebed (perm = intron) & bedtools complement -i transcript -g chromsizes -> makebed (perm = intergenic)
#filter self->original->type = cds, start, stop, makebed (pperm) -> bedtools sort, merge (coding)
#filter self->original->type = utr || type=exon &&id /NR/, makebed (perm) -> bedtools sort, merge = noncoding

#start with file1 = exon bed, file2 = transcript bed (temporary)
my $bedtools_sort = Bio::Tools::Run::BEDTools->new( -command => 'sort' );
my $bedtools_merge = Bio::Tools::Run::BEDTools->new( -command => 'merge' );
my $bedtools_subtract = Bio::Tools::Run::BEDTools->new( -command => 'subtract' );
my $bedtools_complement = Bio::Tools::Run::BEDTools->new( -command => 'complement' );
my $bedtools_coverage = Bio::Tools::Run::BEDTools->new( -command => 'coverage' );

#generate introns
#my $result_file = $bedtools_subtract->run( -bgv1 => 'transcripts.bed', -bgv2 => 'exons.bed' );

#generate intergenic
#my $result_file = $bedtools_complement->run( -bgv => 'transcripts.bed', -genome => 'hg38.chrom.sizes')

$result_file = $bedtools_sort->run( -bgv => $stdin);
$result_file = $bedtools_merge->run( -bgv => $stdin);

}
=end comment
=cut
