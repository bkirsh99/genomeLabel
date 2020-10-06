#!/usr/bin/perl

use JSON;
use strict;
use warnings;
use Data::Dumper;
use 5.010;
use 5.005;
use Config::Simple;
use File::Sort qw( sort_file );
use Storable;
use lib '.';
use vars qw($chr $start $end);
use fetchFeature qw($chr $start $end);
use bedFeature;
#collect arguments sent by ./run.pl
my $history = shift or die;
my $chr = shift or die;
my $start = shift or die;
my $end = shift or die;
my $biotype = shift or die;

#retrieve old history and previously downloaded data
my $retrieve = retrieve($history); #stores previous inputs (AS AN ARRAY OF HASHES FOR EACH BIOTYPE => CHROMOSOME => COORDINATES) and fetchFeature objects (AS A HASH FOR EACH DATABASE) - ARRAY OF SIZE 2
my $pinput = @$retrieve[0]; #previous inputs (HASH)
my $phash = @$retrieve[1]; #previous bedFeature objects (HASH)

#add new history and new downloaded data
my $ninput = appendInput();
my $nhash = appendHoB();
my @store = ($ninput, $nhash);

#overwrite 'history.txt' file with new plus old stuff
store(\@store, $history);

#download data
sub collect_All{
	my %HobedFeatures;	
	print "\nCollecting data - this might take a while...\n\n";
	
	my $chromSizes = fetchFeature->new({
	name => 'chromSizes',
	chrom => $chr,
        chromStart => $start,
        chromEnd => $end,
	url => 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes',
	fixed => 1,
	});

	my $liftOver = fetchFeature->new({
	name => 'liftOver',
	chrom => $chr,
        chromStart => $start,
        chromEnd => $end,
	url => 'http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz',
	fixed => 1,
	});

	my $ncbiRefSeq = fetchFeature->new({
        name => 'ncbiRefSeq',
	chrom => $chr,
        chromStart => $start,
        chromEnd => $end,
        url => 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz',
	});

	my $RepeatMasker = fetchFeature->new({
        name => 'RepeatMasker',
	chrom => $chr,
        chromStart => $start,
        chromEnd => $end,
        trackName => 'rmsk',
	});

	my $ReMap = fetchFeature->new({
        name => 'ReMap',
	chrom => $chr,
        chromStart => $start,
        chromEnd => $end,
        trackName => 'ReMap',
        hub => 'http://remap.univ-amu.fr/storage/public/hubReMap2020UCSC/hub.txt',
	});

	my $FANTOM5 = fetchFeature->new({
        name => 'FANTOM5',
	chrom => $chr,
        chromStart => $start,
        chromEnd => $end,
        trackName => 'robustPeaks',
        hub => 'https://fantom.gsc.riken.jp/5/datahub/hub.txt',
	});

	my $SegWay = fetchFeature->new({
        name => 'SegWay',
	chrom => $chr,
        chromStart => $start,
        chromEnd => $end,
        trackName => join('',$biotype, '_Segway_segmentation'),
        hub => 'http://ftp.ebi.ac.uk/pub/databases/ensembl/encode/integration_data_jan2011/hub.txt',
        genome => 'hg19',
	});	 

	push my @fetchFeatures, ($ncbiRefSeq, $RepeatMasker, $ReMap, $FANTOM5, $SegWay);

	foreach my $feature (@fetchFeatures){
		my @AobedFeatures;
		foreach ( @{ $feature->{array} } ) {
			push @AobedFeatures, bedFeature->new($_);
		}
		$HobedFeatures{ $feature->{name} } = [ @AobedFeatures ];
	}
	
	return \%HobedFeatures;
}

sub appendInput{
	my @coords = [$start , $end];
	foreach my $btype (keys %$pinput) {
		unless( defined $$pinput{$biotype} ){
			$$pinput{$biotype} = { $chr => [@coords] };
		} else {
			foreach my $ch ( keys %{ $$pinput{$biotype} } ){
				unless( defined $$pinput{$biotype}{$chr} ){
					$$pinput{$biotype}{$chr} = [@coords];
				} else {
					push @{ $$pinput{$biotype}{$chr} }, @coords;
				}
			}
		}
	}
	return $pinput;
}

sub appendHoB{
	my $new = collect_All();

	foreach my $source (keys %$phash){
		if (@{ $$new{$source} }){
			foreach (@{ $$new{$source} }){
				push @{ $$phash{$source} }, $_;
			}
		}
	}

	return $phash;
}

#helper subroutines (for debugging maybe)
sub collect_by_Source{
        my ($source, %hash) = @_;
        foreach my $key (keys %hash){
                if ($key eq $source){
			return @{ $hash{$key} };
        	}
	}
}

sub print_allSources{
	my %HobedFeatures = @_;
	foreach my $key (keys %HobedFeatures){
		print "$key\n";
	}
}


sub peek{
	my %hash = @_;
        foreach my $key (keys %hash){
		print Dumper $hash{$key}[0]; 
        }
}

sub sampleoriginal{
        my %hash = @_;
        #return $hash{ncbiRefSeq}[0];
	#return $hash{ReMap}[0];
	#return $hash{RepeatMasker}[0];
	#return $hash{FANTOM5}[0];
	#return $hash{SegWay}[0];
}

