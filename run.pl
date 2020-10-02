#!/usr/bin/perl

use strict;
use warnings;
use 5.010;
use 5.005;
use Getopt::Long;
use Pod::Usage;
use File::Which;
use Storable;
use lib '.';

use Data::Dumper;

#[0] is the spelling used in the SegWay hub url, [1] is the spelling used in the ReMap biotype field
use constant BIOTYPES => {
        GM12878 => ['GM12878', 'GM12878'],
        H1HESC => ['H1hESC', 'WA01'],
        HELAS3 => ['HeLaS3', 'HeLa-S3'],
        HEPG2 => ['HepG2', 'Hep-G2'],
        HUVEC => ['HUVEC', 'Huvec-C'],
        K562 => ['K562', 'K-562'],
};

#get [OPTIONS] from command-line: ./run.pl [chr:start-stop] --biotype [BIOTYPE] --path [/PATH/] <command(s)> <filter(s)>
GetOptions(
  'path=s' => \ (my $path), #required by liftOver
  'biotypes=s' => \ (my $biotype), #required by SegWay and ReMap
  'regulator:s' => \ (my $regulator),
  'repeat:s' => \ (my $repeat),
  'felement:s' => \ (my $felement),
  'makeLabels' => \ (my $makeLabels), #default
  'makeTracks' => \ (my $makeTracks),
  'makeHub' => \ (my $makeHub),
  'getStats:s' => \ (my $getStats),
  'help' => \ (my $help = 0 ),
) or die "Invalid options passed to $0\nPlease refer to the help documentation\n\n";

pod2usage(2) if $help; #help message

die "ERROR: $0 requires user input <chr:start-stop> --path <path> --biotype <biotype>\n" unless @ARGV and $path and $biotype; #minimum requirements: region, path, and biotype

#check for bedtools & bedops installation
my ($bedtools, $bedops) = ( which('bedtools'), which('bedops') );

if (!$bedtools){
   print STDERR "bedtools was not found in your path:\n\n$ENV{PATH}\n\nPlease install bedtools and add it to your path:\n\n";
   print STDERR "git clone https://github.com/arq5x/bedtools2.git\ncd bedtools2\nmake clean\nmake all\n\n";
   exit(1);
}

if (!$bedops){
   print STDERR "bedops was not found in your path:\n\n$ENV{PATH}\n\nPlease install bedops and add it to your path:\n\n";
   print STDERR "https://bedops.readthedocs.io/en/latest/\n";
   exit(1);
}
#verify validity of required arguments
my ($region, $chr, $start, $end) = getRegion(@ARGV);
my ($biotypeS, $biotypeR) = getBiotype(BIOTYPES);
#make input hash for checkHistory()

my $input = { $biotype => { $chr => [$start,$end] } };

foreach ($regulator, $repeat, $felement, $getStats){
	$_ = [split(/,/, $_)] if defined $_;
}

my @filter = [$biotypeR, $regulator, $repeat, $felement];

###DEBUG###
#	print Dumper $region, $chr, $start, $end, $biotypeS, $biotypeR;
#	print Dumper $biotype, $path, $regulator, $repeat, $felement, $makeLabels, $makeTracks, $makeHub, $getStats, $help; 
	print Dumper @filter;

my $history = "history.txt";
#updateHistory();

if ($makeTracks){
	#system "./makeTracks.pl $history $path $region $chr $start $end $biotype @filter"; 
}

#verify and extract region
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

#verify and extract biotype
sub getBiotype{
        my $ref = shift;

        #get rid of non alphanumeric characters
        $biotype =~ s/\W//g;
        #check biotype validity
        foreach ( sort keys %{ $ref } ){
                if ($biotype =~ /^$_$/i){
                        $biotype = $_;
                        return ($ref->{$_}[0], $ref->{$_}[1]);
                }
        }
        die "ERROR: Invalid --biotype. Use --help for help.\n";
}

#check for previous downloads and add data, if needed
sub updateHistory{

        unless (!-f $history){
		my $retrieve = retrieve($history); #stores previous inputs (ARRAY OF TWO HASHES) 
		my @pinput = @$retrieve[0]; #previous inputs
		my $init = getNewRegion(\@pinput); #Hash of two keys: first one is general (may contain multiple arrays), second one is cell-specific, i.e. Segway (may contain multiple arrays)
		print Dumper $init;
		if ( scalar(@{ $$init{'all'} }) == scalar(@{ $$init{'biotype'} }) ){
			my ($loop_index, $s1, $e1, $s2, $e2);
			for($loop_index = 0; $loop_index <= $#{ $$init{'all'} }; $loop_index++) {
			print "compare regions\n";
				unless (@{ $$init{'all'}[$loop_index] }) { 
					unless (@{ $$init{'biotype'}[$loop_index] }) {
						print "nothing to do - both empty\n";
					} else {
					$s1 = $$init{'biotype'}[$loop_index][1];
					$e1 = $$init{'biotype'}[$loop_index][2];
					} 
				} else {
					unless (@{ $$init{'biotype'}[$loop_index] }) {
					$s1 = $$init{'all'}[$loop_index][1];
					$e1 = $$init{'all'}[$loop_index][2];
					} else {
						if ($$init{'all'}[$loop_index][1] == $$init{'biotype'}[$loop_index][1]) { 
							$s1 = $$init{'all'}[$loop_index][1];
						} else {
							$s1 = $$init{'all'}[$loop_index][1];
							$s2 = $$init{'biotype'}[$loop_index][1];
						}
						if ($$init{'all'}[$loop_index][2] == $$init{'biotype'}[$loop_index][2]) {
							$e1 = $$init{'all'}[$loop_index][2];
						} else {
							$e1 = $$init{'all'}[$loop_index][2];
                                                        $e2 = $$init{'biotype'}[$loop_index][2];
						}
					}	
				}
			#print "do not exist\n" unless $s2, $e2;
			#system "./init.pl $history $chr $s1 $e1 $biotypeS"
			#system "./init.pl $history $chr $s2 $e2 $biotypeS" if $s2, $e2;
			#print Dumper $s1, $e1;
			}
		} else {
			print Dumper $$init{'biotype'}[0], $$init{'biotype'}[1], $$init{'biotype'}[2], $$init{'biotype'}[3];
		}
        } else {
                #system "./init.pl $history $chr $start $end $biotypeS";
                print "does not\n";
        }

}

my (%eval, $start2, $end2, @new, @newbiotype, $new1, $new2, %seen);
##CHANGE!!!!!!!!
sub getNewRegion{
        #my $old = shift;
	my $old = [ { 'K562' => { 'chrX' => [ '15200000','15600000' ] } } ,
	{ 'K562' => { 'chrX' => [ '1520000','15800000' ] } } ,
	{ 'K562' => { 'chrY' => [ '15200000','15600000' ] } } ,
        { 'HUVEC' => { 'chrX' => [ '15200000','15600000' ] } } ,
	{ 'HUVEC' => { 'chrX' => [ '15000000','15800000' ] } } ,
	{ 'HUVEC' => { 'chrX' => [ '13000000','15800000' ] } } ,
	];

	foreach my $pinput (@$old){ #old input hash
		foreach my $biotyp (keys %$pinput){
			unless (@newbiotype){
				push @newbiotype, $biotype unless defined $seen{$biotype};
				push @newbiotype, $chr unless defined $seen{$chr};
				push @newbiotype, $start unless defined $seen{$start};
				push @newbiotype, $end unless defined $seen{$end};
				%seen = map {$_ => 1} @newbiotype;
			}
			foreach my $chromosome (keys %{ $$pinput{$biotyp} } ){
				if ($chr eq $chromosome) {
					my @coordinates = @{ $$pinput{$biotyp}{$chromosome} };
					@new = adjustCoord( @coordinates );
					@newbiotype = adjustCoord( @coordinates ) if ($biotyp eq $biotype);
				}
			}
		}
	}
	$eval{'all'} = \@new;
	$eval{'biotype'} = \@newbiotype;
	return \%eval;
}

sub adjustCoord{
	my @coordinates = @_;
	my @new;

	if ($start < $coordinates[0] || $end > $coordinates[1]){
		$end2 = $start < $coordinates[0] ? $coordinates[0] - 1 : undef;
                $start2 = $end > $coordinates[1] ? $coordinates[1] + 1 : undef;
                $new1 = defined $start2 ? [ $chr, $start2, $end ] : [];
                $new2 = defined $end2 ? [ $chr, $start, $end2 ] : [];
               @new = ($new1, $new2);
	}
	return @new;
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
