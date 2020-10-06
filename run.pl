#!/usr/bin/perl

use strict;
use warnings;
use 5.010;
use 5.005;
use Getopt::Long;
use Pod::Usage;
use File::Find;
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

pod2usage(-verbose => 3) if $help; #help message

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

#make input hash to add to history
#my $input = { $biotype => { $chr => [$start,$end] } };

foreach ($regulator, $repeat, $felement, $getStats){
	$_ = [split(/,/, $_)] if defined $_;
}

my $filter = [$biotypeR, $regulator, $repeat, $felement];

#make unique "session id" to identify directory in "./out"
my @strarray;

if (@$filter){
	foreach (@$filter){
        	if (ref($_) eq 'ARRAY'){
                	my $arrstr = join('_', @$_) if defined $_;
                	push @strarray, $arrstr;
        	}
	}
}
unshift @strarray, $biotype;

my $filterstr = join ('.', @strarray) if @strarray;

#establish directory names to store files and create, if needed
my $out = join('@',$region,$biotype);
if (!-d $out){
        mkdir $out or die "Failed to create output directory: $!\n";
}
my $trackdir = join('/',$out,'tracks');
if (!-d $trackdir){
        mkdir $trackdir or die "Failed to create track file directory: $!\n";
}
my $labeldir = join('/',$out,'labels'); #this is where labels (i.e. unfiltered files) will be stored
my $filterdir = join('/',$trackdir,$filterstr); #this is where tracks for this session will be stored (filter-specific)

my $history = "history.txt";

#> STEP 1: Check history to verify if data exists for this region/biotype
# - Call init.pl via updateHistory() to either retrieve data from "history.txt" or extract and append new data from databases
 
updateHistory();
#CHECK:
my $retrieve = retrieve($history);
my @input = @$retrieve[0];
my @hash = @$retrieve[1];
print Dumper \@input, \@hash;


#> STEP 2: Fulfill user's desires now that we have data for the input region and chosen biotype: 
if ( $getStats || (!$getStats && !$makeTracks) ){
       	mkdir $labeldir or die "Failed to create region/biotype-specific label file directory: $!\n" unless (-d $labeldir);
	#system "./makeLabels.pl $history $path $biotypeR $labeldir" if (emptyDir($labeldir)); 
	#system "./getStats.pl $labeldir \@options" unless !$getStats;
}

if ($makeTracks){
	mkdir $filterdir or die "Failed to create filter-specific track file directory: $!\n" unless (-d $filterdir);
        #system "./makeTracks.pl $history $path $filterdir $filter" if (emptyDir($filterdir));
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
		if ( scalar(@{ $$init{'all'} }) == scalar(@{ $$init{'biotype'} }) ){
			my ($loop_index, $s1, $e1, $s2, $e2);
			for($loop_index = 0; $loop_index <= $#{ $$init{'all'} }; $loop_index++) { #compare regions
				unless (@{ $$init{'all'}[$loop_index] }) { 
					unless (@{ $$init{'biotype'}[$loop_index] }) { #both empty, no data needed
						return;
					} else {
					$s1 = $$init{'biotype'}[$loop_index][1];
					$e1 = $$init{'biotype'}[$loop_index][2];
					} 
				} else {
					unless (@{ $$init{'biotype'}[$loop_index] }) {
					$s1 = $$init{'all'}[$loop_index][1];
					$e1 = $$init{'all'}[$loop_index][2];
					} else {
						my ($s2, $e2);
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
			system "./init.pl $history $chr $s1 $e1 $biotypeS";
			system "./init.pl $history $chr $s2 $e2 $biotypeS" if $s2, $e2;
			}
		} else { #Segway and rest require different dowload commands
			system "./init.pl $history $$init{'biotype'}[1] $$init{'biotype'}[2] $$init{'biotype'}[3] $biotypeS 0"; #ARGV = 6 (download only Segway)
			foreach (@{ $$init{'all'} }){
				system "./init.pl $history $$_[0] $$_[1] $$_[2]"; #ARGV = 4 (download all except Segway)
			}  
		}
        } else { #no history file at all, must initialize everything
                system "./init.pl $history $chr $start $end $biotypeS"; #ARGV = 5 (download all)
        }

	appendInput();
}

sub appendInput{
	my $retrieve = retrieve($history);
	my $pinput = @$retrieve[0];
	my $nhash = @$retrieve[1];

        my @coords = [$start , $end];
        foreach my $btype (keys %$pinput) {
                unless( defined $$pinput{$biotype} ){
                        $$pinput{$biotype} = { $chr => [@coords] };                                                                                                                                 } else {
                        foreach my $ch ( keys %{ $$pinput{$biotype} } ){                                                                                                                                            unless( defined $$pinput{$biotype}{$chr} ){
                                        $$pinput{$biotype}{$chr} = [@coords];
                                } else {
                                        push @{ $$pinput{$biotype}{$chr} }, @coords;
                                }
                        }
                }
        }
#	print Dumper $pinput;
	my @store = ($pinput, $nhash);
	#overwrite 'history.txt' file with new plus old stuff
	store(\@store, $history);
}

sub emptyDir{
        my $dirname = shift;
        opendir(my $dh, $dirname) or die "Not a directory";
        return scalar(grep { $_ ne "." && $_ ne ".." } readdir($dh)) == 0;
}

my (%eval, $start2, $end2, @new, @newbiotype, $new1, $new2, %seen);

sub getNewRegion{
	my $old = shift;
	#my $old = [ { 'HUVEC' => { 'chrX' => [ [ 15200000 , 15800000 ] , [ 15000000 , 15800000 ] ] } } ];
	
	foreach my $pinput (@$old){ #old input hash where primary key = biotype
                foreach my $biotyp (keys %$pinput){
                        unless (@newbiotype){
				push @newbiotype, $biotype unless defined $seen{$biotype};
				push @newbiotype, $chr unless defined $seen{$chr};
                                push @newbiotype, $start unless defined $seen{$start};
				push @newbiotype, $end unless defined $seen{$end};
                                %seen = map {$_ => 1} @newbiotype;
                        }
                        foreach my $chromosome (keys %{ $$pinput{$biotyp} } ){ #secondary key = chromosome
                                if ($chr eq $chromosome) {
					foreach my $coordinates ( @{ $$pinput{$biotyp}{$chromosome} } ){
                                        	@new = adjustCoord( @$coordinates );
                                        	@newbiotype = adjustCoord( @$coordinates ) if ($biotyp eq $biotype);
					}
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


=end comment
=cut

__END__

=head1 NAME

genomeLabel - Label genomic features in a region.

=head1 SYNOPSIS

	USAGE: ./run.pl [chr:start-stop] --biotype [BIOTYPE] --path [/PATH/] <command(s)> <filter(s)>
	
	The program accepts any of the following six human cell, tissue, or DNA sample used by the ENCODE Consortium: GM12878, H1-hESC, K562, HeLa-S3, HepG2, and HUVEC.
	A description of these samples is available at L<http://genome.ucsc.edu/cgi-bin/hgEncodeVocab?term=GM12878,K562,H1-hESC,HeLa-S3,HepG2,HUVEC>.
	This option is case- and delimiter-insensitive, and inputs such as H1HESC, helas3, and HEp-g2 are all valid.

	The minimum requirements to run the script are the region, biotype, and path to liftOver.
     
        Additional filter arguments can be passed alongside different commands for added specificity. These correspond to "--repeat," "--regulator," and "--felement," 
	as well as any "--getStats" option. They must be comma-separated and either enclosed within double-quotes, where whitespace is allowed, or listed subsequently
	without whitespace in between commas.
     
	ACCEPTABLE:       --regulator "GATA1,CTCF,     SOX2"
                          --regulator GATA1,CTCF,SOX2
        UNACCEPTABLE:     --regulator GATA1, CTCF, SOX2
     
        If <command(s)> is omitted, the default behaviour of the program is to generate a set of "label files." These can be useful to extract additional statistics or 
	create more complex tracks according to user needs. For simple applications, <command(s)> can be used to output select "track files," track hubs, and summary
	statistics. 

=head1 OPTIONS

=over 4

=item B<chr:start-stop>

	Input region in the format chr:start-end.

=item B<--biotype>

	One of the six human cell types used by the ENCODE Consortium: GM12878, H1-hESC, K562, HeLa-S3, HepG2, HUVEC. A description of
	these cell types is available at L<http://genome.ucsc.edu/cgi-bin/hgEncodeVocab?term=GM12878,K562,H1-hESC,HeLa-S3,HepG2,HUVEC>

=item B<--path>

	Directory path to liftOver binary file in the format /path/to/utility/ (i.e. enclosed by "/" and excluding the name of the utility itself).

=item B<--regulator>

	A list of one or more of the 960 DNA-binding proteins, including transcription factors (TFs), transcription co-activators (TCFs) and 
	chromatin-remodeling factors (CRFs) used by the ReMap Atlas. Access to ReMap documentation is available at L<http://remap.univ-amu.fr>

=item B<--repeat>

	A list of one or more of the families, classes, or names of repeats used by RepeatMasker. Access to RepeatMasker documentation is available
	at L<http://www.repeatmasker.org/webrepeatmaskerhelp.html>

=item B<--felement>

	A list of one or more of the genomic states used by Segway. Access to Segway documentation is available under "Segway Segmentations" at 
	L<http://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=913156841_uydtGuw88KR9Xqpgn3fXaMtXmsVQ&c=chr15&g=hub_4607_genomeSegmentation>

=item B<--makeLabels>

	Default behaviour. Creates unfiltered "labels" that are reported in different files: exon.bed (coding/noncoding exons), intron.bed, transcript.bed,
intergenic.bed, promoter.bed, repeat.bed, cr_module.bed, and f_element.bed. Used to compute statistical data.

=item B<--makeTracks>

	Creates Genome Browser tracks, which are filtered "labels" that are reported in different files: genome.bed (coding/noncoding exons, introns, intergenic
	regions), promoter.bed, repeat.bed, cr_module.bed, and f_element.bed.

=item B<--getStats>

	Generates summary statistics.

=item B<--help>

	Prints this message and exits successfully.

=back

=head1 DESCRIPTION

B<genomeLabel> will annotate the genomic features within the given region and work with the contents thereof.

=cut
