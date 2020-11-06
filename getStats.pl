#!/usr/bin/perl

use JSON;
use strict;
use warnings;
use Data::Dumper;
use 5.010;
use 5.005;
use Storable;
use Readonly;
use Text::Table;
use List::Util qw(max sum);
use Statistics::Basic qw(:all);
use lib '.';
use bedFeature;

my $dir = 'chrX:15200000-15600000';
my $biotype = 'K562';
my ($chr,$start,$end) = split /[:-]+/, $dir;
#print Dumper $chr,$start,$end;

my $statdir = join('/', $dir, 'stat_files'); 
if (!-d $statdir){
        mkdir $statdir or die "Failed to create stats directory: $!\n";
}
my $stats = join('/',$dir,'stats.txt');

# > Comppute statistics:
unless (-e $stats){ 
	my $r = create_empty_file($stats);
        die $r if $r;
	#stats();
} else {
	#%fileStats = %{ retrieve($stats) };
}
#my $statsfiles = select_files($dir);
#print Dumper $statsfiles;
my %elements;
stats($dir);

sub create_empty_file {       
	eval {
	open my $fh, '>', $_[0]
        or die "Cannot create $_[0]: $!\n";
        close $fh or die "Cannot close $_[0]: $!\n";
        };
        return $@;
}

sub stats{
	my $dir = shift;
	my $catfile = join('/',$dir,join('_',$biotype,'cat.bed'));

	open(FH, '<', $catfile) or die "Can't open '$catfile': $!";
        while(<FH>){
                my $string = $_;
                chomp $string;
                my @fields = split('\t', $string);
                my ($chrom, $start, $end, $name) = @fields[0,1,2,3];
		my @id = split(':', $name);
		my $categorydir = join ('/',$statdir, $id[0]);
		if (!-d $categorydir){
 			mkdir $categorydir or die "Failed to create directory: $!\n";
		}
		for (my $index = 1; $index <= $#id; $index++){
                       	my $groupdir = $#id == 1 ? $categorydir : join ('/', $categorydir, "group" . "$index");
			mkdir $groupdir or die "Failed to create directory: $!\n" unless (-d $groupdir);
			my $class = $id[$index];
			my $classfh = join ('/', $groupdir, $class);
			#my $classdir = join ('/', $groupdir, $class);
			#mkdir $classdir or die "Failed to create directory: $!\n" unless (-d $classdir);
 			open(OUT,'+>>',$classfh);
			print OUT $string . "\n";
		}
 
		#push @{ $elements{$category} }, \@fields; 
	}
	close(FH);

	#print Dumper \%elements;

}
