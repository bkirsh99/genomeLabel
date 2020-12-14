package fetchFeature;
use Bio::ToolBox::parser::gff;
use JSON;
use strict;
use warnings;
use File::Temp;
use Data::Dumper;
use IO::Compress::Gzip qw(gzip $GzipError);
my $tmp = File::Temp->new();

sub new {
    my ($class, $args) = @_;
    my $self = bless {
        name => $args->{name},
        url => $args->{url} || undef, 
        hub => $args->{hub} || undef,
        trackName => $args->{trackName} || undef,
        chrom => $args->{chrom},
        chromStart => $args->{chromStart},
        chromEnd => $args->{chromEnd},
        genome => $args->{genome} || 'hg38',
        fixed => $args->{fixed} || undef,
	}, $class;
	
	unless ($self->{fixed}){
		$self->make_AoH;
	} else {
		$self->load;
	}

   return $self;
}

sub make_AoH{

        my $self = shift;
        my @fileDump;

        if (defined $self->{trackName}){
                my $data;
                {
                        local $/ = undef;
                        open(FH, '<', $self->load) or die $!;
                        $data = <FH>;
                        close(FH);
                }
		print "\tParsing JSON object...\n\n";
                my $text = decode_json( $data );
                @fileDump = @{ $text->{ $self->{trackName} } };
        } elsif (defined $self->{url}) {
		my $gzip = $self->load;
		unless (`gunzip -c $gzip | wc -c` == 0){ #make sure awk doesn't return empty file
			my $parser = Bio::ToolBox::parser::gff->new(
                       		file => $gzip,
                       		#do_gene => 0, 
                       		#do_exon => 0, 
                       		do_cds => 1,
                       		do_utr => 1,
                       		do_codon => 1,
               		) or die "Unable to open \$parser->{file}, please ensure '.gff' or '.gtf' extension.\n";

			print "\tParsing GFF annotation file...\n\n";

               		while (my $feature = $parser->next_feature() ) { #each $feature is a SeqFeature object
                       		my %hash = gtfCols(@$feature);
                       		push @fileDump, \%hash;
                	}
		}
	}
        $self->{array} = \@fileDump;
}

sub load{
        my $self = shift;
        my (@command, $file);

        if(defined $self->{url}){
                $file = (split '/', $self->{url})[-1];
		if (!-f $file){
			print "Downloading $file...\n\n";
			system "wget -O- -q $self->{url} > $file";
		}
		unless ($self->{fixed}){
			$tmp .= ".gtf.gz";
			@command = ("Parsing data from $file...\n" , 
			"gunzip -c $file | awk -F \"\\t\" 'BEGIN{OFS=\"\\t\";} ((\$9~/NM/ || \$9~/NR/) && \$1==\"$self->{chrom}\" && ((\$5 >= $self->{chromStart} && \$5 <= $self->{chromEnd}) || (\$4 >= $self->{chromStart} && \$4 <= $self->{chromEnd})))' | gzip > $tmp");
		}
	} else {
		if (defined $self->{hub}){
			@command = ("Requesting track data from $self->{trackName} to UCSC Public Hubs REST API...\n" ,
               		"wget -O- -q 'http://api.genome.ucsc.edu/getData/track?hubUrl=$self->{hub};genome=$self->{genome};track=$self->{trackName};chrom=$self->{chrom};start=$self->{chromStart};end=$self->{chromEnd}' > $tmp");
        	} else {
			@command = ("Requesting track data from $self->{trackName} to UCSC Genome Browser REST API...\n" ,
        		"wget -O- -q 'http://api.genome.ucsc.edu/getData/track?genome=$self->{genome};track=$self->{trackName};chrom=$self->{chrom};start=$self->{chromStart};end=$self->{chromEnd}' > $tmp");
		}
	}
	
	if (@command){
		print $command[0];
 		system($command[1]);
	}
	return $tmp;
}

sub gtfCols{
        my @colVals = @_;
        my %arrayVals = ( SEQID => 0,
	        START => 1,
        	STOP  => 2,
        	STRND => 3,
        	NAME  => 4,
        	ID    => 5,
        	TYPE  => 6,
        	SRC   => 7,
        	SCORE => 8,
        	PHASE => 9,
        	ATTRB => 10,
        	SUBF => 11 );

        my %hash;
        foreach my $key (keys %arrayVals){
                $hash{$key} = @colVals[ $arrayVals{$key} ];
        }
        return %hash;
}

1;
