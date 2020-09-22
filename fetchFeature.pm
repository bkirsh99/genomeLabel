package fetchFeature;
use JSON;
use strict;
use warnings;
use File::Temp;
use Data::Dumper;

my $tmp = File::Temp->new();

BEGIN {
    use Exporter ();

    @fetchFeature::ISA         = qw(Exporter);
    @fetchFeature::EXPORT      = qw();
    @fetchFeature::EXPORT_OK   = qw($region $chr $start $end);

  }
  use vars qw($region $chr $start $end);

sub new {
    my ($class, $args) = @_;
    my $self = bless {
        name => $args->{name},
        url => $args->{url} || undef,
        hub => $args->{hub} || undef,
        trackName => $args->{trackName} || undef,
        chrom => $chr,
        chromStart => $start,
        chromEnd => $end,
        genome => $args->{genome} || 'hg38',
        fixed => $args->{fixed} || 0,
        }, $class;

        #$self->load unless $self->{fixed};
        $self->load;
        $self->make_AoH unless $self->{fixed};

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
                my $text = decode_json( $data );
                @fileDump = @{ $text->{ $self->{trackName} } };
        } else {
                unless ($self->{fixed}) {
                        my $parser = Bio::ToolBox::parser::gff->new(
                                file => $self->load,
                                #do_gene => 0,
                                #do_exon => 0,
                                do_cds => 1,
                                do_utr => 1,
                                do_codon => 1,
                        ) or die "Unable to open \$parser->{file}, please ensure '.gff' or '.gtf' extension.\n";

                        while (my $feature = $parser->next_feature() ) {
                        #each $feature is a SeqFeature object
                                my %hash = gtfCols(@$feature);
                                push @fileDump, \%hash;
                        }
                }
        }
        $self->{array} = \@fileDump;
}

sub load{
        my $self = shift;
        my $command;

        if(defined $self->{url}){
                my $file = (split '/', $self->{url})[-1];
                system "wget -O- -q $self->{url} > $file" if (!-f $file);

                unless ($self->{fixed}){
                        $tmp .= ".gtf.gz";
                        $command = "gunzip -c $file | awk -F \"\\t\" 'BEGIN{OFS=\"\\t\";} ((\$9~/NM/ || \$9~/NR/) && \$1==\"$self->{chrom}\" && ((\$5 >= $self->{chromStart} && \$5 <= $self->{chromEnd}) || (\$4 >= $self->{chromStart} && \$4 <= $self->{chromEnd})))' | gzip > $tmp";
                }

        } else {
                if(defined $self->{hub}){
                        $command = "wget -O- -q 'http://api.genome.ucsc.edu/getData/track?hubUrl=$self->{hub};genome=$self->{genome};track=$self->{trackName};chrom=$self->{chrom};start=$self->{chromStart};end=$self->{chromEnd}' > $tmp";
                }
                else {
                         $command = "wget -O- -q 'http://api.genome.ucsc.edu/getData/track?genome=$self->{genome};track=$self->{trackName};chrom=$self->{chrom};start=$self->{chromStart};end=$self->{chromEnd}' > $tmp";
                }
        }

        if (defined $command){
                system($command);
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
