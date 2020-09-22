package bedFeature;
use strict;
use warnings;
use Data::Dumper;
use Object::ID;

use lib '.';
use fetchFeature;

sub new {
    my $class = shift;
    my $args = shift;
    my $options = shift;

    my $self = bless {}, $class;

        $self->{original} = $args;
        $self->encode;
        $self->{chrom} = $args->{chrom} || $args->{CHROM} || $args->{SEQID} || $args->{genoName} || undef;
        $self->{chromStart} = $args->{chromStart} || $args->{CHROMSTART} || $args->{START} || $args->{genoStart} || undef;
        $self->{chromEnd} = $args->{chromEnd} || $args->{CHROMEND} || $args->{STOP} || $args->{genoEnd} || undef;
        $self->strand($args->{strand} || $args->{STRAND} || $args->{STRND}) if exists $args->{strand} || $args->{STRAND} || $args->{STRND};
        $self->{score} = $args->{score} || $args->{SCORE} || $args->{swScore} || undef;
        $self->{itemRgb} = $args->{reserved} || undef;

        if (defined $args->{SRC} && $args->{SRC} eq 'ncbiRefSeq'){
                $self->{chromStart} = $self->{chromStart} -1;
        }

        if (defined $args->{thickStart}){
                $self->{thickStart} = $args->{thickStart};
        } else {
                $self->{thickStart} = $self->{chromStart};
        }
        if (defined $args->{thickEnd}){
                $self->{thickEnd} = $args->{thickEnd};
        } else {
                $self->{thickEnd} = $self->{chromEnd};
        }

        if(@_){
                my $source = shift;
                $self->filter($options, @_);
        } else {
                $self->filter($options);
        }

    return $self;
}

sub encode{
        my $self = shift;
        my $json = JSON->new->allow_nonref;
        my $encoded = $json->encode($self->{original});
        $self->{name} = $encoded;
}

sub decode{
        my $self = shift;
        my $json = JSON->new->allow_nonref;
        return $json->decode($self->{name});
}

sub filter{
        my $self = shift;
        my $options = shift;

        if(shift){
                my $source = shift;
                print Dumper $options, $source;
        }

=begin comment
        if (defined $self->{original}->{SRC} && $self->{original}->{SRC} eq 'ncbiRefSeq'){
                $self->parse_ncbiRefSeq;
        }
        elsif (defined $self->{original}->{Biotypes}){
                if(defined $options->{biotype} && $self->{original}->{Biotypes} eq $options->{biotype}){
                        if (defined $options->{regulator}){
                                foreach my $element (@{ $options->{regulator} }){
                                        if ($self->{original}->{TF} eq $element){
                                                #$self->{name} = $self->{original}->{name};
                                                $self->{itemRgb} = $self->{original}->{reserved};
                                                return $self;
                                        }
                                }
                        }
                }
        }

        elsif (defined $self->{original}->{sourceScores}){
                #$self->{name} = $self->{original}->{name};
        }
	
        elsif (defined $self->{original}->{repName}){
                if(defined $options->{repeat}){
                        foreach my $key (keys %{ $options->{repeat} }){
                                if(defined $options->{repeat}->{$key}){
                                        foreach my $element (@{ $options->{repeat}->{$key} }){
                                                if ($self->{original}->{repName} eq $element || $self->{original}->{repClass} eq $element || $self->{original}->{repFamily} eq $element){
                                                        #$self->{name} = $element;
                                                        return $self;
                                                }
                                        }
                                }
                        }
                }
        }
=end comment
=cut
}


sub strand {
        my $self = shift;
        if (@_) {
                my $str = $_[0];
                if ($str eq '+' || $str eq '1') {
                        $self->{strand} = '+';
                }
                elsif ($str eq '-' || $str eq '-1') {
                        $self->{strand} = '-';
                } else {
                        $self->{strand} = '.';
                }
        }
        return defined $self->{strand} ? $self->{strand} : '.';
}

sub bed_string {
        my $self = shift;
        my $string = join("\t", (
                $self->{chrom},
                $self->{chromStart},
                $self->{chromEnd},
                defined $self->{name} ? $self->{name} : 'fill',
                defined $self->{score} ? $self->{score} : 0,
                defined $self->{strand} ? $self->{strand} : '.',
                $self->{thickStart},
                $self->{thickEnd},
                defined $self->{itemRgb} ? $self->{itemRgb} : 'undef',
        ) );
        return "$string\n";
}

sub from_bed_string {
        my $self = shift;
        my $string = shift;
        my $json = JSON->new->allow_nonref;
        # check string
        chomp $string;

        my @fields = split('\t', $string);

        $self->{original} = $json->decode($fields[3]);
        $self->{chrom} = $fields[0];
        $self->{chromStart} = $fields[1];
        $self->{chromEnd} = $fields[2];
        $self->{score} = $fields[4];
        $self->{strand} = $fields[5];
        $self->{thickStart} = $fields[6];
        $self->{thickEnd} = $fields[7];
        $self->{itemRgb} = $fields[8];

        return $self;
}

sub get_original_entry{
        my $self = shift;
        return Dumper ? $self->{original} : $self->{original};
}

sub print_original_entry{
        my $self = shift;
        print Dumper $self->{original};
}


1;
