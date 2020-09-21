package bedFeature;
use strict;
use warnings;
use Data::Dumper;
use Object::ID;
#use Bio::Tools::Run::BEDTools;

use lib '.';
use fetchFeature;

use constant LABELS => {
	EXON => ['exon.bed', 'Exonic Region', '0,255,0', 'Green'],
	INTRON => ['intron.bed', 'Intronic Region', '255,0,0', 'Red'],
	INTERGENIC => ['intergenic.bed', 'Intergenic Region', '0,0,255', 'Blue'],
	CODING => ['coding.bed', 'Coding Exonic Region', '0,102,0', 'Dark Green'],
	NONCODING => ['noncoding.bed', 'Non-Coding Exonic Region', '153,255,153', 'Light Green'],
	LINE => ['LINE Element', '255,255,0', ''],
	SINE => ['SINE Element', '255,204,153', ''],
	ALU => ['Alu Element', '204,102,0', ''],
	PROMOTER => ['Promoter Element', '127,0,255', ''],
	FUNC_ELEMENT => ['Functional Element', '102,178,255', ''],
	CIS_REG_MODULE => ['Cis-Regulatory Module', '255,0,127', ''],
};  

sub new {
    my $class = shift;
    my $args = shift;
    my $options = shift;

    my $self = bless {}, $class;

        $self->{original} = $args;
	$self->{name} = "$args";
        $self->{chrom} = $args->{chrom} || $args->{CHROM} || $args->{SEQID} || $args->{genoName} || undef;
        $self->{chromStart} = $args->{chromStart} || $args->{CHROMSTART} || $args->{START} || $args->{genoStart} || undef;
        $self->{chromEnd} = $args->{chromEnd} || $args->{CHROMEND} || $args->{STOP} || $args->{genoEnd} || undef;
        $self->strand($args->{strand} || $args->{STRAND} || $args->{STRND}) if exists $args->{strand} || $args->{STRAND} || $args->{STRND};
        $self->{score} = $args->{score} || $args->{SCORE} || $args->{swScore} || undef;

        #$args->{name} ||= $args->{NAME} ||= $args->{TYPE} ||= $args->{repName} ||= $args->{repClass} ||= $args->{repFamily} ||= $args->{Biotypes} ||= $args->{TF};

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
		#print Dumper @_, $source;
		$self->filter($options, @_);
	} else {
		$self->filter($options);
	}

	#print Dumper $self;
    return $self;
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
#=begin comment
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
#=end comment
#=cut
#=begin comment
	elsif (defined $self->{original}->{sourceScores}){
		#$self->{name} = $self->{original}->{name};
	}
#=end comment 
#=cut 
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

sub parse_ncbiRefSeq{
	my $self = shift;

#	if ($self->{original}->{TYPE} eq 'exon'){
                        #$self->{name} = LABELS->{EXON}[0];
                        #make exon bed
#                }
#               if ($self->{original}->{TYPE} eq 'transcript'){
#make transcript bed
#call bedtools here to sort exon, transcript, and make intron b$#
		#}
#               if ($self->{original}->{TYPE} eq 'CDS' || $self->{original}->{TYPE} eq 'start_codon' || $self->{original}->{TYPE} eq 'stop_codon'){
#make coding bed
		#}
#               if ($self->{original}->{TYPE} eq 'UTR' || ( $self->{original}->{TYPE} eq 'exon' && $self->{original}->{ATTR}->{Parent} eq '/NR/') ){
#make noncoding bed
		#}

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
                defined $self->{itemRgb} ? $self->{itemRgb} : 'not',
        ) );
        return "$string\n";
}

sub get_original_entry{
        my $self = shift;
        return Dumper ? $self->{original} : $self->{original};
}

sub print_original_entry{
        my $self = shift;
        print Dumper $self->{original};
}

sub gtf_Correct{
        my $self = shift;

        if ($self->{original}->{SRC} eq 'ncbiRefSeq'){
                $self->{chromStart} = $self->{chromStart} -1;
        }
}


1;
