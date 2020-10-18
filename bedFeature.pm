package bedFeature;
use strict;
use warnings;
use Data::Dumper;
use Object::ID;
use MIME::Base64;
use Encode;

use lib '.';
use fetchFeature;

use constant LABELS => {
        EXON => ['exon.bed'],
        INTRON => ['intron.bed'],
        INTERGENIC => ['intergenic.bed'],
        CODING => ['coding.bed'],
        NONCODING => ['noncoding.bed'],
        REPEAT => ['repeat.bed'],
        PROMOTER => ['promoter.bed'],
        FUNC_ELEMENT => ['f-element.bed'],
        CIS_REG_MODULE => ['cr-module.bed'],
        TRANSCRIPT => ['transcript.bed'],
};

sub new {
    my $class = shift;
    my $args = shift;
    my $options = shift;

    my $self = bless {}, $class;

        $self->{original} = $args;
	$self->encodeName; #add name field
        $self->{chrom} = $args->{chrom} || $args->{CHROM} || $args->{SEQID} || $args->{genoName} || undef;
        $self->{chromStart} = $args->{chromStart} || $args->{CHROMSTART} || $args->{START} || $args->{genoStart} || undef;
        $self->{chromEnd} = $args->{chromEnd} || $args->{CHROMEND} || $args->{STOP} || $args->{genoEnd} || undef;

        if (defined $args->{SRC} && $args->{SRC} eq 'ncbiRefSeq'){
        	$self->parse_ncbiRefSeq;
        } else {
		$self->strand($args->{strand} || $args->{STRAND} || $args->{STRND}) if exists $args->{strand} || $args->{STRAND} || $args->{STRND};
        	$self->{score} = $args->{score} || $args->{SCORE} || $args->{swScore} || undef;
        	$self->{itemRgb} = $args->{reserved} || undef;

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
	}	

    return $self;
}

sub encodeName{
	my $self = shift;
	my $json_string = JSON->new->allow_nonref->encode($self->{original});
	$self->{name} = $json_string;
}

sub decodeName{
	my $self = shift;
	unless ($self->{name} eq 'intergenic'){
		return JSON->new->allow_nonref->decode($self->{name});
	} else {
		return \{};
	}
}

sub parse_ncbiRefSeq{
	my $self = shift;
	--$self->{chromStart};
	--$self->{chromEnd};
}

sub renameFeature{
	my $self = shift;
	my $file = shift;
	my $flag = shift;
	
	$file = substr($file,rindex($file,'_')+1);

#	print Dumper $file;
#print Dumper LABELS->{CIS_REG_MODULE}[0];
	if ($file eq LABELS->{EXON}[0] || $file eq LABELS->{INTRON}[0] || $file eq LABELS->{INTERGENIC}[0] || $file eq LABELS->{TRANSCRIPT}[0] || $file eq LABELS->{CODING}[0] || $file eq LABELS->{NONCODING}[0]){
		$self->{name} = $self->{original} eq 'BEDTools Complement' ? 'intergenic' : $self->{original}->{ATTRB}->{gene_name};
	} elsif ($file eq LABELS->{REPEAT}[0]){
    		$self->{name} = join(':', $self->{original}->{repClass}, $self->{original}->{repFamily}, $self->{original}->{repName});
        } elsif ($file eq LABELS->{PROMOTER}[0]){
                $self->{name} = $self->{original}->{name};
	} elsif ($file eq LABELS->{CIS_REG_MODULE}[0]){
		$self->{name} = $self->{original}->{TF};
        } elsif ($file eq LABELS->{FUNC_ELEMENT}[0]){
                $self->{name} = $self->{original}->{name};
        }

	if ($flag == 0){
		unless ($file eq LABELS->{TRANSCRIPT}[0]){
			$file = substr($file,0,length($file)-4);
		} else {
			$file = 'intron';
		}
		$self->{name} = join('.', $file, $self->{name});
	}

	return $self;
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
                $self->{name},
                defined $self->{score} ? $self->{score} : 0,
                defined $self->{strand} ? $self->{strand} : '.',
                defined $self->{thickStart} ? $self->{thickStart} : '-',
                defined $self->{thickEnd} ? $self->{thickEnd} : '-',
                defined $self->{itemRgb} ? $self->{itemRgb} : '-',
        ) );
        return "$string\n";
}

#not being used - maybe consider if using another encoding
sub alphanum{
	my $self = shift;
	my $str = $self->{name};
	$str =~ tr/\n/\@/;
	return $str;
}

sub from_bed_string {
        my $self = shift;
	my $string = shift;
	my @feats;
        # check string
	chomp $string;

        my @fields = split('\t', $string);
	my @id = ( $fields[3] =~ /(^{.*?}.*?)}/ ) ? ( $fields[3] =~ /(^{.*?}.*?)}/ ) : $fields[3];
        foreach (@id){
		unless ($_ eq 'intergenic'){
			$_ .= "\}" unless substr($_, -1) eq "\}";
                	my $hash;
                	eval { $hash = JSON->new->allow_nonref->decode( $_ ) };
                	if ($@){
                		$_ .= "\}";
                	}
			$self->{original} = JSON->new->allow_nonref->decode( $_ );
		} else {
			$self->{original} = "BEDTools Complement";
		}
                $self->{name} = $_;
		$self->{chrom} = $fields[0];
        	$self->{chromStart} = $fields[1];
        	$self->{chromEnd} = $fields[2];
		$self->{score} = $fields[4];
        	$self->{strand} = $fields[5];
		$self->{thickStart} = $fields[6];
		$self->{thickEnd} = $fields[7];
		$self->{itemRgb} = $fields[8];

       		push @feats, $self;
	}
	return @feats;
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
