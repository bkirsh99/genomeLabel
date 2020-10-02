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
	unless ($self->{name} eq 'NA'){
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

	if ($file =~ /exon|intron|intergenic|transcript|coding/){
        	$self->{name} = $self->{original} eq 'BEDTools Complement' ? 'Intergenic' : $self->{original}->{ATTRB}->{gene_name};
        } elsif ($file =~ /repeat/){
                $self->{name} = join(".", $self->{original}->{repClass}, $self->{original}->{repFamily}, $self->{original}->{repName});
        } elsif ($file =~ /promoter/){
                $self->{name} = $self->{original}->{name};
	} elsif ($file =~ /cr_module/){
		$self->{name} = $self->{original}->{TF};
        } elsif ($file =~ /f_element/){
                $self->{name} = $self->{original}->{name};
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
                } #else {
		#	$self->{strand} = '.';
		#}
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
		unless ($_ eq 'NA'){
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
