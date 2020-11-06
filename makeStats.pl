#!/usr/bin/perl

use JSON;
use strict;
use warnings;
use Data::Dumper;
use File::Path;
use 5.010;
use 5.005;
use Storable;
use lib '.';
use bedFeature;

use constant LABELS => {
        EXON => ['exon.bed'],
        INTRON => ['intron.bed'],
        INTERGENIC => ['intergenic.bed'],
        CODING => ['coding-exon.bed'], 
        NONCODING => ['noncoding-exon.bed'],
        REPEAT => ['repeat.bed'],
        PROMOTER => ['promoter.bed'],
        FUNC_ELEMENT => ['f-element.bed'],
        CIS_REG_MODULE => ['cr-module.bed'],
        TRANSCRIPT => ['transcript.bed'],
};

my $history = shift;
my $path = shift;
my $biotype = shift;
my $dir = shift;
my $biotypeR = shift;
my ($chr,$start,$end) = split /[:-]+/, $dir;

my $tmpdir = join("/", $dir, "tmp");
if (!-d $tmpdir){
        mkdir $tmpdir or die "Failed to create output directory: $!\n";
}
my $catfile = join("/", $dir, join("_", $biotype, "cat.bed"));

##DEBUGGING:
#print Dumper $chr, $start, $end, $biotype;

#append <chr:start-stop@biotype>/labels/ to file names - must be done before declaring %file in order to update file names
my $var = addDir($tmpdir, LABELS);

#declare arrays for each label (to be filled with bedFeature objects)
my (@exon, @transcript, @intron, @intergenic, @coding, @noncoding, @promoter, @repeat, @f_element, @cr_module);
my %file = ( $var->{EXON}[0] => \@exon, $var->{TRANSCRIPT}[0] => \@transcript,
                $var->{INTRON}[0] => \@intron, $var->{INTERGENIC}[0] => \@intergenic,
                $var->{CODING}[0] => \@coding, $var->{NONCODING}[0] => \@noncoding,
                $var->{PROMOTER}[0] => \@promoter, $var->{REPEAT}[0] => \@repeat,
                $var->{FUNC_ELEMENT}[0] => \@f_element, $var->{CIS_REG_MODULE}[0] => \@cr_module, );

#print Dumper \%file;
my $retrieve = retrieve($history);
my $hash = @$retrieve[1]; #all previously downloaded bedFeature objects (HASH)
newLabels(%$hash);
deleteEmpty();

#append <chr:start-stop@biotype>/labels/ to file names
sub addDir{
	my $dir = shift;
        my $ref = shift;
        my $new;

        while (my ($key, $value) = each ( %{ $ref } )){
                unless ($key =~ /FUNC_ELEMENT|CIS_REG_MODULE/){
			$new->{$key}[0] = join('/', $dir, join('_',"raw",$ref->{$key}[0]));
		}
	}
	$new->{'FUNC_ELEMENT'}[0] = join('/', $dir, join('_',"filtered",$biotype,$ref->{'FUNC_ELEMENT'}[0]));
	$new->{'CIS_REG_MODULE'}[0] = join('/', $dir, join('_',"filtered",$biotype,$ref->{'CIS_REG_MODULE'}[0]));
        return $new;
}

#sort bedFeature objects for bed file conversion
sub newLabels{
        my %hash = @_; #hash of bedFeature arrays, where key = database

        foreach my $key (keys %hash){
                my @array = @{ $hash{$key} };
                foreach (@array){
                        my $hash = $_->decodeName; #$_->{original}
			if ($key eq 'ncbiRefSeq'){
                               	if ($hash->{TYPE} eq 'exon'){
                                       	push @exon, $_;
					push @noncoding, $_ unless $hash->{ATTRB}->{Parent} =~ /NM/;
                               	} elsif ($hash->{TYPE} eq 'transcript'){
                                       	push @transcript, $_;
                               	} elsif ($hash->{TYPE} =~ /CDS|codon/i){
                                       	push @coding, $_;
                               	} elsif ($hash->{TYPE} =~ /UTR/){
                                       	push @noncoding, $_;
                               	}
                        } elsif ($key eq 'FANTOM5'){
                                push @promoter, $_;
                        } elsif ($key eq 'RepeatMasker'){
                                push @repeat, $_;
                        } elsif ($key eq 'SegWay'){
                                push @f_element, $_;
                        } elsif ($key eq 'ReMap'){
				push @cr_module, $_;
                        }
                }
	}
	while (my ($key, $value) = each (%file)){
		to_Bed($key, $value) unless (-f $key);
        }
        manipulate_Bed(); #make corrections
}

#convert bedFeature objects to bed files (with renamed "name" field)
sub to_Bed{
        my $file = shift;
        my $array = shift;

        open(FH, '>', $file) or die "Could not open file '$file'\n";
        foreach (@{ $array }){
		my $renamed = $_->renameFeature($file,0);
		print FH $renamed->bed_string;
        }
        close(FH);
}

#create additional ncbiRefSeq elements and apply liftOver to SegWay
sub manipulate_Bed{
        my $path = join('', $path, 'liftOver');

	my $command = "awk -F \"\\t\" 'BEGIN{OFS=\"\\t\";}; {if(\$2<$start) \$2=$start}; {if(\$3<$start) \$3=$start}; {if(\$7<$start) \$7=$start}; {if(\$8<$start) \$8=$start}; {if(\$2>$end) \$2=$end}; {if(\$3>$end) \$3=$end}; {if(\$7>$end) \$7=$end}; {if(\$8>$end) \$8=$end}; \$1==\"$chr\" && \$3>$start && \$2<$end'";

        system "bedtools sort -i $var->{EXON}[0] | awk 'BEGIN{OFS=\"\\t\";} {print \$1,\$2,\$3,\$4,0,\".\",\$2,\$3}' | uniq | gzip > $var->{EXON}[0].gz";

        system "bedtools sort -i $var->{TRANSCRIPT}[0] | bedtools subtract -a stdin -b $var->{EXON}[0] | awk 'BEGIN{OFS=\"\\t\";} {print \$1,\$2,\$3,\$4,0,\".\",\$2,\$3}' | gzip > $var->{INTRON}[0].gz";

        system "bedtools sort -i $var->{TRANSCRIPT}[0] | bedtools complement -i stdin -g hg38.chrom.sizes | awk 'BEGIN{OFS=\"\\t\";} {print \$1,\$2,\$3,\"intergenic\",0,\".\",\$2,\$3}' | gzip > $var->{INTERGENIC}[0].gz";

        system "bedtools sort -i $var->{CODING}[0] | awk 'BEGIN{OFS=\"\\t\";} {print \$1,\$2,\$3,\$4,0,\".\",\$2,\$3}' | gzip > $var->{CODING}[0].gz";

        system "bedtools sort -i $var->{NONCODING}[0] | bedtools subtract -a stdin -b $var->{CODING}[0] | sort -k1,1 -k2,2n | awk 'BEGIN{OFS=\"\\t\";} {print \$1,\$2,\$3,\$4,0,\".\",\$2,\$3}' | uniq | gzip > $var->{NONCODING}[0].gz";

        system "$path $var->{FUNC_ELEMENT}[0] hg19ToHg38.over.chain.gz tmp1 tmp2 && mv tmp1 $var->{FUNC_ELEMENT}[0] && rm tmp2";

	system "awk 'BEGIN{OFS=\"\\t\";} {gsub(/intron/,\"transcript\",\$4)}; {print}' $var->{TRANSCRIPT}[0] | gzip > $var->{TRANSCRIPT}[0].gz";
	
	system "sed -i /hg/d $var->{PROMOTER}[0]";

        foreach (glob "$tmpdir/*"){
                unless ($_ =~ /.gz/ || $_ =~ /$var->{PROMOTER}[0]|$var->{REPEAT}[0]|$var->{FUNC_ELEMENT}[0]|$var->{CIS_REG_MODULE}[0]/){
			unlink $_ or die "Can't delete $_: $!\n";
			next;
		}
		if ($_ =~ /.gz/){ 
			system "gunzip -f $_";
			$_ =~ s{.gz}{};
		}
		system "awk -F \"\\t\" 'BEGIN{OFS=\"\\t\";} {print \$1,\$2,\$3,\$4 >> \"$catfile\"}' $_" unless $_ =~ /$var->{EXON}[0]|$var->{INTERGENIC}[0]/;
	}

	rmtree $tmpdir;
	system "bedtools sort -i $catfile > tmp && mv tmp $catfile";
}

sub deleteEmpty{
        foreach (glob "$dir/*"){
                if (`cat $_ | wc -c` == 0){
                        print "$_ is empty... Deleting file.\n";
                        unlink $_ or die "Can't delete $_: $!\n";
                }
	}
}
