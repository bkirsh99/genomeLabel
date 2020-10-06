#!/usr/bin/perl

use JSON;
use strict;
use warnings;
use Data::Dumper;
use 5.010;
use 5.005;
use Storable;
use lib '.';
use bedFeature;

use constant LABELS => {
        EXON => ['exon.bed'],
        INTRON => ['intron.bed'],
        INTERGENIC => ['intergenic.bed'],
        CODING => ['coding.bed'], 
        NONCODING => ['noncoding.bed'],
        REPEAT => ['repeat.bed'],
        PROMOTER => ['promoter.bed'],
        FUNC_ELEMENT => ['f_element.bed'],
        CIS_REG_MODULE => ['cr_module.bed'],
        TRANSCRIPT => ['transcript.bed'],
};

my $history = shift;
my $path = shift;
my $biotype = shift;
my $dir = shift;
my ($chr,$start,$end) = split /[:-]+/, substr($dir,0,index($dir,"@"));

##DEBUGGING:
print Dumper $chr, $start, $end, $biotype;
#print Dumper $regulator, $repeat, $felement;

#append <chr:start-stop@biotype>/labels/ to file names - must be done before declaring %file in order to update file names
my $var = addDir($dir, LABELS);

#declare arrays for each label (to be filled with bedFeature objects)
my (@exon, @transcript, @intron, @intergenic, @coding, @noncoding, @promoter, @repeat, @f_element, @cr_module);
my %file = ( $var->{EXON}[0] => \@exon, $var->{TRANSCRIPT}[0] => \@transcript,
                $var->{INTRON}[0] => \@intron, $var->{INTERGENIC}[0] => \@intergenic,
                $var->{CODING}[0] => \@coding, $var->{NONCODING}[0] => \@noncoding,
                $var->{PROMOTER}[0] => \@promoter, $var->{REPEAT}[0] => \@repeat,
                $var->{FUNC_ELEMENT}[0] => \@f_element, $var->{CIS_REG_MODULE}[0] => \@cr_module, );


my $retrieve = retrieve($history);
my $hash = @$retrieve[1]; #all previously downloaded bedFeature objects (HASH)
newLabels(%$hash);

#append <chr:start-stop@biotype>/labels/ to file names
sub addDir{
	my $dir = shift;
        my $ref = shift;
        my $new;

        while (my ($key, $value) = each ( %{ $ref } )){
                $new->{$key}[0] = join('/', $dir, $ref->{$key}[0]);
        }
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
                               	} elsif ($hash->{TYPE} eq 'transcript'){
                                       	push @transcript, $_;
                               	} elsif ($hash->{TYPE} =~ /CDS|codon/i){
                                       	push @coding, $_;
                               	} elsif ($hash->{TYPE} =~ /UTR/ || ( $hash->{TYPE} eq 'exon' && $hash->{ATTR}->{Parent} =~ /NR/) ){
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
	      to_Bed($key, $value);
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

#convert bed files to bedFeature objects (with original "name" field)
sub from_Bed{
        my $file = shift;
	my $flag = shift;
	my @array;
       
	open(FH, '<', $file) or die "Could not open file '$file'\n";
        while(<FH>){
		my $feature = bedFeature->new();
		foreach($feature->from_bed_string($_)){
			if ($flag == 0){
				push @array, $_->renameFeature($file,0);
			}
        	}
	}
        close(FH);
	return \@array;
}

#create additional ncbiRefSeq elements and apply liftOver to SegWay
sub manipulate_Bed{
        my $path = join('', $path, 'liftOver');
        my $command = "awk -F \"\\t\" 'BEGIN{OFS=\"\\t\";}; {if(\$2<$start) \$2=\$7=$start}; {if(\$3>$end) \$3=\$8=$end}; \$1==\"$chr\" && \$3>=$start && \$2<=$end'";

        system "bedtools sort -i $var->{EXON}[0] | bedtools merge -i stdin -c 4 -o distinct | awk 'BEGIN{OFS=\"\\t\";} {print \$1,\$2,\$3,\$4,0,\".\",\$2,\$3}' | uniq | gzip > $var->{EXON}[0].gz";

        system "bedtools sort -i $var->{TRANSCRIPT}[0] | bedtools merge -i stdin -c 4 -o distinct | bedtools subtract -a stdin -b $var->{EXON}[0] | awk 'BEGIN{OFS=\"\\t\";} {print \$1,\$2,\$3,\$4,0,\".\",\$2,\$3}' | gzip > $var->{INTRON}[0].gz";

        system "bedtools sort -i $var->{TRANSCRIPT}[0] | bedtools complement -i stdin -g hg38.chrom.sizes | awk 'BEGIN{OFS=\"\\t\";} {print \$1,\$2,\$3,\"intergenic\",0,\".\",\$2,\$3}' | gzip > $var->{INTERGENIC}[0].gz";

        system "bedtools sort -i $var->{CODING}[0] | bedtools merge -i stdin -c 4 -o distinct | awk 'BEGIN{OFS=\"\\t\";} {print \$1,\$2,\$3,\$4,0,\".\",\$2,\$3}' | gzip > $var->{CODING}[0].gz";

        system "bedtools sort -i $var->{NONCODING}[0] | bedtools subtract -a stdin -b $var->{CODING}[0] | sort -k1,1 -k2,2n | bedtools merge -i stdin -c 4 -o distinct | awk 'BEGIN{OFS=\"\\t\";} {print \$1,\$2,\$3,\$4,0,\".\",\$2,\$3}' | uniq | gzip > $var->{NONCODING}[0].gz";

        system "$path $var->{FUNC_ELEMENT}[0] hg19ToHg38.over.chain.gz tmp1 tmp2 && mv tmp1 $var->{FUNC_ELEMENT}[0] && rm tmp2";

        foreach (glob "$dir/*"){
                if ($_ =~ /.gz/){
                        system "gunzip -f $_";
                        $_ =~ s{.gz}{};
                }
                system "cat $_ | $command > tmp && mv tmp $_"; #sanity check
        }

	system "cat $var->{CODING}[0] $var->{NONCODING}[0] > $var->{EXON}[0] | rm $var->{CODING}[0] $var->{NONCODING}[0]"; 

}

