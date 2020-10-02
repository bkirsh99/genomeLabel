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
        EXON => ['exon.bed', '0,255,0'], #green
        INTRON => ['intron.bed', '255,0,0'], #red
        INTERGENIC => ['intergenic.bed', '0,0,255'], #blue
        CODING => ['coding.bed', '0,102,0'], #dark-green
        NONCODING => ['noncoding.bed', '153,255,153'], #light-green
        REPEAT => ['repeat.bed', '255,0,127'], #pink
        PROMOTER => ['promoter.bed','127,0,255'], #purple
        FUNC_ELEMENT => ['f_element.bed'], #default (Segway)
        CIS_REG_MODULE => ['cr_module.bed'], #default (ReMap)
        TRANSCRIPT => ['transcript.bed'], #not kept, only used to extract intronic/intergenic regions
};

my $history = shift;
my $path = shift;
my $region = shift;
my $chr = shift;
my $start = shift;
my $end = shift;
my $biotype = shift;
#my @filter = shift; #in order: $biotypeR, $regulator, $repeat, $felement
my @filter = ('Huvec-C', ['CTCF','GATA1','GATA2'], ['LINE'], undef);
	my $biotypeR = $filter[0];
	my $regulator = $filter[1] if defined $filter[1];
	my $repeat = $filter[2] if defined $filter[2];
	my $felement = $filter[3] if defined $filter[3];

my @strarray;
foreach (@filter){
	if (ref($_) eq 'ARRAY'){
		my $arrstr = join '_', @$_ if defined $_;
		push @strarray, $arrstr;
	}
}
unshift @strarray, $biotype;
my $filterstr = join '.', @strarray;

##DEBUGGING:
#print Dumper $region, $chr, $start, $end, $biotype;
#print Dumper $regulator, $repeat, $felement;

#establish directory names to store files and create, if needed
my $out = join('@',$region,$biotype);
if (!-d $out){
	mkdir $out or die "Failed to create output directory: $!\n";
}
my $trackdir = join('/',$out,'tracks');

if (!-d $trackdir){
	mkdir $trackdir or die "Failed to create track file directory: $!\n";
}

my $filterdir = join('/',$trackdir,$filterstr); #this is where tracks will be stored (filter-specifiv)
if (!-d $filterdir){
	mkdir $filterdir or die "Failed to create filter-specific track file directory: $!\n";
}

#append <chr:start-stop@biotype>/tracks/<filter>/ to file name - must be done before declaring %file in order to update file names
my $var = addDir(LABELS);

#declare arrays for each label (to be filled with bedFeature objects)
my (@exon, @transcript, @intron, @intergenic, @coding, @noncoding, @promoter, @repeat, @f_element, @cr_module);
my %file = ( $var->{EXON}[0] => \@exon, $var->{TRANSCRIPT}[0] => \@transcript,
                $var->{INTRON}[0] => \@intron, $var->{INTERGENIC}[0] => \@intergenic,
                $var->{CODING}[0] => \@coding, $var->{NONCODING}[0] => \@noncoding,
                $var->{PROMOTER}[0] => \@promoter, $var->{REPEAT}[0] => \@repeat,
                $var->{FUNC_ELEMENT}[0] => \@f_element, $var->{CIS_REG_MODULE}[0] => \@cr_module, );

my $labeldir = join('/',$out,'labels');
if (!-d $labeldir){
	my $retrieve = retrieve($history);
	my $hash = @$retrieve[1]; #all previously downloaded bedFeature objects (HASH)
	newTracks(%$hash);
} else {
	my $hash = parseLabelsFile();
	oldTracks($hash); #all previously downloaded bedFeature objects (converted from BED FILES where name = orginal entry) 
}
genomeTrack(); #merge coding/noncoding exon, intron, and intergenic and delete transcript


#append <chr:start-stop@biotype>/tracks/ to file names
sub addDir{
        my $ref = shift;
        my $new;

        while (my ($key, $value) = each ( %{ $ref } )){
                $new->{$key}[0] = join('/', $filterdir, $ref->{$key}[0]);
        }
        return $new;
}

#filter and sort bedFeature objects for bed file conversion
sub newTracks{
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
                                $_->{itemRgb} = LABELS->{PROMOTER}[1]; #add colour
                                push @promoter, $_;
                        } elsif ($key eq 'RepeatMasker'){
				foreach my $repfilter (@$repeat){
					if ($repfilter =~ /$hash->{repClass}|$hash->{repFamily}|$hash->{repName}/i){
                                		$_->{itemRgb} = LABELS->{REPEAT}[1]; #add colour
                                		push @repeat, $_;
					}
				}
                        } elsif ($key eq 'SegWay'){
				foreach my $fefilter (@$felement){
                                        if ($fefilter =~ /$hash->{name}/i){
                                		push @f_element, $_;
					}
				}
                        } elsif ($key eq 'ReMap'){
				foreach my $regfilter (@$regulator){
                                       	if ($regfilter =~ /$hash->{TF}/i && $hash->{Biotypes} =~ /$biotypeR/i){
						push @cr_module, $_;
					}
				}
                        }
                }
	}
	while (my ($key, $value) = each (%file)){
	      to_Bed($key, $value);
        }
        manipulate_Bed(); #make corrections
}

sub oldTracks{
	my $hash = shift;
	
	foreach my $file (keys %$hash){
		print Dumper $file;
		to_Bed($file, $$hash{$file});
	}
	manipulate_Bed();
}

sub parseLabelsFile{
	my @files = glob "$labeldir/*.bed";
        my %hash;
	
	for(my $i = 0; $i <= $#files; $i++){
		my $file = join('/', $filterdir, $files[$i]=~ m/([^\/]+)$/);
		my $array = from_Bed($files[$i],0);
		$hash{$file} = $array;
	}

        return \%hash;
}

#not being used
sub parseLabelsSource{
	my @files = glob "$labeldir/*.bed";
	my (%hash, $key, @ncbi);
	
	for(my $i = 0; $i < $#files; $i++){
		my $array = from_Bed($files[$i],0);
		if ($files[$i] =~ /exon|intron|intergenic|transcript|coding/){
			$key = "ncbiRefSeq";
			foreach (@$array){
				push @ncbi, $_;
			}
		} elsif ($files[$i] =~ /repeat/){
			$key = "RepeatMasker";
		} elsif ($files[$i] =~ /promoter/){
			$key = "FANTOM5";
		} elsif ($files[$i] =~ /cr_module/){
			$key = "ReMap";
		} elsif ($files[$i] =~ /f_element/){
			$key = "SegWay";
		}
		$hash{$key} = $array unless $key eq "ncbiRefSeq";
	}
	
	$hash{"ncbiRefSeq"} = \@ncbi;
	return \%hash;
}

#convert bedFeature objects to bed files (with renamed "name" field)
sub to_Bed{
        my $file = shift;
        my $array = shift;

        open(FH, '>', $file) or die "Could not open file '$file'\n";
        foreach (@{ $array }){
		my $renamed = $_->renameFeature($file);
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
				push @array, $_->renameFeature($file);
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

        system "bedtools sort -i $var->{EXON}[0] | bedtools merge -i stdin -c 4 -o distinct | awk 'BEGIN{OFS=\"\\t\";} {print \$1,\$2,\$3,\$4,0,\".\",\$2,\$3,\"0,255,0\"}' | uniq | gzip > $var->{EXON}[0].gz";

        system "bedtools sort -i $var->{TRANSCRIPT}[0] | bedtools merge -i stdin -c 4 -o distinct | bedtools subtract -a stdin -b $var->{EXON}[0] | awk 'BEGIN{OFS=\"\\t\";} {print \$1,\$2,\$3,\$4,0,\".\",\$2,\$3,\"255,0,0\"}' | gzip > $var->{INTRON}[0].gz";

        system "bedtools sort -i $var->{TRANSCRIPT}[0] | bedtools complement -i stdin -g hg38.chrom.sizes | awk 'BEGIN{OFS=\"\\t\";} {print \$1,\$2,\$3,\"NA\",0,\".\",\$2,\$3,\"0,0,255\"}' | gzip > $var->{INTERGENIC}[0].gz";

        system "bedtools sort -i $var->{CODING}[0] | bedtools merge -i stdin -c 4 -o distinct | awk 'BEGIN{OFS=\"\\t\";} {print \$1,\$2,\$3,\$4,0,\".\",\$2,\$3,\"0,102,0\"}' | gzip > $var->{CODING}[0].gz";

        system "bedtools sort -i $var->{NONCODING}[0] | bedtools subtract -a stdin -b $var->{CODING}[0] | sort -k1,1 -k2,2n | bedtools merge -i stdin -c 4 -o distinct | awk 'BEGIN{OFS=\"\\t\";} {print \$1,\$2,\$3,\$4,0,\".\",\$2,\$3,\"153,255,153\"}' | uniq | gzip > $var->{NONCODING}[0].gz";

        system "$path $var->{FUNC_ELEMENT}[0] hg19ToHg38.over.chain.gz tmp1 tmp2 && mv tmp1 $var->{FUNC_ELEMENT}[0] && rm tmp2";

        foreach (glob "$filterdir/*"){
                if ($_ =~ /.gz/){
                        system "gunzip -f $_";
                        $_ =~ s{.gz}{};
                }
                system "cat $_ | $command > tmp && mv tmp $_"; #sanity check
        }

}

sub genomeTrack{
	my $file = join("/",$filterdir,"genome.bed");
	
	if (!-f $file){
		system "cat $var->{INTRON}[0] $var->{INTERGENIC}[0] $var->{CODING}[0] $var->{NONCODING}[0] > $file && rm $var->{INTRON}[0] $var->{INTERGENIC}[0] $var->{CODING}[0] $var->{NONCODING}[0] $var->{TRANSCRIPT}[0] $var->{EXON}[0]";
	}
}

