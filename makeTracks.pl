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
        FUNC_ELEMENT => ['f-element.bed'], #default (Segway)
        CIS_REG_MODULE => ['cr-module.bed'], #default (ReMap)
        TRANSCRIPT => ['transcript.bed'], #not kept, only used to extract intronic/intergenic regions
};

my $history = shift;
my $path = shift;
my $biotype = shift;
my $dir = shift;
my $biotypeR = shift;
my $repfilter =shift;
my $felefilter = shift;
my $crmodfilter = shift;
my ($chr,$start,$end) = split /[:-]+/, $dir;

my @repfilter = $repfilter eq "raw" ? [] : split(/_/,$repfilter);
my @regfilter = $crmodfilter eq "raw" ? [] : split(/_/,$crmodfilter);
my @fefilter = $felefilter eq "raw" ? [] : split(/_/, $felefilter);
@repfilter = @repfilter[2..$#repfilter];
@regfilter = @regfilter[2..$#regfilter];
@fefilter = @fefilter[2..$#fefilter];

my $regulator = \@regfilter;
my $repeat = \@repfilter;
my $felement = \@fefilter;
#print Dumper $regulator, $repeat, $felement;
#append <chr:start-stop@biotype>/tracks/<filter>/ to file names - must be done before declaring %file in order to update file names
my $var = addDir($dir, LABELS);

#declare arrays for each label (to be filled with bedFeature objects)
my (@exon, @transcript, @intron, @intergenic, @coding, @noncoding, @promoter, @repeat, @f_element, @cr_module);
my %file = ( $var->{EXON}[0] => \@exon, $var->{TRANSCRIPT}[0] => \@transcript,
                $var->{INTRON}[0] => \@intron, $var->{INTERGENIC}[0] => \@intergenic,
                $var->{CODING}[0] => \@coding, $var->{NONCODING}[0] => \@noncoding,
                $var->{PROMOTER}[0] => \@promoter, $var->{REPEAT}[0] => \@repeat,
                $var->{FUNC_ELEMENT}[0] => \@f_element, $var->{CIS_REG_MODULE}[0] => \@cr_module, );

my $genomefile = join("/",$dir,"raw_genome.bed");
#if ( (!-d $labeldir || emptyDir($labeldir) ) && emptyDir($dir) ){ #no labels or filtered files yet, create from scratch
	my $retrieve = retrieve($history);
	my $hash = @$retrieve[1]; #all previously downloaded bedFeature objects (HASH)
	newTracks(%$hash);
#} elsif ( -d $labeldir && !emptyDir($labeldir) ) { #labels but no filtered files, create from labels
#	my $hash = parseLabelsSource();
#	oldTracks($hash); 
#}
#genomeTrack($genomefile); #merge coding/noncoding exon, intron, and intergenic and delete transcript
deleteEmpty(); #remove empty files (possibly no filter match

#append filters to file names
sub addDir{
        my $dir = shift;
        my $ref = shift;
        my $new;

        while (my ($key, $value) = each ( %{ $ref } )){
                unless ($key =~ /REPEAT|FUNC_ELEMENT|CIS_REG_MODULE/){
                        $new->{$key}[0] = join('/', $dir, join('_',"raw",$ref->{$key}[0]));
                }
        }
        $new->{'REPEAT'}[0] = join('/', $dir, join('_',$repfilter,$ref->{'REPEAT'}[0]));
        $new->{'FUNC_ELEMENT'}[0] = join('/', $dir, join('_',$felefilter,$ref->{'FUNC_ELEMENT'}[0]));
        $new->{'CIS_REG_MODULE'}[0] = join('/', $dir, join('_',$crmodfilter,$ref->{'CIS_REG_MODULE'}[0]));
        return $new;
}

sub emptyDir{
	my $dirname = shift;
    	opendir(my $dh, $dirname) or die "Not a directory";
    	return scalar(grep { $_ ne "." && $_ ne ".." } readdir($dh)) == 0;
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
				unless (@$repeat){
					$_->{itemRgb} = LABELS->{REPEAT}[1]; #add colour
                                        push @repeat, $_;
				} else {
					foreach my $repfilter (@$repeat){
						if ($repfilter =~ /$hash->{repClass}|$hash->{repFamily}|$hash->{repName}/i){
                                			$_->{itemRgb} = LABELS->{REPEAT}[1]; #add colour
                                			push @repeat, $_;
						}
					}
				}
                        } elsif ($key eq 'SegWay'){
				unless (@$felement){
                                        push @f_element, $_;
              			} else {
					foreach my $fefilter (@$felement){
                                        	if ($fefilter =~ /$hash->{name}/i){
                                			push @f_element, $_;
						}
					}
				}
                        } elsif ($key eq 'ReMap'){
				if($hash->{Biotypes} =~ /$biotypeR/i){
					unless (@$regulator){
						push @cr_module, $_;
					} else {
						foreach my $regfilter (@$regulator){
	                                       		if ($regfilter =~ /$hash->{TF}/i){
								push @cr_module, $_;
							}
						}
					}
				}
                        }
                }
	}
	while (my ($key, $value) = each (%file)){
	      to_Bed($key, $value) unless (-f $key);
        }
        manipulate_Bed(); #make corrections
}
=begin comment
sub oldTracks{
	my $hash = shift;
	
}

sub parseLabelsFile{
	my @files = glob "$labeldir/*.bed";
        my %hash;
	
	for(my $i = 0; $i <= $#files; $i++){
		my $file = join('/', $dir, $files[$i]=~ m/([^\/]+)$/);
		my $array = from_Bed($files[$i],0);
		$hash{$file} = $array;
	}

        return \%hash;
}

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
=end comment
=cut

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
		#my $feature = bedFeature->new();
		#foreach($feature->from_bed_string($_)){
		#	if ($flag == 0){
				push @array, $_;#->renameFeature($file,0);
		#	}
        	#}
	}
        close(FH);
	return \@array;
}

#create additional ncbiRefSeq elements and apply liftOver to SegWay
sub manipulate_Bed{
        my $path = join('', $path, 'liftOver');
        my $command = "awk -F \"\\t\" 'BEGIN{OFS=\"\\t\";}; {if(\$2<$start) \$2=$start}; {if(\$3<$start) \$3=$start}; {if(\$7<$start) \$7=$start}; {if(\$8<$start) \$8=$start}; {if(\$2>$end) \$2=$end}; {if(\$3>$end) \$3=$end}; {if(\$7>$end) \$7=$end}; {if(\$8>$end) \$8=$end}; \$1==\"$chr\" && \$3>$start && \$2<$end'";

        system "bedtools sort -i $var->{EXON}[0] | bedtools merge -i stdin -c 4 -o distinct | awk 'BEGIN{OFS=\"\\t\";} {print \$1,\$2,\$3,\$4,0,\".\",\$2,\$3,\"0,255,0\"}' | uniq | gzip > $var->{EXON}[0].tmp.gz" unless (-f $var->{EXON}[0] || -f $genomefile);

        system "bedtools sort -i $var->{TRANSCRIPT}[0] | bedtools merge -i stdin -c 4 -o distinct | bedtools subtract -a stdin -b $var->{EXON}[0] | awk 'BEGIN{OFS=\"\\t\";} {print \$1,\$2,\$3,\$4,0,\".\",\$2,\$3,\"255,0,0\"}' | gzip > $var->{INTRON}[0].tmp.gz" unless (-f $var->{INTRON}[0] || -f $genomefile);

        system "bedtools sort -i $var->{TRANSCRIPT}[0] | bedtools complement -i stdin -g hg38.chrom.sizes | awk 'BEGIN{OFS=\"\\t\";} {print \$1,\$2,\$3,\"intergenic\",0,\".\",\$2,\$3,\"0,0,255\"}' | gzip > $var->{INTERGENIC}[0].tmp.gz" unless (-f $var->{INTERGENIC}[0] || -f $genomefile);

        system "bedtools sort -i $var->{CODING}[0] | bedtools merge -i stdin -c 4 -o distinct | awk 'BEGIN{OFS=\"\\t\";} {print \$1,\$2,\$3,\$4,0,\".\",\$2,\$3,\"0,102,0\"}' | gzip > $var->{CODING}[0].tmp.gz" if (-f "$var->{EXON}[0].tmp.gz");

        system "bedtools sort -i $var->{NONCODING}[0] | bedtools subtract -a stdin -b $var->{CODING}[0] | sort -k1,1 -k2,2n | bedtools merge -i stdin -c 4 -o distinct | awk 'BEGIN{OFS=\"\\t\";} {print \$1,\$2,\$3,\$4,0,\".\",\$2,\$3,\"153,255,153\"}' | uniq | gzip > $var->{NONCODING}[0].tmp.gz" if (-f "$var->{EXON}[0].tmp.gz");

        system "$path $var->{FUNC_ELEMENT}[0] hg19ToHg38.over.chain.gz tmp1 tmp2 && mv tmp1 $var->{FUNC_ELEMENT}[0] && rm tmp2" unless (-f $var->{FUNC_ELEMENT}[0]);

        foreach (glob "$dir/*"){
                if ($_ =~ /.gz/){
                       	system "gunzip -f $_";
                       	$_ =~ s{.gz}{};
                }
		unless (-f $genomefile){
                	if ($_ =~ /tmp/){
				system "cat $_ >> $genomefile" unless ($_ =~ /$var->{EXON}[0]|$var->{TRANSCRIPT}[0]/); 
				unlink $_ or die "Can't delete $_: $!\n";
			} else {
				foreach ($var->{INTRON}[0], $var->{INTERGENIC}[0],$var->{EXON}[0]){
					system "cat $_ >> $genomefile";
				}
			}
			system "bedtools sort -i $genomefile > tmp && mv tmp $genomefile";
		}
		system "cat $_ | $command > tmp && mv tmp $_" if ($_ =~ /$genomefile|$var->{PROMOTER}[0]|$var->{REPEAT}[0]|$var->{CIS_REG_MODULE}[0]|$var->{FUNC_ELEMENT}[0]/); #sanity check
		unlink $_ or die "Can't delete $_: $!\n" if ($_ =~ /$var->{CODING}[0]|$var->{NONCODING}[0]/);
	}

}

=begin comment
sub genomeTrack{
	my $file = shift;
	my @genomefiles = ($var->{INTRON}[0], $var->{INTERGENIC}[0], $var->{CODING}[0], $var->{NONCODING}[0]);	

	if (!-f $file){
		foreach (glob "$dir/*"){
                	unless ($_ =~ /$var->{PROMOTER}[0]|$var->{REPEAT}[0]|$var->{CIS_REG_MODULE}[0]|$var->{FUNC_ELEMENT}[0]/){
				unless ($_ =~ /$var->{EXON}[0]|$var->{TRANSCRIPT}[0]/){
					system "cat $_ >> $file";
				}
				unlink $_ or die "Can't delete $_: $!\n";
			}
                }
		system "bedtools sort -i $file > tmp && mv tmp $file";
	}
	
	#foreach (@genomefiles){
	#	unlink $_ if (-f $_ && $_ ne $var->{INTRON}[0]);
	#}
}
=end comment
=cut

sub deleteEmpty{
	foreach (glob "$dir/*"){
		if (`cat $_ | wc -c` == 0){
                        print "$_ is empty... Deleting file.\n";
                        unlink $_ or die "Can't delete $_: $!\n";
                }
	}
}

