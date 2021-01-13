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
use File::Path qw(make_path);
use IO::Compress::Bzip2 qw(bzip2 $Bzip2Error);
use File::Basename;

use constant LABELS => {
        EXON => ['exon.bed', '0,255,0', {group1=>"transcript",group2=>"gene"}], #green
        INTRON => ['intron.bed', '255,0,0', {group1=>"transcript",group2=>"gene"}], #red
        INTERGENIC => ['intergenic.bed', '0,0,255'], #blue
        CODING => ['coding-exon.bed', '0,102,0', {group1=>"transcript",group2=>"gene"}], #dark-green 
        NONCODING => ['noncoding-exon.bed', '153,255,153', {group1=>"transcript",group2=>"gene"}], #light-green
        REPEAT => ['repeat.bed', '255,0,127', {group1=>"class",group2=>"family",group3=>"name"}], #pink
        PROMOTER => ['promoter.bed', '127,0,255', {group1=>"promoter"}], #purple
        FUNC_ELEMENT => ['f-element.bed', '', {group1=>"state"}], #default (Segway)
        CIS_REG_MODULE => ['cr-module.bed', '', {group1=>"factor"}], #default (ReMap)
        TRANSCRIPT => ['transcript.bed', '', {group1=>"transcript",group2=>"gene"}], #not kept, only used to extract intronic/intergenic regions
};

my $history = shift;
my $path = shift;
my $biotype = shift;
my $dir = shift;
my $biotypeR = shift;
my $repfilter = shift // "raw";
my $felefilter = shift // "raw";
my $crmodfilter = shift // "raw";

my ($chr,$start,$end) = split /[:-]+/, $dir;
my @repfilter = $repfilter eq "raw" ? [] : split(/_/,$repfilter);
my @regfilter = $crmodfilter eq "raw" ? [] : split(/_/,$crmodfilter);
my @fefilter = $felefilter eq "raw" ? [] : split(/_/, $felefilter);
@repfilter = @repfilter[1..$#repfilter];
@regfilter = @regfilter[2..$#regfilter];
@fefilter = @fefilter[2..$#fefilter];

my $regulator = \@regfilter;
my $repeat = \@repfilter;
my $felement = \@fefilter;

#append <chr:start-stop@biotype>/labels/ to file names - must be done before declaring %file in order to update file names
my $var = addDir($dir, LABELS);

#declare arrays for each label (to be filled with bedFeature objects)
my (@exon, @transcript, @intron, @intergenic, @coding, @noncoding, @promoter, @repeat, @f_element, @cr_module);
my %files = ( $var->{EXON}[0] => \@exon, $var->{TRANSCRIPT}[0] => \@transcript,
                $var->{INTRON}[0] => \@intron, $var->{INTERGENIC}[0] => \@intergenic,
                $var->{CODING}[0] => \@coding, $var->{NONCODING}[0] => \@noncoding,
                $var->{PROMOTER}[0] => \@promoter, $var->{REPEAT}[0] => \@repeat,
                $var->{FUNC_ELEMENT}[0] => \@f_element, $var->{CIS_REG_MODULE}[0] => \@cr_module, );

my $genomefile = join("/", $dir, "raw_genome.bed");
my $catfile = join("/", $dir, join("_", $biotype, "cat.bed"));
my @keepfiles = ( $genomefile, $catfile, $var->{PROMOTER}[0], $var->{REPEAT}[0], $var->{FUNC_ELEMENT}[0], $var->{CIS_REG_MODULE}[0] );
my @removefiles = ( $var->{EXON}[0], $var->{TRANSCRIPT}[0], $var->{INTRON}[0], $var->{INTERGENIC}[0], $var->{CODING}[0], $var->{NONCODING}[0] );

my @checkfiles;
foreach (@keepfiles){
	push @checkfiles, $_ unless (-f $_);
}

if(@checkfiles){ #there is some new file to create
	my $retrieve = retrieve($history);
    my $hash = @$retrieve[1]; #all previously downloaded bedFeature objects
	#print Dumper $hash;
    writeFiles($hash, \@checkfiles);
}
#deleteEmpty();

my $giggledir = join('/', $dir, "giggle_files");
if (!-d $giggledir){
	mkdir $giggledir or die "Failed to create GIGGLE file directory: $!\n";
	splitCat($giggledir);
}

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
	$new->{'REPEAT'}[0] = join('/', $dir, join('_',$repfilter,$ref->{'REPEAT'}[0]));
	$new->{'FUNC_ELEMENT'}[0] = join('/', $dir, join('_',$felefilter,$ref->{'FUNC_ELEMENT'}[0]));
	$new->{'CIS_REG_MODULE'}[0] = join('/', $dir, join('_',$crmodfilter,$ref->{'CIS_REG_MODULE'}[0]));
    return $new;
}

#sort bedFeature objects for bed file conversion
sub writeFiles{
	my $href = shift;
	my $arref = shift;

	my %hash = %$href;#@_; #hash of bedFeature arrays, where key = database

    foreach my $key (keys %hash){
    	my @array = @{ $hash{$key} };
        foreach (@array){
        	my $hash = $_->decodeName; #$_->{original}
			if ($key eq 'ncbiRefSeq' && (grep /$genomefile|$catfile/, @$arref)){
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
			} elsif ($key eq 'FANTOM5' && (grep /$var->{PROMOTER}[0]|$catfile/, @$arref)){
				$_->{itemRgb} = LABELS->{PROMOTER}[1]; #add colour
            	push @promoter, $_;
            } elsif ($key eq 'RepeatMasker' && (grep /$var->{REPEAT}[0]|$catfile/, @$arref)){
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
			} elsif ($key eq 'SegWay' && (grep /$var->{FUNC_ELEMENT}[0]|$catfile/, @$arref)){
	            unless (@$felement){
	            	push @f_element, $_;
	            } else {
	            	foreach my $fefilter (@$felement){
	            		if ($fefilter =~ /$hash->{name}/i){
	            			push @f_element, $_;
	            		}
	            	}
	            } 
			} elsif ($key eq 'ReMap' && (grep /$var->{CIS_REG_MODULE}[0]|$catfile/, @$arref)){
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
	if (emptyDir($dir)){
		while (my ($key, $value) = each (%files)){
			to_Bed($key, $value) unless (-e $key);
        }
        manipulate_new($dir);
	} else {
    	while (my ($key, $value) = each (%files)){
        	to_Bed($key, $value) unless (-e $key);
        }
        manipulate_filter($dir) if (@repeat);
	}
}

#convert bedFeature objects to bed files (with renamed "name" field)
sub to_Bed{
        my $file = shift;
        my $array = shift;

		if (@$array){
        	open(FH, '>', $file) or die "Could not open file '$file'\n";
        	foreach (@{ $array }){
			my $renamed = $_->renameFeature($file,0);
			print FH $renamed->bed_string;
        	}
        	close(FH);
        }
}

sub manipulate_filter{
	my $dir = shift;
	my $path = join('', $path, 'liftOver');
	
	system "$path $var->{FUNC_ELEMENT}[0] hg19ToHg38.over.chain.gz tmp1 tmp2 && mv tmp1 $var->{FUNC_ELEMENT}[0] && rm tmp2" if -e $var->{FUNC_ELEMENT}[0];
#	deleteEmpty($dir);
}
		
#create initial files and apply liftOver
sub manipulate_new{
	my $dir = shift;
    my $path = join('', $path, 'liftOver');
	my $command = "awk -F \"\\t\" 'BEGIN{OFS=\"\\t\";}; {if(\$2<$start) \$2=$start}; {if(\$3<$start) \$3=$start}; {if(\$7<$start) \$7=$start}; {if(\$8<$start) \$8=$start}; {if(\$2>$end) \$2=$end}; {if(\$3>$end) \$3=$end}; {if(\$7>$end) \$7=$end}; {if(\$8>$end) \$8=$end}; \$1==\"$chr\" && \$3>$start && \$2<$end'";
	my @genomefiles = ( $var->{EXON}[0], $var->{INTRON}[0], $var->{INTERGENIC}[0] );
	my @catfiles = ( $genomefile, $var->{TRANSCRIPT}[0], $var->{PROMOTER}[0], $var->{REPEAT}[0], $var->{FUNC_ELEMENT}[0], $var->{CIS_REG_MODULE}[0] ); 

	system "bedtools sort -i $var->{EXON}[0] | awk 'BEGIN{OFS=\"\\t\";} {print \$1,\$2,\$3,\$4,0,\".\",\$2,\$3,\"@{(LABELS->{EXON})}[1]\"}' | uniq | gzip > $var->{EXON}[0].gz" if -e $var->{EXON}[0];
    system "bedtools sort -i $var->{TRANSCRIPT}[0] | bedtools subtract -a stdin -b $var->{EXON}[0] | awk 'BEGIN{OFS=\"\\t\";} {print \$1,\$2,\$3,\$4,0,\".\",\$2,\$3,\"@{(LABELS->{INTRON})}[1]\"}' | gzip > $var->{INTRON}[0].gz" if -e $var->{TRANSCRIPT}[0];
    system "bedtools sort -i $var->{TRANSCRIPT}[0] | bedtools complement -i stdin -g hg38.chrom.sizes | awk 'BEGIN{OFS=\"\\t\";} {print \$1,\$2,\$3,\"intergenic\",0,\".\",\$2,\$3,\"@{(LABELS->{INTERGENIC})}[1]\"}' | gzip > $var->{INTERGENIC}[0].gz" if -e $var->{TRANSCRIPT}[0];
    system "bedtools sort -i $var->{CODING}[0] | awk 'BEGIN{OFS=\"\\t\";} {print \$1,\$2,\$3,\$4,0,\".\",\$2,\$3,\"@{(LABELS->{CODING})}[1]\"}' | gzip > $var->{CODING}[0].gz" if -e $var->{CODING}[0];
    system "bedtools sort -i $var->{NONCODING}[0] | bedtools subtract -a stdin -b $var->{CODING}[0] | sort -k1,1 -k2,2n | awk 'BEGIN{OFS=\"\\t\";} {print \$1,\$2,\$3,\$4,0,\".\",\$2,\$3,\"@{(LABELS->{NONCODING})}[1]\"}' | uniq | gzip > $var->{NONCODING}[0].gz" if -e $var->{NONCODING}[0];
    system "$path $var->{FUNC_ELEMENT}[0] hg19ToHg38.over.chain.gz tmp1 tmp2 && mv tmp1 $var->{FUNC_ELEMENT}[0] && rm tmp2" if -e $var->{FUNC_ELEMENT}[0];
	system "awk 'BEGIN{OFS=\"\\t\";} {gsub(/intron/,\"transcript\",\$4)}; {print}' $var->{TRANSCRIPT}[0] | gzip > $var->{TRANSCRIPT}[0].gz" if -e $var->{TRANSCRIPT}[0];
	system "sed -i /hg/d $var->{PROMOTER}[0]" if -e $var->{PROMOTER}[0];
	system "cat $var->{CODING}[0] $var->{NONCODING}[0] > $var->{EXON}[0] && rm *$var->{CODING}[0]* && rm *$var->{NONCODING}[0]*" if (-e $var->{CODING}[0] && -e $var->{NONCODING}[0]);

	foreach (glob "$dir/*.bed"){
		unlink $_ or die "Can't delete $_: $!\n" unless $_ =~ /$var->{PROMOTER}[0]|$var->{REPEAT}[0]|$var->{FUNC_ELEMENT}[0]|$var->{CIS_REG_MODULE}[0]/;
	}
	foreach (glob "$dir/*.{bed,bed.gz}"){
		if ($_ =~ /.gz/){ 
			system "gunzip -f $_";
			$_ =~ s{.gz}{};
			system "cat $_ | $command > tmp && mv tmp $_";
		}
	}
	
	unless (-e $genomefile){
		foreach (@genomefiles){
			system "cat $_ >> $genomefile";
		}
	}
	unless (-e $catfile){
		foreach (@catfiles){
			system "awk -F \"\\t\" 'BEGIN{OFS=\"\\t\";} {print \$1,\$2,\$3,\$4 >> \"$catfile\"}' $_" ;
		}
	}
	foreach my $f (glob "$dir/*.bed"){
		unless (grep {$_ eq $f} @removefiles){
			system "bedtools sort -i $f > tmp && mv tmp $f";
		} else {
			system "rm $f";
		}
	}
}

#empty helper functions
sub emptyDir{
	my $dirname = shift;
    opendir(my $dh, $dirname) or die "Not a directory";
    return scalar(grep { $_ ne "." && $_ ne ".." } readdir($dh)) == 0;
}

sub deleteEmpty{
	my $dirname = shift;
	foreach (glob "$dirname/*"){
    	unless (-d $_){
    		if (`cat $_ | wc -c` == 0){
        		print "$_ is empty... Deleting file.\n";
            	unlink $_ or die "Can't delete $_: $!\n";
			}
		}
	}
}

sub create_empty_file {
	eval {
		open my $fh, '>', $_[0] or die "Cannot create $_[0]: $!\n";
		close $fh or die "Cannot close $_[0]: $!\n";
	};
     return $@;
}

#create giggle directory from cat file 
sub splitCat{
	my $dir = shift;
	my $tmp = File::Temp->newdir();

	if(emptyDir($dir)){
		print "Creating, sorting, and bgzipping giggle files\n";
    	open(FH, '<', $catfile) or die "Can't open '$catfile': $!";
    	while(<FH>){
       		my $string = $_;
        	chomp $string;
        	my @fields = split('\t', $string);
        	my ($chrom, $start, $end, $name) = @fields[0,1,2,3];
        	my @id = split(':', $name);
        	my $categorydir = join ('/',$tmp, $id[0]);
        	if (!-d $categorydir){
          		mkdir $categorydir or die "Failed to create directory: $!\n";
			}
			for (my $index = 1; $index <= $#id; $index++){
				my $groupdir = join ('/', $categorydir, "group" . "$index");
				mkdir $groupdir or die "Failed to create directory: $!\n" unless (-d $groupdir);
            	my $class = $id[$index];
            	my $fh = join ('/', $groupdir, "${class}.bed");
            	open(OUT,'>>',$fh);
            	print OUT $string . "\n";
			}
		}
    	close(FH);
    	renameCat($tmp, $dir, LABELS);
    }
}

sub renameCat{
	my $tmp = shift;
	my $dir = shift;
	my $ref = shift;
	my @subdirs = glob("${tmp}/*");
		
	while (my ($key, $value) = each ( %{ $ref } )){
		my $fname = substr($ref->{$key}[0],0,index($ref->{$key}[0],'.'));
		my $hash = $ref->{$key}[2];
		foreach my $directory (@subdirs){
			my $element = (split('/', $directory))[-1];
			my @groups = glob("${directory}/*");
			if ($fname eq $element){
				foreach(@groups){
					my $group = (split('/', $_))[-1];
					if(defined $hash->{$group}){
						my $new = join('/',$dir,$element,$hash->{$group});
						make_path($new) unless (-d $new);
						system "cp -r ${_}/* $new";
					}
				}
			}
		}
	}
}
