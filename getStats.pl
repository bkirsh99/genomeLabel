#!/usr/bin/perl

use JSON;
use strict;
use warnings;
use Data::Dumper;
use 5.010;
use 5.005;
use Storable;
use Readonly;
use Text::Table;
use List::Util qw(max);
use lib '.';
use bedFeature;

my $dir = 'chrX:15200000-15800000';
my $biotype = 'K562';
my ($chr,$start,$end) = split /[:-]+/, $dir;
#print Dumper $chr,$start,$end;

#my @options = shift;

##DEBUGGING:
#print Dumper $trackdir;
#print Dumper $regulator, $repeat, $felement;
#my $out = substr($dir,0,index($dir,"/"));
my $stats = join('/',$dir,'stats.txt');
my %fileStats;
my %overlapStats;
my @categories = ('promoter','transcript','coding','noncoding','intron');

# > Comppute statistics:
unless (-e $stats){ 
	my $r = create_empty_file($stats);
        die $r if $r;
	#stats();
} else {
	#%fileStats = %{ retrieve($stats) };
}
#my $statsfiles = select_files($dir);
#print Dumper $statsfiles;
gene_table($dir);
#my ($table1,$index) = region_tbl(\%fileStats);

# > Print table(s):
#print_main($table1,$index);
#print_cat($table1,$index); #all main and subs in one table
#print_indv($table1,$index); #all main and subs in separate tables
#my @select = ("promoter", "transcript", "intron");
#print_indv($table1,$index,\@select);
#print_indv($table1,$index);
#print Dumper(\%fileStats;
#gene_summary(\%fileStats); #gene stats

sub create_empty_file {       
	eval {
	open my $fh, '>', $_[0]
        or die "Cannot create $_[0]: $!\n";
        close $fh or die "Cannot close $_[0]: $!\n";
        };
        return $@;
}

sub select_files{
	my $dir = shift;
	my @files = glob "$dir/*.bed";
	my @selectfiles;
        for (0..$#files){
		my $fh = substr($files[$_] , index($files[$_],'/') + 1);
        	my $count = ($fh =~ tr/_//);
		if ($count == 1 || ($count ==2 && substr($fh,8) =~ /$biotype/)){ #unfiltered 
			push @selectfiles, $fh;
		}       
        }
	return \@selectfiles;
}

sub gene_table{
	my $dir = shift;
	my $genomefile = join('/',$dir,'raw_genome.bed');
	my $transcriptfile = join('/',$dir,'raw_transcript.bed');
	my %genehash;
	my @genearray;

	open(FH, '<', $transcriptfile) or die "Can't open '$transcriptfile': $!";
	while(<FH>){
		my $string = $_;
		chomp $string;
		my @fields = split('\t', $string);
        	my ($chrom, $start, $end, $name) = @fields[0,1,2,3];
		my ($category, $transcript, $gene) = split(':', $name);
		push @{ $genehash{$gene} }, { transcript => $transcript, chromosome => $chrom, start => $start,
			end => $end, gene_length => ($end - $start), status => ($transcript =~ /NM/ ? "protein-coding" : "non-coding") };
	}
	close(FH);

	open(FH, '<', $genomefile) or die "Can't open '$genomefile': $!";
        while(<FH>){
                my $string = $_;
                chomp $string;
                my @fields = split('\t', $string);
                my ($chrom, $start, $end, $name) = @fields[0,1,2,3];
                my @elements = split(',',$name);
                foreach my $ele (@elements){
                        my ($category, $transcript, $gene) = split(':', $ele);
                	print Dumper $category,$transcript,$gene;
			#foreach (keys %genehash){
				my ($match) = grep { $transcript eq $_->{transcript} } @{ $genehash{$gene} } if $transcript;
				my $size = $end - $start; 
				++$$match{$category};
				$$match{'intron_length'} += $size if $category eq 'intron';
				$$match{'exon_length'} += $size if $category =~ /exon/;
				$$match{'transcript_length'} += $size if $category eq 'coding-exon';
				#print "element is $chrom, $start, $end, $ele\n";
				#print Dumper $match;
				#print "match!: $gene equals $_ and $transcript is $match\n" if ($transcript, $gene);
			#}
		}
	}
	close(FH);
	print Dumper \%genehash;
}
 
sub gene_summary{
        my $ref = shift;
        my @abs = ('absolute','absCount','absBps','absPct');
        my @rel = ('relative','relCount','relBps','relPct');
        my @genes = ();
        my @cols = ();
        my @rows = ();
        my $load = [];
        my %hash;

        push @cols, +{title => "Gene", align => "center", align_title => "center", align_title_lines => "center"},
		+{title => "promoter", align => "center", align_title => "center"}, +{title => "transcript", align_title => "center", align => "center"}, 
		+{title => "coding", align_title => "center", align => "center"}, +{title => "noncoding", align_title => "center", align => "center"},
		+{title => "intron", align_title => "center", align => "center"};

        my $genelist = $ref->{'transcript'}->{$rel[0]};
        foreach ( sort keys %{$genelist} ){
                push @genes, $_;
	}        

        my $tb = Text::Table->new( ( map { +( ( ref($_) ? $_ : "$_" ) , \'|' ) } @cols ) );

        foreach my $ele1 (keys %{$ref}){ #major: intergenic, repeat, exon, intron, cr, fe, promoter
                if ($ele1 =~ /exon|intron|transcript/){
                        my @rules = ();
                        my $abs1 = $ref->{$ele1}->{$abs[0]};
                        my $rel1 = $ref->{$ele1}->{$rel[0]};
                        foreach my $ele2 (keys %$rel1){ #sub1: coding/noncoding exon or gene (intron or transcript)
                                if (grep {$ele2 =~ /$_/} @genes){ #is a gene
                                        $load = [ "$rel1->{$ele2}->{$rel[1]}\n$rel1->{$ele2}->{$rel[2]}\n$rel1->{$ele2}->{$rel[3]}" ];
                                        $hash{$ele2}{$ele1} = $load;
                                } else { #is coding/noncoding, must go down a level further
                                        my $rel2 = $rel1->{$ele2}->{$rel[0]};
                                        foreach my $ele3 (keys %$rel2){ #sub2: repeat family, gene (exon)
                                                if (grep {$ele3 =~ /$_/} @genes){
                                                        $load = [ "$rel2->{$ele3}->{$rel[1]}\n$rel2->{$ele3}->{$rel[2]}\n$rel2->{$ele3}->{$rel[3]}" ];
                                                        $hash{$ele3}{$ele2} = $load;
                                                }
                                        }
                                }
                        }
                } elsif ($ele1 =~ /promoter/){
                        my $rel1 = $ref->{$ele1}->{$rel[0]};
                        foreach my $ele2 (sort keys %$rel1){ #sub1: gene
                                foreach my $gene (keys %hash){
                                        $load = "$ele2 ($rel1->{$ele2}->{'coords'})";
                                        push (@{ $hash{$gene}{$ele1} }, $load) if substr($ele2,index($ele2,"@")+1) eq $gene;
                                }
                        }
                }

        }

        foreach my $g (sort keys %hash){
                my @row = ();
		push @row, $g;
		foreach (@categories){
                        push @{ $hash{$g}{$_} }, "NA" unless defined $hash{$g}{$_};
                        if ($_ =~ /promoter/ && defined $hash{$g}{$_}){
				@{ $hash{$g}{$_} } = join("\n" , @{ $hash{$g}{$_} });
                        }
			push @row, pop @{ $hash{$g}{$_} };
		}
		push @rows, \@row;
        }

        my %indexes = (
                0 => 1, # Before header row
                1 => 1,  # After header row
		2 => 1   # After units row
        );

        # Calculate body rows for which rules should be printed
        my $count = 2;
        foreach my $row (@rows) {
		my $newlines = max map { tr/\n// } @$row;
		$count += $newlines + 1;
                $indexes{$count} = 1;
        }
	$tb->load( ["Name","FANTOM5 Coordinates","count|avgBps|avgPct", "count|totBps|totPct", "count|totBps|totPct", "count|totBps|totPct", "count|totBps|totPct"] );
        $tb->load(@rows);

        my $rule = $tb->rule('-', '+');
	foreach my $j (0 .. $tb->height) {
		print $rule if exists $indexes{$j};
                print $tb->table($j);
	}

}

sub print_indv2{
	my $tb = $_[0];
	my $idx = $_[1];
	my %index = %{ $idx};

	print $tb->title();
	print $tb->rule('=' , '+');
	my $i = 0;
	OUTER: while ($i <= $tb->body_height) {
		unless (defined $index{$i}){
			print $tb->body($i);
		} elsif ($index{$i} == 0){
			if (defined $_[2]){
				INNER: foreach my $ele (@{ $_[2] }){
						unless ($tb->body($i-1) =~ /$ele/){
							++$i until ($tb->body($i - 1) =~ /$ele/ ||  $i == $tb->body_height);
						}
						#print $tb->body($i-1);
						print $tb->body_rule('=','+');
						shift @{ $_[2] };
						++$i;
						next OUTER;
					}
			} else {
					print $tb->body_rule('=','+');
			}
		} elsif ($index{$i} == 1){
			print $tb->body_rule('.', '+');
		} elsif ($index{$i} == 2){
			print $tb->body_rule('-','+');
		} elsif ($index{$i} == 3){
			print $tb->body_rule('.', '+');
		} elsif ($index{$i} == 4 && $i != $tb->body_height-1){
				print "\n\n";
                                print $tb->rule('=' , '+');
			        print $tb->title();
                                print $tb->rule('=' , '+');
		} else {
			print $tb->body_rule('+') unless defined $index{$i-1} == 3;
		}
		++$i;
	}
}

sub print_indv{
        my $tb = shift;
        my $idx = shift;
        my %index = %{ $idx};

        print $tb->title();
        print $tb->rule('=' , '+');
	foreach my $i (0 .. $tb->body_height) {
		unless (defined $index{$i}){
			print $tb->body($i);
		} else {
			$index{$i} == 0 ? print $tb->body_rule('=','+') : $index{$i} == 1 ? print $tb->body_rule('.', '+') : $index{$i} == 2 ? print $tb->body_rule('-','+') : $index{$i} == 3 ? print $tb->body_rule('.','+') : print $tb->body_rule('+') unless defined $index{$i-1} == 3;
			if ($index{$i} == 4 && $i != $tb->body_height-1){
				print "\n\n";
				print $tb->rule('=' , '+');
				print $tb->title();
        			print $tb->rule('=' , '+');
			}
		}
	}
}


sub print_cat{
	my $tb = shift;
	my $idx = shift; 
	my %index = %{ $idx};

	print $tb->title();
        print $tb->rule('=' , '+');
        foreach my $i (0 .. $tb->body_height) {
                if (defined $index{$i}){
                        $index{$i} == 0 ? print $tb->body_rule('=','+') : $index{$i} == 1 ? print $tb->body_rule('.','+') : $index{$i} == 2 ? print $tb->body_rule('-','+') : $index{$i} == 3 ? print $tb->body_rule('.','+') : print $tb->body_rule('+') unless defined $index{$i-1} == 3;                
		} else {
                        print $tb->body($i);
                }
        }
}

sub stats{
	my @files = glob "$dir/*.bed";
	for (0..$#files){
		coverage($files[$_]);
		promoterStats('promoter', $files[$_]) if $files[$_] =~ /promoter/;
	} 
	groupRepeat('repeat');
	groupTranscript('transcript');
	groupExon('coding','noncoding');
	groupIntron('intron');
	store(\%fileStats,$stats);
}

sub promoterStats{
	my $category = shift;
	my $file = shift;
	my %coords;

	open(FH, '<', $file) or die "Can't open '$file': $!";
        while(<FH>){
		my $string = $_;
		chomp $string;
		my @fields = split('\t', $string);
		my ($chrom, $start, $end) = @fields[0,1,2];
		my $coord = join(":",$chrom,join("-",$start,$end));
		my $name = index($fields[3],".") > 0 ? substr($fields[3],index($fields[3],".")+1) : $fields[3];
		$coords{$name} = $coord;
	}
	close(FH);
	foreach (keys %coords ){
		$fileStats{$category}->{'relative'}{$_}{'coords'} = $coords{$_} if defined $fileStats{$category}->{'relative'}->{$_};
		$overlapStats{$category}->{'promoter'}->{$_}->{'ref'}->{'coords'} = $coords{$_} if defined $overlapStats{$category}->{'promoter'}->{$_};
	}
}

sub promoterStats2{
	my $category = shift;
        my $refcategory = shift;
        my $refsubcategory = shift;
	my $file = shift;
        my %coords;

        open(FH, '<', $file) or die "Can't open '$file': $!";
        while(<FH>){
                my $string = $_;
                chomp $string;
                my @fields = split('\t', $string);
                my ($chrom, $start, $end) = @fields[0,1,2];
                my $coord = join(":",$chrom,join("-",$start,$end));
                my $name = index($fields[3],".") > 0 ? substr($fields[3],index($fields[3],".")+1) : $fields[3];
                $coords{$name} = $coord;
        }
        close(FH);
        foreach (keys %coords ){
                $overlapStats{$category}{$refcategory}{$refsubcategory}->{'map'}->{'promoter'}->{$_}->{'coords'} = $coords{$_} if defined $overlapStats{$category}{$refcategory}{$refsubcategory}->{'map'}->{'promoter'}->{$_}; 
        }
}

sub region_tbl{
	my $ref = shift;
	my @abs = ('absolute','absCount','absBps','absPct');
	my @rel = ('relative','relCount','relBps','relPct');
	my @cols = ();
	my @data = ();
	my $load = [];
	my %index;

	push @cols, +{title => "Category", align => "left", align_title => "center", align_title_lines => "center"}, 
		+{title => "cvgCount", align => "center", align_title => "center", align_title_lines => "center"}, 
		+{title => "cvgBps", align => "center", align_title => "center", align_title_lines => "center"}, +{title => "cvgPct", align => "center", align_title => "center", align_title_lines => "center"};

	my $tb = Text::Table->new( ( map { +( ( ref($_) ? $_ : "$_" ) , \'|' ) } @cols ) );
	
	foreach my $ele1 (keys %{$ref}){ #major: intergenic, repeat, exon, intron, cr, fe, promoter
		my @rules = ();
		my $abs1 = $ref->{$ele1}->{$abs[0]};
		my $rel1 = $ref->{$ele1}->{$rel[0]};
		$load = [ "$ele1", $abs1->{$abs[1]}, $abs1->{$abs[2]}, $abs1->{$abs[3]} ];
		push @data, $load;
		push @data, 0;
		foreach my $ele2 (keys %$rel1){ #sub1: repeat class, coding/noncoding exon , genomic state, TF, gene (intron, promoter, transcript)
			if (defined $rel1->{$ele2}->{$abs[0]}){
				my $abs2 = $rel1->{$ele2}->{$abs[0]};
				my $rel2 = $rel1->{$ele2}->{$rel[0]};
				$load = [ "  $ele2", $abs2->{$abs[1]}, $abs2->{$abs[2]}, $abs2->{$abs[3]} ];
				push @data, $load unless $ele2 eq $ele1;
				push @data, 1 unless $ele2 eq $ele1;
				foreach my $ele3 (keys %$rel2){ #sub2: repeat family, gene (exon)
					if (defined $rel2->{$ele3}->{$abs[0]}){
						my $abs3 = $rel2->{$ele3}->{$abs[0]};
						my $rel3 = $rel2->{$ele3}->{$rel[0]};
						$load = [ "    $ele3", $abs3->{$abs[1]}, $abs3->{$abs[2]}, $abs3->{$abs[3]} ]; 
						push @data, $load unless $ele3 eq $ele2;
						push @data, 2 unless $ele3 eq $ele2;
						foreach my $ele4 (keys %$rel3){ #sub3: repeat name
							$load = [ "      $ele4", $rel3->{$ele4}->{$rel[1]}, $rel3->{$ele4}->{$rel[2]}, $rel3->{$ele4}->{$rel[3]} ];
							push @data, $load unless $ele4 eq $ele3;
						
						}
						push @data, 3;
					} else {
						$load = [ "    $ele3", $rel2->{$ele3}->{$rel[1]}, $rel2->{$ele3}->{$rel[2]}, $rel2->{$ele3}->{$rel[3]} ];
						push @data, $load unless $ele3 eq $ele2;
					}
				}
			} else {
				$load = [ "  $ele2", $rel1->{$ele2}->{$rel[1]}, $rel1->{$ele2}->{$rel[2]}, $rel1->{$ele2}->{$rel[3]} ];
				push @data, $load unless $ele1 eq $ele2;
			}
		}
		push @data, 4;
	}
	my $count = 0;
	#Indexes: 0 = after major, 1 = after sub1, 2 = after sub2, 3 = after sub3, 4 = end of major
	foreach my $row (@data) {
		$index{$count} = $row unless (ref $row eq 'ARRAY');
		$tb->load($row);
		++$count;
	}

	return($tb, \%index);
}	

sub coverage{
        my $file = shift;
	my $out = `echo "$chr\t$start\t$end\tfoo-1" | bedmap --echo-ref-row-id --echo-ref-size --bases-uniq --bases-uniq-f --echo-map-id --echo-overlap-size --count - $file`;
        my ($count, %eleStats);
	my $checkSum = 0;
	
        my @fields = split(/\|/, $out); #split the bedmap outputs
	my $size = $fields[1]; #total size (bps) of input region
	my $cvg = $fields[2]; #unique bps in input region covered by elements of $file
	my $pct = sprintf("%.2f", ($fields[3] *100)); #percent coverage
	my $uniq = $fields[6]; #unique element count

	my ($category, $subcategory);
	my @categories;
	my @element = split(";",$fields[4]) if defined $fields[4];
	my @elementcvg = split(";",$fields[5]) if defined $fields[5]; 
	foreach (@element){
		$category = index($_,".") > 0 ? substr( $_,0,index( $_,".")) : $_;
		push @categories, $category;
		$eleStats{$category}->{'absolute'}->{'absBps'} = $cvg;
       		$eleStats{$category}->{'absolute'}->{'absPct'} = $pct;
		$eleStats{$category}->{'absolute'}->{'absCount'} = $uniq;
		my $subcategory = index($_,".") > 0 ? substr( $_,index( $_,".")+1) : $_;
		$eleStats{$category}->{'relative'}->{$subcategory}->{'relBps'} = defined $eleStats{$category}->{'relative'}->{$subcategory} ? ($eleStats{$category}->{'relative'}->{$subcategory}->{'relBps'} + shift @elementcvg) : shift @elementcvg;
		$eleStats{$category}->{'relative'}->{$subcategory}->{'relCount'}++;
	}
	foreach $category (@categories){
		foreach (keys %{ $eleStats{$category}->{'relative'} }){
			$eleStats{$category}->{'relative'}->{$_}->{'relPct'} = sprintf("%.2f", (($eleStats{$category}->{'relative'}->{$_}->{'relBps'} / $size) * 100));
		}
	$fileStats{$category} = $eleStats{$category};
	}

}

sub groupRepeat{ #breakdown into class -> family -> name
	my $category = shift;
	my (%h1, %h2, %h3, %seen1, %seen2, %seen3, @arr1, @arr2);

	foreach (keys %{ $fileStats{$category}->{'relative'} }){
		my ($class, $family, $name) = split(":",$_);
		push @{ $h2{$family} }, $name unless defined $seen3{$name};
		push @{ $h1{$class} }, $family unless defined $seen2{$family};
		$seen3{$name} = $fileStats{$category}->{'relative'}{$_};
                ++$seen2{$family};
		++$seen1{$class};  
	}
	
	my $hash;
	foreach my $class (keys %seen1){
		my $families = $h1{$class};
		foreach my $family (@$families){
			my $names = $h2{$family};
			foreach my $name (@$names){
				$$hash{$class}{$family}{$name} = $seen3{$name} if (defined $seen3{$name});# && defined $h2{$family} && defined $h1{$class});
			}
		}
	}
	delete($fileStats{$category}->{'relative'});
	$fileStats{$category}->{'relative'} = $hash;
	
	foreach my $class (keys %{ $fileStats{$category}->{'relative'} }){
                foreach my $family (keys %{ $fileStats{$category}->{'relative'}{$class} }){
			foreach my $name (keys %{ $fileStats{$category}->{'relative'}{$class}{$family} }){ 
                        	$fileStats{$category}->{'relative'}{$class}->{'relative'}->{$family}->{'relative'}->{$name} = $fileStats{$category}->{'relative'}{$class}{$family}{$name};
                        	$fileStats{$category}->{'relative'}{$class}->{'relative'}->{$family}->{'absolute'}->{'absBps'} += $fileStats{$category}->{'relative'}{$class}{$family}{$name}->{'relBps'};
				$fileStats{$category}->{'relative'}{$class}->{'relative'}->{$family}->{'absolute'}->{'absPct'} += $fileStats{$category}->{'relative'}{$class}{$family}{$name}{'relPct'};
				$fileStats{$category}->{'relative'}{$class}->{'relative'}->{$family}->{'absolute'}->{'absCount'} += $fileStats{$category}->{'relative'}{$class}{$family}{$name}{'relCount'};
				delete($fileStats{$category}->{'relative'}{$class}{$family}{$name});
			}
			$fileStats{$category}->{'relative'}{$class}->{'absolute'}->{'absBps'} += $fileStats{$category}->{'relative'}{$class}->{'relative'}->{$family}->{'absolute'}->{'absBps'};
			$fileStats{$category}->{'relative'}{$class}->{'absolute'}->{'absPct'} += $fileStats{$category}->{'relative'}{$class}->{'relative'}->{$family}->{'absolute'}->{'absPct'};
			$fileStats{$category}->{'relative'}{$class}->{'absolute'}->{'absCount'} += $fileStats{$category}->{'relative'}{$class}->{'relative'}->{$family}->{'absolute'}->{'absCount'};
			delete($fileStats{$category}->{'relative'}{$class}{$family});
        	}
	}
}

sub groupTranscript{ #output average size and pct coverage based on count
	my $category = shift;

	foreach my $gene (keys %{ $fileStats{$category}->{'relative'} }){
		foreach (keys %{ $fileStats{$category}->{'relative'}{$gene} }){
			$fileStats{$category}->{'relative'}{$gene}{$_} = sprintf("%.2f", ($fileStats{$category}->{'relative'}{$gene}{$_} / $fileStats{$category}->{'relative'}{$gene}->{'relCount'}) ) unless $_ eq 'relCount';
		}
	}
} 

sub groupExon{ #group coding and noncoding stats into exon breakdown
        my @categories = @_;

        foreach my $subcategory1 (@categories){ #different coding status
		foreach my $stat (keys %{ $fileStats{$subcategory1}->{'absolute'} }){
			$fileStats{'exon'}->{'absolute'}{$stat} = $fileStats{$subcategory1}->{'absolute'}{$stat};
		}
	        $fileStats{'exon'}->{'relative'}->{$subcategory1} = $fileStats{$subcategory1}->{'relative'};
	}
	foreach my $subcategory1 (keys %{ $fileStats{'exon'}->{'relative'} }){
		foreach my $subcategory2 (keys %{ $fileStats{'exon'}->{'relative'}{$subcategory1} }){ #different genes
			$fileStats{'exon'}->{'relative'}{$subcategory1}->{'relative'}->{$subcategory2} = $fileStats{'exon'}->{'relative'}{$subcategory1}{$subcategory2};
			$fileStats{'exon'}->{'relative'}{$subcategory1}->{'absolute'}->{'absBps'} += $fileStats{'exon'}->{'relative'}{$subcategory1}{$subcategory2}->{'relBps'};
			$fileStats{'exon'}->{'relative'}{$subcategory1}->{'absolute'}->{'absPct'} += $fileStats{'exon'}->{'relative'}{$subcategory1}{$subcategory2}->{'relPct'};
			$fileStats{'exon'}->{'relative'}{$subcategory1}->{'absolute'}->{'absCount'} += $fileStats{'exon'}->{'relative'}{$subcategory1}{$subcategory2}->{'relCount'};
			delete($fileStats{'exon'}->{'relative'}{$subcategory1}{$subcategory2});
		}
	}
	foreach my $subcategory1 (keys %{ $fileStats{'exon'}->{'relative'} }){
                $fileStats{'exon'}->{'absolute'}{'absCount'} += $fileStats{'exon'}->{'relative'}{$subcategory1}->{'absolute'}->{'absCount'};
	}
	foreach (@categories){
		delete($fileStats{$_});
	}
}

sub groupIntron{
	my $category = shift;
	my %add;

	foreach my $subcategory (keys %{ $fileStats{$category}->{'relative'} } ){ #gene or multiple genes
		my @genes = index($subcategory,",") > 0 ? split(",", $subcategory) : ();
		if (@genes){
			foreach (@genes){
				my $extra = index($_,".") > 0 && substr($_,0,index($_,".")) eq $category ? substr($_,index( $_,".")+1) : $_; #for exceptions in gene names with a "dot" (e.g. GS1-594A7.3)
				foreach my $key (keys %{ $fileStats{$category}->{'relative'}{$subcategory} }){
					$fileStats{$category}->{'relative'}{$extra}{$key} += $fileStats{$category}->{'relative'}{$subcategory}{$key}; 
				}
			}
			delete($fileStats{$category}->{'relative'}{$subcategory});
		}
	}
}

my $f2 = join("@","chrX:15200000-15800000","GM12878/labels/exon.bed");
my $f1 = join("@","chrX:15200000-15800000","GM12878/labels/promoter.bed");
my $f3 = join("@","chrX:15200000-15800000","GM12878/labels/transcript.bed");
#intersect($f1,$f2);
#intersect($f1,$f3);
sub intersect{
        my $ref = shift;
        my $map = shift;
        my $name = join("~",$ref,$map);
        my $out = `bedmap --echo-ref-row-id --echo-ref-size --bases-uniq --bases-uniq-f --echo-map-id --echo-overlap-size $ref $map`;
        my @lines = split "\n", $out; #split output into different lines
	my %eleStats = ();
	my ($category, $subcategory);
	
        foreach (@lines){
                my $string = $_;
                my @fields = split(/\|/, $string); #split the bedmap outputs 
                my ($refid, $refsize, $refbps, $refpct, $mapele, $mapbps) = @fields;
                $refid = substr($refid, 3); #remove leading "id-"
		$refid = readLine($ref, $refid);
		my $refcat = index($refid,".") > 0 ? substr($refid,0,index($refid,".")) : $refid;
		my $refsubcat = index($refid,".") > 0 ? substr($refid,index($refid,".")+1) : $refid;
		$refpct = sprintf("%.2f", $refpct * 100);
		unless ($refpct == 0){
                	++$eleStats{$refcat}{$refsubcat}->{'ref'}->{'count'};
			$eleStats{$refcat}{$refsubcat}->{'ref'}->{'totBps'} += $refsize;
			$eleStats{$refcat}{$refsubcat}->{'ref'}->{'cvgBps'} += $refbps;
			$eleStats{$refcat}{$refsubcat}->{'ref'}->{'cvgPct'} += $refpct;
			my @mapelements = split(";",$mapele) if defined $mapele;
			my @mapbases = split(";",$mapbps) if defined $mapbps;
			foreach (@mapelements){
				$category = index($_,".") > 0 ? substr( $_,0,index( $_,".")) : $_;
				$subcategory = index($_,".")> 0 ? substr( $_,index( $_,".")+1) : $_;
				$eleStats{$refcat}{$refsubcat}->{'map'}->{$category}{$subcategory}->{'cvgBps'} = defined $eleStats{$refcat}{$refsubcat}->{'map'}->{$category}{$subcategory}->{'cvgBps'} ? ($eleStats{$refcat}{$refsubcat}->{'map'}->{$category}{$subcategory}->{'cvgBps'} + shift @mapbases) : shift @mapbases; 
			}
		}
	}
	$overlapStats{$name} = \%eleStats;
	my %refmap;
 	#print Dumper \%overlapStats;
	foreach my $refs (sort keys %eleStats){
		promoterStats($name,$ref) if $refs eq 'promoter';
		#groupRepeat($name) if $refs eq 'repeat';
		my $refhash = $eleStats{$refs};
		foreach my $keys (sort keys %$refhash){
			my $maphash = $$refhash{$keys}{'map'};
			foreach (sort keys %$maphash){
				promoterStats2($name,$refs,$keys,$map) if $_ eq 'promoter';
			}
		}
	}
	#print Dumper \%refmap;
	#promoterStats(,$ref) if $refs eq 'promoter';
	#print Dumper \%overlapStats;
	printOverlap(\%eleStats);
}

sub readLine{
        my $filename = shift;
        my $line_number = shift;
        my $line = '';
        my @fields = ();

        open(FH, "<", $filename) or die "Can't open $filename for reading: $!\n";
        while (<FH>) {
                $line = $_;
                last if ($. == $line_number);
        }
	close(FH);
        if ($line eq '') {
		die "Didn't find line $line_number in $filename\n";
        } else {
                chomp $line;
		@fields = split('\t', $line);
                return $fields[3];
        }
}
=begin comment
sub printBreakdown{
        my $hash = shift;
        my $flag = shift;
        my %index;
        my ($start, $end) = (0,0);

        my @cols = ();
        push @cols,  +{title => "Count", align => "center"}, +{title => "totalBps", align => "center"},
    +{title => "pctBps", align => "center"};

        my $tb        = Text::Table->new( "Category", \'||',
    ( map { +( ( ref($_) ? $_ : "$_" ) , \'|' ) } @cols ) );

        foreach my $category ( sort keys %{ $hash } ) {
        $start = $start + $end;
        $tb->load( [
                   majorElement(LABELS, $category),
                   $$hash{$category}->{'absolute'}->{'absCount'},
                   $$hash{$category}->{'absolute'}->{'absBps'},
                   sprintf("%.2f", $$hash{$category}->{'absolute'}->{'absPct'}),
        ] );

        print $tb->title();
        print $tb->rule('=' , '+');
        foreach my $category ( sort keys %index ) {
		print $tb->body;
                print $tb->rule('-' , '+');
        }
}
=end comment
=cut

sub printOverlap{
        my $hash = shift;
        my %index;
        my %seen;
        my ($start, $end) = (0,0);
	my @rows = ();

        my @cols = ();
        push @cols,  +{title => "Ref", align => "center"}, +{title => "refStats", align => "center"};

	foreach my $refs (sort keys %$hash){
		my $refhash = $$hash{$refs};
		foreach my $keys (sort keys %$refhash){
			my $maphash = $$refhash{$keys}{'map'};
			foreach my $maps (sort keys %$maphash){
				my $mapelehash = $$maphash{$maps};
				foreach (sort keys %$mapelehash){
					push (@cols, +{title => $_, align => "center"} ) unless defined $seen{$_};
                        		$seen{$_} = 1;
				}
                        }
                }
	}

        my $tb = Text::Table->new( ( map { +( ( ref($_) ? $_ : "$_" ) , \'|' ) } @cols ) );
	foreach my $refs (sort keys %$hash){
		my $refhash = $$hash{$refs};
		foreach my $keys (sort keys %$refhash){
			my $refstats = "Count: $$refhash{$keys}{'ref'}{'count'}\ntotBps: $$refhash{$keys}{'ref'}{'totBps'}\ncvgBps: $$refhash{$keys}{'ref'}{'cvgBps'}\n cvgPct: $$refhash{$keys}{'ref'}{'cvgPct'}";
			$refstats .= "\n($$refhash{$keys}{'ref'}{'coords'})" if defined $$refhash{$keys}{'ref'}{'coords'};
			my @ref = ($keys, $refstats);
			my $maphash = $$refhash{$keys}{'map'};
			foreach my $maps (sort keys %$maphash){
				my $mapelehash = $$maphash{$maps};
				foreach (sort keys %seen){
					my $mapstats = "cvgBps: $$mapelehash{$_}{'cvgBps'}" if defined $$mapelehash{$_}{'cvgBps'};
					$mapstats .= "\n($$mapelehash{$_}{'coords'})" if defined $$mapelehash{$_}{'coords'};
					push (@ref, $mapstats ? $mapstats : '-');
				}
			}
			$tb->load(\@ref);
			push @rows, \@ref;
		}
	}

	my %indexes = (
                0 => 1, # Before header row
                1 => 1,  # After header row
        );

        # Calculate body rows for which rules should be printed
        my $count = 1;
        foreach my $row (@rows) {
                my $newlines = max map { tr/\n// } @$row;
                $count += $newlines + 1;
                $indexes{$count} = 1;
        }

        my $rule = $tb->rule('-', '+');
        foreach my $j (0 .. $tb->height) {
                print $rule if exists $indexes{$j};
                print $tb->table($j);
        }

=begin comment         
       foreach my $ref (keys %{ $$hash{$name} } ){
                        my ($mapBps, $mapPct);
                        my $totalref = $$hash{$name}->{$ref}->{'totalBps'};
                        my $overlapref = $$hash{$name}->{$ref}->{'overlapBps'};
                        my $pctref = $overlapref != 0 ? sprintf("%.2f", ( ($overlapref / $totalref) * 100) ) : 0;
                        my @ref = ($ref, $overlapref, $totalref, $pctref);
                        foreach (sort keys %seen){
                                $mapBps = defined $$hash{$name}->{$ref}->{$_} ? $$hash{$name}->{$ref}->{$_} : '-';
                                $mapPct = defined $$hash{$name}->{$ref}->{$_} ? sprintf("%.2f", (( $mapBps / $totalref ) * 100 )) : undef;
                                push (@ref, defined $mapPct ? "$mapBps ($mapPct%)" : '-');
                        }
                       #splice @ref, 1, 1;
                        $tb->load( \@ref );
                }
        print $tb->title();
        print $tb->rule('=' , '+');
        print $tb->body();
        }
}
=end comment
=cut
}
