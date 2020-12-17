use strict;
use warnings;
use Data::Dumper;

sub gff2genemodel{

	my $ingff=shift;
	open GFF, "<$ingff" or die "Cannot open $ingff";
	
	my %gene=();
	my %mrna=();
	my %transcript=();
	my %exon=();
	my %intron=();
	my %cds=();
	my %start_codon=();
	my %stop_codon=();
	my %utr_3=();
	my %utr_5=();

	my %count_type=("transcript"=>0,"exon"=>0,"intron"=>0,"CDS"=>0,"start_codon"=>0,"stop_codon"=>0,"3UTR"=>0,"5UTR"=>0);

	my %gff_feature=("transcript"=>[\%transcript,0],"exon"=>[\%exon,0],"intron"=>[\%intron,0],
		"CDS"=>[\%cds,0],"start_codon"=>[\%start_codon,0],"stop_codon"=>[\%stop_codon,0],
		"3UTR"=>[\%utr_3,0],"5UTR"=>[\%utr_5,0]);
	my $dup=0;

	my @lines=<GFF>;
	foreach my $line(@lines){
		if(($line=~/^#/g)|(length($line) == 1)){
			next;
		}

		chomp $line;

		my @col=split /\t/,$line;

		if(!exists $count_type{$col[2]}){
			next;
		}

		my @attributes=split /;/,$col[8];
		my %att;
		foreach my $attribute(@attributes){
			$attribute =~ s/^\s+|\s+$//g;
			my @kv=split /\s+/,$attribute;
			$att{$kv[0]}=$kv[1];
		}

		if(!$att{ID}){
			my $id;
			if($col[2] eq 'transcript'){
				$id=$att{transcript_id};
			} elsif($col[2] eq 'gene'){
				$id=$att{gene_id};
			} else {
				$id=$col[2].$count_type{$col[2]};
			}
			$att{ID}=$id;	
		}

		foreach(keys %gff_feature){
			if($col[2] eq $_){
				my $hash=$gff_feature{$_}[0];
				$hash->{$att{ID}}->{type}=$col[2];
				$hash->{$att{ID}}->{seqid}=$col[0];
				$hash->{$att{ID}}->{source}=$col[1];
				$hash->{$att{ID}}->{start}=$col[3];
				$hash->{$att{ID}}->{end}=$col[4];
				$hash->{$att{ID}}->{strand}=$col[6];
				$hash->{$att{ID}}->{score}=$col[5];
				$hash->{$att{ID}}->{phase}=$col[7];
				$hash->{$att{ID}}->{attributes}=\%att;
				$count_type{$col[2]}++;
			}
		}

	}

	foreach my $feature(keys %gff_feature){
		unless($feature eq 'transcript'){
			my $hash=$gff_feature{$feature}[0];
			foreach my $id(keys %{$hash}){
				my $parent = $hash->{$id}->{attributes}->{transcript_id};
				if(exists $transcript{$parent}){
					$transcript{$parent}->{$feature}->{$id}=$hash->{$id};
				} else {
					print STDERR "Warning: Parent($parent) for $id not found. Skipped.";
				}
			}
		}
	}

	foreach my $id(keys %transcript){
		my $parent=$transcript{$id}->{attributes}->{gene_id};
			$gene{$parent}->{transcript}->{$id}=$transcript{$id};
	}
	
	close GFF;
	return \%gene;
}


sub printGene{ # gene based coordinates
	my $gene=shift;
	my $diff=shift;
	my @valid_childs=("transcript","exon","intron","CDS","start_codon","stop_codon","3UTR","5UTR");

	if(defined $gene->{transcript}){
		foreach my $mrna(keys %{$gene->{transcript}}){
			my $hash = $gene->{transcript}{$mrna};
			printGene($hash,$diff);
		}
	} else {
		my $att=$gene->{attributes};
		my $attributes="";
		foreach my $k(keys %$att){
			$attributes=$attributes."; ".$k." ".$$att{$k} unless $k eq "ID";
		}
		$attributes=~s/^;//;
		$attributes .= ";";

		if(!$diff){
			$diff=1;
		}
		my $gstart=$gene->{start}-$diff;
		my $gend=$gene->{end}-$diff;

		my $line=join("\t",$gene->{seqid},$gene->{source},$gene->{type},$gstart,$gend,$gene->{score},$gene->{strand},$gene->{phase},$attributes);
		print $line."\n";

		foreach my $k(keys %{$gene}){
			foreach my $valid(@valid_childs){
				if($k eq $valid){
					my $childs=$gene->{$k};
					foreach my $child_id(keys %$childs){
						my $child=$childs->{$child_id};
						printGene($child,$diff);
					}
				}
			}
		}
	}
}


# geneid2mrnaid
# mrnaid2geneid
# return a hash where key is mRNA ID y value is the corresponding gene ID.
sub mrnaid2geneid{
	my $genes=shift;
	my %eq;
	foreach my $gene(keys %$genes){
		my $mrnas=$genes->{$gene}->{mRNA};
		foreach my $mrna(keys %$mrnas){
			$eq{$mrna}=$gene;
			
		}
	}
	return \%eq;
}



# filter transcripts on a gene structure. if returns 0, no transcript associated to the gene
# options are, "longest", that keep only the longest transcript; and "shortest", that keep only shortest transcript. if transcripts are equally length, keep just one. the last that it finds
# usage: $filtered_gene_structure=filterGeneTranscripts($gene_structure, 'longest')
sub filterGeneTranscripts{
	my $id=shift;
	my $gene=shift;
	my $opt=shift; # puede ser longest o shortest

	if(!exists $gene->{transcript}){
		print STDERR "Gene $id not contain transcript childs. Skipped\n";
		return 0;
	}

	my $mrnas=$gene->{transcript};

	my $last_len=0;
	my $valid_id;
	foreach my $id(keys %$mrnas){
		my $len=$$mrnas{$id}->{end} - $$mrnas{$id}->{start};
		if($opt eq 'longest'){
			if($len>=$last_len){
				$valid_id=$id;
			}
			else{
				print STDERR "Warning: $id length ($len) < 0?\n";
			}
		}
		if($opt eq 'shortest'){
			if($last_len<=0){
				$last_len=$len;
				$valid_id=$id;
			}
			else{
				if($len<$last_len){
					$last_len=$len;
					$valid_id=$id;
				}
			}
		}
		if(($opt ne 'longest')&&($opt ne 'shortest')){
			die "$opt is not valid option for filterGeneTranscripts";
		}
	}
	my $filtered_gene=$gene;
	my $valid_mrna=$gene->{transcript}->{$valid_id};
	$gene->{transcript}=();
	$gene->{transcript}->{$valid_id}=$valid_mrna;

	return $gene;
}

return 1;
