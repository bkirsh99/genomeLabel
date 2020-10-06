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

my $dir = join("@",'chrX:15200000-15800000','HUVEC/labels');#=shift;

my ($chr,$start,$end) = split /[:-]+/, substr($dir,0,index($dir,"@"));
#print Dumper $chr,$start,$end;

#my @options = shift;

##DEBUGGING:
#print Dumper $trackdir;
#print Dumper $regulator, $repeat, $felement;

my $stats = join('/',split($dir,0,index($dir,"/")),'stats.txt');
my %fileStats;
#append <chr:start-stop@biotype>/tracks/<filter>/ to file name - must be done before declaring %file in order to update file names
#my $var = adddir(LABELS);

#declare arrays for each label (to be filled with bedFeature objects)
#my (@exon, @transcript, @intron, @intergenic, @coding, @noncoding, @promoter, @repeat, @f_element, @cr_module);
#my %file = ( $var->{EXON}[0] => \@exon, $var->{TRANSCRIPT}[0] => \@transcript,
 #               $var->{INTRON}[0] => \@intron, $var->{INTERGENIC}[0] => \@intergenic,
 #               $var->{CODING}[0] => \@coding, $var->{NONCODING}[0] => \@noncoding,
 #               $var->{PROMOTER}[0] => \@promoter, $var->{REPEAT}[0] => \@repeat,
 #               $var->{FUNC_ELEMENT}[0] => \@f_element, $var->{CIS_REG_MODULE}[0] => \@cr_module, );

###############################

#1. Unzip zipped files;
#unzip();
#2. Compute statistics:
#a)

#tracklabel();
#genometrack();
test2();

sub test2{
	my @files = glob "$dir/*.bed";
	#print Dumper @files;
	for (0..$#files){
#		my $test = join("@",'chrX:15200000-15800000','HUVEC/labels/exon.bed');
#        	coverage($test);
		coverage($files[$_]);
		groupspecial($files[$_]) if ($files[$_] =~ /repeat|transcript|exon/);
#		print Dumper \%fileStats;
		#absoluteCvg($files[$_]);
        	#relativeCvg($files[$_]);
        }
}

sub coverage{
        my $file = shift;
	my $out = `echo "$chr\t$start\t$end\tfoo-1" | bedmap --echo-ref-row-id --echo-ref-size --bases-uniq --bases-uniq-f --echo-map-id --echo-overlap-size - $file`;
        my ($count, %eleStats);
	my $checkSum = 0;
	
        my @fields = split(/\|/, $out); #split the bedmap outputs
	my $size = $fields[1]; #total size (bps) of input region
	my $cvg = $fields[2]; #unique bps in input region covered by elements of $file
	my $pct = $fields[3] * 100; #percent coverage

	$eleStats{$file}->{'absolute'}->{'absBps'} = $cvg;
	$eleStats{$file}->{'absolute'}->{'absPct'} = $pct;

	my @element = split(";",$fields[4]) if defined $fields[4];
	my @elementcvg = split(";",$fields[5]) if defined $fields[5]; 
	foreach (@element){
		#my @subcategories = index($_,",") > 0 ? split(",", $_) : $_;
		#foreach my $subcategory (0..$#subcategories){
			my $category = index($_,".") > 0 ? substr( $_,0,index( $_,".")) : $_;
			my $subcategory = index($_,".") > 0 ? substr( $_,index( $_,".")+1) : $_;
			$eleStats{$category}->{'relative'}->{$subcategory}->{'relBps'} = defined $eleStats{$category}->{'relative'}->{$subcategory} ? ($eleStats{$category}->{'relative'}->{$subcategory}->{'relBps'} + shift @elementcvg) : shift @elementcvg;
			$eleStats{$category}->{'relative'}->{$subcategory}->{'relCount'}++;
		#}
	}
	foreach (keys %{ $eleStats{$file}->{'relative'} }){
		$eleStats{$file}->{'relative'}->{$_}->{'relPct'} = sprintf("%.2f", (($eleStats{$file}->{'relative'}->{$_}->{'relBps'} / $size) * 100));
	}
	
	$fileStats{$file} = $eleStats{$file};

}

sub groupspecial{
	my $file = shift;
	my (%h1, %h2, %h3, %seen1, %seen2, %seen3, @arr1, @arr2);

	foreach (keys %{ $fileStats{$file}->{'relative'} }){
		if($file =~ /repeat/){ #naming: repClass:repFamily:repName
			my ($class, $family, $name) = split(":",$_);
			push @{ $h2{$family} }, $name unless defined $seen3{$name};
			push @{ $h1{$class} }, $family unless defined $seen2{$family};
			++$seen3{$name};
                        ++$seen2{$family};
			++$seen1{$class};  
		}
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
	print Dumper $hash;
} 

=begin comment
sub test{
        my @files = glob "$out/*.bed";
        for (0..$#files){
        absoluteCvg($files[$_]);
        #relativeCvg($files[$_]);
        }
        #my $file2= join("/",$out,"cr_module.bed");
        #my $file = join("/",$out,"f_element.bed");

        #absoluteCvg($file);
        #relativeCvg($file);
        #overlapCvg($file, $file2);
        #check($file);
        #print Dumper %fileStats;
        #print Dumper %groupStats;
        printMain(\%fileStats);
        #printBreakdown(\%fileStats, 0); #0 is a flag for "individual" data
        #printOverlap(\%fileStats,\%groupStats, 1); #1 is a flag for "group" data
}

sub check{
        my $file = shift;
        my $sumCount = 0;
        my $sumBps = 0;
        my $sumPct = 0;
        foreach ( sort keys %{ $fileStats{$file}->{'relative'} } ){
                my $rel = $fileStats{$file}->{'relative'}->{$_};
                $sumCount += $rel->{'relCount'};
                $sumBps += $rel->{'relBps'};
                $sumPct += $rel->{'relPct'};
        }

        print Dumper $fileStats{$file}->{'absolute'}, $sumCount, $sumBps, $sumPct;

#       for (0..$#files){

#        }
}


sub majorElement{
        my $ref = shift;
        my $file = shift;

        foreach ( sort keys %{ $ref } ){
                if ($file eq $ref->{$_}[0]){
                        return $ref->{$_}[1];
                }
        }
}

sub absoluteCvg{
        my $file = shift;
        my %eleStats = ();
        $eleStats{$file}->{'absCount'} = `echo "$chr\t$start\t$end\tfoo-1" | bedtools coverage -a stdin -b $file | cut -f 5`;
        $eleStats{$file}->{'absBps'} = `echo "$chr\t$start\t$end\tfoo-1" | bedtools coverage -a stdin -b $file | cut -f 6`;
        $eleStats{$file}->{'absPct'} = `echo "$chr\t$start\t$end\tfoo-1" | bedtools coverage -a stdin -b $file | cut -f 8` * 100;

        $fileStats{$file}->{'absolute'} = $eleStats{$file};

}

sub relativeCvg{
        my $file = shift;
        my %eleStats = ();
        my $total = $fileStats{$file}->{'absolute'}->{'absBps'};

        open(FH, '<', $file) or die "Can't open '$file': $!";
        while(<FH>){
                my $string = $_;
                chomp $string;
                my @fields = split('\t', $string);
                my $start = $fields[1];
                my $end = $fields[2];
                my $size = $end - $start;
                my @id = ( $fields[3] =~ /(^{.*?}.*?)}/ );
                foreach (@id){
                        $_ .= "\}" unless substr($_, -1) eq "\}";
                        my $hash;
                        eval { $hash = JSON->new->allow_nonref->decode( $_ ) };
                        if ($@){
                                $_ .= "\}";
                        }
                        parse( $file, JSON->new->allow_nonref->decode( $_ ) );

                        ++$eleStats{$_}->{'relCount'};
                        $eleStats{$_}->{'relBps'} += $size;
                }
        }
        close(FH);

        foreach (sort keys %eleStats){
                $eleStats{$_}->{'relPct'} = ( $eleStats{$_}->{'relBps'} / $total ) * 100;
                $fileStats{$file}->{'relative'}->{$_} = $eleStats{$_};
        }

}

sub overlapCvg{
        my $ref = shift;
        my $map = shift;
        my $name = join("~",$ref,$map);
        my $out = `bedmap --echo-ref-row-id --echo-ref-size --bases-uniq --bases-uniq-f --echo-map-id-uniq --echo-overlap-size $ref $map`;
        my @lines = split "\n", $out; #split output into different lines
        my %overlapStats = ();

        foreach (@lines){
                my $string = $_;
                my @fields = split(/\|/, $string); #split the bedmap outputs (ref row-id, ref size, ref bps coverage, ref pct coverage, map elements (;), map bps (;)
                my ($refid, $refsize, $refbps, $refpct, $mapele, $mapbps) = @fields;
                $refid = substr($refid, 3); #remove leading "id-"
                $refid = pop @{ getOriginal ($ref, ( readLine($ref, $refid) ) ) }; #ref name (according to original database label)
                $mapele = getOriginal($map, $mapele); #map name(s) (according to original database label)
                my @mapbparray = split(";",$mapbps) if defined $mapbps;
                my %eleStats = ();
                $eleStats{$refid}->{'size'} = $refsize;
                $eleStats{$refid}->{'overlapBps'} = $refbps;
                $eleStats{$refid}->{'overlapPct'} = $refpct * 100;
                foreach ( @{ $mapele } ){
                        $eleStats{$refid}->{'overlapEle'}->{$_} = shift @mapbparray;
                }
                foreach (sort keys %eleStats){
                        push @{ $overlapStats{$name}->{$_} }, $eleStats{$_};
                }
        }
        group($overlapStats{$name}, $name);
}

sub group{
        my $hash = shift;
        my $name = shift;

        foreach my $ref (sort keys %$hash){
                my %refStats = ();
                my $size = 0;
                my $overlap = 0;
                foreach my $refhash ( @{ $$hash{$ref} } ){
                        my %mapStats = ();
                        $size = $size + $refhash->{'size'};
                        $overlap = $overlap + $refhash->{'overlapBps'};
                        foreach my $mapEle ( sort keys %{ $refhash->{'overlapEle'} } ){
                                $refStats{$ref}->{$mapEle} += $refhash->{'overlapEle'}->{$mapEle};
                        }
                }
                $refStats{$ref}->{'totalBps'} = $size;
                $refStats{$ref}->{'overlapBps'} = $overlap;
                foreach (sort keys %refStats){
                        $groupStats{$name}->{$_} = $refStats{$_};
                }
        }
}

sub parse{
        my $file = shift;
        my $original = shift;
        my $ele;

        if ($file =~ /transcript|exon|intron|coding/) {
                $ele = $original->{ATTRB}->{'gene_name'};
        #} elsif ($file =~ /repeat/) {
        } elsif ($file =~ /f_element/) {
                $ele = $original->{name};
        } elsif ($file =~ /cr_module/) {
                $ele = $original->{TF};
        }

        return defined $ele ? $ele : '';
}


sub getOriginal{
        my $file = shift;
        my @array = ();
        my @AoA = ();

        unless (@_){
                open(FH, '<', $file) or die "Can't open '$file': $!";
                while(<FH>){
                        my $string = $_;
                        chomp $string;
                        my @fields = split('\t', $string);
                        push @AoA, getOriginal($file, $fields[3]);
                }
                close(FH);
                return \@AoA;

        } else {

        my $field = shift;
        my @original = split(';', $field) if defined $field;
        foreach (@original){
                push @array, parse( $file, JSON->new->allow_nonref->decode($_) );
        }
        return \@array;
        }
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

sub printMain{
        my $hash = shift;                                                                                                                                                                                                                                                                                    my $flag = shift;
        my $size = $end - $start;
        my @cols = ();
        push @cols,  +{title => "Count", align => "center"}, +{title => "totalBps", align => "center"},
    +{title => "pctBps", align => "center"};

        my $tb        = Text::Table->new( "Category", \'||',
    ( map { +( ( ref($_) ? $_ : "$_" ) , \'|' ) } @cols ) );

        foreach my $category ( sort keys %{ $hash } ) {
        $tb->load( [
                   majorElement(LABELS, $category),
                   $$hash{$category}->{'absolute'}->{'absCount'},
                   $$hash{$category}->{'absolute'}->{'absBps'},
                   sprintf("%.2f", $$hash{$category}->{'absolute'}->{'absPct'}),
        ] );
        }
        print "Region: $region (totalBps: $size bps)\n";
        print $tb->title();
        print $tb->rule('=' , '+');
        print $tb->body();
}

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
                foreach my $subcategory (sort keys %{ $$hash{$category}->{'relative'} } ){
                ++$end;
        $tb->load( [
                   $subcategory,
                   $$hash{$category}->{'relative'}->{$subcategory}->{'relCount'},
                   $$hash{$category}->{'relative'}->{$subcategory}->{'relBps'},
                   sprintf("%.2f", $$hash{$category}->{'relative'}->{$subcategory}->{'relPct'}),
        ] )
                }
        $index{$category}->{'start'} = $start;
        $index{$category}->{'end'} = $end;
        }

        print $tb->title();
        print $tb->rule('=' , '+');
        foreach my $category ( sort keys %index ) {
                print $tb->body($index{$category}->{'start'});
                print $tb->rule('-' , '+');
                print $tb->body($index{$category}->{'start'} + 1, $index{$category}->{'end'});
        }
}

sub printOverlap{
        my $href = shift;
        my $hash = shift;
        my $flag = shift;
        my %index;
        my %seen;
        my ($start, $end) = (0,0);


        my @cols = ();
        push @cols,  +{title => "Ref/Map", align => "center"}, +{title => "totalBps", align => "center"},
             +{title => "overlapBps", align => "center"}, +{title => "pctCoverage", align => "center"};

        foreach my $name ( sort keys %{ $hash } ){
                while (my ($ref, $map) = each %{ $$hash{$name} }) {
                        foreach my $mapele (sort keys %{ $map } ){
                                push (@cols, +{title => $mapele, align => "center"} ) unless defined $seen{$mapele} or $mapele =~ /Bps/;
                                $seen{$mapele} = 1;
                        }
                }

        my $tb = Text::Table->new( ( map { +( ( ref($_) ? $_ : "$_" ) , \'|' ) } @cols ) );
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

#=begin comment
sub printOverlap{
        my $href = shift;
        my $hash = shift;
        my $flag = shift;
        my %index;
        my %seen;
        my ($start, $end) = (0,0);


#print Dumper $href, $hash;
        my @cols = ();
        push @cols,  +{title => "Ref/Map", align => "center"};

        foreach my $name ( sort keys %{ $hash } ){
                while (my ($ref, $map) = each %{ $$hash{$name} }) {
                        foreach my $mapele (sort keys %{ $map } ){
                                push (@cols, +{title => $mapele, align => "center"} ) unless defined $seen{$mapele};#grep{$_ eq $mapele} @seen;
                                $seen{$mapele} = 1;
                        }
                }

        my $tb = Text::Table->new( ( map { +( ( ref($_) ? $_ : "$_" ) , \'|' ) } @cols ) );
                foreach my $ref (keys %{ $$hash{$name} } ){
                        my @files = split("~", $name);
                        my $refName = $files[0];
                        my $refBps = $$href{$refName}->{'relative'}->{$ref}->{'relBps'};
                        my ($mapBps, $mapPct);
                        #my @ref = ($ref);
                        my @ref = ("$ref ($refBps)"); #if we want total bps of ref
                        foreach (sort keys %seen){
                                $mapBps = defined $$hash{$name}->{$ref}->{$_} ? $$hash{$name}->{$ref}->{$_} : '-';
                                $mapPct = defined $$hash{$name}->{$ref}->{$_} ? sprintf("%.2f", (( $mapBps / $refBps ) * 100 )) : undef;
                                push (@ref, defined $mapPct ? "$mapBps ($mapPct%)" : '-');
                        }
                        splice @ref, 1, 1;
                        $tb->load( \@ref );
                }
        print $tb->title();
        print $tb->rule('=' , '+');
        print $tb->body();
        }
}
#=end comment
#=cut
sub printPretty1{
        my $hash = shift;
        my $flag = shift;

        if ($flag == 0){
        print "Category                  Count        bpTotal     bpPct\n";
        print "========                  =====        =======     =====\n";

                foreach my $category ( sort keys %{ $hash } ) {
        print ""
        . sprintf(
                   "%-20s   %-10d   %-10d   %0.2f",
                   majorElement(LABELS, $category),
                   $$hash{$category}->{'absolute'}->{'absCount'},
                   $$hash{$category}->{'absolute'}->{'absBps'},
                   $$hash{$category}->{'absolute'}->{'absPct'},
        )
        . "% \n";

        print "-----------------------------------------------------------\n";
        print "Subcategory                  Count        bpTotal     bpPct\n";
        print "===========                  =====        =======     =====\n";

                        foreach my $subcategory (sort keys %{ $$hash{$category}->{'relative'} } ){

        print ""
        . sprintf(
                   "%-20s   %-10d   %-10d   %0.2f",
                   $subcategory,,
                   $$hash{$category}->{'relative'}->{$subcategory}->{'relCount'},
                   $$hash{$category}->{'relative'}->{$subcategory}->{'relBps'},
                   $$hash{$category}->{'relative'}->{$subcategory}->{'relPct'},
        )
        . "% \n";

                        }
                }


        } else {
                foreach my $name ( sort keys %{ $hash } ){
        print "Reference                  Map        bpTotal\n";
        print "=========                  ===        =======\n";
                        while (my ($ref, $map) = each %{ $$hash{$name} }) {
        print ""
        . sprintf(
                   "%-20s",
                   $ref,
        )
        . "\n";
                                foreach my $mapele (sort keys %{ $map } ){
        print "\t"
        . sprintf(
                   "%-20s   \t%-10d",
                   $mapele,
                   $$hash{$name}->{$ref}->{$mapele},
        )
        . "\n";
                                }
                        }
        print "-----------------------------------------------------------\n";
                }
        }

}


#append <chr:start-stop@biotype>/tracks/ to file names
sub adddir{
        my $ref = shift;
        my $new;

        while (my ($key, $value) = each ( %{ $ref } )){
                $new->{$key}[0] = join('/', $filterdir, $ref->{$key}[0]);
        }
        return $new;
}

sub majorElement{
        my $ref = shift;
        my $file = shift;

        foreach ( sort keys %{ $ref } ){
                if ($file eq $ref->{$_}[0]){
                        return $_;
                }
        }
}

=end comment
=cut
