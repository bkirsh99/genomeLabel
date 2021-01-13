## Usage
This document outlines how to call the scripts created for genomeLabel. Previous versions of the final scripts are also available and can be identified by an underscore.  

|Script|Sample Call|
|:---:|:---:|
|rename_v1.py|<pre>wget http://fantom.gsc.riken.jp/5/datafiles/latest/extra/Enhancers/Human.sample_name2library_id.txt<br>cat Human.sample_name2library_id.txt \| sed -e "s/, */-/g" \| sed -e "s/ /_/g" \| sed -e "s/(/-/g" \| sed -e "s/)/-/g" \| sed -e "s/:/-/g" \| sed -e "s/'//g" \| sed -e "s/\^//g" \| sed -e "s/\///g" >  Human.sample_name2library_id.sanitized.txt<br>wget https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest/extra/CAGE_peaks_expression/hg38_fair+new_CAGE_peaks_phase1and2_counts_ann.osc.txt.gz<br>zcat hg38_fair+new_CAGE_peaks_phase1and2_counts_ann.osc.txt.gz	> hg38_fair+new_CAGE_peaks_phase1and2_counts_ann.osc.txt<br>wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz<br>mkdir split<br>ulimit -Sn 16384<br>python ./rename_v1.py Human.sample_name2library_id.sanitized.txt hg38_fair+new_CAGE_peaks_phase1and2_counts_ann.osc.txt split /home/liftOver|
|rename_v2.py|python ./rename_v2.py Human.sample_name2library_id.sanitized.txt hg38_fair+new_CAGE_peaks_phase1and2_counts_ann.osc.txt split /home/liftOver filter.txt|
|rename_v3.py|python ./rename_v3.py filter.txt split /home/liftOver|
|rename_v4.py|python ./rename_v3.py filter.txt split|
|plot_v1.py|./plot_v1.py -i "giggle_results/*" -o plot.pdf --stat odds|
|plot_v2.py|./plot_v2.py -i "giggle_results/*" -o plot.pdf --stat odds|
|plot_v3.py|./plot_v3.py -i "giggle_results/*" -o plot.pdf -m names.txt --stat odds|
|verify.py|<pre>mkdir split_sort<br>bash $GIGGLE_ROOT/scripts/sort_bed "split/*.bed" split_sort/ > /dev/null<br>./verify.py -i "split_sort/*" -q "giggle_results/*" -o verify.txt -m names.txt|

