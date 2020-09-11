# genomeLabel
## Data Overview
Overview of the contents and structure of each data set. 

Contents
======
- [NCBI RefSeq - "hg38.ncbiRefSeq.gtf.gz"](#ncbi-refseq)
- [FANTOM5 Enhancers - "F5.hg38.enhancers.bed.gz"](#fantom5)
- [FANTOM5 Promoters - "hg38_fair+new_CAGE_peaks_phase1and2.bed.gz"](#fantom5)
- [RepeatMasker - "rmsk.txt.gz"](#repeatmasker)
- [Segway Encyclopedia - "segway_encyclopedia.bed.gz"](#segway)
- [Segway - ""](#segway)
- [ReMap 2020 CRMs - "remap2020_crm_macs2_hg38_v1_0.bed.gz"](#remap)
- [ReMap 2020 - ""](#remap)
- [UCSC Chromosome Sizes - "hg38.chrom.sizes"](#ucsc)
- [UCSC liftOver - "hg19ToHg38.over.chain.gz"](#ucsc)

NCBI RefSeq
======
```bash
chr1    ncbiRefSeq      transcript      11874   14409   .       +       .       gene_id "DDX11L1"; transcript_id "NR_046018.2";  gene_name "DDX11L1";
chr1    ncbiRefSeq      exon    11874   12227   .       +       .       gene_id "DDX11L1"; transcript_id "NR_046018.2"; exon_number "1"; exon_id "NR_046018.2.1"; gene_name "DDX11L1";
chr1    ncbiRefSeq      exon    12613   12721   .       +       .       gene_id "DDX11L1"; transcript_id "NR_046018.2"; exon_number "2"; exon_id "NR_046018.2.2"; gene_name "DDX11L1";
chr1    ncbiRefSeq      exon    13221   14409   .       +       .       gene_id "DDX11L1"; transcript_id "NR_046018.2"; exon_number "3"; exon_id "NR_046018.2.3"; gene_name "DDX11L1";
chr1    ncbiRefSeq      transcript      14362   29370   .       -       .       gene_id "WASH7P"; transcript_id "NR_024540.1";  gene_name "WASH7P";
chr1    ncbiRefSeq      exon    14362   14829   .       -       .       gene_id "WASH7P"; transcript_id "NR_024540.1"; exon_number "1"; exon_id "NR_024540.1.1"; gene_name "WASH7P";
chr1    ncbiRefSeq      exon    14970   15038   .       -       .       gene_id "WASH7P"; transcript_id "NR_024540.1"; exon_number "2"; exon_id "NR_024540.1.2"; gene_name "WASH7P";
chr1    ncbiRefSeq      exon    15796   15947   .       -       .       gene_id "WASH7P"; transcript_id "NR_024540.1"; exon_number "3"; exon_id "NR_024540.1.3"; gene_name "WASH7P";
chr1    ncbiRefSeq      exon    16607   16765   .       -       .       gene_id "WASH7P"; transcript_id "NR_024540.1"; exon_number "4"; exon_id "NR_024540.1.4"; gene_name "WASH7P";
chr1    ncbiRefSeq      exon    16858   17055   .       -       .       gene_id "WASH7P"; transcript_id "NR_024540.1"; exon_number "5"; exon_id "NR_024540.1.5"; gene_name "WASH7P";
```

FANTOM5
======
Enhancers
------
```bash
chr10   100006233       100006603       chr10:100006233-100006603       35      .       100006509       100006510      0,0,0                                                       2218,34  0,336
chr10   100008181       100008444       chr10:100008181-100008444       101     .       100008274       100008275      0,0,0                                                       232,108  0,155
chr10   100014348       100014634       chr10:100014348-100014634       3       .       100014546       100014547      0,0,0                                                       2112,1   0,285
chr10   100020065       100020562       chr10:100020065-100020562       211     .       100020344       100020345      0,0,0                                                       2210,148 0,349
chr10   100043485       100043744       chr10:100043485-100043744       3       .       100043656       100043657      0,0,0                                                       286,1    0,258
chr10   100114218       100114567       chr10:100114218-100114567       51      .       100114352       100114353      0,0,0                                                       21,80    0,269
chr10   100148595       100148922       chr10:100148595-100148922       3       .       100148772       100148773      0,0,0                                                       229,1    0,326
chr10   100182422       100182522       chr10:100182422-100182522       3       .       100182466       100182467      0,0,0                                                       21,11    0,89
chr10   100184498       100184704       chr10:100184498-100184704       10      .       100184589       100184590      0,0,0                                                       229,51   0,155
chr10   100184852       100185124       chr10:100184852-100185124       2       .       100185033       100185034      0,0,0                                                       2158,66  0,206
```
