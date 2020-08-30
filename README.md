# genomeLabel
Create genomic labels to help identify candidate cis-regulatory regions.

## Genome Annonation:
Annotating the genome can be broken down into two steps:

**1) Structural Annotation -** Identifying elements in the genome (i.e. exons, introns, UTRs, CDSs, etc.)

**2) Functional Annotation -** Attaching biological information to the elements in the genome (e.g. biochemical product, regulatory role, expression, etc.)

# Data:
This program uses:
- NCBI RefSeq Curated Data (https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz)
- FANTOM5 Promoter Data (https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest/extra/CAGE_peaks/hg38_fair+new_CAGE_peaks_phase1and2.bed.gz)
- RepeatMasker Data (https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.sql)
- UCSC Chromosome Sizes Data (https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes)

This program generates the following output:
## a) Labels:

Label | Definition
:---: | :---:
Exon | Intragenic stretch of DNA sequence, including non-coding untranslated regions, that can synthesize a functional RNA molecule, including mRNAs and ncRNAs.
Intron |  Intragenic stretch of non-coding DNA sequence located between two successive exons.
Intergenic |  Stretch of non-coding DNA sequence located between the two successive genes.
Coding | Exonic stretch of DNA sequence that can synthesize a functional protein.
Noncoding | Exonic stretch of DNA sequence that appears twofold: *(1)* as regulatory untranslated regions of mRNAs or *(2)* as translated regions of ncRNAs.
Promoter | CAGE-defined core promoter regions.
LINE | Non-long terminal repeat (non-LTR) retrotransposon that contains an RNA polymerase II promoter.
SINE | Non-long terminal repeat (non-LTR) retrotransposon that does not encode a functional reverse transcriptase protein and relies on other mobile transposons, especially LINEs.
Alu | Most common SINE element, which is highly conserved and often implicated in disease.


# Scripts:
## 1) *makeLabels.pl*
- Command line argument(s): Input region in the format **chr:start-stop**

## 2) *getStats.pl*
- Command line argument(s): Input region in the format **chr:start-stop**

## 3) *makeTracks.pl*
- Command line argument(s): Input region in the format **chr:start-stop**

## 4) *makeHubs.pl*
- Command line argument(s):
  - Input region in the format **chr:start-stop**
  - Path to the UCSC utility **bedToBigBed** in the format **/path/bedToBigBed**. If you are unsure where this file lies, you can run the following command:
  
  ``` bash
  find / -type f -name "bedToBigBed" 2>/dev/null
  ```



*Inspired by:*
- **davetang/defining_genomic_regions:** https://github.com/davetang/defining_genomic_regions
- **cschlaffner/TrackHubGenerator:**  https://github.com/cschlaffner/TrackHubGenerator
