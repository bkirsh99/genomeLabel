# genomeLabel
## User Guide
Create genomic labels to help identify candidate cis-regulatory regions.

Contents
======
- [What is genomeLabel?](#what-is-genomelabel)
  - [Background](#background)
- [Installation](#installation)
  - [Dependencies](#dependencies)
  - [Quick Start](#quick-start)
  - [genomeLabel](#genomelabel)
- [Running genomeLabel](#running-genomelabel)
- [Output of genomeLabel](#output-of-genomelabel)
  - [Example Directory Tree](#example-directory-tree)
  - [Example Data Tracks](#example-data-tracks)

What is genomeLabel?
======
**genomeLabel** is a command line tool for the fully automatic generation of genomic label data that can be summarized into statistically relevant information and visualized in a genome browser.
The labels annotate the genome both structurally and functionally, revealing insights into potential regulatory roles of certain regions.

genomeLabel is implemented in **perl** and automatically executes tools provided by bedtools and UCSC to annotate the Genome Reference Consortium Human Build 38 (**hg38**) assembly.  

For display by the UCSC Genome Browser, assembly hubs need to be hosted on a publicly accessible web server.

Background:
------
## Genome Annonation:
Annotating the genome can be broken down into three main steps: (1) identifying portions of the genome that do not code for proteins, (2) identifying elements on the genome, a process called *gene prediction*, and (3) attaching biological information to these elements. This yields two types of annotation:

**1) Structural Annotation -** Identifying elements in the genome (i.e. exons, introns, UTRs, CDSs, etc.)

**2) Functional Annotation -** Attaching biological information to the elements in the genome (e.g. biochemical product, regulatory role, expression, etc.)

Installation
======
Dependencies
------
genomeLabel contains different perl scripts that require bedtools and [UCSC utilities](http://hgdownload.soe.ucsc.edu/admin/exe).

Quick Start
------
## Install bedtools:
```bash
git clone https://github.com/arq5x/bedtools2.git
cd bedtools2
make clean && make all
```
## Clone repository:



## Data:
This program uses:
- NCBI RefSeq Curated Data (https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz)
- FANTOM5 Promoter Data (https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest/extra/CAGE_peaks/hg38_fair+new_CAGE_peaks_phase1and2.bed.gz)
- FANTOM5 Enhancer Data (https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest/extra/enhancer/F5.hg38.enhancers.bed.gz) - **Active, in-vivo transcribed enhancers**
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
Enhancer | Regulatory DNA sequence that, when bound by specific proteins (i.e transcription factors), enhance the transcription of an associated gene.

In order of priority,

**EXON** (**CDS** > **exon** > **UTR**) > **INTRON** > **INTERGENIC**

#### Examples:
![Image of Labelling Priority](https://docs.google.com/drawings/d/e/2PACX-1vQ51t4D1h96WMh588J429qSXkb_Fa6Cg_PhF3FHI4t2yPqMk1nzN0g54jFnf6wyD3hjs0qZS0brCaf3/pub?w=960&h=720)

# Structure:

``` bash
.
├── data
│   ├── hg38.chrom.sizes
│   ├── hg38.ncbiRefSeq.gtf.gz
│   ├── hg38_fair+new_CAGE_peaks_phase1and2.bed.gz
│   └── rmsk.txt.gz
├── hubs
│   └── myHub_chrX:15200000-15800000
│       ├── hg38
│       │   ├── codingTrack.bb
│       │   ├── genomicCodingTrack.bb
│       │   ├── genomicTrack.bb
│       │   ├── promoterTrack.bb
│       │   ├── repeatTrack.bb
│       │   └── trackDb.txt
│       ├── genomes.txt
│       └── hub.txt
├── out
│   └── chrX:15200000-15800000
│       ├── dataFiles
│       │   ├── hg38.chrom.sizes
│       │   ├── hg38.ncbiRefSeq.gtf.gz
│       │   ├── hg38_fair+new_CAGE_peaks_phase1and2.bed.gz
│       │   └── rmsk.txt.gz
│       ├── intersectFiles
│       │   ├── Alu.bed.gz
│       │   ├── LINE.bed.gz
│       │   └── SINE.bed.gz
│       ├── labelFiles
│       │   ├── Alu.bed.gz
│       │   ├── LINE.bed.gz
│       │   ├── SINE.bed.gz
│       │   ├── coding.bed.gz
│       │   ├── exon.bed.gz
│       │   ├── intergenic.bed.gz
│       │   ├── intron.bed.gz
│       │   ├── noncoding.bed.gz
│       │   └── promoter.bed.gz
│       ├── trackFiles
│       │   ├── codingTrack.bed
│       │   ├── genomicCodingTrack.bed
│       │   ├── genomicTrack.bed
│       │   ├── promoterTrack.bed
│       │   └── repeatTrack.bed
│       └── summary.txt
├── getStats.pl
├── makeHubs.pl
├── makeLabels.pl
├── makeTracks.pl
└── run.pl
```
- **Directory Content:*

# Scripts:
## 1) *run.pl*
- Command line argument(s):
  - Input region in the format **chr:start-stop**
- Output(s): 

Directory | Contents
:---: | :---:
data | Data used to generate labels and tracks for *all regions* (i.e. raw NCBI RefSeq, FANTOM5, RepeatMasker, and UCSC data)
out | Regional directories created for each *input region* in the format **chr:start-stop**
out/chr:start-stop | Regional subdirectories and summary file created for each *input region*
out/chr:start-stop/dataFiles | Data used to generate labels and tracks for *input region* (i.e. filtered NCBI RefSeq, FANTOM5, RepeatMasker, and UCSC data)
out/chr:start-stop/labelFiles | Labels created for each feature in the *input region* (i.e. exon, intron, coding, etc.)
out/chr:start-stop/trackFiles | Tracks created for different combinations of features in the *input region* (i.e. genomic, repetitive elements, etc.)
out/chr:start-stop/summary.txt | File created for summary statistics in the *input region* (e.g. coverage and average)


## 2) *makeLabels.pl*
- Command line argument(s): Input region in the format **chr:start-stop**

## 3) *getStats.pl*
- Command line argument(s): Input region in the format **chr:start-stop**

## 4) *makeTracks.pl*
- Command line argument(s): Input region in the format **chr:start-stop**

## 5) *makeHubs.pl*
- Command line argument(s):
  - Input region in the format **chr:start-stop**
  - Path to the UCSC utility **bedToBigBed** in the format **/path/bedToBigBed**. If you are unsure where this file lies, you can run the following command:
  
  ``` bash
  find / -type f -name "bedToBigBed" 2>/dev/null
  ```



*Inspired by:*
- **davetang/defining_genomic_regions:** https://github.com/davetang/defining_genomic_regions
- **cschlaffner/TrackHubGenerator:**  https://github.com/cschlaffner/TrackHubGenerator
