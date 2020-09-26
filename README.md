# genomeLabel
## User Guide
Create genomic labels to help identify candidate cis-regulatory regions.

Contents
======
- [What is genomeLabel?](#what-is-genomelabel)
  - [Background](#background)
  - [Data](#data)
  - [Usage](#usage)
- [Installation](#installation)
  - [Dependencies](#dependencies)
  - [Quick Start](#quick-start)
  - [genomeLabel](#genomelabel)
- [Running genomeLabel](#running-genomelabel)
  - [run.pl](#run.pl)
  - [makeHubs.pl](#makehubs.pl)
  - [getStats.pl](#getstats.pl)
- [Output of genomeLabel](#output-of-genomelabel)
  - [Example Directory Tree](#example-directory-tree)
    - [Directory Contents](#directory-contents)
  - [Example Data Tracks](#example-data-tracks)

What is genomeLabel?
======
**genomeLabel** is a command line tool for the fully automatic generation of genomic label data that can be summarized into statistically relevant information and visualized in a genome browser.
The labels annotate the genome both structurally and functionally, providing insight into the potentially regulatory role of certain regions.

The genomeLabel tool is implemented in **perl** and automatically executes commands provided by bedtools and UCSC to annotate the Genome Reference Consortium Human Build 38 **(hg38)** assembly. For that, the **liftOver** UCSC utility must be downloaded to convert Segway's hg19 data.

For display by the UCSC Genome Browser, a script 'makeHubs.pl' can be used. Hoewever, assembly hubs need to be hosted on a publicly accessible web server provided by the user.

Background:
------
## Genome Annonation:
Annotating the genome can be broken down into three main steps: (1) identifying portions of the genome that do not code for proteins, (2) identifying elements on the genome, a process called *gene prediction*, and (3) attaching biological information to these elements. This yields two types of annotation:

**1) Structural Annotation -** Identifying elements in the genome (i.e. exons, introns, UTRs, CDSs, etc.)

**2) Functional Annotation -** Attaching biological information to the elements in the genome (e.g. biochemical product, regulatory role, expression, etc.)

Data:
------
This program takes in annotation data from genomic databases as input to generate custom label tracks.
  
|INPUT | OUTPUT|
|:--: | :--:|
|[NCBI RefSeq genes, curated subset](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz)|exon, intron, coding, noncoding, intergenic|
|[UCSC hg38.chrom.sizes](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes)|intergenic|
|[FANTOM5 TSS (CAGE) peaks](https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest/extra/CAGE_peaks/hg38_fair+new_CAGE_peaks_phase1and2.bed.gz)|promoter|
|[RepeatMasker, soft-masked](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.sql)|repetitive element|
|[Segway genomic "states"](https://noble.gs.washington.edu/proj/encyclopedia/segway_encyclopedia.bed.gz)|functional element|
|[UCSC hg19ToHg38.over.chain.gz](http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz)|functional element|
|[ReMap ChIP-seq peaks](http://remap.univ-amu.fr/storage/remap2020/hg38/MACS2/remap2020_crm_macs2_hg38_v1_0.bed.gz)|cis-regulatory module|

For an overview of the contents and format of each data set, please read up on additional [data](https://github.com/bkirsh99/genomeLabel/blob/master/data.md) information.

Usage:
------
```javascript
USAGE: ./run.pl [chr:start-stop] --biotype [BIOTYPE] --path [/PATH/] <command(s)> <filter(s)>
     
     The program accepts any of the following six human cell types used by the ENCODE Consortium: GM12878, H1-hESC, K562, HeLa-S3, HepG2, and HUVEC.
     The minimum requirements to run the script are the region, biotype, and path to liftOver.
     
     If <command(s)> is omitted, the default behaviour of the program is to generate a set of "label files." These can be useful to extract additional statistics or create more complex tracks according to user needs.
     
     For simple applications, <command(s)> can be used to output select "track files," track hubs, and summary statistics. 
     Filter arguments can be passed alongside different commands for added specificity. 
          

    Commands: (Any number of commands can be used at a single time.)
     --------
      
      --makeLabels          Make "label files" (exon.bed, intron.bed, intergenic.bed, transcript.bed, coding.bed, noncoding.bed, promoter.bed, repeat.bed, 
                              f_element.bed,  and cr_module.bed. Default operation.
      --makeTracks          Make "track files" (genomic.bed, promoter.bed, repeat.bed, f_element.bed, and cr_module.bed). 
      --makeHub             Make a tack data hub that can be imported into the UCSC Genome Browser.
      --getStats            Extract summary statistics.

      --help                Print this message and exit successfully.
    
    
    Filters:
     --------
    
    Filters for makeLabels and makeTracks specify the subclasses of features to be included in label and track files. If both commands are executed in the same call, all filters will apply to each set of output files.
    
      --repeat                 A list of one or more repeats - must follow the naming conventions adopted by RepeatMasker - may be a repName, repClass, repFamily, 
                                 or a combination of those.
      --regulator              A list of one or more transcriptional regulators - must follow the naming conventions adopted by ReMap - may be a transcription factor
                                 (TF), transcription co-activator (TCF), chromatin-remodeling factor (CRF), or a combination of those.
      
      Filters for getStats specify the relativity of statistics to be computed.

      --genes                     A list of the genes and their respective transcripts contained in the input region.
      --coverage-absolute         The fraction of distinct bases from the input region covered by each major feature (exon, intron, intergenic, promoter, repetitive
                                    element, functional element, and cis-regulatory module). 
      --coverage-absolute-ex      The fraction of distinct bases from the input region covered by coding and noncoding exons.
      --coverage-absolute-fe      The fraction of distinct bases from the input region covered by
    
```

|Label | Definition|
|:---: | :---:|
|![#00ff00](https://via.placeholder.com/15/00ff00/000000?text=+) Exon | Intragenic stretch of DNA sequence, including non-coding untranslated regions, that can synthesize a functional RNA molecule, including mRNAs and ncRNAs.|
|![#ff0000](https://via.placeholder.com/15/ff0000/000000?text=+) Intron |  Intragenic stretch of non-coding DNA sequence located between two successive exons.|
|![#0000ff](https://via.placeholder.com/15/0000ff/000000?text=+) Intergenic |  Stretch of non-coding DNA sequence located between the two successive genes.|
|![#006600](https://via.placeholder.com/15/006600/000000?text=+) Coding | Exonic stretch of DNA sequence that can synthesize a functional protein.|
|![#99ff99](https://via.placeholder.com/15/99ff99/000000?text=+) Noncoding | Exonic stretch of DNA sequence that appears twofold: *(1)* as regulatory untranslated regions of mRNAs or *(2)* as translated regions of ncRNAs.|
|![#7f00ff](https://via.placeholder.com/15/7f00ff/000000?text=+) Promoter | CAGE-defined core promoter regions.|
|![#ffff00](https://via.placeholder.com/15/ffff00/000000?text=+) LINE | Non-long terminal repeat (non-LTR) retrotransposon that contains an RNA polymerase II promoter.|
|![#ffcc99](https://via.placeholder.com/15/ffcc99/000000?text=+) SINE | Non-long terminal repeat (non-LTR) retrotransposon that does not encode a functional reverse transcriptase protein and relies on other mobile transposons, especially LINEs.|
|![#cc6600](https://via.placeholder.com/15/cc6600/000000?text=+) Alu | Most common SINE element, which is highly conserved and often implicated in disease.|
|![#7f007f](https://via.placeholder.com/15/7f007f/000000?text=+) Enhancer | Regulatory DNA sequence that, when bound by specific proteins (i.e transcription factors), enhance the transcription of an associated gene.|
|![#66b2ff](https://via.placeholder.com/15/66b2ff/000000?text=+) Functional Element | Putative transcriptional and regulatory DNA sequences identified by [Segway](https://www.biorxiv.org/content/10.1101/086025v3.full), which reflect DNA binding and accessibility  across different cell types (e.g. quiescent, enhancer, promoter, bivalent, etc.).|
|![#ff007f](https://via.placeholder.com/15/ff007f/000000?text=+) Cis-Regulatory Module (CRM) | Putative binding sites for transcriptional regulators identified by [ReMap](https://academic.oup.com/nar/article/48/D1/D180/5608991), which reflect the intersection of overlapping ChIP-seq peaks across different datasets of TRs and biotypes.|

Installation
======
Dependencies:
------
genomeLabel contains different perl scripts that require bedtools, liftOver, and bedToBigBed. The latter two are made available by [UCSC utilities](http://hgdownload.soe.ucsc.edu/admin/exe).

Quick Start:
------
## Install perl and required perl modules
```bash
sudo apt-get install perl
cpan File::Which
cpan File::Basename
cpan Tie::IxHash
```
## Install bedtools
```bash
git clone https://github.com/arq5x/bedtools2.git
cd bedtools2
make clean && make all
```
## Clone repository
```bash
git clone https://github.com/bkirsh99/genomeLabel.git
```
## Navigate to the genomeLabel directory
```bash
cd ${genomeLabel_DIR}
```
genomeLabel:
------
To initialize the program with all the necessary data files, *run.pl* must be the first script to be run. After that, standalone scripts can be used for either the same or different input regions.
To run any script from the genomeLabel directory, simply enter:
```bash
chmod +x ./script.pl
./script.pl [ARGS]
```

Running genomeLabel
======
This program contains 5 perl scripts and must be initialized by running **run.pl**. 
- If a script requires an input region as a command line argument, it must be provided in the format **chr:start-stop**.
- If a script requires the path to a UCSC binary utility (i.e. **liftOver** and **bedToBigBed**), it must be provided in the format **/path/to/utility/** (i.e. enclosed by "/" and excluding the name of the utility itself). If you are unsure where this file lies, you can run the following command:

  ``` bash
  find / -type f -name "bedToBigBed" 2>/dev/null
  ```

|Script|Argument(s)|Output(s)|Requirement(s)|
|:---:|:---:|:---:|:---:|
|*run.pl*|input_region|At minimum all files in **./data** and **./out/chr:start-stop/dataFiles**|NA|
|*makeLabels.pl*|input_region|All files in **./out/chr:start-stop/labelFiles**|*run.pl* to download and filter basic input data|
|*makeTracks.pl*|input_region|All files in **./out/chr:start-stop/trackFiles**|*makeLabels.pl* to combine labels into tracks|
|*getStats.pl*|input_region|All files in **./out/chr:start-stop/intersectFiles**, as well as **./out/chr:start-stop/summary.txt**|*makeTracks.pl* to extract relevant track metrics|
|*makeHubs.pl*|input_region UCSC_path email|All files in **./hubs/myHub_chr:start-stop**|*makeTracks.pl* to compile all tracks into a hub|

run.pl
------
This script can make calls to all other scripts and generate labels, tracks, summary statistics, and data hubs. It prompts the user for each of these subroutines with a (Y/N) question. 
It can be the only script used every time the program is run, or the user may choose to run subroutines individually depending on their desired output. Either way, it is the first script that must be ran for **each new input region** to initialize the genomeLabel tool with appropriately filtered data.

makeHubs.pl
------
This script produces a directory named after the input region in the format **myHub_chr:start-stop**. It follows UCSC's [Public Hub Guidelines](http://genomewiki.ucsc.edu/index.php/Public_Hub_Guidelines) by generating the following:
```bash
├── myHub_chrX:15200000-15800000 - directory containing track hub files
    ├── hg38 - directory containing hg38 (GRCh38) human assembly data
    |   ├── [...].bb - tracks in this directory
    |   └── trackDb.txt - display properties of tracks in this directory
    ├── genomes.txt - list of genome assemblies included in the hub (i.e. hg38 only)
    └── hub.txt - short description of hub properties
```
To use the output of *makeHubs.pl* in the UCSC Genome Browser, copy the complete hub folder (e.g. myHub_chrX:15200000-15800000) to a publicly accessible web server. Then, go to https://genome.ucsc.edu/index.html, click on My Data -> Track Hubs -> My Hubs, and add the link to your publicly available hub.txt file into the URL window.

Alternatively, you can load individual tracks from './out/chr:start-stop/trackFiles' as custom tracks on UCSC or IGV.

getStats.pl
------
The statistics for a region are provided in both totality and relative metrics. The relativity parameters vary by label type:

| |Total Region|Exonic|Intronic|Intergenic|Repetitive|Cis-Regulatory|
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
|Total Region|
|Coding Element|
|Noncoding Element|
|LINE Element|
|SINE Element|
|Alu Element|
|Promoter Element|
|Enhancer Element|
|Functional Element|
|Cis-Regulatory Module|

Output of genomeLabel
======
Example Directory Tree:
------

``` bash
.
├── data
│   ├── F5.hg38.enhancers.bed.gz
│   ├── hg19ToHg38.over.chain.gz
│   ├── hg38.chrom.sizes
│   ├── hg38.ncbiRefSeq.gtf.gz
│   ├── hg38_fair+new_CAGE_peaks_phase1and2.bed.gz
│   ├── lifted.bed
│   ├── remap2020_crm_macs2_hg38_v1_0.bed.gz
│   ├── rmsk.txt.gz
│   ├── segway_encyclopedia.bed.gz
│   └── unlifted.bed
├── hubs
│   └── myHub_chrX:15200000-15800000
│       ├── hg38
│       │   ├── codingTrack.bb
│       │   ├── crmTrack.bb
│       │   ├── enhancerTrack.bb
│       │   ├── functional_elementTrack.bb
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
│       │   ├── F5.hg38.enhancers.bed.gz
│       │   ├── hg38.chrom.sizes
│       │   ├── hg38.ncbiRefSeq.gtf.gz
│       │   ├── hg38_fair+new_CAGE_peaks_phase1and2.bed.gz
│       │   ├── remap2020_crm_macs2_hg38_v1_0.bed.gz
│       │   ├── rmsk.txt.gz
│       │   └── segway_encyclopedia.bed.gz
│       ├── intersectFiles
│       │   ├── Alu.bed.gz
│       │   ├── LINE.bed.gz
│       │   └── SINE.bed.gz
│       ├── labelFiles
│       │   ├── Alu.bed.gz
│       │   ├── LINE.bed.gz
│       │   ├── SINE.bed.gz
│       │   ├── coding.bed.gz
│       │   ├── crm.bed.gz
│       │   ├── enhancer.bed.gz
│       │   ├── exon.bed.gz
│       │   ├── functional_element.bed.gz
│       │   ├── intergenic.bed.gz
│       │   ├── intron.bed.gz
│       │   ├── noncoding.bed.gz
│       │   └── promoter.bed.gz
│       ├── trackFiles
│       │   ├── codingTrack.bed
│       │   ├── crmTrack.bed
│       │   ├── enhancerTrack.bed
│       │   ├── functional_elementTrack.bed
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
### Directory Contents
**********

|Directory | Contents|
|:---: | :---:|
|data | Data used to generate labels and tracks for *all regions* (i.e. raw NCBI RefSeq, FANTOM5, RepeatMasker, and UCSC data)|
|out | Regional directories created for each *input region* in the format **chr:start-stop**|
|out/chr:start-stop | Regional subdirectories and summary file created for each *input region*|
|out/chr:start-stop/dataFiles | Data used to generate labels and tracks for *input region* (i.e. filtered NCBI RefSeq, FANTOM5, RepeatMasker, and UCSC data)|
|out/chr:start-stop/labelFiles | Labels created for each feature in the *input region* (i.e. exon, intron, coding, etc.)|
|out/chr:start-stop/trackFiles | Tracks created for different combinations of features in the *input region* (i.e. genomic, repetitive elements, etc.)|
|out/chr:start-stop/intersectFiles | Files created to compute the overlap between different features in the *input region* (e.g. percent distribution of Alu elements in intronic vs. intergenic regions)|
|out/chr:start-stop/summary.txt | File created for summary statistics in the *input region* (e.g. coverage and average)|
|hubs/myHub_chr:start-stop | Track hub created to store *all* tracks generated for the *input region*|

Example Data Tracks:
------
![Image of Labelling Schema](https://docs.google.com/drawings/d/e/2PACX-1vQ51t4D1h96WMh588J429qSXkb_Fa6Cg_PhF3FHI4t2yPqMk1nzN0g54jFnf6wyD3hjs0qZS0brCaf3/pub?w=960&h=720)

This is an example of the priority levels behind the labelling of the **"Genomic Coding"** track, which is only one of the many tracks generated:
1. **Genomic Track** - Exonic, intronic, and intergenic elements (genomicTrack.bed)
2. **Coding Track** - Coding and noncoding exonic elements (codingTrack.bed)
3. **Genomic Coding Track** - Coding and noncoding exonic, intronic, and intergenic elements (genomicCodingTrack.bed)
4. **Promoter Track** - Promoter elements (promoterTrack.bed)
5. **Enhancer Track** - Enhancer elements (enhancerTrack.bed)
6. **Repeat Track** - Repetitive elements, including LINEs, SINEs, and Alu (repeatTrack.bed)
7. **Functional Element Track** - Putative cell type-agnostic functional elements (funtional_elementTrack.bed)
8. **Cis-Regulatory Module (CRM) Track** - Putative cell type and transcriptional regulator-agnostic binding sites (crmTrack.bed)

**Note:** In order of classification priority, we have: **EXON** (**CDS** > **exon** > **UTR**) > **INTRON** > **INTERGENIC**

![Image of All Tracks](https://docs.google.com/drawings/d/e/2PACX-1vQu6f0fA1D0bx0ZyYR9ErXmIUtKLcLBKDiEyCQNUNW80IOYvIM7_Ods73hpA_wV-2shq4CUGxppKLrZ/pub?w=960&h=720)

*This tool was inspired by:*
- **davetang/defining_genomic_regions:** https://github.com/davetang/defining_genomic_regions
- **cschlaffner/TrackHubGenerator:**  https://github.com/cschlaffner/TrackHubGenerator
