# genomeLabel
Create genomic labels to help identify candidate cis-regulatory regions.

# Data:
This program uses:
- NCBI RefSeq Curated Data (https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz)
- FANTOM5 Promoter Data (https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest/extra/CAGE_peaks/hg38_fair+new_CAGE_peaks_phase1and2.bed.gz)
- RepeatMasker Data (https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.sql)
- UCSC Chromosome Sizes Data (https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes)

# Scripts:
## 1) *makeLabels.pl*
- Command line argument(s): Input region in the format **chr:start-stop**

## 2) *makeTracks.pl*
- Command line argument(s): Input region in the format **chr:start-stop**

## 3) *makeHubs.pl*
- Command line argument(s):
  - Input region in the format **chr:start-stop**
  - Path to the UCSC utility **bedToBigBed** in the format **/path/to/utilities/**. If you are unsure, you can run the following command:
  
  ``` bash
  find / -type f -name "bedToBigBed" 2>/dev/null
  ```



*Inspired by:*
- **davetang/defining_genomic_regions:** https://github.com/davetang/defining_genomic_regions
- **cschlaffner/TrackHubGenerator:**  https://github.com/cschlaffner/TrackHubGenerator
