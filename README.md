# **Inbred Strain** 
## Introduction
This report describes the genome sequence results for the IM and M-AB and a comparison of genome sequences between the two inbred lines. The codes and text files used in this analysis are summarized below.
## Dependencies
### SNPs call　<br>
  - fastq-dump (2.11.3)<br>
  - fastp (0.20.1)<br>
  - bwa-mem (0.7.17)<br>
  - gatk (4.4.0.0)<br>
  - vcftools (0.1.16)<br>
  - samtools (1.13)<br>
  - snpeff (4.3)<br><br>
### Phylogenetic analysis <br>
- python (3.10.9)<br>
  - vcf2phylip (2.0)<br>
- modeltest-ng (0.1.7)<br>
 - raxml-ng (1.2.0)<br><br>
### Genomic characteristics analysis <br>
  - R () <br>
    - ggplot2
    - openxlsx
    - patchwork
    - ggsignif
    - dplyr
    - ggbreak
## SNPs call
`snpscall.sh`completes the process from downloading fastq to SNPs calling. <br>
***Threading should be adjusted arbitrarily, as processing takes time.***<br>
```
bash snpscall.sh
```

## Phylogenetic analysis 
After SNPs Calls are completed, perform the following`phy.sh`.
```
bash phy.sh
```

## Genomic characteristics analysis
After merging the VCF in Phylogenetic analysis (`vcf-merge vcf/fin_vcf/* > vcf/merge/merged.vcf`) and completing the "merged.vcf", the following R can be run to calculate the data needed for the paper.

```
Rscript　anlysis.R
```

## Data availability
### BioProject
 - PRJNA1002090
### SRA
 - M-AB (26th)
   - M-AB1: SRR25514300
   - M-AB2: SRR25514299
   - M-AB3: SRR25514298
 - IM (47th)
   - IM1: SRR25514304
   - IM2: SRR25514303
   - IM3: SRR25514302
 - *AB
   - *AB1: SRR25514325
   - *AB2: SRR25514324
   - *AB3: SRR25514323
 - India
   - India1: SRR25514329
   - India2: SRR25514328
   - India3: SRR25514327

## Citation
Ken-ichiro Sadamitsu, Fabien Velilla, Minori Shinya, Makoto Kashima, Yukiko Imai, Toshihiro Kawasaki, Hiromi Hirata and Noriyoshi., "SakaiEstablishment of a zebrafish inbred strain through full sib-pair mating with capable of genetic manipulation in the embryo."  doi: <br>

<br>
