# **Inbred Strain** 
## Introduction

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
    - 

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


## Citation
Ken-ichiro Sadamitsu, Fabien Velilla, Minori Shinya, Makoto Kashima, Yukiko Imai, Toshihiro Kawasaki, Hiromi Hirata and Noriyoshi., "SakaiEstablishment of a zebrafish inbred strain through full sib-pair mating with capable of genetic manipulation in the embryo."  doi: <br>

<br>
