# **Inbred Strain** 
## Introduction
This report describes the genome sequence results for the IM and M-AB and a comparison of genome sequences between the two inbred lines. The codes and text files used in this analysis are summarized below.
## Dependencies
### Download DNA sequence & GTF<br>
  - seqkit (2.1.0)
### SNPs call　<br>
  - fastq-dump (2.11.3)<br>
  - fastp (0.20.1)<br>
  - bwa-mem (0.7.17)<br>
  - gatk (4.4.0.0)<br>
  - vcftools (0.1.16)<br>
  - samtools (1.13)<br>
  - snpeff (4.3)<br>
### Phylogenetic analysis <br>
- python (3.10.9)<br>
  - vcf2phylip (2.0)<br>
- modeltest-ng (0.1.7)<br>
 - raxml-ng (1.2.0)<br>
### Genomic characteristics analysis <br>
  - R (4.2.3) <br>
   - ggplot2 (3.4.2)
   - openxlsx (4.2.5.2)
   - patchwork (1.1.2)
   - ggsignif (0.6.4)
   - dplyr (1.1.2)
   - ggbreak (0.1.1)
   - stringr (1.5.0)
   - UpSetR (1.4.0)
   - reshape2 (1.4.4)
   - sets (1.0-24)
## Download DNA sequence & GTF 
Reference data is downloaded by the following command.
```
wget https://ftp.ensembl.org/pub/release-101/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna.primary_assembly.fa.gz
```
For fa file as reference, only Chr 1~25 and MT are used.
```
for i in `cat chr_mt.txt`
do
seqkit grep -p ${i} Danio_rerio.GRCz11.dna.primary_assembly.fa >> Danio_rerio.GRCz11.dna.primary_assembly-only-chr.fa
done
```
GTF file
```
wget https://ftp.ensembl.org/pub/release-109/gff3/danio_rerio/Danio_rerio.GRCz11.109.gff3.gz
```

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
