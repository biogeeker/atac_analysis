The workflow is designed for ATAC-seq analysis in article "ATAC-seq reveals the landscape of open chromatin and cis-regulatory elements in the *Phytophthora sojae* genome". The sequencing data can be downloaded from NCBI BioProject accession number: PRJNA761250. Now, let us start the process of analysis.

The first step is to remove low quality reads and adaptors in raw reads:
```
# Here we use trim_galore do this job.
$ trim_galore -q 25 --phred33 --length 35 -e 0.1 --stringency 4 --paired -o clean/ ./S025_jiyi-A_MY_AHVMTLCCXY_S5_L006_R1_001.fastq.gz ./S025_jiyi-A_MY_AHVMTLCCXY_S5_L006_R2_001.fastq.gz
# We can check the clean data using tool - fastqc, if the result of QC passes, we can start mapping the reads to the genome. 
$ fastqc -o ./ clean/*gz
```

The second step is to map clean reads to the reference genome:
```
$ bowtie2 --very-sensitive -X 1000 -x psojae -1 S025_jiyi-A_MY_AHVMTLCCXY_S5_L006_R1_001_val_1.fq.gz -2 S025_jiyi-A_MY_AHVMTLCCXY_S5_L006_R2_001_val_2.fq.gz | samtools sort -@ 6 -O bam -o MY_sorted.bam
```

The third step is carried out to post-alignment QC:
```

```
