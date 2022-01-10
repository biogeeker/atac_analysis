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
# filtering, retain properly paired reads
# -f 3: only include alignments marked with the SAM flag 3, which means "properly paired and mapped"
$ samtools view -bh -f 3 MY_sorted.bam > MY_sorted_filtered.bam

# remove PCR duplicates. Here we use the GATK v4.1.
$ export PATH=/home/wy/disk/aATAC-seq/tools/gatk-4.1.4.1/:$PATH
$ gatk MarkDuplicates --REMOVE_DUPLICATES TRUE -I MY_sorted_filtered.bam -M MY_dups_metrics.txt -O MY_noDups.bam

# plot fragment size figures
$ gatk CollectInsertSizeMetrics -I MY_noDups.bam -O MY_insert_size_metrics.txt -H MY_histogram.pdf -M 0.5
```

The forth step is to call peaks:
```
# We use macs2 to call peaks.
# only call NFR peaks
./split_bam.sh in.bam 147


```
