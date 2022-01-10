The workflow is designed for ATAC-seq analysis in article [ATAC-seq reveals the landscape of open chromatin and cis-regulatory elements in the Phytophthora sojae genome]

remove low quality reads and adaptors
```
$ trim_galore -q 25 --phred33 --length 35 -e 0.1 --stringency 4 --paired -o clean/ ./S025_jiyi-A_MY_AHVMTLCCXY_S5_L006_R1_001.fastq.gz ./S025_jiyi-A_MY_AHVMTLCCXY_S5_L006_R2_001.fastq.gz ; trim_galore -q 25 --phred33 --length 35 -e 0.1 --stringency 4 --paired -o clean/ ./S026_jiyi-A_IF_AHVMTLCCXY_S6_L006_R1_001.fastq.gz ./S026_jiyi-A_IF_AHVMTLCCXY_S6_L006_R2_001.fastq.gz
```
