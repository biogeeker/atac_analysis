# *P. sojae* ATAC-seq analysis pipeline
>The workflow is designed for ATAC-seq analysis in article "*ATAC-seq reveals the landscape of open chromatin and cis-regulatory elements in the *Phytophthora sojae* genome*". The sequencing data can be downloaded from NCBI: *[BioProject accession number PRJNA761250]*.

### Remove low quality reads and adaptors in raw data
```
# Here we use trim_galore do this job.
$ trim_galore -q 25 --phred33 --length 35 -e 0.1 --stringency 4 --paired -o clean/ ./S025_jiyi-A_MY_AHVMTLCCXY_S5_L006_R1_001.fastq.gz ./S025_jiyi-A_MY_AHVMTLCCXY_S5_L006_R2_001.fastq.gz
# We can check the clean data using tool - fastqc, if the result of QC passes, we can start mapping the reads to the genome. 
$ fastqc -o ./ clean/*gz
$ multiqc ./
```

### Map clean reads to the *P. sojae* reference genome
```
# Firstly, build genome index
$ bowtie2-build Psojae11.fasta psojae
# Secondly, start the alignment
$ bowtie2 --very-sensitive -X 1000 -x psojae -1 S025_jiyi-A_MY_AHVMTLCCXY_S5_L006_R1_001_val_1.fq.gz -2 S025_jiyi-A_MY_AHVMTLCCXY_S5_L006_R2_001_val_2.fq.gz | samtools sort -@ 6 -O bam -o MY_sorted.bam
```

### Post-alignment quality control
```
# Filtering, retain properly paired reads
# -f 3: only include alignments marked with the SAM flag 3, which means "properly paired and mapped"
$ samtools view -bh -f 3 MY_sorted.bam > MY_sorted_filtered.bam

# Remove PCR duplicates. Here we use the GATK v4.1. Create environment variable.
$ export PATH=/home/wy/disk/aATAC-seq/tools/gatk-4.1.4.1/:$PATH
$ gatk MarkDuplicates --REMOVE_DUPLICATES TRUE -I MY_sorted_filtered.bam -M MY_dups_metrics.txt -O MY_noDups.bam

# By the way, making index for BAM file
$ samtools index MY_noDups.bam

# plot Nucleosome Fragment Size distribution
$ gatk CollectInsertSizeMetrics -I MY_noDups.bam -O MY_insert_size_metrics.txt -H MY_histogram.pdf -M 0.5
```

### Call peaks

>If you followed original protocol for ATAC-Seq, you should get Paired-End reads. If so, I would suggest you just use --format BAMPE to let MACS2 pileup the whole fragments in general. But if you want to focus on looking for where the ‘cutting sites’ are, then --nomodel --shift -100 --extsize 200 should work. --By Liu

```
# We use macs2 to call peaks
# Only call NFR peaks, NFR < 130 bp, determined by Nuclesome Fragment Size Plot
$ split_bam_fragment.sh MY_noDups.bam 130
$ macs2 callpeak -f BAMPE -g 8.8e7 --keep-dup all -n MY_NFR -t MY_NFR.bam --outdir macs2/MY_NFR 2> macs2.log
# Filtering blacklist regions: this is very important in calling peaks, but non-model speceis don not have marked black list regions well. I speculated that extreme peaks may represent the candidata black list regions. And so, we can mask it using seqtk and re-run masc2 call peaks again.
```

### Detect TF footprints using [RGT-HINT](http://www.regulatory-genomics.org/motif-analysis/introduction/)
>For peak calling, The paper mentioned "all reads aligning to + strand were offset by +4 bp, all reads aligning to the - strand are offset -5 bp". The offsets are only really important when doing TF footprinting. HINT-ATAC can correct Tn5 bias automatically and thus shifting reads unnessarily.
```
# Here we use HINT-ATAC for footprintings and need to modify config file by yourself
$ rgt-hint footprinting --atac-seq --organism=ps11 --paired-end --output-location=my_atac --output-prefix=MY_NFR_Footprint ./MY_NFR.bam ./MY_NFR_MACS2_peaks.narrowPeak
# Generate bias corrected bw file for genome browser visualizing
$ rgt-hint tracks --bc --bigWig --organism=ps11 ./MY_NFR.bam ./MY_NFR_MACS2_peaks.narrowPeak --output-location=my_atac --output-prefix=MY_NFR

# footprints map to JASPAR database and annotate motif ids
$ rgt-motifanalysis matching --organism=ps11 --input-files my_atac/MY_NFR_Footprint.bed --output-location=my_atac
# Also, we can modify customized footprints file to generate footprintings by meme-suite motifs coordinate files. We can visualize the peak-dip-peak of motifs. It is very interesting.
```

### Annotate peaks and enrich motifs
```
# Use R ChIPseeker for annotating peaks 

library("GenomicFeatures")
library("ChIPseeker")
# Create txdb file
TxDb_Ps <- makeTxDbFromGFF("Psojae.gtf") # GTF format
txdb = TxDb_Ps
txs <- transcripts(txdb)

# Read peaks
ATACpeak <- readPeakFile("MY_NFR_MACS2_peaks.narrowPeak")

# Annotate peaks
y <- annotatePeak(ATACpeak, tssRegion=c(-3000, 3000), TxDb=txdb, addFlankGeneInfo=TRUE, flankDistance=4000)

write.table(as.data.frame(y), file="MY_NFR_MACS2_peaks.narrowPeakAnno.xls", sep='\t', quote = F)
# Visualize the result
plotAnnoPie(y)
# Note: the table header position may be incorrect in the first row
```

### Other commands
>Calculate openness "cut sites" using NucleoATAC package
```
$ dnase_cut_counter.py -A MY_NFR_MACS2_peaks.narrowPeak MY_NFR130.bam cut_sites.txt
```

>HMMRATAC is designed to analyze ATAC-seq data, maybe you will try it
```
samtools sort ATACseq.bam -o ATACseq.sorted.bam
samtools index ATACseq.sorted.bam ATACseq.sorted.bam.bai
samtools view -H ATACseq.sorted.bam| perl -ne 'if(/^@SQ.*?SN:(\w+)\s+LN:(\d+)/){print $1,"\t",$2,"\n"}' > genome.info
java -jar HMMRATAC_V1.2.4_exe.jar -b ATACseq.sorted.bam -i ATACseq.sorted.bam.bai -g genome.info
awk -v OFS="\t" '$13>=10 {print}' NAME_peaks.gappedPeak > NAME.filteredPeaks.gappedPeak
awk -v OFS="\t" '$5>=10 {print}' NAME_summits.bed > NAME.filteredSummits.bed
```

