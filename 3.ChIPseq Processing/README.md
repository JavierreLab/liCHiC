# ChIP-seq

ChIP-seq is a method used to analyze protein interactions with DNA. ChIP-seq combines chromatin immunoprecipitation (ChIP) with massively parallel DNA sequencing to identify the binding sites of DNA-associated proteins.

## Dependencies

* [Trim Galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)  
* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)  
* [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)  
* [samtools](https://github.com/samtools/samtools)  
* [sambamba](https://lomereiter.github.io/sambamba/)
* [deepTools](https://deeptools.readthedocs.io/en/develop/index.html)  
* [macs2](https://github.com/macs3-project/MACS/wiki/Install-macs2)  


**Note**: Check all dependencies of each software

## Summary of the workflow

1.  **Trimming and Quality**: Trim Galore and FastQC
2.  **Alignment**: Bowtie2
3.  **Filtering**: samtools and sambamba
4.  **Coverage**: deepTools
5.  **Peak Calling**: Macs2

## 1. Trimming and Quality (Trim Galore and FastQC)

The first step of the pipeline is the trimming of the reads to remove any possible adapters present in the reads, using trim_galore:

```bash
trim_galore $FASTQ1 $FASTQ2 --output_dir $FASTQ_DIR --paired --basename $NAME -a $ADAPTERS -A $ADAPTERS --fastqc_args '--outdir $OUTDIR/quality' --cores $THREADS"
```

With these optiones trim_galore will trim our reads and then perform the quality check in the trimmed fastq file. But it's recomended to perform a quality check before doing the trimming.

## 2. Alignment (bowtie2)

Once the reads are cleaned from any adapters we can align them to the reference genome, we use the same one we used for liCHi-C, GRCh38.

```bash
bowtie2 -x $INDEX -1 $FASTQ1 -2 $FASTQ2 --very-sensitive -k 2 -t -p $THREADS -S $OUTDIR/$NAME.sam
```

We obtained a sam file with all the reads aligned. Then we tranform the sam to bam and sort it, and we generate an index for the bam file, all these using both samtools and sambamba.

```bash
samtools view -@ $THREADS -bS $OUTDIR/$NAME.sam > $OUTDIR/$NAME.bam"  >> $JOB
sambamba sort -t $THREADS -o $OUTDIR/$NAME.sort.bam $OUTDIR/$NAME.bam
sambamba index -t $THREADS -p $OUTDIR/$NAME.bam
```

## 2. Filtering (samtools and sambamba)

The filtering steps is the part were we filter out reads based on different criteria. We follow the criteria established by the ENCODE Consortium in their Data Standards: [link](https://docs.google.com/document/d/1lG_Rd7fnYgRpSIqrIfuVlAz2dW1VaSQThzk836Db99c/edit#)

Briefly the filtering has several steps:

a.  Filtering unmapped, mate_is_unmapped, secondary and duplicates, and for paired-end read the properly paired

```bash
sambamba view -h -t $THREADS -f bam -F $FILT_ARG $BAM > ${PATH_PREFIX}.clean.bam
samtools stats ${PATH_PREFIX}.clean.bam > ${PATH_PREFIX}.2b_filt1.stats
samtools flagstat ${PATH_PREFIX}.clean.bam > ${PATH_PREFIX}.3b_filt1.flagstat
```

b.  Filtering Blacklist together with chr MT

```bash
bedtools intersect -v -abam ${PATH_PREFIX}.clean.bam -b $BLACKLIST > ${PATH_PREFIX}_filtered.bam
samtools stats ${PATH_PREFIX}_filtered.bam > ${PATH_PREFIX}.2c_blacklist.stats
samtools flagstat ${PATH_PREFIX}_filtered.bam > ${PATH_PREFIX}.3c_blacklist.flagstat
```

c.  Fix mate and repeat filtering

```bash
sambamba sort -t $THREADS -n -o ${PATH_PREFIX}_namesort.bam ${PATH_PREFIX}_filtered.bam
samtools fixmate -m ${PATH_PREFIX}_namesort.bam ${PATH_PREFIX}_fixmate.bam
sambamba view -h -t $THREADS -f bam -F $FILT_ARG ${PATH_PREFIX}_fixmate.bam > ${PATH_PREFIX}_fixmate.clean.bam
samtools stats ${PATH_PREFIX}_fixmate.clean.bam > ${PATH_PREFIX}.2d_filt2.stats
samtools flagstat ${PATH_PREFIX}_fixmate.clean.bam > ${PATH_PREFIX}.3d_filt2.flagstat
```

d.  Mark duplicates

```bash
sambamba sort -t $THREADS -o ${PATH_PREFIX}_fixmate.sort.bam ${PATH_PREFIX}_fixmate.clean.bam
samtools markdup -@ $THREADS ${PATH_PREFIX}_fixmate.sort.bam ${PATH_PREFIX}.markdup.bam
```

e.  Report library complexity

```bash
sambamba sort -n -t $THREADS -o ${PATH_PREFIX}.markdup.sort.bam ${PATH_PREFIX}.markdup.bam
$REPORT
```

f.  Remove duplicates

```bash
samtools markdup -r -@ $THREADS ${PATH_PREFIX}_fixmate.sort.bam ${PATH_PREFIX}.markdup.bam 2> ${PATH_PREFIX}.4e_markdup.txt
sambamba view -h -t $THREADS -f bam -F $FILT_ARG ${PATH_PREFIX}.markdup.bam | sambamba sort -t 10 /dev/stdin -o ${PATH_PREFIX}.filt.nodup.bam
samtools stats ${PATH_PREFIX}.filt.nodup.bam > ${PATH_PREFIX}.2e_markdup.stats
samtools flagstat ${PATH_PREFIX}.filt.nodup.bam > ${PATH_PREFIX}.3e_markdup.flagstat
```
