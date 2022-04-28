# PCHi-C (Promoter Capture Hi-C)

PCHi-C is a 3C method to detect the promoter interactome of a given sample. For a complete understanding of both the experimental protocol and the computational pipeline we recommend reading the following papers:

* [Promoter Capture Hi-C: High-resolution, Genome-wide Profiling of Promoter Interactions](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6102006/pdf/jove-136-57320.pdf)

* [HiCUP: pipeline for mapping and processing Hi-C data](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4706059/pdf/f1000research-4-7903.pdf)

* [CHiCAGO: Robust Detection of DNA Looping Interactions in Capture Hi-C data](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0992-2)

* [Detecting chromosomal interactions in Capture Hi-C data with CHiCAGO and companion tools](https://www.nature.com/articles/s41596-021-00567-5)

## Dependencies

* [HiCUP](https://www.bioinformatics.babraham.ac.uk/projects/hicup/) (0.8.2)
* [Chicago R Package](https://bioconductor.org/packages/release/bioc/html/Chicago.html)  
* bowtie2 (2.3.2)
* R (>3.6)

**Note**: Check all dependencies of each software

## Summary of the workflow

1.  **Mapping and filtering**: HiCUP
2.  **Capture Efficiency**: HiCUP miscellaneous script
3.  **Interaction calling**: Chicago

## 1. Mapping and Filtering (HiCUP)

#### Pre-steps

To run HiCUP first we need to prepare several files: genome digest and bowtie2 index. To generate these files we need to download the reference genome from Ensembl repository, in our case we use the release 104 (03/2021) for the version [GRCh38](http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/).

* **Genome Digest**

To generate the genome digest we need to know which restriction enzyme we are using, in our case we use HindIII which target sequence is A^AGCTT. To perform this step we used the ```hicup_digester``` function, which looks for the target sequence in the genome and generate all the restriction fragments for this specific restriction enzyme.

```bash
hicup_digester --genome Human_GRCh37 --re1 A^AGCTT,HindIII *.fa
```

* **Bowtie2 Index**

To run HiCUP we need a genome index and since HiCUP works both with bowtie and bowtie2, we can choose between the two, in our case we use bowtie2 because it's faster and allow parallelization. So to perform the genome indexing we used the ```bowtie2-build``` function, which takes the fasta files of the reference genome and indexes them.

```bash
bowtie2-build 1.fa,2.fa,...,MT.fa Human_GRCh38
```

#### Running HiCUP

Once we have the necessary files we can run HiCUP in our fastq files. To run HiCUP we need a config file, to obtain a template file we can use the function ```hicup --example```, the config file looks like this:

```bash
#Directory to which output files should be
#written (optional parameter)
#Set to current working directory by default
Outdir:

#Number of threads to use
Threads: 1

#Suppress progress updates (0: off, 1: on)
#Quiet:0

#Retain intermediate pipeline files (0: off, 1: on)
Keep:0

#Compress outputfiles (0: off, 1: on)
Zip:1

#Path to the alignment program (Bowtie or Bowtie2)
#Remember to include the executable Bowtie/Bowtie2 filename.
#Note: ensure you specify the correct aligner i.e.
#Bowtie when using Bowtie indices, or Bowtie2 when using Bowtie2 indices.
#In the example below Bowtie2 is specified.
Bowtie2: /usr/local/bowtie2/bowtie2

#Path to the reference genome indices
#Remember to include the basename of the genome indices
Index: /pathzº/Genomes/Human/GRCh37/Human_GRCh37

#Path to the genome digest file produced by hicup_digester
Digest: Digest_Human_genome_HindIII_None_12-32-06_15-04-2021.txt

#FASTQ format (valid formats: 'Sanger', 'Solexa_Illumina_1.0',
#'Illumina_1.3' or 'Illumina_1.5'). If not specified, HiCUP will
#try to determine the format automatically by analysing one of
#the FASTQ files. All input FASTQ will assumed to be in that
#format.
Format: 

#Maximum di-tag length (optional parameter)
Longest: 1000

#Minimum di-tag length (optional parameter)
Shortest: 150

#FASTQ files to be analysed, placing paired files on adjacent lines
sample_one_1.fq.gz
sample_one_2.fq.gz

sample_two_1.fq.gz
sample_two_2.fq.gz
```

Once we have filled the config file with our data, we can run HiCUP using the following command:

```bash
hicup --config config_file.txt
```

This give us a filtered bam file with the valid unique di-tags (paired reads) and a report in html format of all the process.

## 2. Capture Efficiency

Once we have the bam file from HiCUP we have to assess for the capture efficiency, to do this we split the reads that are captured by the capture system and the non-captured reads. To perform this step we use a miscellaneous script inside HiCUP.

To use this script we need the bam file from HiCUP and a bed file with the coordinates of the captured restriction fragments.

For this purpose we have generated a bed file with the coordinates of the captured fragments and also the annotation for each fragment, this extra information will be necessary for the other steps in further analysis.

```bash
perl HICUP/0.8.2/Misc/hicup_capture --baits baits.bed hicup.bam
```

With this step we get a bam file with only those read with at least one end captured.

## 3. Interaction calling (Chicago)

#### Pre-steps

To run Chicago first we need to prepare several files, what they called the Design Directory. To create all these files we can use several scripts (Chicago Tools) available in their [bitbucket repository](https://bitbucket.org/chicagoTeam/chicago/src/master/chicagoTools/), together with well documented information about each file and script usage. Furthemore, we need to reformat the bam file to a specific format used by Chicago as input, called chinput.

* **rmap & baitmap files**

The first 2 files we have to generate are the rmap and baitmap. The rmap contains all the possible restriction fragments and an ID for each of them, to get this information we have to provide the digest file generated by HiCUP in the first step. The baitmap is a subset of the rmap containing only those fragments that are captured (i.e. bait-fragments) and an extra column with the annotation information for each captured fragment, here is where we use again the bed file from the Capture Efficiency step.

```bash
perl create_baitmap_rmap.pl Digest_Human_genome_HindIII_None_12-32-06_15-04-2021.txt baits.bed
```

* **NPerBin, NBaitsPerBin & Proximal Other End files**

To generate the 3 other files of the DesignDir we need to have generated previously the rmap and baitmap and have both of them in the desired DesignDir, then we only need to run the following script.

```bash
python makeDesignFiles.py --designDir=./designDir --minFragLen=150 --maxFragLen=40000
```

* **Chinput**
 
To generate the chinput files from the bam files that we have already generated, we need the rmap and baitmap files, using the following command:

```bash
bam2chicago.sh captured.hicup.bam design.baitmap design.rmap path/chinputs/output_prefix
```

#### Running Chicago

Once we have all the design files and the chinputs we can run the interaction calling with Chicago.

We can use the ```runChicago.R``` script available in their bitbucket or make our own script to run Chicago. If we choose to use their script it will be as follows:

```bash
Rscript runChicago.R --design-dir DESIGN-DIR <input-files> <output-prefix>
Rscript runChicago.R --design-dir ./designDir path/chinputs/output_prefix.chinput sample_chicago
```

This script has a lot of different arguments, if you want to explore them check the bitbucket repository.

In case we want to make our own script the main command we need to use are the following:

```R
## We need the DesignDir, the chinput folder, the output directory and the output prefix

DesignDir <- "./designDir"
DataPath <- "path/chinputs/"
outputDirectory <- "path/chicago_output/"
outprefix <- "sample_chicago"


## Loading library
library(Chicago)


## Determining chinput file path
chinput <- list.files(DataPath,full.names = T, pattern= "chinput")


## Setting experiment with the DesignDir
cd <- setExperiment(designDir = DesignDir)

## Setting the seed to be reproducible (this parameter is not controlled in runChicago.R)
cd@settings$brownianNoise.seed <- 103

## Reading sample
cd <- readSample(file = chinput, cd=cd)

## Run the chicago pipeline
cd <- chicagoPipeline(cd, outprefix = outprefix, printMemory = T)

saveRDS(cd, paste0(outputDirectory,"/",outprefix, ".Rds"))
## Exporting results with 2 cutoffs, 0 and 5
exportResults(cd, file.path(outputDirectory, paste0(outprefix,"_cutoff_0")), cutoff = 0,format = "interBed")
exportResults(cd, file.path(outputDirectory, paste0(outprefix,"_cutoff_5")), cutoff = 5,format = "interBed")
```
