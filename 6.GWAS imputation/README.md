# GWAS Tissue Set Enrichment Analysis liCHi-C  

Enrichment of GWAS signals in liCHi-C interacting regions have to be tested to confirm the tissue and context are relevant to disease  

For a step to step tutorial, read [this](https://github.com/JavierreLab/Gwas_TissueSetEnrich/blob/main/example/README.md)  
For an in depth explanation on usage of the scripts, read [this](https://github.com/JavierreLab/Gwas_TissueSetEnrich/blob/main/scripts/README.md)

#### Dependencies

* R
* liftover
* tabix
* bgzip
* Chicago R package

## Summary of workflow:
1. Format Summary Statistics
2a. Poor Man's Imputation and Compute Posterior Probability
2b. Blockshifter
3. Plot



## Workflow
This is built to work in marenostrum only, not in local, as it requires too much computational energy and would not work in local.  
Easily adaptable to other supercomputers (CSUC) with the same SLURM and GREASY use. 

#### Prerequisits

* Download reference genotypes from 1000Genome project [version GRCh38](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/). Filter by sample ancestry (e.g. EUR) and index filtered reference genotype using tabix for faster reading. Script for this - **0_Download1KG_GRCh38.sh**

* Obtain GWAS datasets of interest. These will be "summary statistics" files  
   * [GWAS Catalog](https://www.ebi.ac.uk/gwas/downloads/summary-statistics) Full Summary Statistics are at FTP sites
   * [UKBB](https://docs.google.com/spreadsheets/d/1kvPoupSzsSFBNSztMzl04xMoSC3Kcx3CrjVf4yBmESU/edit?ts=5b5f17db#gid=178908679)

* Compute peak matrix with information of significant interactions of all cell types of interest as described [here](https://github.com/JavierreLab/liCHiC/tree/main/1.liCHiC%20Processing).
```
makePeakMatrix.R --twopass --notrans <rds_table> <out_prefix>
```
* Optionally, compute your own "approximate LD blocks" BED file based on recombination frequency data. I recommend you use the file **1cM_LDblocks_GRCh38.bed** obtained as explained [here].

## 1. Format Summary Statistics
Summary Statistics could have been generated in many different formats and genome assemblies. For this reason, after downloading each one, best practice is to format the file. Script for this - **1_Format_SummaryStatistics.sh, which calls 1_format_ss.R**. Summary Statistics from GWAS Catalog in the harmonized format don't need to be liftover, they are already in the genome assembly GRCh38. However, the ones from the UKBB are in the GRCh37 and need to be liftover to GRCh38.
```
./1_Format_SummaryStatistics.sh [-q -c -l -h] <Work.Dir> <Sum.Stats>
./1_Format_SummaryStatistics.sh /path/to/work/dir ChronicLymphocytic_UKBB_20001_1055-gwas.imputed_v3.both_sexes.Build37.tsv
```

## 2a.1. Poor Man's Imputation
Poor Man's imputation is performed to associate more variants to the trait of interest, taking into account LD and R-squared thresholds.
For every LD block, variants from the reference panel (1KG) which are not associated to the trait are given the p-value of the variant (already associated to the trait at hand) with which it is in largest linkage disequilibrium. If the minimum R-squared (o.6) is not achieved, the two variants are not similar enough and therefore you cannot predict one pval with the other, so imputation is not possible. Percentage of successful PMI tends to be about 60%. Once PMI can been computed the posterior probability of all the associated variants to the trait at hand may be computed. The posterior probability is computes using wakefield's implementation - explained before. 

## 2a.2. Compute Posterior Probability 
#### Why compute a posterior probability?  
These posterior probabilities (PPi) are calculated for each SNP in a trait's summary statistic to take into account Linkage Disequilibrium (LD refers to the presence of a statistical association between allelic variants within a population due to the history of recombination, mutation, and selection in a genomic region) assuming that there is only one "causal" variant within a region, while the rest have no biological impact.  

To compute a PPi we need the following information:
* **p-Value** The p-value of the variant in the summary statistic tells us whether the variant is or not associated to the particular trait at hand
* **Minor Allele Frequency** MAF of the variant is frequency at which the less common allele of a variable site is found in a population.  
* **Sample Size** Number of individuals which partook in the study
* **Type of Study** Was the study a case/control (Rheumatoid Arthitis) or a quantitative study (Platelet Volume)? Calculation of variance shrinkage changes depending on study type
* **Proportion of Cases** Proportion of individuals which were cases. If study is quantitative, proportion = 1
* **Prior Probability**  
* **LD Blocks** Posterior Probability must be calculated taking into account LD assuming that there is only one "causal" variant within a block  

Steps 1 and 2 are performed with the same script **2a_compute_pmi_GRCh38.R**:
```
Rscript 2a_compute_pmi_GRCh38.R <Chromosome> <Summary.Statistics> <Gwas.Type> <Sample.Size> <Proportion.Cases> <LD.Regions> <Out.Dir> <Tabix.Bin = /apps/HTSLIB/1.8/INTEL/bin/tabix>
Rscript 2a_compute_pmi_GRCh38.R 22 /path/to/ChronicLymphocytic_UKBB_20001_1055.bed.bgz CC 361141 0.000656253 /path/to/1cM_LDblocks_GRCh38.bed /path/to/out/dir/ /path/to/tabix/bin
```

## 2b. BlockShifter  
Compute z-scores for enrichment of variants associated to a particular trait within promoter-interacting regions taking into account correlation structure in both GWAS and PIR datasets. (Examines GWAS signal enrichment at PIRs overcoming LD and interaction fragment correlation). 
```
Rscript 2b_blockshifter_GRCh38.R  <PeakMatrix> <Gwas.PPi> <Test.Set> <Control.Set> <Metric> <Out.Dir>
Rscript 2b_blockshifter_GRCh38.R  path/to/liHiC.peakmatrix.txt /path/to/gwasppi/CLL_ukbb.pmi Ery_500K,MK_500K,Mon_500K,nB_1M,nCD4_500K,nCD8_500Kmix,Ery_cambridge,MK_cambridge,Mon_cambridge,nB_cambridge,nCD4_cambridge,nCD8_cambridge,nB_100K,nB_250K,nB_500K,nB_50K,CLP_WT_merge_45,CMP_WT_merge_45,HSC_WT_merge_15 EP_cambridge ppi /path/to/out/dir
```

## Acknowledgments  
* [Olly Burren ](https://github.com/ollyburren/CHIGP/tree/master/R)  
* My lab mates :)

