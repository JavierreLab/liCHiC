# GWAS Tissue Set Enrichment Analysis liCHi-C  

Enrichment of GWAS signals in liCHi-C interacting regions have to be tested to confirm the tissue and context are relevant to disease  

For a step to step tutorial, read [this](https://github.com/JavierreLab/Gwas_TissueSetEnrich/blob/main/example/README.md)  
For an in depth explanation on usage of the scripts, read [this](https://github.com/JavierreLab/Gwas_TissueSetEnrich/blob/main/scripts/README.md)

## Credits

The method and software was co-developed by [Chris Wallace](http://chr1swallace.github.io/) and [Olly Burren](http://ollyburren.github.io/), and presented in [Javierre et al., 2016](https://doi.org/10.1016/j.cell.2016.09.037).

The scripts in this repository have been adapted from the original ones from [CHIGP](https://github.com/ollyburren/CHIGP).

### Dependencies

* R
* liftover
* tabix
* bgzip
* CHiCAGO R package
* Only for generating the plots: Pyhton3 with pandas (>= 1.2.3) and matplotlib (>=3.5.0)

## Summary of workflow:
1. Format Summary Statistics
2. Poor Man's Imputation and Compute Posterior Probability
3. Blockshifter
4. Plots

## Workflow
This is built to work only in the supercomputer Marenostrum, not in local, as it requires too much computational energy and would not work in local.  
Easily adaptable to other supercomputers (e.g., CSUC) with the same SLURM and GREASY use. 

#### Requirements

* Download reference genotypes from 1000Genome project [version GRCh38](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/). Filter by sample ancestry (e.g. EUR) and index filtered reference genotype using tabix for faster reading. Script for this - **0_Download1KG_GRCh38.sh**

* Obtain GWAS datasets of interest. These will be "summary statistics" files  
   * [GWAS Catalog](https://www.ebi.ac.uk/gwas/downloads/summary-statistics) Full Summary Statistics are at FTP sites
   * [UKBB](https://docs.google.com/spreadsheets/d/1kvPoupSzsSFBNSztMzl04xMoSC3Kcx3CrjVf4yBmESU/edit?ts=5b5f17db#gid=178908679)

* Compute peak matrix with information of significant interactions of all cell types of interest as described [here](https://github.com/JavierreLab/liCHiC/tree/main/1.liCHiC%20Processing).
```
makePeakMatrix.R --twopass --notrans <rds_table> <out_prefix>
```
* Optionally, compute your own "approximate LD blocks" BED file based on recombination frequency data. I recommend you use the file **1cM_LDblocks_GRCh38.bed** obtained as detailed [here](https://carrerasresearch-my.sharepoint.com/:f:/g/personal/plopez_carrerasresearch_org/EgvlXxvJdzlGrjJxwOEf6n4B0eJQNQTw3Zq3Ka02Yw6ewQ?e=qz68eV).

## 1. Format Summary Statistics
Summary Statistics could have been generated in many different formats and genome assemblies. For this reason, after downloading each one, best practice is to format the file. Script for this: **1_Format_SummaryStatistics.sh**, which calls **1_format_ss.R**. 
The Summary Statistics from the GWAS Catalog in the harmonized format are already in the genome assembly GRCh38. The ones from the UKBB are in the GRCh37 and, therefore, they have to be liftover to GRCh38.
```
./1_Format_SummaryStatistics.sh [-q -c -l -h] <Work.Dir> <Sum.Stats>
# Example:
./1_Format_SummaryStatistics.sh /path/to/work/dir ChronicLymphocytic_UKBB_20001_1055-gwas.imputed_v3.both_sexes.Build37.tsv
```

## 2a. Poor Man's Imputation and Compute Posterior Probability
### 2a.1. Poor Man's Imputation
Poor Man's imputation is performed to associate more variants to the trait of interest, taking into account LD and R-squared thresholds.
For every LD block, variants from the reference panel (1KG) which are not associated to the trait are given the p-value of the variant (already associated to the trait at hand) with which it is in largest linkage disequilibrium. If the minimum R-squared (0.6) is not achieved, the two variants are not similar enough and therefore you cannot predict one pval with the other, so imputation is not possible. Percentage of successful PMI tends to be about 60%. 

### 2a.2. Compute Posterior Probability 
Once PMI have been computed, the posterior probability of all the associated variants to the trait at hand may be computed. The posterior probability is computes using wakefield's implementation.

#### Why compute a posterior probability?  
These posterior probabilities (PPi) are calculated for each SNP in a trait's summary statistic to take into account Linkage Disequilibrium (LD refers to the presence of a statistical association between allelic variants within a population due to the history of recombination, mutation, and selection in a genomic region) assuming that there is only one "causal" variant within a region, while the rest have no biological impact.  

To compute a PPi we need the following information:
* **p-Value** The p-value of the variant in the summary statistic tells us whether the variant is or not associated to the particular trait at hand
* **Minor Allele Frequency** MAF of the variant is frequency at which the less common allele of a variable site is found in a population.  
* **Sample Size** Number of individuals which partook in the study
* **Type of Study** Was the study a case/control (e.g., Rheumatoid Arthitis) or a quantitative study (e.g., Platelet Volume)? Calculation of variance shrinkage changes depending on study type
* **Proportion of Cases** Proportion of individuals which were cases. If study is quantitative, proportion = 1
* **Prior Probability**  
* **LD Blocks** Posterior Probability must be calculated taking into account LD assuming that there is only one "causal" variant within a block  

Info on Sample size and proportion of cases may be found at the source of the Summary Statistics.

Steps 1 and 2 are performed with the same script **2a_compute_pmi_GRCh38.R**:
```
Rscript 2a_compute_pmi_GRCh38.R <Chromosome> <Summary.Statistics> <Gwas.Type> <Sample.Size> <Proportion.Cases> <LD.Regions> <Out.Dir> <Tabix.Bin = /apps/HTSLIB/1.8/INTEL/bin/tabix>
# Example:
Rscript 2a_compute_pmi_GRCh38.R 22 /path/to/ChronicLymphocytic_UKBB_20001_1055.bed.bgz CC 361141 0.000656253 /path/to/1cM_LDblocks_GRCh38.bed /path/to/out/dir/ /path/to/tabix/bin
```

## 2b. BlockShifter  
Compute z-scores for enrichment of variants associated to a particular trait within promoter-interacting regions taking into account correlation structure in both GWAS and PIR datasets. (Examines GWAS signal enrichment at PIRs overcoming LD and interaction fragment correlation). 
```
Rscript 2b_blockshifter_GRCh38.R  <PeakMatrix> <Gwas.PPi> <Test.Set> <Control.Set> <Metric> <Out.Dir>
# Example:
Rscript 2b_blockshifter_GRCh38.R  path/to/liHiC.peakmatrix.txt /path/to/gwasppi/CLL_ukbb.pmi Ery_500K,MK_500K,Mon_500K,nB_1M,nCD4_500K,nCD8_500Kmix,Ery_cambridge,MK_cambridge,Mon_cambridge,nB_cambridge,nCD4_cambridge,nCD8_cambridge,nB_100K,nB_250K,nB_500K,nB_50K,CLP_WT_merge_45,CMP_WT_merge_45,HSC_WT_merge_15 EP_cambridge ppi /path/to/out/dir
```

## 3. Plots
### 3a. Zscore plot Myeloid vs. Lymphoid
Use the jupyter notebook **3a_Zscore_plot_MvsL.ipynb** modifying the desired arguments. For using this script is necessary to have installed python3 with the libraries pandas (>= 1.2.3) and matplotlib (>3.5.0).

### 3b. Heatmap Blockshifter Z-score
Use the script **3_plot_heatmap_blockshifter.R** modifying the desired arguments.


## Acknowledgments  
* [Olly Burren ](https://github.com/ollyburren/CHIGP/tree/master/R)  
* My lab mates :)
