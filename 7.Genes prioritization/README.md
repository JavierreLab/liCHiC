# GWAS Genes prioritization
Capture Hi-C Omnibus Gene Score (COGS) algorithm allows the prioritization of relevant genes for each GWAS trait and cell type.

Briefly, this method takes into account linkage disequilibrium, interactions and functional SNP annotation. It uses the posterior probabilities computed in the previous step to compute a gene-score for all the genes involved in liCHi-C significant interactions in at least one cell type. This gene-score is composed of the contributions of 3 components: 
* coding SNPs in the annotated genes as computed by VEP; 
* promoter SNPs, that are those located in promoter regions (that is, in a bait or in the flanking HindIII fragments); 
* SNPs overlapping Promoter-Interacting Regions.

For more details, go to [Javierre et al., 2016](https://doi.org/10.1016/j.cell.2016.09.037).


## Credits
The method and software were co-developed by [Chris Wallace](http://chr1swallace.github.io/) and [Olly Burren](http://ollyburren.github.io/), and presented in [Javierre et al., 2016](https://doi.org/10.1016/j.cell.2016.09.037).

The scripts in this repository have been adapted from the original ones from [CHIGP](https://github.com/ollyburren/CHIGP).


## Dependencies
* R
  * data.table
  * GenomicRanges (bioconductor)
  * reshape2
* Python3 with pandas (>= 1.2.3)
* For Reactome Pathway Analysis:
  * ReactomePA
  * org.Hs.eg.db
  * clusterProfiler
  * httr
  * tidyverse
  * ggplot2

## Required input files
* Coding SNPs file
* Annotated Peakmatrix with information of significant interactions of all cell types of interest
* Annotation files for your capture system


## File preparation
### 1. Generate coding SNPs file using Ensembl Variant Effect Predictor (VEP)
You can use [this file](https://carrerasresearch-my.sharepoint.com/:t:/g/personal/plopez_carrerasresearch_org/EZeqg2w41yJOre19s2pnHZMBoH-Pihz43eYEVdJuf_cBSA?e=65ZsT2). It was generated using [VEP](https://www.ensembl.org/info/docs/tools/vep/index.html) with the genome version 104 and the variants from the EUR cohort from the 1000 Genomes Project (GRCh38) as detailed [here](https://github.com/JavierreLab/liCHiC/tree/main/6.GWAS%20imputation). 

Optionally, you can also generate your own coding SNPs file following these steps:
1. Run VEP for each chromosome. Example command:
```
./vep --cache --dir /path/to/VEP_cache -i /path/to/1000Genomes_GRCh38/EUR.chr22.phase3.GRCh38.vcf.gz -o /output/path/coding_SNPs_VEP_1KG_GRCh38/EUR.chr22.phase3.GRCh38.cSNPs_VEP.txt --coding_only --verbose --force_overwrite --offline
```
2. Format the VEP output format to 3-column format, using the script **1.FormatVEPoutput.py**

### 2. Annotate peakmatrix
The Peakmatrix obtained from [liCHiC pipeline](https://github.com/JavierreLab/liCHiC/tree/main/1.liCHiC%20Processing) contains the interactions classified as significant by CHiCAGO (CHiCAGO score >5) in at least one cell type. It has to be annotated by Gene Name or by Ensg ID, and become a gene-centered matrix. Moreover, it requires some formatting changes. This can be done with these two scripts: 
* **2.Annotate_Peakmatrix_1_GeneNames.R** - as arguments it takes:
  * pk.file - Original peakmatrix
  * annotation - This file will be unique to your capture system. Example: File with bait's annotation by [Gene Names](https://carrerasresearch-my.sharepoint.com/:u:/g/personal/plopez_carrerasresearch_org/EST6A-ksn_lOtzkB6W9y1DIBDPqsHlR9AqFHC-GJ356iZg?e=6fXder) or [Ensembl ID](https://carrerasresearch-my.sharepoint.com/:u:/g/personal/plopez_carrerasresearch_org/EbeWZMBtpHpJlHxDpkWdor4BDVSUa_Z22HT6UfgjUqWNew?e=Ip58Ot). 
  Remember to change the output name inside the script.
* **3.Annotate_Peakmatrix_2.R** - as arguments it takes:
  * pk - annotated peakmatrix from the previous step
  * annotation - this will also be unique to your capture system. [Example](https://carrerasresearch-my.sharepoint.com/:t:/g/personal/plopez_carrerasresearch_org/EemaI8jhNjVCjK_aw7V3z4oBMqYL6FFJgHGL99DuPhHCkg?e=iGejla). 

## Workflow
### 4. Generate Resource Files
Use **4.GenerateResourceFiles.R**. This script generates the files neccesary to quickly compute gene scores.
Arguments:
* prefix
* cSNP_file - [coding SNPs file](https://carrerasresearch-my.sharepoint.com/:t:/g/personal/plopez_carrerasresearch_org/EZeqg2w41yJOre19s2pnHZMBoH-Pihz43eYEVdJuf_cBSA?e=65ZsT2) with all the chromosomes generated in the step 1.
* interaction_file - annotated peakmatrix generated in the previous step.
* pchic.thresh - CHiCAGO score treshold. Recommended to use 5.
* res_frag_file - Restriction map from CHiCAGO (.rmap)
* region_bed_file - LD blocks used in [GWAS Tissue Set Enrichment Analysis](https://github.com/JavierreLab/liCHiC/tree/main/6.GWAS%20imputation).
* out_dir

It will generate 3 R objects that will be used in the next step (frags.by.ld.RData, csnps.by.ld.RData, interactions.RData)

### 5. Compute Gene Score
Use **5.ComputeGeneScore.R**. This script computes basic gene scores by integrating functional and chromatin conformational data for a given GWAS.
* pmi_file - PMI file obtained during poor man's imputation (step before blockshifter)
* out_file - full path and name of output file
* csnps - generateResourceFiles.R RData object of HindIII fragments with metadata of associated ensg genes, ld blocks, type of fragment (promoter/interaction, meaning pr/oe) and fragment id (csnps.by.ld.RData)
* int - generateResourceFiles.R RData object of coding SNPs with metadata of associated ensg genes, ld blocks and HindIII fragment ids (interactions.RData)
* frags - generateResourceFiles.R RData object of ensg genes with their associated PCHi-C interactions and all extra data found in peakmatrix (frags.by.ld.RData)
* target.gene.cSNPs.only - T/F if T, only coding SNPs in target gene are included in analysis
* include.interactions.only - T/F if switch that allows us to include promoter and coding snp component(s) so we can examine contribution of tissue specific interactions to the gene score. If TRUE, it will only use the component of the SNPs located in PIRs. If you want to use the three components, use FALSE.
* decompose  - T/F switch allows us to compute geneScores for sets of tissues but also all individual tissues. Note that in this case if a set has one tissue its set name will be replaced with the tissue name so as to avoid duplication.

This script has to be run separately for each trait. As output, it will generate a table with the gene scores associated to each gene:
disease, ensg, name, biotype, strand, baitCh, all_gene_score, coding_gene_score, cell_type_specific_gene_scores, promoter_gene_score	

### 6. Reactome Pathway Analysis
First, we need a summary table containing the the gene-score associated to each gene and trait, taking as input the tables generated in the previous step. It can be generated with the script **6.1.Pre_Reactome.ipynb**. We will need the output generated in the step *3.1.1 Select only prioritized genes* (GeneScore_allTraits_prioritizedGenes_DatasetName.tsv).

Then, we can run the Reactome Pathway Analysis with the script **6.2.Reactome_Pathway_Analysis.R**.
