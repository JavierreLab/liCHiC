# GWAS Genes prioritization
Capture Hi-C Omnibus Gene Score (COGS) algorithm allows the prioritization of relevant genes for each GWAS trait and cell type

## Credits
The method and software was co-developed by [Chris Wallace](http://chr1swallace.github.io/) and [Olly Burren](http://ollyburren.github.io/), presented in [Javierre et al., 2016](https://doi.org/10.1016/j.cell.2016.09.037).

The scripts in this repository have been adapted from the original ones from [CHIGP](https://github.com/ollyburren/CHIGP).


Briefly, this method considers linkage disequilibrium to estimate the posterior probability of each SNP being casual for each trait. Then these SNPs are used to compute a gene-score for all the genes involved in liCHi-C significant interactions in at least one cell type. This gene-score is composed by 3 components: coding SNPs annotated by VEP (104); SNPs located in promoter regions; and SNPs overlapping Other-Ends.

## Prerequisites
* R (developed with 3.2.2)
* data.table
* snpStats (bioconductor)
GenomicRanges (bioconductor)
reshape2
yaml
data.tree



## Workflow
