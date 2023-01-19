# Histone marks enrichment

This analysis computes enrichment of different histone marks from ChIP-seq in target regions. These target regions could be of any nature, in our case we used the regions from the different networks that we obtained with liCHi-C.


## Summary of the workflow

1.  **Preparing Data**
2.  **Computing Enrichment**: [regioneReloaded](https://bioconductor.org/packages/release/bioc/html/regioneReloaded.html)


## 1. Preparing Data

The first step of the analysis is to prepare the data in the proper format to be analysed.

Basically, like the ChAs analysis, we need 2 type of data, but in this case in the same format. We need 2 set of regions to compute the enrichment of one set over the other. The first set of regions came from our liCHi-C data, these are the nodes of the network of interest, for example, the other-end of the P-OE subnetwork, or only the Promoters of these regions. The other set of regions is the regions obtained from the ChIP-seq analysis.

In the article we computed the enrichment of both the Promoter and the Other Ends, separately, classifying the nodes based on the expression of the genes in the case of Promoter, and the interacting gene in the case of Other Ends.

To prepare the data we need to transform the coordinates of both sets to GenomicRanges objects.

```R

ibed <- "path/to/interactions.ibed"

# we load the ibed file using the load_interactions function from HiCapture
interactions <- HiCaptuRe::load_interactions(ibed)

# we retrieve the coordinates of the Promoters and the Other Ends separately
coordinates_P <- as.data.frame(interactions[interactions$int == "P_OE"])[,c(1,2,3)]
coordinates_OE <- as.data.frame(interactions[interactions$int == "P_OE"])[,c(10,11,12)]

## These regions could be splitted by the expression level of the genes, to keep it simple in this example we'll use the whole set

# then we have to transform them to GenomicRanges
P_GR <- GenomicRanges::makeGRangesFromDataFrame(coordinates_P,seqnames.field = "seqnames1",start.field = "start1",end.field = "end1")
OE_GR <- GenomicRanges::makeGRangesFromDataFrame(coordinates_OE,seqnames.field = "seqnames2",start.field = "start2",end.field = "end2")

# then we create a GRanges list with names
lichic_GRs <- list(Promoters=P_GR, OtherEnds=OE_GR)

# Now we prepare the data of ChIP-seq
# we do the same, read the bed files and turn them into GenomicRanges
features <- c("path/to/H3K4me3.bed","path/to/H3K4me1.bed")
names(features) <- c("H3K4me3","H3K4me1")

# then we have to transform them to GenomicRanges
ftlist <- list()
for (i in 1:length(features) )
{
  df <- read.table(feat_files[i])
  ftlist[[i]] <- makeGRangesFromDataFrame(df,seqnames.field = "V1",start.field = "V2",end.field = "V3")
}
names(ftlist) <- names(features)
```

So we have 2 GenomicRanges list, one with the coordinates from liCHi-C and the other one with the coordinates from ChIP-seq.

## 2. Computing Enrichment

The enrichment is computed using a permutation test that is implemented in the package [regioneR](https://bioconductor.org/packages/3.16/bioc/html/regioneR.html). But for doing it simultaneously for several sets against several sets, regioneReloaded has been developed.


```R
# we load regioneReloaded and the library of our genome
library(regioneReloaded)
library(BSgenome.Hsapiens.UCSC.hg38)

# since the test need randomizations we have to set a seed
set.seed(103)

# now we performe the crosswise permutation test, using 5000 sampling of the data and 100 permutations. 
# The function that we use for randomization if "randomizeRegions" (NOTE: It can take a while, but you 
# can do it in parallel using the argument mc.cores) and the evalutation function is "numOveralps".
lichic_chip <-  crosswisePermTest(
  Alist = lichic_GRs,
  Blist = ftlist,
  min_sampling = 5000,
  ranFUN = "randomizeRegions", 
  evFUN = "numOverlaps", 
  ntimes = 100,
  genome = "hg38",
  mc.cores = 20
)

# Finally we generate the matrix that summarise this analysis
lichic_chip <- makeCrosswiseMatrix(lichic_chip)

# we can use the function to plot implemented in the package or take the matrix and plot it ourselves
plot_enrich <- plotCrosswiseMatrix(lichic_chip)
```

This is the simplest example, if we want to test tens of region sets against tens of ChIP-seq datasets, we have to prepare the corresponding GenomicRanges list and perform the test to all of them at once.

## Acknowledgments  
* [Roberto Malinverni](https://github.com/RMalinverni)  
