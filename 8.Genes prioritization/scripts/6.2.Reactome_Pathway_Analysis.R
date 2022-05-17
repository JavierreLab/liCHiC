## Reactome Pathway Analysis

## Load libraries
#BiocManager::install("ReactomePA")
library(ReactomePA)
library(clusterProfiler)
library(data.table)
library(org.Hs.eg.db)
library(httr)
library(data.table)
library(tidyverse)
library(ggplot2)

## ARGUMENTS
dataset <- "bloodprecursors" #"onlyblood" "onlycambridge" "onlyprecursors" "bloodprecursors"
biotype <- "PC" #PC or All
universe_set <- "baits" #prior or baits
outpath <- "/home/plopez/scratch/method/COGS/6.2.ReactomePA_repeat/"

# args <- commandArgs(trailingOnly = T)
# dataset <- args[1]
# biotype <- args[2] #PC or All
# universe_set <- args[3] #prior or baits

## INPUT FILE
tab <- fread(paste0(outpath,"tables_prioritizedGenes/GeneScore_allTraits_prioritizedGenes_",dataset,".tsv"), header=T)

## (1) UNIVERSE GENE SET
## Select universe gene set and convert to entrezIDs
if (universe_set == "prior"){
  if (biotype == "PC"){
    universe <- tab$ensg[tab$biotype == "protein_coding"] #3189 take protein coding
  } else if (biotype == "All"){
    universe <- tab$ensg # all biotypes of genes
  }
  universe <- universe[!is.na(universe)]
  universe_entrez <- bitr(universe, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  universe_entrez_genes <- universe_entrez$ENTREZID #3122
  
} else if (universe_set == "baits"){
  universe_path = "/home/plopez/scratch/db/COGS/0_AnnotatePeakmatrix/annot.pk.ensg.GRCh38.txt"
  # Universe genes
  df.uni = fread(universe_path,header=TRUE)
  
  if (biotype == "PC"){
    uni.genes <- df.uni$ensg[df.uni$biotype == 'protein_coding'] #21359 only protein coding
  } else if (biotype == "All"){
    # uni.genes <- df.uni$ensg #32647 all biotypes
  }
  
  universe_entrez <- bitr(uni.genes, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db") #17297
  universe_entrez_genes <- as.character(unique(universe_entrez$ENTREZID)) #17261
}


## (2) REACTOME
count = 0
for (i in 6:length(names(tab))) { #change the range when changing the number of traits evaluated 
  
  print(names(tab)[i])

  if (biotype == "PC"){
    ensg <- tab$ensg[(tab$biotype == "protein_coding") & (tab[[i]] >= 0.5)] #only protein coding
  } else if (biotype == "All"){
    ensg <- tab$ensg[(tab[[i]] >= 0.5)] # all biotypes of genes
  }
 
  ensg <- ensg[!is.na(ensg)]
  
  ## Skip traits without genes prioritized
  if (length(ensg) == 0){
    next
  }
  print(length(ensg))
  genes_entrez <- bitr(ensg, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db") #17297
  genes_entrez <- as.character(unique(genes_entrez$ENTREZID))

  reactome_enrichment <- enrichPathway(gene = genes_entrez,
                                       organism = "human",
                                       pAdjustMethod = "fdr",
                                       universe = universe_entrez_genes, #UNIVERSE
                                       pvalueCutoff = 0.05,
                                       readable=TRUE)
  
  ## Save table
  # print(paste0(outpath,dataset,"_universe",biotype,universe_set,"/universe",biotype,universe_set,"_reactome_",names(tab)[i],"_table.tsv"))
  write.table(reactome_enrichment,paste0(outpath,dataset,"_universe",biotype,universe_set,"/universe",biotype,universe_set,"_reactome_",names(tab)[i],"_table.tsv"),
  row.names = FALSE, sep="\t", quote=FALSE)
  
  
  # Create summary table and save plot
  if (dim(reactome_enrichment)[1] !=0){
    print("enriched")

    if (count == 0){
      df <- reactome_enrichment[,]
      df$trait <- names(tab)[i]
      df$num_genes <- length(ensg)
      count = 1
    } else {
      df1 <- reactome_enrichment[,]
      df1$trait <- names(tab)[i]
      df1$num_genes <- length(ensg)
      df <- rbind(df,df1)
    }

    pdf(paste0(outpath,dataset,"_universe",biotype,universe_set,"/universe",biotype,universe_set,"_reactome_",names(tab)[i],"_dotplot.pdf"))
    for (i in 1:1) {
      print(dotplot(reactome_enrichment))
    }
    dev.off()
  }
}

## Store summary table with all the enrichments
df <-df[c("trait","ID","Description","GeneRatio","BgRatio","pvalue","p.adjust","qvalue","geneID","Count","num_genes")]
# write.table(df,paste0(outpath,dataset,"_universe",biotype,universe_set,"_reactome_allTraits_table.tsv"),
#   row.names = FALSE, sep="\t", quote=FALSE)



## (3) PLOT ALL TRAITS
laure <- df
laure$geneR <- as.numeric(unlist(lapply(laure$GeneRatio, function(x) str_split(x,"/")[[1]][1])))/as.numeric(unlist(lapply(laure$GeneRatio, function(x) str_split(x,"/")[[1]][2])))

to.remove = c('INS','HDL','UC')
laure = laure[!(laure$trait %in% to.remove),]

## Change order axis
level_order <- c('SBP','GLC','LDL','TC','TG',
                 'T2D', 'CD',
                 'HL_ukbb',  
                 'NHL_ukbb', 'PLHMN_ukbb','LL_gwas','LL_ukbb',
                 'MC_gwas','MPC', 'PC')
laure$trait <- factor(laure$trait, levels=level_order)
laure <- laure[order(laure$trait), ]


laure$Description <- factor(laure$Description,level=rev(unique(laure$Description)))

ggplot(laure,aes(x=trait,y=Description)) + 
  geom_point(aes(col=p.adjust,size=geneR)) + 
  scale_color_viridis_c(direction = -1) +
  ggpubr::theme_classic2()  +
  labs(title = paste0(dataset,": ",biotype," ",universe_set)) +
  theme(axis.text.x=element_text(size=rel(1.4),angle=90),
        axis.text.y=element_text(size=rel(0.9)),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()
  ) +
  guides(size=guide_legend(title="Gene Ratio"),
         color=guide_colorbar()) 

output_name = '/home/plopez/scratch/method/COGS/6.2.ReactomePA_repeat/plots/'
pdf(paste0(output_name,"Bloodprecursors_reactome_2.pdf"),width = 10,height = 10)
dev.off()



## Add secondary axis

laure$Description <- factor(laure$Description,level=rev(unique(laure$Description)))

laure2 <- as_tibble(laure) %>% group_by(trait) %>% summarise(n_genes=c(unique(num_genes),rep(NA,n()-1)))
laure3 <- laure
laure3$n_genes <- laure2$n_genes

pdf(paste0(output_name,"Bloodprecursors_reactome_2_withcounts.pdf"),width = 11,height = 10)
ggplot(laure3,aes(x=trait,y=Description,label=n_genes)) + 
  geom_point(aes(col=p.adjust,size=geneR)) + 
  scale_color_viridis_c(direction = -1) +
  # scale_colour_gradient(low = "#7d89c5ff", high = "#2c355eff") +
  ggpubr::theme_classic2()  +
  theme(axis.text.x=element_text(size=rel(1.4),angle=90,hjust = 0.95,vjust = 0.5),
        axis.text.y=element_text(size=rel(0.9)),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 0.5)
  ) +
  guides(size=guide_legend(title="Gene Ratio"),
         color=guide_colorbar()) +
  geom_text(y = 58, # Set the position of the text to always be at '14.25'
            hjust = 0.5,
            size = 2,
            color='black') +
  coord_cartesian(# This focuses the x-axis on the range of interest
                  clip = 'off')+
  theme(plot.margin = unit(c(1,3,1,1), "lines"))

dev.off()






