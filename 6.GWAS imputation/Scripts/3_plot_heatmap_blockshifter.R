## Script to plot gwas blockshifter results

## Load libraries
library(data.table)
library(gplots)
library(ggplot2)

## ARGUMENTS
ld = "1"
Res.Dir <- paste0("/home/plopez/scratch/method/Gwas_TissueSetEnrich/merge4GWAS_precursors_",ld,"cM/results/tissue_set_enrichment/")
traits.names <- c(
  'SBP','DBP','FNBD',  'DEP_ukbb',  'AUT_ukbb','LDL','BMI',  'GLC', 'TC','TG','T2D',
  'PBC', 'T1D', 'RA','SLE', 'MS',  'CD', 'CEL',
  'BC_gwas', 'LeuC_gwas', 'CD4regTC_gwas','nCD8_gwas', 
  'HL_ukbb', 'LEU_ukbb','NHL_ukbb','NHL_gwas','PLHMN_ukbb', 'LL_C91_ukbb','CLL_ukbb',   'LL_gwas', 'LL_ukbb', 
  'MonC_gwas','MC_gwas','MWCC_gwas','MCHC','MCH',  'MCV',  'MPC',   'PC',   'RBCC'
)
target <- c("HSC_WT_merge_15","CLP_WT_merge_45","CMP_WT_merge_45",
            "nCD4_500K","nCD4_cambridge",
            "nCD8_500Kmix","nCD8_cambridge",
            "nB_1M","nB_cambridge",
            "Mon_500K","Mon_cambridge",
            "MK_500K","MK_cambridge",
            "Ery_500K","Ery_cambridge"
)
##


files <- list.files(Res.Dir,full.names = T,pattern=".txt")

df <- data.table::fread(files[1])
gwas.name <- gsub(".txt","",basename(files))[1]
df <- df[,z,test]
names(df)[ncol(df)] <- as.character(gwas.name)

for(f in files[2:length(files)]){
  tmp <- data.table::fread(f)
  gwas.name <- gsub(".txt","",basename(f))[1]
  tmp <- tmp[,z,test]
  df <-merge(df,tmp)
  names(df)[ncol(df)] <- gwas.name
  
}


## Select traits and cell types
df <- df[match(target, df$test),]

rwn <-  df$test #cell types
df$test <- NULL # remove column with the cell type name

#plot only the selected traits
df <- subset(df,select=traits.names)
#

cln <- names(df)#traits names
df <- as.matrix(df) 
rownames(df) <- rwn
colnames(df) <- cln

#change the name of the traits
abrev = read.table("/home/plopez/scratch/method/Gwas_TissueSetEnrich/scripts/gwas_datasets_complete_abbreviations_modified.csv",
                   header = FALSE, sep = "\t")
for(i in 1:53) {                                        # Head of for-loop
  colnames(df)[colnames(df) == abrev[i,1]] <- abrev[i,2]  # Code block
}


## PLOT HEATMAP
output_name <- "/home/plopez/Data/liCHiC/Gwas_TissueSetEnrich/merge4GWAS_precursors_1cM/"
png(paste0(output_name,"_Sc_heatmap_definitivecolors3.png"),width = 1000, height = 650)

DF <- melt(df)
names(DF) <- c("cell","trait","value")
# DF$cell <- gsub("_.*","",DF$X2)
# DF$cell <- factor(DF$cell,levels = c("Mon","Ery","MK","CMP","HSC","CLP","nB","nCD8","nCD4"))
DF$cell <- factor(DF$cell,levels = rev(levels(DF$cell)))

minVal=min(DF$value)
maxVal=max(DF$value)
ggplot(DF, aes(x = trait, y = cell)) +
  geom_tile(aes(fill = value), color = "grey") +
  scale_fill_gradientn(colours =  c("#73D055","#95D840","white","#576b9c"),
                       values =scales::rescale(c(maxVal,maxVal/2,0,minVal)),oob = scales::squish) +
  ggpubr::theme_classic2() + 
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=1))


dev.off()

## PDF WITHOUT DENDOGRAM
pdf(paste0(output_name,"_Sc_heatmap_definitivecolors3.pdf"),width = 40,height = 20)

DF <- melt(df)
names(DF) <- c("cell","trait","value")
# DF$cell <- gsub("_.*","",DF$X2)
# DF$cell <- factor(DF$cell,levels = c("Mon","Ery","MK","CMP","HSC","CLP","nB","nCD8","nCD4"))
DF$cell <- factor(DF$cell,levels = rev(levels(DF$cell)))

minVal=min(DF$value)
maxVal=max(DF$value)
ggplot(DF, aes(x = trait, y = cell)) +
  geom_tile(aes(fill = value), color = "grey") +
  scale_fill_gradientn(colours =  c("#73D055","#95D840","white","#576b9c"),
                       values =scales::rescale(c(maxVal,maxVal/2,0,minVal)),oob = scales::squish) +
  ggpubr::theme_classic2() + 
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=1))


dev.off()


## PDF SCALE WITHOUT DENDOGRAM
pdf(paste0(output_name,"_Sc_heatmap_definitivecolors3_SCALE.pdf"),width = 40,height = 20)

DF <- melt(df)
names(DF) <- c("cell","trait","value")
# DF$cell <- gsub("_.*","",DF$X2)
# DF$cell <- factor(DF$cell,levels = c("Mon","Ery","MK","CMP","HSC","CLP","nB","nCD8","nCD4"))
DF$cell <- factor(DF$cell,levels = rev(levels(DF$cell)))

minVal=min(DF$value)
maxVal=max(DF$value)
ggplot(DF, aes(x = trait, y = cell)) +
  geom_tile(aes(fill = value), color = "grey") +
  scale_fill_gradientn(colours =  c("#73D055","#95D840","white","#576b9c"),
                       values =scales::rescale(c(maxVal,maxVal/2,0,minVal)),oob = scales::squish) +
  ggpubr::theme_classic2() + 
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=1),
        legend.key.size = unit(8, 'cm'),legend.position="bottom") +
  guides(fill=guide_colorbar(ticks.colour = NA)) 


dev.off()








