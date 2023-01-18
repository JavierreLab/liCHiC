library(data.table)
# library(JavierreLabR) #version 1.1.4 -> to run in local
library("JavierreLabR", lib.loc="/gpfs/projects/bsc08/shared_projects/IJC_3Dchromatin/programs/Rpackages/")

## Things to change before running the script:
# (1) Arguments
# (3) Output name (line 16) depending if it is GeneName of Ensg

## ARGUMENTS
args<-list(
  pk.file = "/gpfs/scratch/bsc08/bsc08246/method/peakmatrices/titration_cambridgeBR/titration_cambridgeBR.txt", #peakmatrix
  # pk.file = "/home/plopez/scratch/method/peakmatrices/peakmatrix_test.txt",
  # annotation = "/home/plopez/scratch/db/COGS/0_AnnotatePeakmatrix/baits_coordinates_annotated_ensembl_gene_id_GRCm38_104.bed" 
  annotation = "/gpfs/scratch/bsc08/bsc08246/db/COGS/0_AnnotatePeakmatrix/baits_coordinates_annotated_gene_name_GRCm38_104.bed" #bait annotation with gene_names or ensg_name
)
out.name = str_replace(args[['pk.file']],".txt","_GeneNames.txt")
# args <- commandArgs(trailingOnly=T)

# (1) LOAD PEAKMATRIX 
pk <- load_peakmatrix(args[['pk.file']])

# (2) ANNOTATE PEAKMATRIX
annot.pk <- annotate_interactions(pk,args[['annotation']])

# (3) CONVERT TO DATAFRAME
pk.df <- as.data.frame(annot.pk)

#conserve only the desired columns
celltypes = colnames(pk.df)[26:(length(colnames(pk.df))-1)]
conserve = c(c("seqnames1","start1","end1","ID_I","P.id1","seqnames2","start2","end2","ID_II","P.id2","dist"),celltypes)
pk.df <- pk.df[,conserve]

#rename columns
names(pk.df) <- c(c("baitChr","baitStart","baitEnd","baitID","baitName","oeChr","oeStart","oeEnd","oeID","oeName","dist"),celltypes)

#replace non-annotated by NA (in baitName) or by "." (in oeName), to conserve the same notation as in the previous peakmatrices
pk.df$baitName[pk.df$baitName == "non-annotated"] <- "NA"
pk.df$oeName[pk.df$oeName == "non-annotated"] <- "."

# (4) STORE PEAKMATRIX
fwrite(pk.df,file = out.name,col.names = T,row.names = F,quote = F,sep = "\t")
