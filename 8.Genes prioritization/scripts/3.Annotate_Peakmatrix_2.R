# Annotate pchic - center it on ensg 
library(data.table)

## ARGUMENTS: 
## Paths to change manually:
## pk: path to peakmatrix
## annotation: path to bait annotation (in this case, using the version GRCh38)
## STEP 10: output path to to store the annotated peakmatrix

args<-list(
  pk = "/gpfs/scratch/bsc08/bsc08246/method/peakmatrices/titration_cambridgeBR/titration_cambridgeBR_GeneNames.txt",  #peakmatrix 
  annotation = "/gpfs/scratch/bsc08/bsc08246/db/COGS/0_AnnotatePeakmatrix/annot.pk.ensg.GRCh38.txt" #bait annotation with gene names and ensg names
)

# (0) START OF ANNOTATION HERE 
pchic <- fread(args[['pk']])
# pchic <- pchic[,-c(29:30),with=F] #remove last 2 cols (clusterID and clusterPostProb) -> we don't have these columns now, no need to do these steps
pchic[,intID:=paste0(baitID,":",oeID)] #add a new col ('intID': interaction ID with the bait Id + other end Id)

# (1) Remove all interactions with NA in baitname and "." in other end
pchic <- pchic[!(is.na(pchic$baitName) & pchic$oeName=="."),] 

# (2) Keep only bait to bait interaction to duplicate and reverse them 
b2b <- pchic[pchic$oeName!=".",] 

# reverse
setcolorder(b2b, c("oeChr","oeStart","oeEnd","oeID","oeName","baitChr","baitStart","baitEnd","baitID","baitName")) #reorder the first 10 cols
names(b2b)[1:10] <- c("baitChr","baitStart","baitEnd","baitID","baitName","oeChr","oeStart","oeEnd","oeID","oeName") #change the cols names of the first 10 cols

# join with original data table to duplicate
pchic <- rbind(pchic,b2b)

# (3) remove all interactions with NA in baitname and gene annotation in other end - cannot be removed in step one because we have to reverse them first
pchic <- pchic[!is.na(pchic$baitName),]

# (4) remove any REAL duplicates
pchic <- unique(pchic[,-c(ncol(pchic)),with=F]) 

# (5) count number of genes annotated on bait for each interaction
pchic$NumGenes <- stringr::str_count(string = pchic$baitName, ",")+1 #I have changed to ",". Before, it counted ";" but now genes names are separated by commas

# 1171109, should be 1173527 only strand biotype and ensg variables are added in (??) 806318 
sum(pchic$NumGenes) 

# (6) extract name of each gene annotation
name <- do.call(c,strsplit(pchic$baitName,",")) #I have also changed ";" to ","
#table(sort(name) == sort(olly3$name)) # all true

# (7) replicate interaction as many times as number of genes annotated on its baits
pchic<-tidyr::uncount(pchic,NumGenes) #duplicate each row as many times as its NumGenes

# (8) add name of individual gene annotation to each interaction
pchic <-cbind(name,pchic)

# (9) annotate with ensg,biotype and strand
meta<-fread(args[['annotation']])
annot.pk <- unique(dplyr::right_join(meta,pchic)) #join tables by "name" (in meta) and baitID
setcolorder(annot.pk,c("ensg","name","biotype","strand","baitChr","baitStart","baitEnd","baitID","baitName"))

# (10) write annotated peakmatrix 
fwrite(annot.pk,file = paste0("/gpfs/scratch/bsc08/bsc08246/method/peakmatrices/annot_",basename(args[['pk']])),col.names = T,row.names = F,quote = F,sep = "\t")

#check NA: annot.pk[is.na(annot.pk$ensg),]


