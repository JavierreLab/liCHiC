
## 3_interaction_heatmap.R - @ A Nieto 2020 - mn1
## Generating plots with results from autoclass clustering of interactions from peakmatrix of choice - heatmap of interactions, not included in 3_interaction_heatmap.R 
## 
##
## Usage:
## Rscript 3_interaction_heatmap.R <peakmatrix.path> <tab.name> <factor.ds> <cluster.order> <exclude.clusters> 
##
## <peakmatrix.path> - path to peakmatrix 
## <tab.name>        - directory of results , basename to match part of *.tab name --> test-[basename].tab 
## <factor.ds>       - value from 0-1, factor of downsampling to apply to dataset for faster visualisation (default is 0.5)
## <cluster.order>   - in which order to position clusters in visualisation, by "hclust","genomicdist" or "name"
## <exclude.clusters> - optional, vector of clusters to be excluded, given in format "0,1,2"


args <- commandArgs(trailingOnly = T)
options(stringsAsFactors = FALSE)

if(length(args)!=4){
  if(length(args)!=5){
    message("Usage: Rscript 3_interaction_heatmap.R  <peakmatrix.path> <tab.name> <factor.ds> <cluster.order> <exclude.clusters>")
    message("       cluster.order is one of 'hclust', 'genomicdist' or 'name'")
    message("       exclude.clusters is a vector of clusters to be excluded from visualisation, given in format '0,1,2'")
    stop()
  }
}


args <- commandArgs(trailingOnly = T)

## sourcing sven's scripts and my compute_Sc.R script
source("/gpfs/scratch/bsc08/bsc08246/method/autoclass/scripts/autoclass.R")
source("/gpfs/scratch/bsc08/bsc08246/method/autoclass/scripts/plot_CRM_heatmap.R")
source("/gpfs/scratch/bsc08/bsc08246/method/autoclass/scripts/compute_Sc.R")


message("loading libraries...")
suppressMessages(library(data.table))
suppressMessages(library(reshape))
suppressMessages(library(ggplot2))
suppressMessages(library(cluster))
suppressMessages(library(gplots))
suppressMessages(library(GenomicRanges))
suppressMessages(library(GenomicInteractions))
suppressMessages(library(ggplot2))
suppressMessages(library(stats))
suppressMessages(library(tidyverse))
suppressMessages(library(viridis))

## read args ----
peakmatrix.path <- args[1]
tab.name <- args[2]
factor.ds <- as.numeric(args[3])
if(factor.ds > 1 | factor.ds < 0 | is.na(factor.ds)){
  message("factor.ds is not within limits. Must be a value between 0 and 1")
  message("Automatically setting to 0.5")
  factor.ds <- 0.5
}

cluster.order <- args[4] 
if(!(cluster.order %in% c("hclust","genomicdist","name","manual"))){
  message("cluster.order is not recognised. Must be one of 'hclust', 'genomicdist', 'name' or 'manual'")
  message("Automatically setting to 'hclust")
  cluster.order <- "hclust"
}

exclude.clusters <- args[5]
if(!(is.na(exclude.clusters))){
  exclude.clusters<-as.numeric(strsplit(args[5],",")[[1]])
}



## read peakmatrix with interactions that have been clustered ----
pmBOe_4 <- fread(peakmatrix.path)

pks<-pmBOe_4

ncol_sv<-ncol(pks)


## asinh transform the data (no previous scaling line in run_autoclass_MS_Data.r script from sven) ----

pks <- pks[,c(12:ncol_sv),with=F]
pks <- asinh(pks)


## Run make_cqbm function to merge *-1.tab file and peakmatrix data => cqbm table ----
dir <- dirname(tab.name)
output_name <- gsub(".tab","",tab.name)
output_name <- gsub("test-","",output_name)


## adding ID column to peakmatrix ----
pks <- as.data.frame(pks)
ncol <- length(pks[1,])
cat("Adding an id column...\n")
pks[, "autoclassid"] <- 1:length(pks[,1])
ncol <- ncol+1
id <- ncol	

## merge results from autoclass clustering (tab file) and info of interactions (peakmatrix) by ID  ----
cqbm <- make_cqbm(class_tab_fname = tab.name, qbm = pks, qbm_idname = names(pks)[id])



## order clusters based on cluster.order ----
if(cluster.order=="genomicdist"){
  message("you chose to order your clusters based on mean genomic distance")
  
  coordinates <- pmBOe_4[,c(1:3,6:8)]
  coordinates <- tibble::add_column(coordinates,c(1:nrow(coordinates)))
  colnames(coordinates)[7] <- "interactionID"
  
  region1 <- makeGRangesFromDataFrame(coordinates[,1:3], seqnames.field = "baitChr",start.field = "baitStart", end.field = "baitEnd")
  region2 <- makeGRangesFromDataFrame(coordinates[,4:7], seqnames.field = "oeChr",start.field = "oeStart", end.field = "oeEnd")
  
  coords <- GenomicInteractions(region1,region2)
  
  dists <- calculateDistances(coords, method = "midpoint")
  coordinates <- tibble::add_column(coordinates, dists)
  
  clusters <- dplyr::as_tibble(cqbm[,c(1,3)]) 
  colnames(clusters) <- c("ID", "cluster")
  
  inter_lengths <- dplyr::inner_join(coordinates,clusters,by = c("interactionID" = "ID"))
  inter_lengths <- dplyr::group_by(inter_lengths,cluster)
  inter_lengths <- dplyr::arrange(inter_lengths,cluster)
  
  mean_genomic_len <- dplyr::group_map(inter_lengths[,8:9], ~ colMeans(.x, na.rm = TRUE))
  mean_genomic_len <- unlist(mean_genomic_len, use.names = FALSE)
  
  cluster <- as.integer(unique(inter_lengths$cluster))
  
  df0 <- as.data.frame(cbind(cluster,mean_genomic_len))
  colnames(df0) <- c("cluster","mean_distance")
  
  df0 <- dplyr::arrange(df0, df0$mean_distance)
  cluster.order.vector<- rev(df0$cluster)
  rownames(df0) <- df0$cluster
  
  cqbm <- cqbm[order(match(cqbm$V3,cluster.order.vector)),]
  
  
}else if(cluster.order=="name"){ #creo que no funciona, porque usa df0
  
  message("you chose to order your clusters based on name")
  
  cluster.order.vector <- c(1:nrow(df0)-1)
  cqbm <- cqbm[order(match(cqbm$V3,cluster.order.vector)),]
  
}else if(cluster.order=="manual"){
  
  message("you chose to order your clusters and columns manually. If you want to change their order, please modify the lines 156 and 158 in the script.")
  
  cluster.order.vector <- c('0','1','20','8','12','25','14','9','7','13','18','26','23','5','6','11','19','31','28','17','32','22','21','24','4','29','10','33','27','30','34','16','15') #cluster order
  cqbm <- cqbm[order(match(cqbm$V3,cluster.order.vector)),]
  cqbm <- cqbm[c('V1','V2','V3','Mon_500K','Ery_500K','MK_500K',"CMP_WT_merge_45","HSC_WT_merge_15","CLP_WT_merge_45",'nB_1M','nCD8_500Kmix','nCD4_500K')] #column order
  
}else if(cluster.order=="hclust"){
  
  message("you chose to order your clusters based on hclust of cell-type & cluster specificity scores")
  clusters <- sort(as.integer(unique(cqbm$V3)))
  
  
  i <- 0
  for(cl in clusters){                                                                      
    message(paste("Cluster",cl,"of a total of",length(clusters),"being processed..."))      
    
    rows <- which(cqbm$V3==cl, arr.ind = T)
    if(length(rows)>1){
      
      s <- specificity_score(cqbm, subset=rows, cellcols=c(4:ncol(cqbm)),name=cl)                 
      
      if(i==0){                                                                               
        df0 <- s                                                                               
      }                                                                                       
      else{                                                                                   
        df0 <- rbind(df0,s)                                                                     
      }                                                                                       
      
      i <- i+1      
    }
  }                   
  
  
  # at this point we have a matrix of cell_types by clusters showing:                        
  # Specificity Score of each cell type for each cluster
  
  for(col in 1:ncol(df0)){
    df0[,col] <- as.double(df0[,col])
    
  }
  
  df0.mat <- as.matrix(df0[,-c(ncol(df0))])
  
  hc <- hclust(dist(df0.mat),method = "complete")
  dd <- as.dendrogram(hc)
  
  cluster.order.vector<-labels(dd)
  cqbm <- cqbm[order(match(cqbm$V3,cluster.order.vector)),]
  
}

colnames(cqbm) <- gsub("_WT","",colnames(cqbm))

## downsample dataset for visualisation ----
message("factor to downsample by is: ",factor.ds)

## This is done because without subsetting, plotting of heatmap in too computationally costly. By randomly sampling the clusters, we hope to obtain the most accurate representation
## of the real heatmap

for(cl in unique(cqbm$V3)){
  print(cl)
  sbst <- cqbm[cqbm$V3==cl, ]
  
  num_sbst <- round(nrow(sbst)*factor.ds)
  set.seed(103) # should be 103, changed to see what happens with visualisation
  
  if(nrow(sbst)>num_sbst){
  sbst <- sbst[sample(1:nrow(sbst), num_sbst), ]
  }

  if(cl==unique(cqbm$V3)[1]){
    
    newtbl <- sbst
    
  }else{ newtbl <- rbind(newtbl,sbst) }

}

## plotting heatmap of interactions ----
if(cluster.order=="hclust"){
  pdf(paste0(output_name,"_interactions_heatmap.pdf"))
      plot_CRM_heatmap(cqbm = newtbl,tfrange = 4:ncol(newtbl),
                       margin=c(8,2),
                       mycol = rainbow(length(unique(cqbm$V3))),
                       clustercols=F,excludeclasses = exclude.clusters,
                       separateRows=NULL)
  dev.off()
  
  pdf(paste0(output_name,"_interactions_heatmap_average.pdf"))
  plot_CRM_heatmap(cqbm = newtbl,tfrange = 4:ncol(newtbl),
                   margin=c(8,2),
                   mycol = rainbow(length(unique(cqbm$V3))),
                   clustercols=F,excludeclasses = exclude.clusters,clustmethod = "average",
                   separateRows=NULL)
  dev.off()
  
}else if (cluster.order=="manual"){
  ## store PDF
  pdf(paste0(output_name,"_interactions_heatmap_manualorder.pdf"))
  plot_CRM_heatmap(cqbm = newtbl,tfrange = 4:ncol(newtbl),
                   margin=c(8,2),
                   mycol = rainbow(length(unique(cqbm$V3))),
                   clustercols=F,excludeclasses = exclude.clusters,
                   separateRows=NULL, 
                   hmcol="viridis",
                   colclustmethod="complete")
  dev.off()
  
  pdf(paste0(output_name,"_interactions_heatmap_manualorder_clusterCol.pdf"))
  plot_CRM_heatmap(cqbm = newtbl,tfrange = 4:ncol(newtbl),
                   margin=c(8,2),
                   mycol = rainbow(length(unique(cqbm$V3))),
                   clustercols=T,excludeclasses = exclude.clusters,
                   separateRows=NULL, 
                   hmcol="viridis",
                   colclustmethod="complete")
  dev.off()
  
  
}else{
  pdf(paste0(output_name,"_interactions_heatmap.pdf"))
  plot_CRM_heatmap(cqbm = newtbl,tfrange = 4:ncol(newtbl),
                   margin=c(8,2),
                   mycol = rainbow(length(unique(cqbm$V3))),
                   clustercols=F,excludeclasses = exclude.clusters)
  dev.off()
  
}

