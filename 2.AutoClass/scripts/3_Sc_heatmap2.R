
## 3_Sc_heatmap.R - @ A Nieto 2020 - mn1
## Generating plots with results from autoclass clustering of interactions from peakmatrix of choice
## 
## Usage:
## Rscript 3_Sc_heatmap.R <peakmatrix.path> <tab.name> <cluster.order>
## 
## <peakmatrix.path> - path to peakmatrix 
## <tab.name>        - full path to *.tab, basename to match part of *.tab name --> test-[basename].tab 
## <cluster.order>   - optional - in which order to position clusters in visualisation, by "hclust","genomicdist" or "name", default being hclust



## MODIFY ARGUMENTS
#args <- commandArgs(trailingOnly = T)
args <- c('/home/plopez/scratch/method/peakmatrices/merge_original_autoclass.txt',
          '/home/plopez/scratch/method/autoclass/merge_original_autoclass_job1_2/test-merge_original_autoclass_job1_2.tab',
          'name')
options(stringsAsFactors = FALSE)

if(length(args)!=2){
  if(length(args)!=3){
    message("Usage: Rscript 3_Sc_heatmap.R  <peakmatrix.path> <tab.name> <cluster.order>")
    message("       cluster.order is one of 'hclust', 'genomicdist' or 'name', an optional argument")
    stop()
  }
}


## sourcing sven's scripts with functions needed to plot ----
source("/home/plopez/scratch/method/autoclass/scripts/autoclass_paula.R")
source("/home/plopez/scratch/method/autoclass/scripts/plot_CRM_heatmap.R")
source("/home/plopez/scratch/method/autoclass/scripts/compute_Sc.R")

message("loading libraries...")
suppressMessages(library(data.table))
suppressMessages(library(reshape))
suppressMessages(library(ggplot2))
suppressMessages(library(cluster))
suppressMessages(library(gplots))
suppressMessages(library(GenomicRanges))
suppressMessages(library(GenomicInteractions))
suppressMessages(library(stats))
suppressMessages(library(tidyverse))


## read peakmatrix with interactions that have been clustered (autoclass) ----
peakmatrix.path <- args[1]
pmBOe_4 <- fread(peakmatrix.path)

pks<-pmBOe_4


ncol_sv<-ncol(pks)


## asinh transform the data (no previous scaling line in run_autoclass_MS_Data.r script from sven) ----

pks <- pks[,c(12:ncol_sv),with=F]
pks <- asinh(pks)


## Run make_cqbm function to merge *-1.tab file and peakmatrix data => cqbm table ----
tab.name <- args[2]
dir <- dirname(tab.name)
output_name <- gsub(".tab","",tab.name)
output_name <- gsub("test-","",output_name)

## adding id column to peakmatrix
pks <- as.data.frame(pks)
ncol <- length(pks[1,])
cat("Adding an id column...\n")
pks[, "autoclassid"] <- 1:length(pks[,1])
ncol <- ncol+1
id <- ncol	

## merge autoclass clustering results with peakmatrix, which contains information of each interaction ----
# it contains the interactions that survive the filtering of cluster_prob > 0.5
cqbm <- make_cqbm(class_tab_fname = tab.name, qbm = pks, qbm_idname = names(pks)[id])



## calculate cluster specificity scores ----
## uses compute_Sc.R function sourced - Andrea N 2020 

clusters <- sort(as.integer(unique(cqbm$V3)))


i <- 0
message("calculating cluster specificity scores..\n")
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

df0$cluster <- as.integer(rownames(df0))


cluster.order <- args[3]
if(is.na(cluster.order)){
  message("no cluster order defined, using default hclust order")
  cluster.order <- "hclust"
}else if(!(cluster.order %in% c("hclust","genomicdist","name"))){
  message("Error! given cluster.order not recognised. Must be one of three choices:")
  message("    - hclust : order set by hierarchical clustering of clusters based on their mean cell type specificity scores")
  message("    - genomicdist  : order set by mean genomic distance of the interactions in each cluster")
  message("    - name   : order set by name of cluster (0-22 if 23 custers)")
  stop()
}


## calculate mean genomic distances of clusters ----
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
  
  df <- as.data.frame(cbind(cluster,mean_genomic_len))
  colnames(df) <- c("cluster","mean_distance")
  
  df1 <- merge(df0,df)
  df1 <- dplyr::arrange(df1, df1$mean_distance)
  
  cluster.order.vector <- df1$cluster 
  
  df1$cluster <- factor(df1$cluster, levels=unique(df1$cluster))
  
  x1 <- ggplot(df1, aes(cluster,mean_distance)) + geom_col()
  pdf(paste0(output_name,"_distance_plot.pdf"))
    x1
  dev.off()
  
  df0 <- df0[match(cluster.order.vector,df0$cluster),]
  
}else if(cluster.order=="name"){
  
  # cluster.order.vector <- c(1:nrow(df0)-1)
  cluster.order.vector <- c(0,1,
                            20,8,12,25,14,
                            9,
                            7,13,18,26,23,
                            5,6,11,19,31,
                            28,17,
                            32,22,
                            21,24,
                            4,29,10,
                            33,27,30,34,16,
                            15
                            )
  df0 <- df0[match(cluster.order.vector,df0$cluster),]
  
}

##
## plot Sc heatmap ---
##

message("plot Sc heatmap")

#set rownames to cluster names
rownames(df0) <- df0$cluster

# set column names and order
colnames(df0) <- gsub("_WT","",colnames(df0))

keep <- !colnames(df0)=="cluster"
heatmap_matrix <- as.matrix(df0[,keep])

## CHANGE COLUMNS ORDER

## ordering cell types based on hierarchical clustering
col.hc <- hclust(dist(t(heatmap_matrix)))
col.dd <- as.dendrogram(col.hc)

## here you can change the weight of each column (following the order of colnames(df0))
## <0 -> go to the left, and >0 -> go to the right
row.dd.reordered <- reorder(col.dd,c(1,1,-1000,-100,1,1,-200,1,1)) 

if(!exists("row.dd.reordered")){
  row.dd.reordered<-T
}

## COLOR PALETTE
my_palette <- colorRampPalette(c("#73D055","#95D840","white","#39568C"))(n = 299)


## PNG WITHOUT DENDOGRAM
png(paste0(output_name,"_Sc_heatmap_names_laureviridis.png"),width = 1000, height = 1000)
DF <- melt(heatmap_matrix)

DF$cell <- gsub("_.*","",DF$X2)
DF$cell <- factor(DF$cell,levels = c("Mon","Ery","MK","CMP","HSC","CLP","nB","nCD8","nCD4"))
DF$Y <- factor(DF$X1,levels = rev(cluster.order.vector))

minVal=min(DF$value)
maxVal=max(DF$value)
ggplot(DF, aes(x = cell, y = Y)) +
  geom_tile(aes(fill = value), color = "grey") +
  scale_fill_gradientn(colours =  c("#73D055","#95D840","white","#576b9c"),
                       values =scales::rescale(c(maxVal,maxVal/2,0,minVal)),oob = scales::squish) +
  ggpubr::theme_classic2()
dev.off()

## PDF WITHOUT DENDOGRAM SCALE
pdf(paste0(output_name,"_Sc_heatmap_names_laureviridis_scale2.pdf"),width = 20,height = 20)
DF <- melt(heatmap_matrix)

DF$cell <- gsub("_.*","",DF$X2)
DF$cell <- factor(DF$cell,levels = c("Mon","Ery","MK","CMP","HSC","CLP","nB","nCD8","nCD4"))
DF$Y <- factor(DF$X1,levels = rev(cluster.order.vector))

minVal=min(DF$value)
maxVal=max(DF$value)
ggplot(DF, aes(x = cell, y = Y)) +
  geom_tile(aes(fill = value), color = "grey") +
  scale_fill_gradientn(colours =  c("#73D055","#95D840","white","#576b9c"),
                       values =scales::rescale(c(maxVal,maxVal/2,0,minVal)),oob = scales::squish) +
  theme(legend.key.size = unit(8, 'cm'),legend.position="bottom") +
  guides(fill=guide_colorbar(ticks.colour = NA))
  # ggpubr::theme_classic2()

dev.off()


if(cluster.order=="hclust"){
  message("Clustering order is set using Hierarchical Clustering method")
  pdf(paste0(output_name,"_Sc_heatmap_hclust_colsordered.pdf"),width = 15,height = 10)
  gplots::heatmap.2(heatmap_matrix,col=my_palette, margins = c(8,7),density.info="none", trace="none",
                    cexRow=0.8,cexCol=1,sepwidth=c(0.005,0.005),colsep=1:ncol(df0),rowsep=1:nrow(df0),sepcolor="grey",#breaks=colors,
                    Rowv =TRUE,Colv = row.dd.reordered,symm=F,symkey=F,symbreaks=F,scale="none")
  dev.off()
}else if(cluster.order=="genomicdist"){
  message("Clustering order is set based on mean genomic distance of clusters")
  pdf(paste0(output_name,"_Sc_heatmap_gendist.pdf"),width = 15,height = 10)
  gplots::heatmap.2(heatmap_matrix,col=my_palette, margins = c(10,5),density.info="none", trace="none",cexRow=1.5,cexCol=1.5,sepwidth=c(0.005,0.005),colsep=1:ncol(df0),rowsep=1:nrow(df0),sepcolor="grey",#breaks=colors,
                    Rowv =F,Colv = row.dd.reordered,symm=F,symkey=F,symbreaks=,scale="none",main="Sc Heatmap",ylab = "ordered by genomic distance")
  dev.off()
}else{
  message("Clustering order is set based on names of clusters")
  pdf(paste0(output_name,"_Sc_heatmap_names_laureviridis.pdf"),width = 20,height = 20)
  gplots::heatmap.2(heatmap_matrix,col=my_palette, margins = c(10,5),density.info="none", trace="none",
                    cexRow=1.5,cexCol=1.5,sepwidth=c(0.005,0.005),colsep=1:ncol(df0),rowsep=1:nrow(df0),sepcolor="grey",#breaks=colors,
                    Rowv =F,Colv = row.dd.reordered,symm=F,symkey=F,symbreaks=T,scale="none",main="Sc Heatmap",ylab = "ordered by name")
  dev.off()
  ## png
  png(paste0(output_name,"_Sc_heatmap_names_laureviridis.png"),width = 4000,height = 4000, res= 400)
  gplots::heatmap.2(heatmap_matrix,col=my_palette, margins = c(10,5),density.info="none", trace="none",
                    cexRow=1.5,cexCol=1.5,sepwidth=c(0.005,0.005),colsep=1:ncol(df0),rowsep=1:nrow(df0),sepcolor="grey",#breaks=colors,
                    Rowv =F,Colv = row.dd.reordered,symm=F,symkey=F,symbreaks=T,scale="none",main="Sc Heatmap",ylab = "ordered by name")
  dev.off()
}
