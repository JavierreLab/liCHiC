                                                                                                       
## 1_format_ss.R @ A. Nieto 2021 & P. Lopez
##
## This script will format your summary statistics to BED5 format
## If coordinates are in build37, liftover to build38 will be performed -lifted file will be found in same directory as original ss
## If no coordinates are in file, they will be annotated
##
## Input : GWAS catalog f or UKBB summary statistic file
##         If format is not correct (because ss not from GWAS catalog not from UKBB, format manually before running)
##
## Output: BED file with 5 columns: chr | start | end | rsid | pval
##


options(scipen=999)                                                                                                                                                                    
suppressMessages(library(data.table))
suppressMessages(library(stringr))
suppressMessages(library(parallel))
suppressMessages(library(future))

# ARGS: SS file, liftover, genomic coordinates, out_dir
args <- commandArgs(trailingOnly=T)
# args <- c("/home/plopez/Scripts/GWAS/test/data/summary_statistics/original/RheumatoidArthritis-chr22.Build37.f.tsv","F","F","/home/plopez/Scripts/GWAS/test/data/summary_statistics/formatted/")
# args <- c("/home/plopez/Data/GWAS/summary_statistics_GRCh38/original/UKBB/LymphoidLeukaemia_UKBB_20001_1048-gwas.imputed_v3.both_sexes.Build37.tsv","T","F","/home/plopez/Data/GWAS/summary_statistics_GRCh38/test_format_ss/")

ss <- args[1]


Do.Liftover <- as.logical(args[2])
Add.GenCoords <- as.logical(args[3]) #Don't used anymore
Out.Dir <- args[4]

# checking arguments
if(length(args) != 4){
  
  # stop execution, arguments not provided correctly
  message("Usage:")
  message("Rscript 1_format_ss.R <Summary.Statistic> <Do.Liftover> <Add.GenCoords> <Out.Dir>")
  stop()
}

trait <- gsub("-.*","",basename(ss)) #obtain trait from the basename of the SS

message("Reading in SS...")
Sum.Stat <- fread(ss)


# identifying format of file
if(sum(c("variant","minor_allele","minor_AF","n_complete_samples","AC","ytx","beta","se","tstat","pval")
   %in% colnames(Sum.Stat))==10){
  
  message("ukbb format identified")
  Format <- "ukbb"
  
}else if(sum(c("chromosome","base_pair_location","variant_id","effect_allele","p_value","beta","standard_error")
            %in% colnames(Sum.Stat))==7){
  
  message("gwas catalog format identified")
  Format <- "catalog"
  
}else{
  
  message("no format identified. Please enter a valid format. Aborting mission!")
  stop()
  
}


################ FUNCTIONS ################################################################
###########################################################################################

# function which performs liftover from GRCh37 to GRCh38 
liftover.function <- function(sum.stat,out.dir){
  message("Chosen to perform liftover of chromosomal coordinates from Build37 to Build38")
  
  # create a temporary directory
  tmp.dir <- paste0(Out.Dir,"tmp/")
  system(paste0("mkdir ",tmp.dir))
  
  # preparing table to be lifted
  to.liftover <- data.table(chr = paste0("chr",sum.stat$chromosome), start = sum.stat$base_pair_location,end = as.integer(as.numeric(sum.stat$base_pair_location)+1), id = sum.stat$variant_id) #data.table with chr | start | end | id (bed file)
  to.liftover.file <- paste0(tmp.dir,"to_liftover.", basename(ss))
  fwrite(to.liftover,to.liftover.file,col.names = F,row.names = F,quote = F,sep = "\t")
  
  
  lifted.file <- gsub("uild37","uild38.lifted",ss)
  unmapped.file <- paste0(tmp.dir,"unmapped.", basename(ss))
  message("Executing liftover")
  
  # liftover 
  cmd <- paste0("/gpfs/scratch/bsc08/bsc08246/method/Gwas_TissueSetEnrich/scripts/LiftOver/liftOver ", to.liftover.file, " /gpfs/scratch/bsc08/bsc08246/method/Gwas_TissueSetEnrich/scripts/LiftOver/hg19ToHg38.over.chain.gz ", lifted.file, " ", unmapped.file)
  system(cmd)
  
  lifted.table <- fread(lifted.file,fill = T,na.strings = "")
  
  system(paste0("rm -r ",tmp.dir))
  
  return(lifted.table)
}

######## Functions to add coordinates: DON'T USED ANYMORE in the version GRCh37 to GRCh38
# auxiliary function to add coordinates to a gwas catalog formatted file without coordinates
#input: ss and 1KG.reference (chrom, position, rsID)
# rsid2coords <- function(chr,sum.stat){                                                                                                                                                   
#   
#   ref.table <-fread(paste0("/home/plopez/Scripts/GWAS/1000Genomes/rsid2coords/1KG.reference.chr",chr,".txt"))                                         
#   message("Adding coordinates for chromosome ",chr)
#   
#   
#   temp.ss <- merge(sum.stat,ref.table[,1:3,with=F],by.x="variant_id",by.y="rsid",all.x = T)                                                                                            
#   
#   temp.ss <- temp.ss[!is.na(temp.ss$CHROM),]
#   
#   temp.ss$chromosome[!is.na(temp.ss$CHROM)] <- temp.ss$CHROM[!is.na(temp.ss$CHROM)]                                                                                                    
#   temp.ss$base_pair_location[!is.na(temp.ss$POS)] <- temp.ss$POS[!is.na(temp.ss$POS)]                                                                                                  
#   
#   temp.ss$CHROM <- NULL                                                                                                                                                                
#   temp.ss$POS <- NULL       
#   
#   
#   
#   return(temp.ss)                                                                                                                                                                      
#   
# }

# # function to add coordinates to a gwas catalog file
# add.coords.function <- function(sum.stat){
#   res <- mclapply(1:22,rsid2coords,sum.stat=sum.stat,mc.cores = availableCores()-1)
#   res.table <- do.call("rbind",res)
#   return(res.table)
# }


# function to format a ukbb summary statistics file to BED5 format
format.ukbb.function <- function(add.coords,liftover,sum.stat,out.dir=Out.Dir){
  
  message("starting formatting...")
  ## ukbb format
  
  # extracting genomic coordinates from variant name
  sum.stat <- sum.stat[, c("chromosome", "base_pair_location","allele1","allele2") := tstrsplit(variant, ":", fixed=TRUE)]
  sum.stat <- sum.stat[,.(chromosome,base_pair_location,variant,allele1,allele2,pval,beta,se)] #remove unnecessary columns
  colnames(sum.stat) <- c("chromosome","base_pair_location","variant_id","other_allele","effect_allele","p_value","beta","standard_error") #rename columns
  
  if(liftover==TRUE){
    message("Performing liftover")
    # removing any variant without coordinates
    sum.stat <- sum.stat[!is.na(sum.stat$base_pair_location),]
    sum.stat <- sum.stat[!is.na(sum.stat$variant_id),]
    sum.stat <- sum.stat[!is.na(sum.stat$p_value),]
    
    # liftover
    lifted.table<-liftover.function(sum.stat=sum.stat,out.dir=out.dir) #table with 4 cols: chr | start | end | rsid
    names(lifted.table) <- c("chr","start","end","rsid")
    
    message("Liftover was ", signif(nrow(lifted.table)/nrow(sum.stat)*100, digits = 4), "% successful")
    
    # formatting 
    bedfile <- merge(lifted.table, sum.stat[,variant_id,p_value], by.x = "rsid", by.y="variant_id")
    setcolorder(bedfile, c("chr","start","end","rsid","p_value"))
    colnames(bedfile) <- c("chr","start","end","rsid","pval")
    
  }
  
  return(bedfile)
}


# function to format a gwas catalog summary statistics file to BED5 format
format.catalog.function <- function(add.coords,liftover,sum.stat,out.dir=Out.Dir){
  ## gwas catalog f format  
  if(liftover==TRUE){
    message("Performing liftover")
    # removing any variant without coordinates
    sum.stat <- sum.stat[!is.na(sum.stat$base_pair_location),]
    sum.stat <- sum.stat[!is.na(sum.stat$variant_id),]
    sum.stat <- sum.stat[!is.na(sum.stat$p_value),]
    
    # liftover
    lifted.table<-liftover.function(sum.stat=sum.stat,out.dir=out.dir)
    names(lifted.table) <- c("chr","start","end","rsid")
    
    message("Liftover was ", signif(nrow(lifted.table)/nrow(sum.stat)*100, digits = 4), "% successful")
    
    # formatting 
    bedfile <- merge(lifted.table, sum.stat[,variant_id,p_value], by.x = "rsid", by.y="variant_id")
    setcolorder(bedfile, c("chr","start","end","rsid","p_value"))
    colnames(bedfile) <- c("chr","start","end","rsid","pval")
    
  }else{
    
    if(add.coords==TRUE){
      message("Adding coordinates")
      # adding coordinates
      res.table <- add.coords.function(sum.stat=sum.stat)
      
      message("Adding genomic coordinates was ", signif((nrow(res.table)/nrow(sum.stat))*100,digits = 3),"% successful")
      sum.stat<-res.table
    }
    
    message("Formatting")
    ## formatting after adding genomic coordinates or after nothing (!liftover & !add.coords)
    
    # extracting columns with information on chr | pos | rsid | pval
    idx <- which(colnames(sum.stat)==c("chromosome") | colnames(sum.stat)==c("base_pair_location") | colnames(sum.stat)==c("variant_id") | colnames(sum.stat)==c("p_value"))
    bedfile <- sum.stat[,idx,with=F]
    
    # ensuring order of columns is as wanted
    setcolorder(bedfile,c("chromosome","base_pair_location","variant_id","p_value"))
    
    # adding start-stop columns instead of only position column - needed for BED5 format  
    bedfile <- bedfile[,c(1,2,2,3,4)]
    
    # renaming columns
    colnames(bedfile) <- c("chr","start","end","rsid","pval")
    bedfile$end <- as.integer(bedfile$end)+1  
    bedfile <- bedfile[!(is.na(bedfile$end)),]
    
    bedfile$chr <- paste0("chr",bedfile$chr)
    bedfile <- bedfile[order(chr,start),]
  }
  return(bedfile)
}

###########################################################################################
###########################################################################################


message(Add.GenCoords, "= add gen coords")
message(Do.Liftover, "= do liftover")

# calling functions 
if(Format=="ukbb"){
  
  Bed.Res <- format.ukbb.function(add.coords=Add.GenCoords,liftover=Do.Liftover,sum.stat=Sum.Stat)
  
}else if(Format=="catalog"){
  
  Bed.Res <- format.catalog.function(add.coords=Add.GenCoords,liftover=Do.Liftover,sum.stat=Sum.Stat)
  
}else{
  message("Wrong format. Aborting mission!")
}


Out.File <- paste0(Out.Dir,gsub("-.*",".bed",basename(ss)))
message("Writing file to ",Out.File)

fwrite(x = Bed.Res, file = Out.File, col.names = F, row.names = F, quote = F, sep = "\t",na = "NA")


#borrar
# add.coords <- Add.GenCoords
# liftover <- Do.Liftover
# sum.stat <- Sum.Stat

