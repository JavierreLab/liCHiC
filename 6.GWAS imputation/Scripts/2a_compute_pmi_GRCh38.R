## 2a_compute_pmi.R @ A Nieto 2021 - mn1
##
## This script calculate posterior probability (PPi) for a particular summary statistics (Trait) and performs
## Poor Man's Imputation (PMI), to associate more variants to the trait of interest. Explanation on how, below.
## 
##
## The summary statistics (SS) file must follow these pre-requisits:
## - Genomic Coordinates in GRCh37 - if not, lifOver should be used
## - Format should be | CHROM | START | STOP | rsid | pval | => 1_Format_SummayStatistics.sh (BED5)
## - chromosome X should be named "chrX" not just "X"
## - BED5file has to be bgzipped and tabix indexed:
## > sort -k1,1 -k2,2n BED5file 
## > bgzip -c BED5file
## > tabix -p bed BED5file
##
##
## Reference Panel pre-requisites:
## - ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502
## - 1000 Genomes release 20130502
## - EUR subset => 0_Download1KG.sh
## - If you have downloaded your own copy of the reference panel, you have to change the path to your VCF files in the script (line 152) -> DONE
##
##
## LD Regions pre-requisites:
## - must have column names as chr, start, end
## - must be consecutive but non-overlapping regions 
##
##
#### Usage:
#### Rscript ComputePMI.R <Chromosome> <Summary.Statistics> <Gwas.Type> <Sample.Size> <Proportion.Cases> <Out.Dir> <Tabix.Bin>
####
#### <Chromosome>         - Chromosome (1-22)
#### <Summary.Statistics> - SS with appropriate format (bed file, bgzipped, tabix indexed)
#### <Gwas.Type>          - either "QUANT" = quantitative, e.g. height, or "CC" = Cases/Control, e.g. Rheumatoid Arthritis
#### <Sample.Size>        - Number of samples included in the study
#### <Proportion.Cases>   - What proportion of the sample size were Cases ? If Gwas.Type=QUANT, Proportion.Cases=1
#### <LD.Regions>         - BED file with approximate LD blocks. Format chr | start | end, consecutive and non-overlapping regions 
#### <Out.Dir>            - Directory where to write the output table. filename will be Basename_SumStat.pmi
#### <Tabix.Bin>          - Path to tabix bin (optional) - default = /apps/HTSLIB/1.8/INTEL/bin/tabix
####
#### > Rscript 2a_compute_pmi.R 22 /path/to/summary/stat/RA.bed.bgz CC 25708 0.215 /path/to/out/directory /path/to/tabix/bin
##
## I have used "SNP" and "variant" interchangeably throughout comments
##


suppressMessages(library(snpStats))
suppressMessages(library(data.table))
suppressMessages(library(parallel))
suppressMessages(library(future))


args <- commandArgs(trailingOnly=T)


## Are all arguments input correctly?
if(length(args) != 7){
  
  # argument for tabix bin added, not using default
  if(length(args) == 8){
    message("You are changing the path to tabix_bin to ", args[8], " from default /apps/HTSLIB/1.8/INTEL/bin/tabix")
    Tabix.Bin <- args[8]
  }
  
  # stop execution, arguments not provided correctly
  message("Usage:")
  message("Rscript 2a_compute_pmi.R <Chromosome> <Summary.Statistics> <Gwas.Type> <Sample.Size> <Proportion.Cases> <LD.Regions> <Out.Dir> <Tabix.Bin = /apps/HTSLIB/1.8/INTEL/bin/tabix>")
  stop()
  
}else # if number of arguments = 6, default tabix bin must be provided
  Tabix.Bin <- "/apps/HTSLIB/1.8/INTEL/bin/tabix"

## Message to show user the number of CPUs to compute with 
message("Number of available CPUs to work with - ",  availableCores())



## Read in arguments
Chr <- as.numeric(args[1])
Gwas.Tbx <- args[2]
Gwas.Type <- args[3]
Sample.Size <- as.numeric(args[4])
Prop.Cases <- as.numeric(args[5])
LD.Regions <- args[6]
Out.Dir <- args[7]

Trait <- gsub("\\..*","",basename(Gwas.Tbx))
Out.File <- paste0(Trait , "_chr", Chr, ".pmi")

message("Treating ", Trait, " trait as ", Gwas.Type)
message("This trait has a sample size of ", Sample.Size, " and a proportion of ", Prop.Cases, " cases")
message("This is computation for chromosome ", Chr)
message("The results will be written in to ", Out.Dir, " as ", Out.File)

## Set prior probability
Prior.Prob <- 0.0001


####################################################################
## Functions used to compute posterior probabilities of variants  ##
####################################################################

## 1)
logsum <- function(x) {
  my.max <- max(x) # take out the maximum value in log form)
  my.res <- my.max + log(sum(exp(x - my.max )))
  return(my.res)
}

## 2) compute variance shrinkage for case control study
Var.data.cc <- function(f, N, s) {
  1 / (2 * N * f * (1 - f) * s * (1 - s))
}


## 3) compute variance shrinkage for quantitative trait study
Var.data <- function(f, N) {
  1 / (2 * N * f * (1 - f))
}

## 4) compute approx bayes factors and resultant posterior probabilities
## based on the assumption of one causal variant in a region/LD block
approx.bf.p <- function(p, f, type, N, s, pi_i, suffix=NULL) {
  if(type=="QUANT") {
    sd.prior <- 0.15
    V <- Var.data(f, N)
  } else {
    sd.prior <- 0.2
    V <- Var.data.cc(f, N, s)
  }
  z <- qnorm(0.5 * p, lower.tail = FALSE)
  ## Shrinkage factor: ratio of the prior variance to the total variance
  r <- sd.prior^2 / (sd.prior^2 + V)
  ## Approximate BF  # I want ln scale to compare in log natural scale with LR diff
  lABF = 0.5 * (log(1-r) + (r * z^2))
  sBF <- logsum(lABF + log(pi_i))
  ppi<- exp(lABF + log(pi_i))/(exp(sBF) + 1)
  ret <- data.frame(V,z,r,lABF,ppi)
  if(!is.null(suffix))
    colnames(ret) <- paste(colnames(ret), suffix, sep=".")
  return(ret)  
}

#######################################
## Function to Read Reference Panel  ##
#######################################
  
## Constructs an SnpMatrix of variants in a genomic region chrA:start-end on chromosome Chr A
vcf2sm <- function(region,Chr){
  
  # contructing path to 1KG VCF file of specific chromosome
  path.to.vcf <- paste0("/gpfs/scratch/bsc08/bsc08246/db/1000Genomes_GRCh38/Filtered_EUR/EUR.chr",Chr,".phase3.GRCh38.vcf.gz")

  # get header - obtain names of samples
  my.pipe<-pipe(paste0("tabix -H ",path.to.vcf, " ",region))
  header<-tail(scan(my.pipe,what=character(),sep="\n",quiet=TRUE),n=1)
  close(my.pipe)
  cnames<-unlist(strsplit(header,"\t"))
  
  # read in all variants of reference panel within region with "PASS" 
  tmp<-as.data.frame(fread(cmd=paste0(Tabix.Bin, " ", path.to.vcf," ",region," | grep PASS"),sep="\t",header=FALSE,stringsAsFactors=FALSE))
  colnames(tmp)<-cnames
  
  # remove med indels 
  idx<-which(nchar(tmp$ALT)>=8 | nchar(tmp$REF) >=8)
  
  if(length(idx)>0)
    tmp<-tmp[-idx,]
  
  # gt = genotypes
  gt<-tmp[,10:ncol(tmp)]
  if(nrow(gt)==0)
    return(NA)
  
  # info of SNPs
  info<-tmp[,1:9]
  
  # contructing SnpMatrix
  sm<-apply(gt,1,function(x) sub("0\\|0","1",x))
  sm<-apply(sm,1,function(x) sub("(0\\|1)|(1\\|0)","2",x))
  sm<-apply(sm,1,function(x) sub("1\\|1","3",x))
  # set anything else to a missing value
  sm<-t(apply(sm,1,function(x) as.raw(sub("[0-9]\\|[0-9]","0",x))))
  colnames(sm)<-make.names(info$ID,unique=TRUE)
  rownames(sm)<-colnames(gt)
  sm<-new("SnpMatrix", sm)
  
  return(list(map=info,gt=sm,r=region))
}


######################################################################
## Function to Merge Reference Panel and Summary Statistic Variants ##
######################################################################

## Filter variants from reference panel and merge data with trait-associated variants for a region chrA:start-end on Chr A
get1KGInfo<-function(region,Chr,Gwas.Tbx,verbose=TRUE,MAF.thresh=0.01,HWE.thresh=25,CALLRATE.thresh=0.95){
  if(verbose)
    message(paste0("Processing region ",region))
  
  # get 1KG SNPs for this region
  smt<-vcf2sm(region=region, Chr=Chr)
  
  # if not variants found in region, return NA
  if(any(is.na(smt)))
    return(NA)
  
  # extracting info on MAF, CALLRATE, HWE of variants
  sm.gt<-col.summary(smt$gt)
  
  
  # filtering by HWE, MAF and CALLRATE thresholds 
  if(!missing(HWE.thresh) & !missing(MAF.thresh) & !missing(CALLRATE.thresh))
    idx<-unique(c(with(sm.gt,which(MAF < MAF.thresh | z.HWE^2 > HWE.thresh | Call.rate < CALLRATE.thresh))))
  
  # filtering any duplicated variants, based on position only as we are working on one chromosome only
  idx<-unique(c(idx,which(duplicated(smt$map$POS))))
  
  # applying the filtering 
  if(length(idx)>0){
    sm<-list(map=smt$map[-idx,],gt=smt$gt[,-idx])
  }else{
    sm<-list(map=smt$map,gt=smt$gt)
  }
  
  # creating pipe to read in gwas data using tabix
  message("Reading GWAS file")
  tabix_pipe<-pipe(paste0(Tabix.Bin," ", Gwas.Tbx," chr",region))
  
  # reading variants in SS only within defined region
  # return NA if no variants found
  tmp.bed<-tryCatch(read.delim(tabix_pipe,sep="\t",header=FALSE,stringsAsFactors=FALSE),
                    error=function(e){
                      message(paste(region,'no variants found'))
                      return(NA)
                    })
  # return NA if no variants found
  if(is.na(tmp.bed)){
    return(NA)
  }
  if(any(is.na(tmp.bed[,1:3])))
    return(NA)
  
  # adding column names to SS table 
  names(tmp.bed)<-c('chr','start','end','rs','pval')
  
  # filtering any duplicated variants in SS, based on position only as we are working on one chromosome only
  tmp.bed<-tmp.bed[!duplicated(tmp.bed$start),]
  
  # add pvals to the reference panel variants based on GWAS - merge by genomic position
  # all.x=TRUE, keep whole table of reference panel, even if after merge we still have not associated a pval to it - will be used for PMI
  sm$map<-merge(sm$map,tmp.bed[,c('start','pval')],by.x='POS',by.y='start',all.x=TRUE)
  
  # checking that at least one variant has been associated a pval
  idx<-which(!is.na(sm$map$pval))
  
  # for sparser genotyping some regions only have no variants we cannot
  # impute these so bail !
  if(length(idx)==0)
    return(NA)
  
  # to return the sm object, reference panel and SS has to have at least one variant in region, even after filtering
  # also, after merge, at least one variant has to have a pval associated to it
  return(sm)
}

###################################################################
## Perform Poor Man's Imputation and Compute Posterior Prability ##
###################################################################

## Performs PMI taking into account LD and R-squared thresholds. 
## Computes PPi using wakefield's implementation
computePMI<-function(region,Chr,Gwas.Tbx,verbose=T,MAF.thresh=0.01,HWE.thresh=25,CALLRATE.thresh=0.95,R2.thresh=0.6){
  
  # obtain table of variants, some of them with associated pvals (trait association)
  sm<-get1KGInfo(region,Chr,Gwas.Tbx,verbose=T,MAF.thresh=MAF.thresh,HWE.thresh=HWE.thresh,CALLRATE.thresh=CALLRATE.thresh)
  
  # if something wrong with sm, bail
  if(any(is.na(sm)))
    return(sm)
  
  # identifying which variants have associated pval
  idx<-which(!is.na(sm$map$pval))
  
  # for sparser genotyping some regions only have no variants we cannot
  # impute these so bail !
  if(length(idx)==0)
    return(NA)
  
  # at this point we have all info needed to compute PMI
  message("starting PMI")
  
  # extract genotype info
  gt<-sm$gt
  
  # defining variant by their chromosomal position
  colnames(gt)<-sm$map$POS
  
  # calculate measures of linkage disequilibrium between pairs of SNPs (see ?snpStats::ld)
  # SNPs with pval VS SNPs without pval
  tld<-ld(gt[,idx],gt[,-idx],stats="R.squared")
  
  
  # select the max row for each column
  
  # initialise fb matrix (each row represents a possible variant to impute, 4 columns)
  fb<-matrix(NA,nrow=ncol(tld),ncol=4)
  
  # column 1 - "name" of variant to impute, i.e. their chromosomal position
  fb[,1]<-as.numeric(colnames(tld))
  if(nrow(tld)==1){
    # in this case they must match rownames as only 1 SNP with associated pval!
    fb[,2]<-as.numeric(rownames(tld))
    fb[,3]<-as.numeric(tld)
  }else{
    # more than one variant has an associated pval
    
    # column 2 - "name" of variant to impute with
    # choose, for each row/SNP with unknown pval, which SNP with known pval it is in highest LD with (i.e. with max R-squared) 
    fb[,2]<-as.numeric(rownames(tld)[apply(tld,2,which.max)]) 
    
    # column 3 - extract R-squared statistic of LD between SNPs from column 1 & 2
    fb[,3]<-as.numeric(apply(tld,2,max)) 
  }
  
  # any possible imputation below cut off drop on the floor
  # if minimum R-squared is not achieved, the two variants are not similar enough and therefore you cannot predict pval with the other, so imputation is not possible
  fb<-matrix(fb[fb[,3]>=R2.thresh,],ncol=4)
  
  # remove any NaNs
  fb<-fb[fb[,3]!="NaN",]
  
  
  if(!is.matrix(fb)){
    fb<-rbind(fb,rep(NA,4))
    message("Not matrix adding NA's")
  }
  
  # column 4 - adding newly associated pval
  # imputing variants with pvals from SNPs in high enough LD 
  fb[,4]<-sm$map[match(fb[,2],sm$map$POS),]$pval
  
  
  # finally fill in pvals into full sm matrix
  sm$map$pval<-apply(cbind(sm$map$pval,fb[match(sm$map$POS,fb[,1]),4]),1,function(x){ 
    tmp<-x[!is.na(x)]
    ifelse(length(tmp)==1,tmp,NA)
  })	
  
  # add genomic position of which variant was used to impute 
  sm$map$imp.snp.pos<-apply(cbind(sm$map$POS,fb[match(sm$map$POS,fb[,1]),2]),1,function(x){
    tmp<-x[!is.na(x)]
    ifelse(length(tmp)>1,tmp[2],NA)
  })
  
  # add information on R-squared statistics of LD measures
  sm$map$imp.r2<-apply(cbind(sm$map$POS,fb[match(sm$map$POS,fb[,1]),3]),1,function(x){	
    tmp<-x[!is.na(x)]
    ifelse(length(tmp)>1,tmp[2],NA)
  })
  
  ## CURRENTLY PMI METHOD EFFECTIVELY REMOVES SNPS THAT ARE IN ORIGINAL GWAS BUT DON'T  (*)
  ## TAG ANY OTHER SNPS (WITH THRESHOLDS). THIS CAUSES DOWN SAMPLING WHICH MIGHT BE 
  ## A REQUIRED BEHAVIOUR. 
  
  # obtain MAF values for variants
  sm$map$maf<-col.summary(sm$gt)$MAF
  
  # what is the percentage imputation
  kidx<-which(!is.na(sm$map$pval))
  
  # keep is a table with the SNPs that have associated pvals (imputed and original from SS taking (*) into account)
  keep<-sm$map[kidx,]
  
  # message to show how well we did 
  perc<-(nrow(keep)/nrow(sm$map)) * 100
  message(paste(floor(perc),'% imputed for',nrow(keep),'/',nrow(sm$map),region))
  
  # preparing final data table to be written out
  outf<-data.table(keep)
  
  # change colname for chromosome
  setnames(outf,'#CHROM','chr')
  
  # computing PPi using Wakefield implementation
  outf$ppi<-with(outf,signif(approx.bf.p(as.numeric(pval),maf,Gwas.Type,Sample.Size,Prop.Cases,Prior.Prob)$ppi,digits=3))
  
  
  # create output file
  ret<-outf[,.(POS,chr,ID,pval,maf,ppi,imp.snp.pos,imp.r2)]
  ret$start<-ret$POS # start = POS
  ret$POS <- ret$POS+1 # end = POS+1
  setnames(ret,c('end','chr','rsid','pval','maf','ppi','imp.snp.pos','imp.r2','start'))
  setcolorder(ret,c('chr','start','end','rsid','maf','pval','ppi','imp.snp.pos','imp.r2'))
  ret
}

computePPI<-function(region,Chr,Gwas.Tbx,verbose=TRUE){
  sm<-get1KGInfo(region,Chr,Gwas.Tbx,TRUE)

  if(any(is.na(sm)))
    return(sm)
  idx<-which(!is.na(sm$map$pval))
  ## for sparser genotyping some regions only have no variants we cannot
  ## impute these so bail !
  if(length(idx)==0)
    return(NA)
  ## next filter to remove those w/o pvals
  sm$map<-sm$map[idx,]
  sm$gt<-sm$gt[,idx]
  sm$map$maf<-col.summary(sm$gt)$MAF
  ## next compute ppi using implementation of Wakfield.
  sm$map$ppi<-signif(approx.bf.p(as.numeric(sm$map$pval),sm$map$maf,gwas_type,n_samples,prop_cases,pi_i)$ppi,digits = 3)
  ret<-data.table(sm$map)
  setnames(ret,c('start','chr','rsid','ref','alt','qual','filter','info','format','pval','maf','ppi'))
  ret<-ret[,.(start,chr,rsid,pval,maf,ppi)]
  ret$end<-ret$start + 1
  setcolorder(ret,c('chr','start','end','rsid','maf','pval','ppi'))
  message(paste('Got',nrow(ret)))
  ret
}


# Reading Approx. LD blocks (regions) file 
regions <- fread(LD.Regions)
regions <- paste0(regions$chr,":",regions$start,"-",regions$end)

# keeping regions only for chromosome Chr
regions <- regions[grepl(paste0("^",Chr,":"),regions)]

message("Ready to process ", length(regions)," number of regions for chromosome ", Chr, ", trait ", Trait)


# Perform PMI and compute PPi for each region
# Performed in parallel, depending on number of available cores.
pmi <- mclapply(regions,computePMI,Chr=Chr,Gwas.Tbx=Gwas.Tbx,MAF.thresh=0.01,HWE.thresh=25,CALLRATE.thresh=0.95,R2.thresh=0.6,mc.cores=availableCores())
message("pmi step complete")
pmi
results <- pmi

# Joining tables for all regions - i.e. construct one table for the whole chromosome Chr
# Remove any regions which have no variants to include in results table
no.res<-sapply(results,is.data.table)
no.cover<-sum(!no.res)

message(paste(no.cover,ifelse(no.cover==1,"region is","regions are"),"not covered"))
if(length(no.res)>0)
  results<-results[no.res]

results<-results[sapply(results,is.data.table)]
results<-rbindlist(results)

# write results in chosen output directory 
write.table(results,file=paste0(Out.Dir, Out.File),row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)
message("Results written to ", Out.Dir, Out.File)




