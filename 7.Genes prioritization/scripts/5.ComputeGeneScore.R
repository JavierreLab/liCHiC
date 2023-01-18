## computeGeneScore.R 

##' FUNCTION: This script computes basic gene scores by integrating functional and chromatin conformational data
##'           for a given GWAS


##' EXPLANATION OF ARGUMENTS:
##' pmi_file                   - PMI file obtained during poor man's imputation (step before blockshifter)
##' out_file                   - full path and name of output file
##' csnps                      - generateResourceFiles.R RData object of HindIII fragments with metadata of associated ensg genes, ld blocks, type of fragment (promoter/interaction, meaning pr/oe) and fragment id 
##' int                        - generateResourceFiles.R RData object of coding SNPs with metadata of associated ensg genes, ld blocks and HindIII fragment ids
##' frags                      - generateResourceFiles.R RData object of ensg genes with their associated PCHi-C interactions and all extra data found in peakmatrix
##' target.gene.cSNPs.only     - T/F if T, only coding SNPs in target gene are included in analysis
##' include.interactions.only  - T/F if switch that allows us to include promoter and coding snp component(s) so we can examine contribution of tissue specific interactions to the gene score *
##' decompose                  - T/F switch allows us to compute geneScores for sets of tissues but also all indivdual tissues 
##'                                  note that in this case if a set has one tissue its set name will be replaced with the tissue name so as to avoid duplication.

##' RESULTS:
##' table with gene scores associated to each gene for a specific disease
##' disease,ensg,name,biotype,strand,baitCh,all_gene_score,coding_gene_score,cell_type_specific_gene_scores,promoter_gene_score		

## * If include.interactions.only == TRUE -> it will only use the component of the SNPs in PIRs. If we want to use the three components, use FALSE	


suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(GenomicRanges)))
suppressWarnings(suppressMessages(library(reshape2)))

## Give as argument the trait abbreviation
input <- commandArgs(trailingOnly = T)
trait = input[1]
# trait = "RA"

# args<-list(
#   pmi_file = file.path(paste0("/home/plopez/scratch/method/Gwas_TissueSetEnrich/posterior_probability/1cM/",trait,"/",trait,".pmi")),
#   out_file = file.path(paste0("/home/plopez/scratch/method/COGS/results/1cM/",trait,"_geneScores.tab")),
#   csnps = file.path("/home/plopez/scratch/method/COGS/4.GenerateResourceFiles/1cM/liCHic_GRCh38_csnps.by.ld.RData"),
#   int = file.path("/home/plopez/scratch/method/COGS/4.GenerateResourceFiles/1cM/liCHic_GRCh38_interactions.RData"),
#   frags = file.path("/home/plopez/scratch/method/COGS/4.GenerateResourceFiles/1cM/liCHic_GRCh38_frags.by.ld.RData"),
#   ## if this is set to true then only coding SNPs in target gene are 
#   ## included. If set to false then coding SNPs in interactions in genes
#   ## other than target are counted.
#   target.gene.cSNPs.only=FALSE,
#   include.interactions.only = FALSE,
#   decompose = FALSE
# )

args<-list(
  pmi_file = file.path(paste0("/gpfs/scratch/bsc08/bsc08246/method/Gwas_TissueSetEnrich/posterior_probability/1cM/",trait,"/",trait,".pmi")),
  out_file = file.path(paste0("/gpfs/scratch/bsc08/bsc08246/method/COGS/results/repeat_1cM_titration_cambridgeBR/",trait,"_geneScores_0.1cM_onlycambridge.tab")),
  csnps = file.path("/gpfs/scratch/bsc08/bsc08246/method/COGS/4.GenerateResourceFiles/repeat_1cM_titration_cambridgeBR/titration_cambridgeBRcsnps.by.ld.RData"),
  int = file.path("/gpfs/scratch/bsc08/bsc08246/method/COGS/4.GenerateResourceFiles/repeat_1cM_titration_cambridgeBR/titration_cambridgeBRinteractions.RData"),
  frags = file.path("/gpfs/scratch/bsc08/bsc08246/method/COGS/4.GenerateResourceFiles/repeat_1cM_titration_cambridgeBR/titration_cambridgeBRfrags.by.ld.RData"),
  ## if this is set to true then only coding SNPs in target gene are
  ## included. If set to false then coding SNPs in interactions in genes
  ## other than target are counted.
  target.gene.cSNPs.only=FALSE,
  include.interactions.only = FALSE,
  decompose = FALSE
)

print(args)
# if(!interactive()){
#   args<-getArgs(verbose=TRUE,numeric=c('target.gene.cSNPs.only'))
# }



if(sum(names(args)=='sets')==0)
  args[['sets']]<-''

if(is.null(args[['include.interactions.only']]))
  args[['include.interactions.only']]<- FALSE

if(is.null(args[['decompose']]))
  args[['decompose']]<- FALSE


####### LOADING DATA

#Load the annotated peakmatrix 
ints<-get(load(args[['int']]))
## grab tissue names
tmp.tnames<-names(ints)[16:length(names(ints))]

#Load h3 fragments (all promoter fragments + adjacent fragments + interaction fragments; all annotated with ENSG and LDblockID)
frags.gr<-get(load(args[['frags']]))

#Load coding SNPS (chr, pos, ENSG, LDblockId, fragmentID)
cs.gr<-get(load(args[['csnps']]))

options(stringsAsFactors=FALSE) #https://community.rstudio.com/t/what-does-stringsasfactors-in-r-mean/35626

#Load posterior probabilities in a GRanges object
pmi.file=args[['pmi_file']]
disease<-sub("\\.pmi$","",basename(pmi.file))
test.pmi<-fread(pmi.file,header=TRUE)
#setnames(test.pmi,c('chr','start','end','rs','r2','i.pos','maf','pval','ppi','delta'))
pmi.gr<-with(test.pmi,GRanges(seqnames=Rle(chr),ranges=IRanges(start=end,end=end),ppi=ppi,rs=rsid)) #!!!!

# IMPORTANT: when loading pmi.gr, we have stablished the start coord as the end (so, it's the SNP coord +1 )
# This is correct because when doing the overlap with cs.gr, it has two coordenates and it overlaps with
# the end coord (start+1)
################# 

####### REMOVE MHC REGION
## always remove MHC as this buggers things up
# They removed a region of 10Mb. Now, when using GRCh38, we don't need to modify its coordenates, as MHC are still in this 10Mb region
# https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh38.p13
mhc.gr<-GRanges(seqnames=Rle(6),ranges=IRanges(start=25e6,end=35e6))
ol<-as.matrix(findOverlaps(pmi.gr,mhc.gr))
if(length(ol)>0)
  pmi.gr<-pmi.gr[-ol[,1],]
###################################

# overlap SNPs in the posterior prob file with coding SNPs
cs.pmi.gr<-mergeByOverlaps(pmi.gr,cs.gr)

# #######
# test<-data.table(as.data.frame(cs.gr))
# p = 
# test.gr<-with(test,GRanges(seqnames=Rle(seqnames),ranges=IRanges(start=start,end=start),ensg=ensg,ld.id=ld.id,frag.id=frag.id))
# cs.pmi.gr.test<-mergeByOverlaps(pmi.gr,test.gr)
# #######

## add gene and frag info to pmi.gr
tmp<-cs.pmi.gr$pmi.gr #this is a GRanges object inside cs.pmi.gr
mcols(tmp)<-cbind(mcols(tmp),cs.pmi.gr[,c('ensg','ld.id','frag.id')])

## in this the coding SNPs that spoilt the party potentially are those
## where frag.id is !=0 -> These are SNPs in regions without fragments (checked rmap)
cs.pmi.gr<-tmp #here the SNP coord is still +1
to.adj<-subset(cs.pmi.gr,frag.id!=0)  #remove snps without fragment ID -> it would be the second component

## here we want to make sure we don't double count.
to.adj$ppi<- -to.adj$ppi #negative
#collapse snps that fall in the same ensg, ld.id and frag.id, and sum their ppi (sppi)
to.adj<-with(to.adj,data.table(ppi=ppi,ensg=ensg,ld.id=ld.id,frag.id=frag.id)) #convert to a table
to.adj<-to.adj[,list(sppi=sum(ppi)),by="ensg,ld.id,frag.id"]

#collapse snps with the same ensg and ldbloc and sum their ppi. Here we don't remove snps without fragID -< first component
cs<-with(cs.pmi.gr,data.table(ppi=ppi,ensg=ensg,ld.id=ld.id,frag.id=frag.id)) 
cs<-cs[,list(sppi=sum(ppi)),by="ensg,ld.id"]



##### remove all coding SNPs this is an option
if(args[['target.gene.cSNPs.only']]){
  ol<-as.matrix(findOverlaps(pmi.gr,cs.gr)) #grab index of cSNPs
  pmi.nocs.gr<-pmi.gr[-ol[,1],] #remove them
  pmi.frags<-mergeByOverlaps(pmi.nocs.gr,frags.gr)
}else{
  pmi.frags<-mergeByOverlaps(pmi.gr,frags.gr)
}

pmi.frags<-with(pmi.frags,data.table(rs=rs,ppi=ppi,ensg=ensg,ld.id,type=type,frag.id=id))


# collapse snps by fragID, LDblock, gene and type (promoter or interaction)
noncoding<-pmi.frags[,list(sppi=sum(ppi)),by="frag.id,ld.id,ensg,type"]    
#conserve only promoter
noncoding.prom<-subset(noncoding,type=="promoter")
## this does not change between selections
noncoding.prom<-noncoding.prom[,list(sppi=sum(sppi)),by="frag.id,ensg,ld.id"]  

#store only interactions fragments (these are other ends)
noncoding.interactions<-subset(noncoding,type=="interaction")
noncoding.interactions$sid<-with(noncoding.interactions,paste(ensg,frag.id,sep=":"))
ints$sid<-with(ints,paste(ensg,oeID,sep=":")) #generate ID (ensg:fragID)
setkey(noncoding.interactions,sid)

## we can generate one of these for each set of tissues
if(!file.exists(args[['sets']])){
  tissues<-split(tmp.tnames,tmp.tnames)
  #tissues[['all']]<-names(ints)[16:32]
  #tissues<-c(tissues,'all')
  ## decompose switch allows us to compute geneScores for sets of tissues but also 
  ## all indivdual tissues note that in this case if a set has one tissue its set name 
  ## will be replaced with the tissue name so as to avoid duplication.
}else if (args[['decompose']]){
  sets<-get(load(args[['sets']]))
  ## remove single tissue sets
  sets<-sets[sapply(sets,length)>1]
  ## allow ease of selecting a sets 
  names(sets)<-paste('set',names(sets),sep=".")
  tissues<-c(sets,split(tmp.tnames,tmp.tnames))
}else{
  tissues<-get(load(args[['sets']]))
  names(tissues)<-paste('set',names(tissues),sep=".")
}

## esto no lo entiendo: creo que simplemente cambia el nombre de la columna, para indicar que solo estamos mirando interaccions (3rd component)
## si no tuvieramos esta opción marcada, estaríamos mirando los 3 componentes, no solo el tercero
if(args[['include.interactions.only']])
  names(tissues)<-paste(names(tissues),'interactions_only',sep="_") #if FALSE, it doesn't do anything

#add a tissue entry with all the tissues
tissues[['all']]<-tmp.tnames

#add ID (ensg:fragID:Ldblock) to the list of fragments with coding SNPs
to.adj$uid<-with(to.adj,paste(ensg,frag.id,ld.id,sep=":"))

gint<-do.call("rbind",lapply(seq_along(tissues),function(i){
  t<-names(tissues)[i]
  print(paste("Processing",t))
  
  ## grab gene and oeID back
  sids<-ints[which(rowSums(ints[,tissues[[i]],with=FALSE])>0),]$sid #grab interactions significant for the tissue
  idx<-which(noncoding.interactions$sid %in% sids)
  
  #if(t!='all'){
  #	idx<-which(noncoding.interactions$sid %in% ints[ints[[t]],]$sid)
  #}else{
  #	idx<-1:nrow(noncoding.interactions)
  #}
  
  snoncoding.interactions<-noncoding.interactions[idx,list(sppi=sum(sppi)),by="frag.id,ensg,ld.id"] #grab other ends of the significant interactions of the trait
  
  ## allow a switch that allows us to promoter component so we can examine 
  ## contribution of tissue specific interactions to the gene score
  if(t == 'all' | !args[['include.interactions.only']]){
    # Here we are including the SECOND component (promoter SNPs)
    all.genes<-rbind(noncoding.prom,snoncoding.interactions) #keep the interactions of the tissue + promoter fragments with SNPS (ncSNps if target==TRUE)
  }else{
    all.genes<-snoncoding.interactions
  }
  
  all.genes$uid<-with(all.genes,paste(ensg,frag.id,ld.id,sep=":"))
  #to.adj has a list of fragments with cSNPs, and the sppi is the sum of the ppi of the coding snps
  #in the same gene, LDblock and fragment
  # It is the SECOND component (SNPs im promoter regions)
  setcolorder(to.adj,names(all.genes)) 
  
  #all.genes have all the SNPs in the PIRs, and the sppi is the sum of the ppi of the noncoding SNPs
  #in the same gene (gene of the interaction -the promoter- not the other end), LDblock and fragment
  # It is the third component (SNPs in PIRs)
  all.genes$uid<-with(all.genes,paste(ensg,frag.id,ld.id,sep=":"))
  
  if(!args[['target.gene.cSNPs.only']]){
    #select those fragments with cSNPs, and add them to all.genes
    to.adj.m<-subset(to.adj,uid %in% all.genes$uid) 
    all.genes<-rbind(all.genes,to.adj.m)
  }
  
  #Take the repeated rows of gene-Ldblock-fragment and sum their ppi
  all.genes<-all.genes[,list(sppi=sum(sppi,na.rm=TRUE)),by="frag.id,ld.id,ensg"]
  
  #Collapse the fragments in the same gene and LDblock and sum their ppi
  all.genes<-all.genes[,list(sppi=sum(sppi,na.rm=TRUE)),by="ld.id,ensg"]
  setcolorder(cs,names(all.genes))
  
  ## allow a switch that allows us to coding snp component so we can examine 
  ## contribution of tissue specific interactions to the gene score
  if(t == 'all' | !args[['include.interactions.only']]){
    total<-rbind(cs,all.genes) #join all.genes and all the cSNPs (collapsed by ensg-LDblock) -> cs is the FIRST component (cSNPs in the annotated gene)
  }else{
    total<-all.genes
  }
  
  # BLOCKSCORE: sum the previous two components with the first one to obtain the blockscore
  ld.score<-total[,list(block_score=sum(sppi,na.rm=TRUE)),by="ensg,ld.id"]
  
  # GENESCORE
  gs<-ld.score[,list(gene_score=1-prod(1-block_score)),by="ensg"]
  gs$tissue<-t
  gs
}))

## now add in scores for coding snps and promoter.snps
## We are computing the FIRST and SECOND components separataly
## They are independent of the tissues

## promoters
noncoding.prom$uid<-with(noncoding.prom,paste(ensg,frag.id,ld.id,sep=":")) #it has the SNPs in pmi.gr, both coding and non-coding
if(!args[['target.gene.cSNPs.only']]){
  to.adj.m<-subset(to.adj,uid %in% noncoding.prom$uid) #select only coding SNPs ¿?
  noncoding.prom.cor<-rbind(noncoding.prom,to.adj.m)
}else{
  noncoding.prom.cor<-noncoding.prom
}

noncoding.prom.cor<-noncoding.prom.cor[,list(sppi=sum(sppi,na.rm=TRUE)),by="frag.id,ensg,ld.id"]
noncoding.prom.cor<-noncoding.prom.cor[,list(block_score=sum(sppi,na.rm=TRUE)),by="ensg,ld.id"]
noncoding.prom.cor<-noncoding.prom.cor[,list(gene_score=1-prod(1-block_score)),by="ensg"]
noncoding.prom.cor$tissue<-'promoter'
fi<-rbind(gint,noncoding.prom.cor) #add it to the general table

##coding snps
coding<-cs[,list(gene_score=1-prod(1-sppi)),by="ensg"]
coding$tissue<-'coding'
fi<-rbind(fi,coding)  #add it to the general table

foo<-melt(fi,id.vars=c("ensg","tissue")) #reshape the table
results<-data.table(dcast(foo,ensg~tissue+variable),key="ensg")

## add in details (information about genes)
details<-ints[,.(ensg,name,biotype,strand,baitChr)]
setkey(details,ensg)
details<-unique(details)
r<-results[details]

missing<-subset(r,is.na(all_gene_score)) # rows with all NA
merged<-subset(r,!is.na(all_gene_score)) # the rest
merged$disease<-disease

setcolorder(merged,c('disease',names(details),names(results)[names(results)!="ensg"]))
write.table(merged,file=args[['out_file']],sep="\t",row.names=FALSE,quote=FALSE)
message(paste("Written",args[['out_file']]))

