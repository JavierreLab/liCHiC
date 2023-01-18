## generateResourceFiles.R

## REQUIREMENTS


## FUNCTION: This script generates the files neccesary to quickly compute gene scores 

## EXAMPLE USAGE:

# Rscript generateResourceFiles.R --prefix="test_"
# --cSNP_file="../DATA/support/cSNPs_w_ENSG.e75_chr22.bed" \\ 
# --interaction_file="../DATA/chic/misfud_et_al.pm.chr22.tab" --pchic.thresh=5 \\
# --res_frag_file='../DATA/support/Digest_Human_HindIII_chr22.bed'
# --region_bed_file='../DATA/support/0.1cM.regions.b37_chr22.bed' \\ 
# out_dir='../DATA/RDATA/'

suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(GenomicRanges)))

## ARGUMENTS
args<-list(
  prefix = 'titration_cambridgeBR',
  cSNP_file = file.path("/gpfs/scratch/bsc08/bsc08246/db/COGS/coding_SNPs_VEP_1KG_GRCh38_formatted2/EUR.phase3.GRCh38.cSNPs_VEP.f.txt"),
  interaction_file = file.path("/gpfs/scratch/bsc08/bsc08246/method/peakmatrices/annot_titration_cambridgeBR_GeneNames.txt"), #use the annotated peakmatrix
  pchic.thresh = 5,
  res_frag_file = file.path("/gpfs/projects/bsc08/shared_projects/IJC_3Dchromatin/genomes/human/GRCh38/current/designDir0393561/Digest_Human_GRCh38_104_HindIII_None_15-21-17_11-11-2021.txt.rmap"), #use the rmap*
  region_bed_file = file.path("/gpfs/scratch/bsc08/bsc08246/db/LDblocks/1cM_LDblocks_GRCh38.bed"), #try with 0.1cM and with 1cM!
  out_dir = file.path("/gpfs/scratch/bsc08/bsc08246/method/COGS/4.GenerateResourceFiles/repeat_1cM_titration_cambridgeBR/")
)
print(args)
#*don't generate the res_frag_file directly from the digestion output of hicup because the ID's are not the correct ones

cs<-fread(args[['cSNP_file']],header=TRUE)

cs$start <- cs$pos
cs$end <- cs$start+1
cs$pos<-NULL# remove pos column
setcolorder(cs,c('ensg','chr','start','end')) #reorder the order of the columns
#cs$start<-cs$start+1 #I removed this

##load in the PCHIC data
int<-fread(args[['interaction_file']],header=TRUE)

##convert chicago scores to bools
for(n in names(int)[16:length(names(int))]){
	int[[n]]<-int[[n]]>args[['pchic.thresh']]
}

## load in h3 frags for doing promoter stuff
h3<-fread(args[['res_frag_file']],header=FALSE)
setnames(h3,c('chr','start','end','id'))
## get rid of y chromosome
h3<-subset(h3,chr!='Y')
#Create a GRanges object with the h3 fragments
h.gr<-with(h3,GRanges(seqnames=Rle(chr),ranges=IRanges(start=start,end=end),id=id))

#Create a GRanges object with all the other ends in the peakmatrix (no repeats)
oe.gr<-with(int,GRanges(seqnames=Rle(oeChr),ranges=IRanges(start=oeStart,end=oeEnd),id=oeID))
oe.gr<-oe.gr[!duplicated(oe.gr$id),]

#Create a GRanges object with all the baits in the peakmatrix (promoters, id=ensg:baitId) (no repeats)
proms.gr<-with(int,GRanges(seqnames=Rle(baitChr),ranges=IRanges(start=baitStart,end=baitEnd),ensg=ensg,id=paste(ensg,baitID,sep=":")))
proms.gr<-proms.gr[!duplicated(proms.gr$id),]

#Overlap between hind3 fragments and promoters and keep the baits Ids --> fm is a list of the HindII fragments' IDs that are captured (baits)
fm<-unique(subsetByOverlaps(h.gr,proms.gr)$id)
fm<-unique(c(fm-1,fm,fm+1))  #add contiguous baits Ids
## grab them back 
fm.gr<-h.gr[h.gr$id %in% fm,] #fm.gr is a GRanges object with all the promoter HindIII + continguous fragments

#merge by overlaps h3 fragments and promoters fragments from the peakmatrix
#maxgap=1L indicates that it allows a gap of 1bp between two regions to consider that they are overlapping.
#maxgap=0L indicates that the two regions are adjacent (no base pair between them). In this case, the result would be the same
#maxgap=-1L indicates that they have to be TRULLY overlapping
# So, what we are doing in "m" is to put together each promoter fragment and their adjacent fragments
m<-mergeByOverlaps(proms.gr,h.gr,maxgap=1L)
#?Check that all the seqnames (seqnames are the chromosomes, ex: 1 2 23 24) are the same (in proms.gr and in h.gr), if not: remove
rem.idx<-which(as.character(seqnames(m$proms.gr))!=as.character(seqnames(m$h.gr)))
if(length(rem.idx)>0){
	m<-m[-rem.idx,]
}
#change name of the last col of m: id --> id.1
names(m)[5] <- "id.1"	

#create data.table with the merge of promoters and h3 fragments (m) -> ensg | frag.id | orig.id (orig.id is the original fragment ID where the promoter is, the other two fragments are the adjacent ones)
g2f.prom<-with(m,data.table(ensg=ensg,frag.id=id.1,orig.id=sub("ENSG[0-9]+\\:(.*)","\\1",id)))

#create GRanges object with the promoter fragments (and their adjacents fragments) in the peakmatrix
pf.gr<-m$h.gr
pf.gr<-pf.gr[!duplicated(pf.gr$id),] #same as fm.gr 55756

af.gr<-c(pf.gr,oe.gr)  #join promoters with other ends
af.gr<-af.gr[!duplicated(af.gr$id),]


############### LD
## next we compute LD overlap
ld.blocks<-fread(args[['region_bed_file']],header=TRUE)
setnames(ld.blocks,c('chr','start','end'))
## bed files are zero based - need to be consistent throughout ¿??????????????? our files start with position 1
# ld.blocks$start<-ld.blocks$start+1
ld.gr<-with(ld.blocks,GRanges(seqnames=Rle(chr),ranges=IRanges(start=start,end=end),ld.id=1:nrow(ld.blocks))) #creatre GRanges object, add LD_ID (row number)

#Overlap between af.gr (all fragments: other ends + promoters) and ld.gr (LD blocks)
wld<-mergeByOverlaps(af.gr,ld.gr)

##frags overlapping boundary (one h3 fragment that overlap with >1 LD blocks: it is repeated >1 in wld)
bidx<-which(wld$id %in% wld[duplicated(wld$id),]$id)

sings<-wld[-bidx,] #the same as wld but without the h3 fragments that overlap with >1 LD blocks
reps<-wld[bidx,] #h3 fragments that overlap with >1 LD blocks 

upper<-which(with(reps,end(af.gr)>=end(ld.gr))) #h3 fragment ends further than the LD block
lower<-which(with(reps,start(af.gr)<=start(ld.gr))) #h3 fragment start before the LD block -> for each fragment (same af.gr), one row would be in upper and the other in lower

#Change the end coordinate from af.gr upper fragments to the end coordinate of the LD block
end(reps[upper,]$af.gr)<-end(reps[upper,]$ld.gr)
#Change the start coordinate from af.gr lower fragments to the start coordinate of the LD block
start(reps[lower,]$af.gr)<-start(reps[lower,]$ld.gr)

## put back together
wldf<-rbind(sings,reps)
frag.gr<-wldf$af.gr  #create GRanges object (id is the fragment bait id)
frag.gr$ld.id<-wldf$ld.id #add LD id

##flatten and turn into data.table
frag<-data.table(as.data.frame(frag.gr))

##add in promoter gene assignments (fragy)
fragy<-merge(frag,g2f.prom,by.x="id",by.y="frag.id")
fragy$type<-'promoter' #fragy contains all the promoter fragments and adjacent fragments, with its corresponding LDblock ID and ENSG

##add in oe assignments (fragz)
g2f.int<-int[,.(ensg,oeID)]  #table with all the other ends ID and the ENSG from the bait (not the one of the oe)
fragz<-merge(frag,g2f.int,by.x="id",by.y="oeID")
fragz$orig.id=fragz$id
fragz$type<-'interaction'

#Put together promoters and other end (including assignments) -> total of frags == fragy + fragz
frags<-rbind(fragy,fragz)

#uid = ensg:baitID:Ld_block
frags$uid<-paste(frags$ensg,frags$id,frags$ld.id,sep=":")

#order by UID
setkey(frags,uid)

#remove duplicates
frags<-unique(frags)
frags.gr<-with(frags,GRanges(seqnames=Rle(seqnames),ranges=IRanges(start=start,end=end),ensg=ensg,ld.id=ld.id,type=type,id=id))


#############################
## now sort out cSNPs

cs<-subset(cs,chr!="Y")
cs$uid<-with(cs,paste(chr,start,sep=":")) #uid -> chr:start
setkey(cs,uid)
cs.unique<-unique(cs) #remove duplicates
# cs.gr<-with(cs,GRanges(seqnames=Rle(chr),ranges=IRanges(start=start,end=end),cid=uid))
cs.gr<-with(cs.unique,GRanges(seqnames=Rle(chr),ranges=IRanges(start=start,end=end),cid=uid))

##assign LD blocks
# ol<-as.matrix(findOverlaps(cs.gr,ld.gr)) #????

csdf<-mergeByOverlaps(cs.gr,ld.gr) #merge cs and LDblocks
csdf$cs.gr$ld.id<-csdf$ld.id #add LDblocks to the GRanges object csdf$cs.gr
cs.gr<-csdf$cs.gr

lu<-data.table(cid=cs.gr$cid,ld.id=cs.gr$ld.id) #table with SNP ID and LD_block ID
fin<-merge(cs,lu,by.x="uid",by.y="cid") #merge cs and lu by uid and cid
cs.gr<-with(fin,GRanges(seqnames=Rle(chr),ranges=IRanges(start=start,end=end),ensg=ensg,ld.id=ld.id)) #Granges object of cSNPs and LDblocks

## next we need to assign these to frag.id's in our chic data
cs.gr$frag.id<-integer(length=length(cs.gr)) #add a new column for frag.id
ol<-as.matrix(findOverlaps(cs.gr,af.gr)) #find overlap between coding snps and all fragments in the pk
cs.gr[ol[,1],]$frag.id<-af.gr[ol[,2],]$id #assign fragment ID



## SAVE REQUIRED R OBJECTS

# Granges object centered in fragments 
# (each frag is classified in interaction or promoter,
# with the LDblock and the ENSG of the gene in the promoter or related to the interaction)
frag.file<-paste0(args[['out_dir']],'/',args[['prefix']],'frags.by.ld.RData')
save(frags.gr,file=frag.file) 
message(paste('Written',frag.file))

# GRanges object centered at the cSNPs, with LDblock and fragID
cs.file<-paste0(args[['out_dir']],'/',args[['prefix']],'csnps.by.ld.RData')
save(cs.gr,file=cs.file)
message(paste('Written',cs.file))

# Peakmatrix TRUE/FALSE for interactions
int.file<-paste0(args[['out_dir']],'/',args[['prefix']],'interactions.RData')
save(int,file=int.file)
message(paste('Written',int.file))
