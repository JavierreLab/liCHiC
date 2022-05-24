## 2_autoclass2tab_javierrelab.R - @ A. Nieto 2020 - mn1
## this script takes as input the modified version of the file "*.case-data-1" from output of autoclass
## filtering interactions by their probability of belonging to cluster assigned, pmin=0.5
##
## Usage: 
## 2_autoclass2tab_javierrelab.R <tab.file> <pmin>
##
## tab.file : full path to *.tab file with the following format: $1=Case $2=Class  $3=Prob
## pmin     : probability values cutoff. Any interactions clustered with a prob < 0.5 is removed from the dataset  
##
##


args <- commandArgs(trailingOnly = T)



if(length(args)!=2){
    message("Usage: Rscript 2_autoclass2tab_javierrelab.R <tab.file> <pmin>")
    stop()
}

tab.file <- args[1]
pmin <- args[2]

## running modification step
dir <- dirname(tab.file)
cmd <- paste0("awk '/^#Case#/ { found = 1 } found { print }' ", dir, "test.case-data-1 | awk '{ print $1 \"\t\" $2 \"\t\" $3}' > ", tab.file) 
system(cmd)


## reading *.tab file, with the following format: $1=Case $2=Class  $3=Prob
df <- data.table::fread(tab.file, header=T)
df <- df[df$Prob >= pmin]


## rearrange columns and rows
df <- as.data.frame(cbind(df$`#Case#`,df$Prob,df$Class))
df <- dplyr::arrange(df, df$V3)


## write table
data.table::fwrite(df, tab.file, col.names=F, row.names=F, quote=F, sep="\t")

