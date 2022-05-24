##########################################################################################################################################################################################################
# Javierre Lab @ A Nieto-Aliseda - mn1
# Oct 2020
#  
# This script prepares the files needed as input for autoclass clustering method
# Also prepares the shell script with all the commands to run all prepared "versions" of the autoclass clustering with different samples of 30000 interactions in parallel
# This includes: db2, hd2, model, s-params, r-params
# Explanation of files found in README in my drive "clustering_interactions"
# Uses scripts from sven stewiz (babraham institute) 
#
#
# Usage:
#
# Rscript 1_launch_autoclass.R <work.dir> <project.name> <peakmatrix.path> <seeds> <prop.train>
#
# <work.dir>         - directory where we want to write the prep-files into
# <project.name>    - name of project, which will be the same for all seeds (usually "train")
# <peakmatrix.path> - path to peakmatrix with chicago score information
# <seeds>           - when sampling 30,000 interactions for "training" subset, seed must be used for reproducibility. Must be a comma-separated list of maimum 11 seeds (eg: 103,81,5)
# <prop.train>      - proportion of interactions we want to used for "training" subset, rather than the default 30,000 interactions
# 
##########################################################################################################################################################################################################		


## TO CHANGE MANUALLY:
## - path to autoclass 
## - path to autoclass2tab
## - path to where directories of projects will be written
## - path to where shell file will be written - shell file which will launch autoclass for the designated projects

args <- commandArgs(trailingOnly = T)

if(length(args)!=4){
  if(length(args)!=5){
    message("Usage: Rscript 1_launch_autoclass.R  <work.dir> <project.name> <peakmatrix.path> <seeds> <prop.train>")
    message("       seeds must be a comma-separated list of maximum 11 (eg: 103,81,5)")
    message("       prop.train is a number from 0-1, optional argument")
    stop()
  }
}


## SOURCE SCRIPTS, LOAD LIBRARIES  ----

## sourcing sven's script
source("/gpfs/projects/bsc08/shared_projects/IJC_3Dchromatin/programs/autoclass/autoclass.R")

## define paths where different scripts/software are found (autclass software, autoclass2tab script)
autoclasspath="/gpfs/projects/bsc08/shared_projects/IJC_3Dchromatin/programs/autoclass/autoclass/autoclass"
autoclass2tabpath="/gpfs/projects/bsc08/shared_projects/IJC_3Dchromatin/programs/autoclass/2_autoclass2tab_javierrelab.R"
scriptspath="/gpfs/projects/bsc08/shared_projects/IJC_3Dchromatin/programs/autoclass/"

## loading libraries
suppressMessages(library(data.table))
suppressMessages(library(reshape))
suppressMessages(library(ggplot2))                                                                          
suppressMessages(library(cluster))                                                                          
suppressMessages(library(gplots))                                                                           


## READING ARGUMENTS ----
#setting directory we want to work in
work.dir <- args[1]
setwd(work.dir)

system("mkdir jobs")
system("mkdir logs")

project.name <- args[2]
peakmatrix.path <- args[3]
seeds <- args[4]
prop.train <- args[5]

seeds <- as.integer((strsplit(x = seeds,split = ","))[[1]])
## checking how many seeds have been input
## number of seeds == number of subprojects that will be run in parallel, maximum is 11 for optimal time-resource efficiency
if(length(seeds)>11) {
  message("More that 11 seeds have been given. Please only give a maximum of 11")
  message("Aborting mission! Using only first 11 provided")
  seeds <- seeds[1:11]
}

if(is.na(prop.train)){
  prop.train <- 0.0411  # default 0.0411 why? this is the ratio of 30000/728839 (num_interactions in babraham peakmatrix)
}

## LOAD PEAKMATRIX ----

## reading peak matrix
pmBOe_4 <- fread(peakmatrix.path)
pks<-pmBOe_4

## defining number of columns
ncol_sv<-ncol(pks)

## subsetting columns of chicago scores, removing annotation of interactions
pksc_samp <- as.data.frame(pks[, 12:ncol_sv, with = F])

## asinh-transforming the scores
pksc_samp <- asinh(pksc_samp)


## PREPARE PROJECTS ----

## files to run autoclass will be prepared using list of seeds to sample the 30,000 interactions
## files will be prepared for clustering of 3 different samples of 30,000 interactions for <seeds> = 103,81,5
## directory name will be <project.name>_<seed>

## for all seeds we want we:
##	1) set the seed
##	2) sample 30000 interactions
##	3) define project name
##	4) define directory of where project will be set up
##	5) prepare files db2, hd2, model, s-params, r-params using the function autoclass() sourced from autoclass.R
##	6) write command lines into shell script to launch whole pipeline automatically. For all seeds, command lines will be written in same shell script to launch them in parallel using greasy

v <- c()
for(s in seeds){
  #setwd to where project directories will be written	
  setwd(work.dir)

  set.seed(s)
  
  ints2sample <- as.numeric(prop.train)*nrow(pksc_samp)
  sample <- pksc_samp[sample(1:nrow(pksc_samp), ints2sample), ]
  
  # project
  proj <- paste0(project.name,"_",s)
  dir <- paste0(getwd(),"/",proj,"/")
  autoclass(as.data.frame(sample), project=proj, dir=dir, runautoclass=F, runautoclass2tab=F, interactive=F, scale=F, delog=F, rel_error=0.1, single_normal_cn=1:ncol(sample))
  
  
  # creating test.db2
  setwd(dir)
  test_proj <- paste0("test")
  test_dir <- paste0(getwd(),"/",test_proj,"/")
  autoclass(as.data.frame(pksc_samp), project=test_proj, dir=test_dir, runautoclass=F, runautoclass2tab=F, interactive=F, scale=F, delog=F, rel_error=0.1, single_normal_cn=1:ncol(pksc_samp))
  system("mv test.db2 ../")
  setwd(dir)
  system("rm -r test/")
  
  # write sh file
  cinit <- paste0("cd ", dir)
  c1 <- paste0("[# -1 #] ",autoclasspath, " -search ", dir, proj, ".db2 ", dir, proj,".hd2 ", dir, proj,".model ", dir, proj,".s-params < ", scriptspath, "auxfile")
  c2 <- paste0("[# -1 #] ",autoclasspath, " -reports ", dir, proj,".results-bin ", dir, proj,".search ", dir, proj,".r-params")
  c3 <- paste0("[# -1 #] ",autoclasspath, " -predict ", dir, "test.db2 ", dir, proj,".results-bin ", dir, proj,".search ", dir, proj,".r-params")
  c4 <- paste0("[# -1 #] /apps/R/3.6.1/INTEL/bin/Rscript ", autoclass2tabpath , " ", dir, "test-", proj, ".tab 0.5")
  c5 <- paste0("[# -1 #] /apps/R/3.6.1/INTEL/bin/Rscript ", scriptspath ,"3_Sc_heatmap.R " ,peakmatrix.path ," ",dir,"test-",proj,".tab")
  
  v <- c(v,cinit,c1,c2,c3,c4,c5)
  
  
  # write log file
  message("writing log file")
  c0 <- paste0("project name      : ", project.name)
  c1 <- paste0("working directory : ", work.dir)
  c2 <- paste0("seed              : ",s)
  c3 <- paste0("peakmatrix        : ",peakmatrix.path)
  c4 <- paste0("proportion of interactions used in training set: ", prop.train)
  logfile <- c(c0,c1,c2,c3,c4)
  write(logfile, file = paste0(dir,"parameters_", proj,".log"))
  
}

write(v, paste0(work.dir,"jobs/", project.name, ".sh"))   



## write greasy file (configuration file to launch job)

c0 <- "#!/bin/bash"
c1 <- paste0("#SBATCH --job-name=", "\'autoclass_", project.name, "\'")
c2 <- paste0("#SBATCH --workdir=", "\'", work.dir, "logs/\'")
c3 <- paste0("#SBATCH --output=", "\'", work.dir, "logs/", project.name, "_%j.out\'")
c4 <- paste0("#SBATCH --error=", "\'", work.dir, "logs/", project.name, "_%j.err\'")  
c5 <- paste0("#SBATCH --cpus-per-task=4")  
c6 <- paste0("#SBATCH --ntasks=", length(seeds))  
c7 <- paste0("#SBATCH --time=3:00:00")
c8 <- paste0("#SBATCH --mail-type=all")
c9 <- ""
c10 <- ""
c11 <- paste0("/apps/GREASY/latest/INTEL/IMPI/bin/greasy ", work.dir,"jobs/" , project.name, ".sh")

greasy_file <- c(c0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11)
write(greasy_file,paste0(work.dir,"jobs/",project.name, ".greasy"))

# launch job
system(command = paste0("sbatch ",work.dir,"jobs/",project.name, ".greasy"))
message("MAXIMUM 11 TASKS!/SEEDS!")  








