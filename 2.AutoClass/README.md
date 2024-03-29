# Run AutoClass  

AutoClass is an unsupervised Bayesian classification system that seeks a maximum posterior probability classification.  
It can be found in the following link: https://ti.arc.nasa.gov/tech/rse/synthesis-projects-applications/autoclass/autoclass-c/  

## Credits

In this pipeline we will be using the R wrapper of the software, which was created at the Babraham Institute (Sven Sewitz). The whole pipeline is an adaptation of the scripts shared by Sven with Javierre Lab.

The calculation of Cluster Specificity Score is based on the [scripts generated by Steven Hill](https://github.com/Steven-M-Hill/PCHiC-specificity-score-analysis).


## What do we want to cluster?  

We have what we call a peakmatrix. This is a number matrix, with rows representing PCHi-C interactions and columns representing cell type (+ metadata describing the interaction i.e. coordinates of both ends of this interaction in the first 11 columns). In this matrix we have the interactions that have been found statistically significant by CHiCAGO in at least one of the cell types present in the matrix. An interaction is found significant when the chicago score >= 5. 
In the example folder we may find a peak matrix of dimensions 728,838 rows (interactions) by 28 columns (11 being metadata and the remaining 17 corresponding to cell types. This peakmatrix is originally from **Javierre et al. Cell 2016**. 

For a step to step tutorial, read [this](https://github.com/JavierreLab/liCHiC/tree/main/2.AutoClass/example/README.md)  
For an in depth explanation of the usage of the scripts, read [this](https://github.com/JavierreLab/liCHiC/tree/main/2.AutoClass/scripts/README.md)
  
#### Dependencies
* Autoclass Software - [link](https://ti.arc.nasa.gov/tech/rse/synthesis-projects-applications/autoclass/autoclass-c/)
* R
  * data.table 
  * GenomicRanges
  * GenomicInteractions
  * gplots
  * reshape
  * ggplot2
  * stats
  * gtools
  * RColorBrewer

    
## Summary of workflow:  

1.  Preparing the files required by AutoClass to run  
2.  Running AutoClass  
3.  Managing AutoClass' output files  
4.  Generating visualisation of results  

## Workflow  
This is built to work in marenostrum, but can work in local as well. If you want it to work in local, there are many script sourcing dependencies that must be changed manually.  

### Launching the pipeline

In reality, running this pipeline consists of solely one line of code: 

```bash
./1_launch_autoclass.R <workdir> <project_name> <peakmatrix_path> <seeds> <prop_train>
```
But what does this actually do? 

* a - creates a new folder in working directory **workdir** with name "project\_seed/" for all seeds
* b - creates **jobs** and **logs** folders if they do not already exist in **workdir**
* c - sets proportion of training dataset to 4% of all interactions if this has not been predefined as the argument <prop_train>. This default percentage is based on the 30,000 interactions that were sampled from the 728,838 interactions in the peakmatrix published by **Javierre et al. Cell 2016**  
* d - loads the peakmatrix and performs an asinh transformation on the data
* e - for each seed:
  * sets the seed
  * samples <prop_train> of interactions from dataset
  * creates the project directory
  * writes input files required by the autoclass software into the project directory (\*see detail below)
  * writes a log file with the parameters used in this script
* f - writes a job file (\*.sh) and its corresponding configuration file (\*.greasy) into the jobs directory, which contain the tasks to be run in the following steps of the pipeline. (\*see detail below)
* g - launches \*.greasy file to be run in marenostrum  

#### Input files required by AutoClass   
AutoClass requires various files as input in order to run:   
  - file.db2: The actual data values are in this data file (i.e. the chicago scores)  (training set)
  - file.hd2: This is a header file, describing the specific data format and attribute definitions  
  - file.model: This file contains a description of the model that will be used for classification   
  - file.s-params: This is the search parameters file. Here, we set the parameters by which the search for clusters is performed  
  - file.r-params: This is the report parameters file. Here, we set the parameters by which the reports are written once the search has been performed  
  - test.db2: The test set of the chicago score to extrapolate clustering done on training set  
  
These files are prepared by running a function the function named “autoclass”, which is found within the script *autoclass.R*, provided by Sven Sewitz (Babraham Institute). In essence, this function is a wrapper for the actual software.

This function uses a directory as “base” where templates for these input files are found. This directory is defined within the autoclass.R script. In the template files:  
  - template.db2, template.hd2, template.model: 
These files are empty  
  
  - template.s-params: 
Customised search parameters are added here. In this case, so that the results can be reproducible, the following parameters are used
start_fn_type = "block"  
randomize_random_p = false  
rel_error = 0.1  
min_report_period = 30000  
force_new_search_p = true  
break_on_warnings_p = false
max_duration = 9720 (2h42)
Otherwise, autoclass includes pseudo-randomness in search (seed set based on UTC)
  
 - template.r-params:  
Customised report parameters are added here. In this case, we have  
xref_class_report_att_list = 0, 1, 2   
report_mode = "data"   
comment_data_headers_p = true   
  
### Tasks of following steps in pipeline
By running *1_launch_autoclass.R*, you are not only preparing the input files for AutoClass, but launching the whole pipeline. These are the next steps that are performed in the pipeline:

* setwd to directory of /path/to/project_seed/
* run AutoClass C in -search mode. This means training set is clustered  
* run AutoClass C in -reports mode. This writes reports of the clustering done in previous steps  
* run AutoClass C in -predict mode. This predicts clustering of test.db2 based on clustering of training set  

After these steps we will have obtained the "optimal" clustering of the interactions from out peakmatrix, given by AutoClass.

* parse result files to obtain a simplified table of the results  

For an easy interpretation the clustering results, this step in the pipeline aims to re-format the data into a simple 3 column table (\*.tab). The file **test.case-data-1** is parsed to extract  information on which interactions belong to which cluster and with which probability, creating a table of 3 columns: 
column 1 containing ID of the interaction (IDs are defined in step 2 when autoclass function is run);    
column 2 contains the probability with which the interaction is classified into its cluster;    
column 3 contains the ID of the cluster the interaction has been classified into   

In the process of re-formatting, we also have to remove interactions from our data which were not classified as being part of a cluster with a high enough probability. The default in the complete pipeline us a minimum probability of 0.5 (pmin=0.5) . However, this can be changed manually if necessary. See the [Usage](https://github.com/JavierreLab/RunAutoClass/blob/main/scripts/README.md) page for more information.
If you observe that your clustered datasets doesn't include all the interactions in your original dataset (peakmatrix), this may be the reason. Most clusterings will mantain around 95% of the original data. However, we have encountered cases whereby in this filtering step, whole clusters of interactions were removed from the data.

* Run *3_Sc_heatmap.R* script which generates heatmap of cluster and cell type specific specificity scores of the dataset.   
  
  These steps are run in parallel for the whole <project_name> with its different seeds. In other words, each <project_name> will have y number of subprojects, where y corresponds to the number of training sets generated with different seeds.  
  These subprojects will work through the pipeline all in parallel to one another. A maximum of 11 subprojects must be run in parallel for optimum time-resource efficiency. This is defined by the configuration file. If more than 11 are given as input, the pipeline will be stopped. 
   
### Visualisation of results  

#### Specificity Score Heatmap
*3_Sc_heatmap2.R* will generate a PDF file with a heatmap of cluster and cell type specific specificity scores of the dataset
  
#####  What is cluster cell type specificity score?  
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5123897/  

#####  Definition of Specificity Scores  
Consider a set of cell types I. Let xi denote the measured value of a quantitative property (such as CHiCAGO interaction score or gene expression level) for cell type i ∈ I. Then, the specificity score sc for a given cell type c ∈ I is a weighted mean of the differences xc – xi for i ≠ c,

          `sc=1∑i≠cdc,i∑i≠cdc,i(xc−xi)`  

where the weights dc,i are distances between cell type c and cell types i, calculated using the complete dataset (e.g., CHiCAGO interaction scores for all interactions or expression values for all genes; distances calculated using Euclidean distance metric). The distance weights are introduced to account for imbalances in the distances between cell types. For example, among the cell types considered here are three types of macrophages that are likely to have very similar profiles of the measured property compared with other analyzed cell types (and so the distances between macrophage samples will also be smaller than between macrophages and other cell types). The distance weights focus the calculation of sc on cell types that are relatively more distant from cell type c. In this example therefore, they will result in the calculation of sc for each type of macrophage placing relatively little weight on the other types of macrophages. Without this weighting, specificity scores for macrophages would be smaller on average simply because macrophages are over-represented among the cell types considered.  

#####  Calculation of Cluster Specificity Scores  
For a given Autoclass cluster, a specificity score Sc was calculated for each cell type c using the equation above, with xi defined as the mean asinh-transformed CHiCAGO score for cell type i (mean calculated across all interactions in the given cluster). The distance weights weights dc,i were calculated based on the full set of CHiCAGO interaction scores.  

####  Workflow to generate heatmap  

1) Annotating the peakmatrix:  
test.db2 and test-project\_name\_seed.tab are merged by the interaction ID. This annotates the interactions by classifying them into clusters and giving us probability values for their classification. This also filters the interactions, as those without a “match” in the tab from step 4 will be discarded in the merging process.  
  
2) Calculating cluster Sc:  
The cluster Sc is calculated using a function defined in the script *calculate_specificity_score.R* (@A. Nieto 2020).  
However, a function was found published with the Cell paper https://github.com/Steven-M-Hill/PCHiC-specificity-score-analysis which does the same computations.  
They both results in a data matrix of dimensions NxM, N being the number of clusters (rows), M being the number of cell types (columns).  

3) Ordering the clusters:   
This will be based on mean genomic distance, name of however the hierchical clustering orders them.

4) Plot generation:  
Plot will be saving in current working directory with same basename as that of the project_name.  
  
  
#### Interactions Heatmap
 *3_interaction_heatmap2.R* will generate a PDF file with a heatmap of interactions as rows, clustered together as defined by autoclass and columns as cell types, where the colours define the asinh transformed chicago score of the interactions for each cell type. 
 
 ####  Workflow to generate heatmap  
 
 1) Annotating the peakmatrix:  
 test.db2 and test-project\_name\_seed.tab are merged by the interaction ID. This annotates the interactions by classifying them into clusters and giving us probability values for their classification. This also filters the interactions, as those without a “match” in the tab from step 4 will be discarded in the merging process   
 
 2) Randomly downsample dataset of interactions:
 This is done for the solely for visualisation purposes. Why? because without subsetting, plotting of heatmap is too computationally costly. By randomly sampling the clusters, we hope to obtain the most accurate representation of the real heatmap. The clusters are randomly subsetted individually so as to mantain size proportions  
    
 3) Ordering the clusters:   
 This will be based on mean genomic distance, name of however the hierchical clustering orders them  
 
 4) Excluding clusters:
 This is optional. You may exlude a large (uninformative) cluster from the visualisation to make it more appealing  
 
 5) Plot generation:  
 Sourced function from one of sven's scripts called *plot_CRM_heatmap.R*. In summary,this script plots individual heatmaps for each cluster of interactions, where the interactions are ordered by herarchical clustering. All individual clusterings (one for each cluster/class defined by autoclass C) are then "pasted" together  
 
##### Interpret output file  
The function outputs a PDF 5 pages long. The 1st gives us the limits of where one cluster ends and another begins; the 2nd gives us the main figure - the heatmap. Here, each line/row represents an interaction, which are clustered given autoclass C’s classification and ordered by the cluster’s mean genomic distances. Each column represents a cell type; pg 3&4 are the dendrograms of the hierarchical clusterings of rows and columns, if decided to do so. If set to false, these pages will remain blank; and the last page gives us the legend of the heatmap.  
  
I have not get found the "optimal" way of running this script in the cluster. For this reason, this final function is separate from the rest of the pipeline - up to now I have been running this in local. However, in marenostrum this can be run too, but perhaps the proportion of interactions visualised will be very low. The reason for this downsampling is that when each individual heatmap is created, the function hclust() is run to order the interactions by hierarchical clustering. This function crashes when too many interactions are found in a cluster - due to memory errors. 

  
  
*ISSUES*
- set seplwd=0 to full heatmap.3 function 
- make a more optimal way of writing sh script
- Add option of how to order clusters, by clustering? by mean genomic distance? for full heatmap - classOrder=c(a,b,c,d,e) argument
- tutorial + example datasets
- add logfile for autoclass parameters/arguments used  


  
  

