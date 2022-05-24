# USAGE

If you run any of the following scripts without arguments, a usage message will print on page

## 1_launch_autoclass.R
```bash
Rscript 1_launch_autoclass.R <work.dir> <project.name> <peakmatrix.path> <seeds> <prop.train>
```  
Argument        | Input
----------------|-------------------------------------------------------------------------------------------------------------
work.dir        | Full path to directory new project will be created 
project.name    | Name of new project  
peakmatrix.path | Full path to peakmatrix
seeds           | Comma-separated list of seeds to set when randomly sampling the proportion of interactions chosen to use as training set  
prop.train      | Optional - proportion of interactions chosed to use as training set for clustering. If not set, default is 0.0411, which is the ratio of 30000/728839 (num_interactions in babraham peakmatrix)  
  
  
## 2_autoclass2tab_javierrelab.R  

```bash
Rscript 2_autoclass2tab_javierrelab.R <tab.file> <pmin>
```
Argument   | Input
-----------|----------------------------------------------------------------------------------------------------------------
tab.file   | full path to \*.tab file with the following format: $1=Case $2=Class  $3=Prob (non-existing at this point)
pmin       | probability values cutoff. Any interactions clustered with a prob < 0.5 is removed from the dataset  
  
**Output**    
\*.tab file with the following format: $1=Case $2=Class $3=Prob (each row is one interaction). File will be found as stated in tab.file
  
## 3_Sc_heatmap.R
```bash
Rscript 3_Sc_heatmap.R <peakmatrix.path> <tab.name> <cluster.order>
```
Argument         | Input
-----------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
peakmatrix.path  | Full path to peakmatrix
tab.name         | Full path to \*.tab file with the following format: $1=Case $2=Class  $3=Prob
cluster.order    | Optional - in which order to position clusters in visualisation, by "hclust","genomicdist" or "name", default being hclust

**Output**    
PDF of heatmap in same directory where tab.name is found, suffix Sc_heatmap.pdf  


## 3_interaction_heatmap.R
```bash
Rscript 3_interaction_heatmap.R <peakmatrix.path> <tab.name> <factor.ds> <cluster.order> <exclude.clusters> 
```
Argument         | Input
-----------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
peakmatrix.path  | Full path to peakmatrix  
tab.name         | Full path to \*.tab file with the following format: $1=Case $2=Class  $3=Prob
factor.ds        | Value from 0-1, factor of downsampling to apply to dataset for faster visualisation (default is 0.5)  
cluster.order    | In which order to position clusters in visualisation, by "hclust","genomicdist" or "name"  exclude.clusters | Optional, vector of clusters to be excluded, given in format "0,1,2"  

**Output**    
PDF of heatmap in same directory where tab.name is found, suffic interactions_heatmap.pdf 

