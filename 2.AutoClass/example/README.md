
# TUTORIAL

For the purpose of the tutorial, all the examples and datasets are only for chromosome 22. This means that the jobs probably don't have to be queued and would run well in command line. However, for the purpose of practice, try it!

Find all necessary scripts in **/gpfs/projects/bsc08/shared_projects/IJC_3Dchromatin/programs/autoclass/**  
Find all necessary data in **/gpfs/projects/bsc08/shared_projects/IJC_3Dchromatin/programs/autoclass/example_github/**  

As the Autoclass C software is already installed, we can start immediately.  

### Things to know
There are some things that you must know before starting to run the pipeline. I will be bulletpointing all those that come to mind here: 

- There is a lot of documentation describing the whole works of the Autoclass software. Please read this in detail if you want to understand how the software the clustering of interactions is done. Otherwise, you'll just have to trust that what I have done is correct! 
- If you don't use the set up built in marenostrum, you will have to change quite a few paths manually in the scripts:
  - /path/to/templates/ in autoclass wrapper script autoclass.R (Mikhail Spivakov 2011)
  - /path/to/autoclass.R (Mikhail Spivakov 2011) in 1_launch_autoclass.R, 3_Sc_heatmap.R and 3_interaction_heatmap.R
  - /path/to/autoclass software in 1_launch_autoclass.R
  - /path/to/2_autoclass2tab_javierrelab.R in 1_launch_autoclass.R
  - /path/to/all/autoclass/scripts/ in 1_launch_autoclass.R
  - /path/to/plot_CRM_heatmap.R (Mikhail Spivakov 2014) in 3_Sc_heatmap.R and 3_interaction_heatmap.R
  - /path/to/calculate_specificity_score.R in 3_Sc_heatmap.R and 3_interaction_heatmap.R
- All original scripts from Sven Sewitz can be found here (??)
  


### 0 - Preparation  
1) Download Autoclass C Software [here](https://ti.arc.nasa.gov/tech/rse/synthesis-projects-applications/autoclass/autoclass-c/)
2) Prepare templates for input files needed to run Autoclass C (see what these files are in introduction readme or in autoclass documentation), or download them from this github [here](https://github.com/JavierreLab/RunAutoClass/tree/main/example/templates)
3)  



### 1 - Launching Autoclass   
```bash
Rscript 1_launch_autoclass.R /gpfs/projects/bsc08/shared_projects/IJC_3Dchromatin/programs/autoclass/example_github/ example /gpfs/projects/bsc08/shared_projects/IJC_3Dchromatin/programs/autoclass/example_github/PCHiC.peakmatrix.cutoff5.chr22.txt 103
```
  
### 2 - Visualisation of interaction heatmap
```bash
Rscript 3_interaction_heatmap.R /gpfs/projects/bsc08/shared_projects/IJC_3Dchromatin/programs/autoclass/example_github/PCHiC.peakmatrix.cutoff5.chr22.txt /gpfs/projects/bsc08/shared_projects/IJC_3Dchromatin/programs/autoclass/example_github/example_103/test-example_103.tab 1 hclust
```
Suggestion: Run this manually and changing parameters manually in script to obtain heatmap wanted.
