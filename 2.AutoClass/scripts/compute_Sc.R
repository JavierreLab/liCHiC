## compute_Sc.R - @ A Nieto-Aliseda - mn1
## script to cell type calculate specificity score given a set of interactions


library(data.table)                                                                                                                                                                    
                                                                                                                                                                                       
                                                                                                                                                                                       
# peakmat - subset of peakmatrix we want specificity scores for each cell type for (i.e. A
# subset - vector of row indices from which to subset peakmat (i.e. a cluster defined by a
# cellcols - columns of peakmatrix where chicago scores for cell types of interest are found
specificity_score <- function(peakmat,subset=NA,cellcols,name=NA){                                                                                                                     
  options(stringsAsFactors = F)                                                                                                                                                                                       
  #extracting the chicago scores from the peak matrix
  cell_types <- colnames(peakmat)[cellcols]                                                                                                                                            

  chi_scores <- as.matrix(peakmat[,cellcols])                                                                                                                                   
  #chi_scores <- asinh(chi_scores)
                                                                                                                                                                                       
  #calculating distances between cell types based on chicago scores 
  distance_matrix <- dist(t(chi_scores), method = "euclidean", upper= TRUE)                                                                                                            
  distance_matrix <- as.matrix(distance_matrix)                                                                                                                                        
                                                                                                                                                                                       
  # sum of distances for each cell type given Sc formula
  sum_distances <- colSums(distance_matrix)                                                                                                                                            
                                                                                                                                                                                       
  #applying 1/x to sum of distances for weighting given the Sc formula
  weighting <- lapply(sum_distances, function(x) 1/x)                                                                                                                                  
                                                                                                                                                                                       
  #mean CHiCAGO score interaction for every cell type
                                                                                                                                                                                       
  #at this point subsetting chicago_scores for specific cluster if necessary
  if(!is.na(subset)){                                                                                                                                                                  
    chi_scores <-chi_scores[subset,]                                                                                                                                      
  }                                                                                                                                                                                  
                                                                                                                                                                                       
  mean_chiscores <- colMeans(chi_scores)                                                                                                        
                                                                                                                                                                                       
  Sc <- rep(NA, length(mean_chiscores))                                                                                                                                                
                                                                                                                                                                                       
  for(x in 1:length(mean_chiscores)){                                                                                                                                                  
    dmultiplier <- distance_matrix[,x] #obtain distance of c-i from distance matrix , need                                                                                             
    dmultiplier <- dmultiplier[dmultiplier != 0] #remove 0 value      
    
    XcXi <- mean_chiscores[x] - mean_chiscores[-c(x)]                                                                                                                                  
                                                                                                                                                                                       
    dXcXi <- XcXi                                                                                                                                                                      
                                                                                                                                                                                       
    for(i in 1:length(XcXi)){                                                                                                                                                          
      dXcXi[i] <- dmultiplier[i]*XcXi[i]                                                                                                                                               
    }                                                                                                                                                                                  
                                                                                                                                                                                       
    sum_dXcXi <- sum(dXcXi)                                                                                                                                                            
                                                                                                                                                                                       
    Sc[x] <- as.double(weighting[x])*sum_dXcXi                                                                                                                                         
                                                                                                                                                                                       
  }                                                                                                                                                                                    
                                                                                                                                                                                       
  Sc <- as.data.frame(rbind(cell_types,Sc))                                                                                                                                            
  colnames(Sc) <- cell_types                                                                                                                                                           
  Sc <- Sc[-c(1),]                                                                                                                                                                     
  Sc[1,] <- as.double(Sc[1,])                                                                                                                                                          
                                                                                                                                                                                       
  if(!is.na(name)){                                                                                                                                                                    
  rownames(Sc) <- name                                                                                                                                                                 
  }                                                                                                                                                                                    
  return(Sc)                                                                                                                                                                           
}                                                                                                                                                                                     
                                                                                                                  
                                                                                                                                                                                       


