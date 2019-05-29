
contrast_order <- function(peak_matrix, clustering_vector)
{
  intensity_matrix <- peak_matrix$intensity
  clustering <- clustering_vector
  number_of_clusters <- max(clustering)
  number_of_ions <- length(peak_matrix$mass)
  
  all_clusters_contrast_matrix <- matrix(nrow = number_of_clusters, ncol=number_of_ions)
  contrast_cube <- array(dim = c(number_of_clusters, number_of_clusters, number_of_ions))
  
  for(clusA in 1:number_of_clusters)
  {
    for(ion in 1:number_of_ions)
    {
      all_clusters_contrast_matrix[clusA,ion] <- mean(intensity_matrix[which(clustering == clusA), ion])/
                                        mean(intensity_matrix[,ion])
    }
  }
  
  for(clusA in 1:number_of_clusters)
  {
    for(clusB in 1:number_of_clusters)
    {
      for(ion in 1:number_of_ions)
      {
        contrast_cube[clusA, clusB, ion] <- mean(intensity_matrix[which(clustering == clusA), ion])/
                                            mean(intensity_matrix[which(clustering == clusB), ion])
      }
    }
  }
  
  contrast_cube[7, 2, 23]
  
  
}

