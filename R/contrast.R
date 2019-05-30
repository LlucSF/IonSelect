
generate_contrast_data_structure <- function(peak_matrix, clustering_vector)
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
  
 return(list(cube = contrast_cube, all = all_clusters_contrast_matrix))
  
}



get_indexes_from_name <- function(name)
{
  rest <- c()
  indexes <- strsplit(name,split =" vs ")
  
  tmp <- strsplit(indexes[[1]][1],split = "Clus_")
  rest[1] <- as.numeric(tmp[[1]][2])
  
  tmp <- strsplit(indexes[[1]][2],split = "Clus_")
  rest[2] <- as.numeric(tmp[[1]][2])
  
  return(rest)
}



order_by_contrast <- function(list, name, contrast_data, pair = T)
{
  list_indexes <- list$Index

  if(pair)
  {
    clusters <- get_indexes_from_name(name)
    contrast_values <- contrast_data$cube[clusters[1],clusters[2],list_indexes]
    contrast_order <- order(contrast_values, decreasing = T)
  }
    else
    {
      contrast_values <- c()
      contrast_order <- c()
      for(ion in 1:nrow(list))
      {
        cluster <- which(list[ion,] != "--")[1]
        contrast_values <- c(contrast_values, contrast_data$all[cluster,as.numeric(as.character(list$Index[ion]))])
      }
      contrast_order <- order(contrast_values, decreasing = T)
    }
  
  list <- list[contrast_order,]
  list$contrast <- contrast_values[contrast_order]
  return(list)
}




order_results_by_contrast <- function(results, peak_matrix, clustering_vector)
{
  contrast_data <- generate_contrast_data_structure(peak_matrix, clustering_vector)
  for(list_name in names(results$ionsFromVolcano))
  {
    if(length(strsplit(list_name, split = " vs ")[[1]]) == 2 && dim(results$ionsFromVolcano[[list_name]])[2] != 1)  
    {
      results$ionsFromVolcano[[list_name]] <- order_by_contrast(results$ionsFromVolcano[[list_name]], list_name, contrast_data, T)
    }
    
    if(length(strsplit(list_name, split = " vs ")[[1]]) == length(unique(clustering_vector)) && dim(results$ionsFromVolcano[[list_name]])[2] != 1)  
    {
      results$ionsFromVolcano[[list_name]] <- order_by_contrast(results$ionsFromVolcano[[list_name]], list_name, contrast_data, F)
    }
  }
  
  for(list_name in names(results$ionsFromZeros))
  {
    if(length(strsplit(list_name, split = " vs ")[[1]]) == 2 && dim(results$ionsFromZeros[[list_name]])[2] != 1)  
    {
      results$ionsFromZeros[[list_name]] <- order_by_contrast(results$ionsFromZeros[[list_name]], list_name, contrast_data, T)
    }
    
    if(length(strsplit(list_name, split = " vs ")[[1]]) == length(unique(clustering_vector)) && dim(results$ionsFromZeros[[list_name]])[2] != 1)  
    {
      results$ionsFromZeros[[list_name]] <- order_by_contrast(results$ionsFromZeros[[list_name]], list_name, contrast_data, F)
    }
  }
  
  
  
  return(results)
}

list_name <- names(results$ionsFromVolcano)[1]
results <- r
results <- order_results_by_contrast(results, peak_matrix, clustering_vector)

