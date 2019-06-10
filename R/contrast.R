#########################################################################
#
#     Copyright (C) 2019 - Lluc Sementé Fernàndez
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.
############################################################################

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
  list_indexes <- as.numeric(as.character(list$Index))

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
  list$Contrast <- contrast_values[contrast_order]
  return(list)
}


order_results_by_contrast <- function(results, peak_matrix, clustering_vector)
{
  # ordering ions from Volcano by contast
  contrast_data <- generate_contrast_data_structure(peak_matrix, clustering_vector)
  lists_to_remove <- c()
  for(list_name in names(results$ionsFromVolcano))
  {
    if(length(strsplit(list_name, split = " vs ")[[1]]) == 2 && dim(results$ionsFromVolcano[[list_name]])[2] != 1)  
    {
      results$ionsFromVolcano[[list_name]] <- order_by_contrast(results$ionsFromVolcano[[list_name]], list_name, contrast_data, T)
    }
    else
    {
      if(length(strsplit(list_name, split = " vs ")[[1]]) == length(unique(clustering_vector)) && dim(results$ionsFromVolcano[[list_name]])[2] != 1)  
      {
        results$ionsFromVolcano[[list_name]] <- order_by_contrast(results$ionsFromVolcano[[list_name]], list_name, contrast_data, F)
      }
      else
      {
        lists_to_remove <- c(lists_to_remove,list_name)
      }
    }
  }
  
  list_indexes <- c()
  for(name in lists_to_remove)
  {
    list_indexes <- c( list_indexes, which(names(results$ionsFromVolcano) == name))
  }
  
  tmp_list <- list()
  cnt <- 1
  old_names <- c()
  for(l in 1:length(results$ionsFromVolcano))
  {
    if(!any(list_indexes==l))
    {
      old_names <- c(old_names,names(results$ionsFromVolcano)[l])
      tmp_list[[cnt]] <- results$ionsFromVolcano[[l]]
      cnt <- cnt + 1
    }
  }
  names(tmp_list) <- old_names
  results$ionsFromVolcano <- tmp_list
 
  # ordering ions from Zeros by contast  
  lists_to_remove <- c()
  for(list_name in names(results$ionsFromZeros))
  {
    if(length(strsplit(list_name, split = " vs ")[[1]]) == 2 && dim(results$ionsFromZeros[[list_name]])[2] != 1)  
    {
      results$ionsFromZeros[[list_name]] <- order_by_contrast(results$ionsFromZeros[[list_name]], list_name, contrast_data, T)
    }
    else
    {
      if(length(strsplit(list_name, split = " vs ")[[1]]) == length(unique(clustering_vector)) && dim(results$ionsFromZeros[[list_name]])[2] != 1)  
      {
        results$ionsFromZeros[[list_name]] <- order_by_contrast(results$ionsFromZeros[[list_name]], list_name, contrast_data, F)
      }
      else
      {
        lists_to_remove <- c(lists_to_remove,list_name)
      }
    }
  }
  
  list_indexes <- c()
  for(name in lists_to_remove)
  {
    list_indexes <- c( list_indexes, which(names(results$ionsFromZeros) == name))
  }
  
  tmp_list <- list()
  cnt <- 1
  old_names <- c()
  for(l in 1:length(results$ionsFromZeros))
  {
    if(!any(list_indexes==l))
    {
      old_names <- c(old_names,names(results$ionsFromZeros)[l])
      tmp_list[[cnt]] <- results$ionsFromZeros[[l]]
      cnt <- cnt + 1
    }
  }
  names(tmp_list) <- old_names
  results$ionsFromZeros <- tmp_list
  
  #Cleaning up matrixes
  for(big_list in 1:2)
  {
    sub_list_end <- (length(results[[big_list]])-1)
    if(sub_list_end == 0)
    {
      sub_list_end == 1
    }
    
    if(sub_list_end > 0)
    {
      for(sub_list in 1:(length(results[[big_list]])-1))
      {
        noisy_ions <-c()
        for(ion in 1:nrow(results[[big_list]][[sub_list]]))
        {
         if(as.character(results[[big_list]][[sub_list]][ion,1]) == as.character(results[[big_list]][[sub_list]][ion,2]))
         {
           noisy_ions <- c(noisy_ions, ion)
         }
        }
        if(!is.null(noisy_ions))
        {
          results[[big_list]][[sub_list]] <- results[[big_list]][[sub_list]][-noisy_ions, ]
        }
      }
    }
  }
  
  return(results)
}



merge_and_reorder_results <- function(results, peak_matrix, clustering_vector)
{
  # New list to store the results
  merged_list <- list()
  new_results <- list()
  len_volcano <- length(names(results$ionsFromVolcano))
  len_zero <- length(names(results$ionsFromZeros))
  
  # Vector with all the list names
  name_list <- unique(names(results$ionsFromVolcano),
                      names(results$ionsFromZeros)
                      )
  last_name <- name_list[length(name_list)]
  
  # For each cluster comparasion
  for(name in name_list)
  {
    flag_V <- any(names(results$ionsFromVolcano) == name) # Exists a list with that name in the volcano?
    flag_Z <- any(names(results$ionsFromZeros) == name) # Exists a list with that name in the zeros?
    if(len_zero > 0 & flag_Z)
    {
      results$ionsFromZeros[[name]]$Source <- rep("zero", times = nrow(results$ionsFromZeros[[name]]))
    }

    if(len_volcano > 0 & flag_V)
    {
      results$ionsFromVolcano[[name]]$Source <- rep("volcano", times = nrow(results$ionsFromVolcano[[name]]))
    }
    
    if(flag_V & flag_Z)
    {
      merged_list[[name]] <- rbind(results$ionsFromVolcano[[name]],results$ionsFromZeros[[name]])
    }
    
    if(flag_V & !flag_Z)
    {
      merged_list[[name]] <- results$ionsFromVolcano[[name]]
    }
    
    if(!flag_V & flag_Z)
    {
      merged_list[[name]] <- results$ionsFromZeros[[name]]
    }
    
    # Reordering the new list.
    if(name != last_name)
    {
      for(row in 1:nrow(merged_list[[name]]))
      {
        if(merged_list[[name]][row,1] == "DownRegulated")
        {
          merged_list[[name]]$Contrast[row] <- 1/merged_list[[name]]$Contrast[row]
        }
      }
    } else
      {
        for(row in 1:nrow(merged_list[[name]]))
        {
          if(any(merged_list[[name]][row,] == "DownRegulated"))
          {
            merged_list[[name]]$Contrast[row] <- 1/merged_list[[name]]$Contrast[row]
          }
        }
      }

    merged_list[[name]] <- merged_list[[name]][order(merged_list[[name]]$Contrast,decreasing = T),]
    rownames(merged_list[[name]]) <- c()
  }
  
  
  
  new_results$ions <- merged_list 
  new_results$data <- results$ionsData
  return(new_results)
}





