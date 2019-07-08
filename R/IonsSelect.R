#########################################################################
#
#     Copyright (C) 2018 - Lluc Sementé Fernàndez
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

#' TestIonSelect
#'
#' @description Performs a test that calculates and select the up-regulated and/or down-regulated
#' ions from a peak matrix between defernt clusters.
#'
#' @param PeakMtx List. An rMSIprocPeakMatrix. Must contain at least the following categories: \itemize{
#' \item PeakMtx$intensity. A matrix containg the intensities of the peaks for each pixel (rows = pixels, cols = ions).  
#' \item PeakMtx$mass. A vector containg the masses of each peak. Must be in the same order with the columns of the intensity marix.   
#' \item PeakMtx$numPixels. Number of pixels (rows in your matrix).   
#' }
#' @param clusters Numeric Vector. A vector containing the cluster index for each pixel coded in numbers from 1 to number of clusters.
#' @param clusterSubset Numeric vector containing the index of the clusters which will be avaluated.
#' @param precentile Numeric Vector. Percentile restriction of the test for each measure (null ion concentration,p-value, fold change). 
#' @param zeroThreshold Numeric. Intensity value below which an ion is considered a null ion in a pixel.
#' @return List containing three elements. The results from the Volcano test, the results from the Zero test and the values of p, FC & Zero scores. 
#' @export
#'

TestIonSelect <- function(PeakMtx, clusters, zeroThreshold = 0, percentile = c(1,10,10),clusterSubset = sort(unique(clusters)))
{
  ### Imput check and format ###
  clusterSubset <- sort(clusterSubset)
  name_clusters <- unique(clusters)
  
  if(length(clusterSubset)<2)
  {
    return( writeLines("Wrong cluster indexes. At least two clusters must be selected"))
  }
  
  for(i in 1:length(clusterSubset))
  {
    if(clusterSubset[i] > length(name_clusters) | (clusterSubset[i] == 0))
    {
    return( writeLines("Wrong cluster indexes. Indexes must be between 1 and number of clusters Ex: c(1,3,4) to select clusters 1,3 and 4 from a clustering with 4 centers"))
    }
  }

  size_clusters <- c()
  for(i in clusterSubset)
  {
    size_clusters <- c(size_clusters, length(which(clusters == i)))
  }
  
  ID_clusters <- list()
  for(i in 1:length(clusterSubset))
  {
    ID_clusters[[i]] <- which(clusters == clusterSubset[i])
  }
  
  intensityData <- matrix(ncol = length(PeakMtx$mass), nrow = length(unlist(ID_clusters)))
  min <- 0
  for(i in 1:length(ID_clusters))
  {
    pixel_range <- (min + 1):(length(ID_clusters[[i]]) + min)
    intensityData[pixel_range,] <- PeakMtx$intensity[ID_clusters[[i]],]
    ID_clusters[[i]] <- pixel_range - 1
    min <- min + length(ID_clusters[[i]])
  }


  ### C call ###
  Test <- IonSelectC(   zPercentil = percentile[1], 
                        pPercentil = percentile[2], 
                        fcPercentil = percentile[3], 
                          numPixels = sum(PeakMtx$numPixels),
                          SP_Pixels = c(PeakMtx$numPixels),
                         # numPixels = nrow(intensityData), 
                         # SP_Pixels = nrow(intensityData),                    
                            numCols = length(PeakMtx$mass),
                           massAxis = PeakMtx$mass, 
                        numSamples = length(PeakMtx$numPixels),
                       nPTestGroups = length(clusterSubset), 
                      R_pTestGroups = (1:length(clusterSubset))-1,
                       ClustersSize = size_clusters, 
                     ClustersPixels = ID_clusters,
                               data = intensityData,
                      zeroThreshold = zeroThreshold
                     )
  

  ### Output format ###
  #dummy variables
  d <- NULL             
  cnt <- 0              #dummy counter
  IonsData   <- list() 
  emptyRows  <- c()     #variable used to discart the empty rows in each cluster comparation
  emptyList  <- c()     #variable used to discard the empty lists in the ions data
  emptyList1 <- c()     #variable used to discard the empty lists in the volcano symbols
  emptyList2 <- c()     #variable used to discard the empty lists in the ions data symbols
  lname <- c()          #variable containing the name of the lists
  cname <- c()          #names of the columns
  
  
  #Discard data from the comparation between same clusters and naming the lists 
  for(j in 1:(2^length(clusterSubset)))
  {
    if(is.null(dim(Test[[1]][[j]])))
    {
     d <- c(d,j)
    }
      else
      {
        cluster_indexes <- which(intToBits(j-1)==1)
        if(length(cluster_indexes) >= 2)
        {
          cname <- paste("Clus_", sort(clusterSubset[cluster_indexes], decreasing = F),sep = "")
          lname <- c(lname,paste("Clus_", sort(clusterSubset[cluster_indexes], decreasing = F), sep = "", collapse = " vs "))
          colnames(Test[[1]][[j]]) <- cname
          colnames(Test[[2]][[j]]) <- cname
        }
      }
  }
  Test[[1]] <- Test[[1]][-d]
  Test[[2]] <- Test[[2]][-d]
  names(Test[[1]]) <- lname
  names(Test[[2]]) <- lname
  

  a <- rep(paste("Clus_",clusterSubset, " vs",sep = ""),times = length(clusterSubset))
  b <- rep(paste("Clus_",clusterSubset,sep = ""),each = length(clusterSubset))
  ionname <- paste(a,b)
  for(i in 1:length(clusterSubset)^2)
  {
    IonsData[[i]] <- as.data.frame(Test[[3]][ ,((3*i-2):(3*i))])
    IonsData[[i]] <- cbind(IonsData[[i]], as.character(format(PeakMtx$mass, nsmall = 3, digits = 3)))
    IonsData[[i]] <- cbind(IonsData[[i]], as.character(format(1:length(PeakMtx$mass), nsmall = 3, digits = 3)))
    colnames(IonsData[[i]]) <- c("Zero", "p value", "FC", "Ion", "Index") 
  }
  names(IonsData) <- ionname
  sameClusterData <- which(unlist(lapply(lapply((strsplit(ionname,split = " vs ")), unique),length))==1)
  Test[[3]] <- IonsData[-sameClusterData]
  
  
  #Encoding the symbol matrixes & deleting rows and null lists
  for(h in 1:2)
  {
    emptyList <- c()
    for(i in 1:length(Test[[h]]))
    {
      emptyRows <- c() #reset
      
      Test[[h]][[i]] <- cbind(Test[[h]][[i]], as.character(format(PeakMtx$mass, nsmall = 3, digits = 3)))
      colnames(Test[[h]][[i]])[ncol(Test[[h]][[i]])] <- "Ion" 
      Test[[h]][[i]] <- cbind(Test[[h]][[i]], as.character(format(1:length(PeakMtx$mass), nsmall = 3, digits = 3)))
      colnames(Test[[h]][[i]])[ncol(Test[[h]][[i]])] <- "Index" 
      
      for(j in 1:length(PeakMtx$mass))
      {
        for(k in 1:(dim(Test[[h]][[i]])[2]-1))
        {
          if((Test[[h]][[i]][j,k] == 3) | (Test[[h]][[i]][j,k] == 2))
          {
            Test[[h]][[i]][j,k] <- "DownRegulated"
          }

          if((Test[[h]][[i]][j,k] == 1) | (Test[[h]][[i]][j,k] == 4))
          {
            Test[[h]][[i]][j,k] <- "UpRegulated"
          }

          if(Test[[h]][[i]][j,k] == 0)
          {
            Test[[h]][[i]][j,k] <- "--"
          }
        }

        #If one row is empty, remove it later
        if(all(Test[[h]][[i]][j,-((-1:0)+(dim(Test[[h]][[i]])[2]))] == "--"))
        {
          emptyRows <- c(emptyRows,j)
        }
      }

      #remove the rows
      if(length(emptyRows) > 0)
      {
        if(length(emptyRows) == length(PeakMtx$mass))
        {
          emptyList <- c(emptyList,i)
        }
          else
          {
            Test[[h]][[i]] <- as.data.frame(Test[[h]][[i]][-emptyRows,],make.names = F)
          }
      }
    }
    
    if(!is.null(emptyList))
    {
      Test[[h]] <- Test[[h]][-emptyList]
      if(h == 1)
      {
        emptyList1 <- emptyList
      }
    }
  }

  # Contrast parameter & clean up
  Test <- order_results_by_contrast(Test, PeakMtx, clusters)
  Test <- merge_and_reorder_results(Test, PeakMtx, clusters)

  return(Test)
}

#' plotIonDatabyComparedClusters
#'
#' @description Plots the data related to the selected ions from the comparasion between the given clusters. 
#'
#' @param testResults The results of the test. 
#' @param clusterIndex Numeric Vector. A vector containing the cluster index for each pixel coded in numbers from 1 to number of clusters.
#' @param max_ions Numeric. Maximum nunber of ions to be displayed.
#' 
#'

plotIonDatabyComparedClusters <- function(testResults, clusterIndex, max_ions)
{
  clusterIndex <- sort(clusterIndex,decreasing = F)
  name <- ""
  for(i in 1:length(clusterIndex))
  {
    if(i == length(clusterIndex))
    {
      name <- paste(name,"Clus_",clusterIndex[i],sep = "")
    } else
      {
        name <- paste(name,"Clus_",clusterIndex[i]," vs ",sep = "")
      }
  }
  
 
  if(any(names(testResults$ions) == name))
  {
    ion_rows <- as.numeric(as.character(testResults$ions[[name]]$Index))
    df <- testResults$data[[name]][ion_rows, ]
    df$Contrast <- testResults$ions[[name]]$Contrast
    df$Source <- testResults$ions[[name]]$Source
    colnames(df) <- c("Z","p","FC","Ion","Index","Contrast","Source")
    if(nrow(df)>=max_ions)
    {
      df <- df[1:max_ions,]
    }
  } else
    {
      writeLines(paste("Wrong list index. List",name,"is not avaliable. Only the following are allowed: \n"))
      return(print(names(testResults[[df_name]])))
    }
  
    options(warn = -1)
  
    df$Ion <- as.character(df$Ion)
    labs <- factor(1:length(df$Ion), labels = df$Ion)
    df$Ion <- as.numeric(labs)

    g <- ggplot2::ggplot(#*data = df, mapping = ggplot2::aes(x = Ion)#
      ) +
         ggplot2::geom_col(data = subset(df, df$Source == "volcano"),
                           mapping = ggplot2::aes(y = -log2(FC),
                                                  fill = -log10(p),
                                                  x = df$Ion[df$Source =="volcano"],
                                                  colour = Source
                                                  )
                          ) +
      ggplot2::geom_col(data = subset(df, df$Source == "zero"),
                        mapping = ggplot2::aes(y = log2(Z),
                                               colour = Source,
                                               x = df$Ion[df$Source =="zero"]
                                               ),
                        fill = "white"
                        ) +
      ggplot2::labs(y = "log2(Fold Change)",
                    x = "m/z",
                    fill = "-log10(p)",
                    title = name) +
      ggplot2::scale_y_continuous(sec.axis = ggplot2::sec_axis(~.*1, name = "log2(Zero score)")) +
      ggplot2::scale_fill_continuous(type = "viridis") + scale_colour_manual(values = c("blue","red")) +
      ggplot2::theme_bw() +
      ggplot2::guides(colour = guide_legend(override.aes=list(shape=21))) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = -60, hjust = 0)) + 
      ggplot2::scale_x_discrete(breaks = 1:max_ions, labels = as.character(labs),position = "bottom", drop = T)
  print(g)
  options(warn = 0)
}



