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
#' @param PeakMtx A matrix of peaks containing mass spectra data
#' @param Probability Probability restriction of the test.
#' @param SegmentedMtx A matrix containing a segmentated peak matrix.
#' @param SgmtsChoose A vector containing which clusters should be analyzed.
#' @return List containing three elements. The results from the Volcano test, the results from the Zero test and the values of p, FC & Zero scores. 
#' @export
#'

  TestIonSelect <- function (PeakMtx, Probability, SegmentedMtx, SgmtsChoose)
  {
    for(i in 1:length(SgmtsChoose))
    {
      if (SgmtsChoose[i]>=length(SgmtsChoose))
      {
      return( writeLines("Wrong Cluster Indexes"))
      }
    }

  numPixels <- 0
  for (samples in 1:length(PeakMtx$numPixels))
  {
    numPixels <- numPixels + PeakMtx$numPixels[samples]
  }

  ClustrIds <- list()
  for (j in 1:length(SegmentedMtx$size))
  {
    cnt <- 0
    ClustrIds[[j]] <- 0
    for (Id in 1:numPixels)
    {
      if (SegmentedMtx$cluster[Id] == j)
      {
        cnt <- cnt+1
        ClustrIds[[j]][cnt] <- Id-1
      }
    }
  }

  Test <- IonSelectC(m_focalProb = Probability, numPixels = numPixels, SP_Pixels = PeakMtx$numPixels,
                     numCols = length(PeakMtx$mass), massAxis = PeakMtx$mass, numSamples = samples,
                     nPTestGroups = length(SgmtsChoose), R_pTestGroups = SgmtsChoose,
                     ClustersSize = SegmentedMtx$size, ClustersPixels = ClustrIds, data = PeakMtx$intensity)

  d <- NULL
  cnt <- 0
  IonsData <- list()
  for(j in 1:(2^length(SgmtsChoose)))
  {
    if(is.null(dim(Test[[1]][[j]])))
    {
     d <- c(d,j)
    }

    if(Test[[3]][,(3*j-1)][1]!=1)
    {
      cnt <- cnt + 1
      IonsData[[cnt]]<-Test[[3]][,((3*j-2):(3*j))]
    }
  }

  Test[[1]] = Test[[1]][-d]
  Test[[2]] = Test[[2]][-d]
  Test[[3]] = IonsData
  
  for (h in 1:2)
  {
    for (i in 1:length(Test[[h]]))
    {
      for (j in 1:length(PeakMtx$mass))
      {
        for (k in 1:dim(Test[[h]][[i]])[2])
        {
          if ((Test[[h]][[i]][j,k] == 3)||(Test[[h]][[i]][j,k] == 4)) 
          {
            Test[[h]][[i]][j,k] <- "DownRegulated"
          }
          
          if ((Test[[h]][[i]][j,k] == 1)||(Test[[h]][[i]][j,k] == 2)) 
          {
            Test[[h]][[i]][j,k] <- "UpRegulated"
          }
          
          if (Test[[h]][[i]][j,k]==0)
          {
            Test[[h]][[i]][j,k] <- "--"
          }
        }
      }
    }
  }

  return(Test)
  }

  
