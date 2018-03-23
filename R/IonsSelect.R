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

#' TestIonsSelect.
#'
#' @description Performs a test that calculates and select the up-regulated and/or down-regulated
#' ions from a peak matrix between defernt clusters.
#'
#' @param PeakMtx A matrix of peaks containg mass spectra data
#' @param Probability Probability restriction of the test.
#' @param Segmentation_matrix A matrix containing a segmentated peak matrix.
#' @param SgmtsChoose A vector containing which clusters should be analaized.
#' @return List containing 3 elements. The results from the Volcano test, the results from the Zero test and the values of p, FC & Zero scores. 
#' @export
#'

  TestIonsSelect <- function (PeakMtx, Probability, SegmentedMtx, SgmtsChoose)
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

  Test <- IonsSelectC(m_focalProb = Probability, numPixels = numPixels, SP_Pixels = PeakMtx$numPixels,
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

  return(Test)
  }


