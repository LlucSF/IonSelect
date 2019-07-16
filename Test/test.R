#########################################################################
#
#     Copyright (C) 2019 - Lluc Sementé Fernàndez <lluc.semente@urv.cat>
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

# CAUTION!
# Be aware that all the paths include in this file should be changed for those on your computer
# and following the structures of directory of your OS.
# CAUTION!

#------------------------------------------------------------------------------------------------#
#----------------------------- INSTALLATION -----------------------------------------------------#
#------------------------------------------------------------------------------------------------#

# In order to install and proced with the test you have two options: Installing all the rMSI enviorment
# or just load the data, using the base R methods, included in the .Rdata files included in this folder. 
# Also you can just load the output of the package without installing it from results.Rdata.



#----------------------#
# 1.  Install the rMSI & rMSIproc R packages from GitHub.This will allow you to visualize ions and
#     clusters easily. To install packages from GitHub you first need the devtools package.
  install.packages("devtools")

#     Now install rMSI & rMSIproc packages. Its recomended to follow the instructions of each package,
#     repository but running the following lines should work in most cases.
  devtools::install_github("prafols/rMSI", ref = "0.8")
  devtools::install_github("prafols/rMSIproc", ref = "0.2")
  
#     This package requires the GSL library. If you are using windows you cna obtain more information in the 
#     Readme of the repository. Also you need to install the following packages to link C++ code to R.  
  install.packages("Rcpp")
  install.packages("RcppGSL")
  
#  Finally, install the rMSIKeyIon package.
  devtools::install_github("LlucSF/rMSIKeyIon")
  
  

#----------------------#  
# 2.  In this case you just need to install the following packages and libraries:
#     This package requires the GSL library. If you are using windows you cna obtain more information in the 
#     Readme of the repository. Also you need to install the following packages to link C++ code to R.  
  install.packages("Rcpp")
  install.packages("RcppGSL")

#  Finally, install the rMSIKeyIon package.
  devtools::install_github("LlucSF/rMSIKeyIon")

#------------------------------------------------------------------------------------------------#
#----------------------------- TEST DATA LOADING ------------------------------------------------#
#------------------------------------------------------------------------------------------------#

# First we need to load the peak matrix containing the information of the brain. You can found it
# in the test folder inside IonSelect GitHub's repo

#   Using rMSIproc
peak_matrix <- rMSIproc::LoadPeakMatrix("~/path/to/file/Brain_peak_matrix.zip")
#   Base load method
load("~/path/to/file/peak_matrix.Rdata")


# Second, create a clustering or load the included in the test folder
new_clustering_vector <- kmeans(peak_matrix$intensity, centers = 7)$cluster 
load("~/path/to/file/clustering.Rdata") # Clustering will be stored in clustering_vecto


# With rMSIproc you can visualize the clustering and ions easily
rMSIproc::plotClusterImage(peakMatrix = peak_matrix, clusters = clustering_vector)
rMSIproc::plotPeakImage(peakMatrix = peak_matrix, column = 123)


#------------------------------------------------------------------------------------------------#
#----------------------------- ALGORITHM RUN ----------------------------------------------------#
#------------------------------------------------------------------------------------------------#

# Time to run the algorithm
results <- rMSIKeyIon::TestIonSelect(PeakMtx = peak_matrix,
                                    clusters = clustering_vector,
                                    percentile = c(1,10,10),
                                    zeroThreshold = 2.5e-4)
# The output is a list of two lists. The first one, called "ions" contains the up/down-regulated ions
# found at each comparison. The second one, called "data" contains all the data for each ion at each
# comparison. 

# If you want to load the results without installing any package just run the following line.
load("~/path/to/file/results.Rdata")






