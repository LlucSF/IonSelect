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

#------------------------------------------------------------------------------------------------#
#----------------------------- INSTALLATION -----------------------------------------------------#
#------------------------------------------------------------------------------------------------#

# In order to proced with the test you need to install the rMSI & rMSIproc R packages from GitHub.
# To install packages from GitHub you first need the devtools package.
install.packages("devtools")

# Now install rMSI & rMSIproc packages. Its recomended to follow the instructions of each package,
# but running the following lines should work in most cases.
devtools::install_github("prafols/rMSI", ref = "0.8")
devtools::install_github("prafols/rMSIproc", ref = "0.2")

# Now install the IonSelect package.
devtools::install_github("LlucSF/IonSelect")

#------------------------------------------------------------------------------------------------#
#----------------------------- DATA LOADING -----------------------------------------------------#
#------------------------------------------------------------------------------------------------#

# First we need to load the peak matrix containing the information of the brain. You can found it
# in the test folder inside IonSelect GitHub's repo


peak_matrix <- rMSIproc::LoadPeakMatrix("~/path/to/file/Brain_peak_matrix.zip")

# Create a clustering or load the included in the test folder

new_clustering_vector <- kmeans(peak_matrix$intensity, centers = 7)$cluster 
load("~/path/to/file/clustering.Rdata") # Clustering will be stored in clustering_vector

rMSIproc::plotClusterImage(peakMatrix = peak_matrix, clusters = clustering_vector)
#------------------------------------------------------------------------------------------------#
#----------------------------- ALGORITHM RUN ----------------------------------------------------#
#------------------------------------------------------------------------------------------------#

# Time to run the algorithm

results <- IonSelect::TestIonSelect(PeakMtx = peak_matrix,
                                    clusters = clustering_vector,
                                    percentile = 1,
                                    zeroThreshold = 2.5e-4)



results <- IonSelect::TestIonSelect(PeakMtx = peak_matrix,
                                    clusters = new_clustering_vector,
                                    percentile = 1,
                                    zeroThreshold = 2.5e-4)










