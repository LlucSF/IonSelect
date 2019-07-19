# rMSIKeyIon
R package for selecting statistically significant ions between different regions in MSI datasets. It uses the classical p-value/ FC method and introduces a third novel method to treat undetected MS peaks more accurately.


## Installation

In order to install the package you have two options: Install all the rMSI environment (steps 1 and 2)
or just load the data (using the base R methods) included in the .Rdata files in the test folder (just step 2). 
Alternatively, you can just load the output of the package without installing it from results.Rdata if you are only interested in checking the demo results

1.Install the rMSI & rMSIproc R packages from GitHub.This will allow you to visualize ions and clusters easily. To install packages from GitHub you first need the devtools package.
```
install.packages("devtools")
```

Now install rMSI & rMSIproc packages. Its recomended to follow the instructions of each package, repository but running the following lines.
```
devtools::install_github("prafols/rMSI", ref = "0.8")
devtools::install_github("prafols/rMSIproc", ref = "0.2")
```

2.Install the following packages and libraries.This package requires the GSL library. If you are using Windows, you can obtain more information on how to install it below. Also you need to install the following packages to link the C++ code to R.  
```
  install.packages("Rcpp")
  install.packages("RcppGSL")
```  

Finally, install the rMSIKeyIon package.
```
  devtools::install_github("LlucSF/rMSIKeyIon")
```
### GSL for windows
This package requires the GSL library. Windows users should download a binary version of it. 
In the following link you can find information about installing it and a download repository. http://joonro.github.io/blog/posts/installing-gsl-and-cythongsl-in-windows.html#. Linux users can obtain GSL using the default package manager of its distribution.

## Testing

First, we need to load the peak matrix containing the information of the MSI data. You can get a sample dataset from the test folder of this github repository.
```
# Using rMSIproc
peak_matrix <- rMSIproc::LoadPeakMatrix("~/path/to/file/Brain_peak_matrix.zip")

#Base load method
load("~/path/to/file/peak_matrix.Rdata")
```

Then, create an image segmentation by clustering or using regions of interest (ROI). Alternatively, you can load the included in the demo data folder.
```
new_clustering_vector <- kmeans(peak_matrix$intensity, centers = 7)$cluster 
load("~/path/to/file/clustering.Rdata") # Clustering will be stored in clustering_vector
```

With rMSIproc you can visualize the clustering and ions easily
```
rMSIproc::plotClusterImage(peakMatrix = peak_matrix, clusters = clustering_vector)
rMSIproc::plotPeakImage(peakMatrix = peak_matrix, column = 123)
```


Finally, it is time to run the rMSIKeyIon algorithm
```
results <- rMSIKeyIon::TestIonSelect(PeakMtx = peak_matrix,
                                    clusters = clustering_vector,
                                    percentile = c(1,10,10),
                                    zeroThreshold = 2.5e-4)
```
rMSIproc outputs a list made up of two R lists. The first one, called "ions" contains the up/down-regulated ions
found at each comparison. The second one, called "data" contains all the data belonging to each ion at each
comparison. 

If you are only interested on visualizing the demo data results without installing any package, yo can just run the following line.
```
load("~/path/to/file/results.Rdata")
```
