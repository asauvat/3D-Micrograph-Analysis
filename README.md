# Introduction to 3D Micrograph Analysis with R and Python
This project provides a tutorial on segmenting 3D micrographs using both **R** and **Python**.

## System Requirements

### R and Python

To use the resources provided in this project, you will need to install the following software:

- [R](https://cran.r-project.org/)
- [RStudio](https://www.rstudio.com/products/rstudio/download/)
- [Java](https://www.java.com/fr/download/manual.jsp)
- [Anaconda](https://www.anaconda.com/download).

Please follow the distributors' instructions for install. 

> **Note** : The code has been tested only on **Ubuntu/Debian** systems, although it should work on other operating systems as well.

### CONDA environment for R reticulate
In a bash terminal, create and activate a new Conda environment, then install the necessary Python packages:

```sh
conda create --name r-reticulate python=3.10
conda activate r-reticulate
conda install -c conda-forge numpy scipy scikit-image tifffile imagecodecs napari pyqt
```

### R Packages
Before running the R scripts, install the following R packages:

```R
install.packages(pkgs=c('rJava','reticulate','magrittr','xml2'), repos = "http://cloud.r-project.org")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("EBImage","RBioFormats"))
```

## Resources provided
- The `MICROGRAPHS` folder contains pre-processed 3D stacks in OME-TIFF format. These stacks stem from the confocal acquisition of MCA205 cells aggregated into spheroids, after fixation and staining using Hoechst 33342. Cell treatments are indicated in image labels (CTR : control; OXA : oxaliplatin, 24h). Some cells express a fluorescent biomarker, visualized in the FITC channel.

- The file `DIU_SPH3D image analysis.R` is an R script that performs segmentation and data interpretation on the provided images.
