## Dimensionality reduction, clustering, and differential expression with single-cell RNA-seq - GrasPods Tutorial

**If you're reading this before the tutorial, please install the R package dependencies first -- if you're working from a completely clean R install this can take a few minutes.**

### Introduction

The purpose of the tutorial is to provide an introduction to dimensionality reduction, visualization, basic clustering and differential expression methods applicable to single-cell RNA-seq (scRNA-seq) data. We'll be using a set of single-cell RNA-seq data for 2700 PBMCs freely available from 10X genomics (https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.0.0/pbmc3k). 

For the most part, the methods described in this tutorial can also be used in the context of RNA-seq analysis for bulk samples. Certain code chunks are optimized for speed to handle the high volume of data from scRNA-seq. 

### Contents

* `data/pbmc3k_seurat_extracted.rda`: R data object containing the PBMC expression matrix (pre-processed from Seurat)
* `Rmd/tutorial.Rmd`: Rmarkdown file containing the tutorial

### Getting set up

This tutorial assumes basic working knowledge of R and RStudio. All steps are documented in `Rmd/tutorial.Rmd`, and more nicely in `Rmd/tutorial.md`. 

If you run into errors while installing packages, please make sure you have the latest version of RStudio installed.

To render the document:

* Clone the repository, i.e. `git clone https://github.com/Irrationone/graspods_scrnaseq_clustering_tutorial.git`or download a copy from under `Releases`
* Navigate to the cloned repository in RStudio
* Open the report file, and set the `Rmd` directory to your working directory
* Run the chunks under **Setup** to install any packages you might be missing
* Render the report with `rmarkdown::render("tutorial.Rmd", "all")`. This will run through 2 rendering steps to generate a table of contents within the markdown file. 

Please create an issue in this repository (or email `alzhang` at `bcgsc.ca`) if you have any questions/comments. 


