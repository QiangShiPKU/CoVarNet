# **Welcome to the CoVarNet!**
CoVarNet is a computational framework aiming to unravel the coordination among multiple cell types by analyzing the covariance in the frequencies of cell types across various samples.

For more details, see our Nature publication: https://www.nature.com/articles/s41586-025-09053-4


## **Installation**
```
devtools::install_github(repo = "https://github.com/QiangShiPKU/CoVarNet")
library(CoVarNet)
```


## **Tutorials**
* [Discovery of cellular modules in scRNA-seq data](https://htmlpreview.github.io/?https://github.com/QiangShiPKU/CoVarNet/blob/main/vignette/tutorial_discovery.html)
* [Recovery of cellular modules in scRNA-seq data and spatial transcriptomics data](https://htmlpreview.github.io/?https://github.com/QiangShiPKU/CoVarNet/blob/main/vignette/tutorial_recovery.html)
* [Trajectory inference for individuals](https://htmlpreview.github.io/?https://github.com/QiangShiPKU/CoVarNet/blob/main/vignette/tutorial_trajectory.html)


## **Requirements**
The R/Python packages listed below are required for running CoVarNet. These versions are used for testing the CoVarNet code. Other versions might work too.
* R (v4.1.2).
* R packages: dplyr(v1.1.4), NMF(v0.30.1), Seurat(v5.1.0), cluster(v2.1.6), sp(2.1-4), spdep(v1.3-5), igraph(v1.6.0), circlize(v0.4.15), ComplexHeatmap (v2.15.4), ggsci(v3.0.3), grid(v4.1.2), psych(v2.4.3), RColorBrewer(v1.1-3), ggplot2(v3.5.0), viridis(v0.6.5), tidytext(v0.4.1), dendextend(v1.17.1), anndata(v0.7.5.6), reticulate(v1.40.0).
* Python (v3.12.2, only for Tutorial 3).
* Python packages: Scanpy (v1.11.0), Palantir (v1.3.3)
