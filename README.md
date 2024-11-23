# **Welcome to the CoVarNet!!**
CoVarNet is a computational framework aiming to unravel the coordination among multiple cell types by analyzing the covariance in the frequencies of cell types across various samples.


## **Installation**
```
devtools::install_github(repo = "https://github.com/QiangShiPKU/CoVarNet")
library(CoVarNet)
```


## **Tutorials**
* [Discovery of cellular modules in scRNA-seq data](./vignette/tutorial_discovery.html)
* [Recovery of cellular modules in scRNA-seq data and spatial transcriptomics data](./vignette/tutorial_recovery.html)


## **Requirements**
The R packages listed below are required for running CoVarNet. The version numbers indicate the package versions used for testing the CoVarNet code. Other R versions might work too.
* R (v4.1.2).
* R packages: dplyr(v1.1.4), NMF(v0.30.1), Seurat(v5.1.0), cluster(v2.1.6), sp(2.1-4), spdep(v1.3-5), igraph(v1.6.0), circlize(v0.4.15), ComplexHeatmap (v2.15.4), ggsci(v3.0.3), grid(v4.1.2), psych(v2.4.3), RColorBrewer(v1.1-3), ggplot2(v3.5.0), viridis(v0.6.5), tidytext(v0.4.1), dendextend(v1.17.1).

