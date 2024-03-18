# **Welcome to the CoVarNet!**
CoVarNet is a method that uses the frequencies of individual cell types/subsets in samples to identify co-varing multicellular networks or named as cellular modules (CMs).

## Installation 
```
devtools::install_github(repo = "https://github.com/QiangShiPKU/CoVarNet")
library(CoVarNet)
library(NMF)
```

## **Quick start**
### 1. Get nodes of each network
```
node <- getNode(mat, K = 12)
```
**Input:** a matrix (mat) indicating the frequencies of individual cell types/subsets (rows) across samples (columns) and a integer (K) indicating the expected number of CMs to be identified. The example dataset can be accessed by
```
mat <- mat_freq_raw
```
**Output:** a data frame listing nodes and weights for individual networks.

### 2. Get edge of networks
```
edge <- getEdge(mat)
```
**Output:** a data frame listing several attributes of network edges.

### 3. Visualize networks
```
creNet(node, mat)
```
**Output:** visualize networks.

## **Requirements**
* R (R version >=3.5.0).
* R libraries: NMF, dplyr, psych, reshape2, igraph, circlize, ggsci
