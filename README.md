
# scPharm
***

## Overview
***

R package for pharmacological cell subpopulations identification, drug prioritization and potetial drug combinations identification from single-cell RNA-seq of tumour tissue.
scPharm is a computational framework tailored for scRNA-seq data, integrating pharmacogenomics profiles to uncover therapeutic heterogeneity within tumours at single-cell resolution. The tool not only prioritizes tailored drugs but also provides insight into combination therapy regimen and drug toxicity in cancers. 

## Installation
***

Within R, set first:
```
install.packages("remotes")
```
To install scPharm then just type:
```
remotes::install_github("WangHYLab/scPharm", dependencies = T)
```

## Data inputs format
***

scPharm use as input single cell data in the form of specific S4 objects.

## Usage
***

#### Function examples:

##### Identify pharmacological cell subpopulations

```
result <- scPharmIdentify(seurat.object, type = "tissue", cancer = "LUAD")
```

##### Compute Dr(drug prioritization)

```
Dr <- scPharmDr(scPharmIdentify.result)
```

##### Identify potential drug combinations

```
combo <- scPharm(scPahrmIdentify.result, scPharmDr.result)
```

##### Predict drug side effects

```
Dse <- scPharmDse(scPharmIdentify.result)
```

## References
***

Tian P, Zheng J, Xu Y, et al. scPharm: identifying pharmacological subpopulations of single cells for precision medicine in cancers. Published online December 12, 2023. doi:10.1101/2023.12.11.571182
        
        
  
***

