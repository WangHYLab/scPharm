
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
data("sysdata", package = "copykat")
scPharmIdentify.result <- scPharmIdentify(object = example_data, type = "tissue", cancer = "BRCA")
```

##### parameter of scPharmIdentify function

object: a seurat object of patient or cell line cells(example data can be download from https://github.com/WangHYLab/scPharm/blob/main/example_data.rds).  
type: the source of cell, cell line or tumor tissue. Can be set to 'cellline' or 'tissue'.  
cancer: the TCGA cancer type of the cells (eg: BRCA). cancer='pan' means calculating in the context of pan-cancer.  
drug: the drug name for identifying pharmacological cell subpopulations. If not specified all drugs from GDSC2 project will be taken.  
nmcs: number of components to compute and store for MCA. Default:50  
nfeatures: number of genes used to make cell ID.  
cores: number of CPU cores to use. This parameter can only be set to 1 on windows platform. Default:1  
features: character vector of feature names to run MCA. If not specified all features will be taken.  
slot: slot of seurat object used to run MCA. Default:data.  
assay: assay of seurat object used to run MCA. Default:RNA.

##### value

A seurat object added pharmacological label to cells about different drugs in GDSC2 project. Label for drugs in seurat object meta data named with 'scPharm_label_DURGID_DRUGNAME'.

##### Compute Dr(drug prioritization)

```
scPharm.Dr <- scPharmDr(scPharmIdentify.result)
```

##### parameter of scPharmDr function

object: a seurat object of patient or cell line cells after running scPharmIdentify function.

##### value

A data frame consisted of six columns. Column 1 and 2 are id and name of drugs in GDSC2 project. Column 3 and 4 are the size of sensitive and resistant subpopulations to specific drug. Column 5 and 6 are score(Dr) and prioritization to different drugs for specific patient.

##### Identify potential drug combinations

```
scPharm.combo <- scPharmCombo(scPahrmIdentify.result, scPharmDr.result)
```

##### parameter of scPharmCombo function

object: a seurat object after running scPharmIdentify function  
score: output of scPharmDr function.  
drug: name of drug as drug 1. Default:NULL  
topN: number of drugs at top of Dr list as drug 1. Default:1

##### value

A data frame consisted of five columns. Column 1 is name of drug 1. Column 2 and 3 are id and name of drug 2. Column 4 and 5 are the value of combiantion effects and the type of combination.

##### Predict drug side effects

```
scPharm.Dse <- scPharmDse(scPharmIdentify.result)
```

##### parameter of scPharmDse function

object: a seurat object of patient cells after running scPharmIdentify function.

##### value

A data frame consisted of three columns. Column 1 and 2 are id and name of drugs in GDSC2 project. Column 3 is the value of side effects to different drugs.

## References
***

Tian P, Zheng J, Xu Y, et al. scPharm: identifying pharmacological subpopulations of single cells for precision medicine in cancers. Published online December 12, 2023. doi:10.1101/2023.12.11.571182
  
***

