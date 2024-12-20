% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scPharmIdentify.R
\encoding{UTF-8}
\name{scPharmIdentify}
\alias{scPharmIdentify}
\title{Identify pharmacological cell subpopulations}
\usage{
scPharmIdentify(
  object,
  type,
  cancer,
  drug = NULL,
  nmcs = 50,
  nfeatures = 200,
  cores = 1,
  features = NULL,
  slot = "data",
  assay = "RNA",
  threshold.s = -1.751302,
  threshold.r = 1.518551,
  bulkdata = bulkdata,
  gdscdata = gdscdata
)
}
\arguments{
\item{object}{Seurat object}

\item{type}{the source of cell, cell line or tumor tissue. Can be set to 'cellline' or 'tissue'.}

\item{cancer}{the TCGA cancer type of cell. A character or vector.(eg: BRCA or c('LUAD', 'LUSC')). cancer='pan' means calculating in the context of pan-cancer.}

\item{drug}{the drug name for identifying. If not specified all drugs from GDSC2 project will be taken.}

\item{nmcs}{number of components to compute and store for MCA. Default:50}

\item{nfeatures}{number of genes used to make cell ID}

\item{cores}{number of CPU cores to use. This parameter can only be set to 1 on windows. Default:1}

\item{features}{character vector of feature names. If not specified all features will be taken.}

\item{slot}{slot of seurat object used to run MCA}

\item{assay}{assay of seurat object used to run MCA}

\item{threshold.s}{the threshold to label sensitive cells. Default: -1.751302}

\item{threshold.r}{the threshold to label resistant cells. Default: 1.518551}

\item{bulkdata}{file for all cancer cell lines, automatically loaded in scPharm.}

\item{gdscdata}{pharmacological file for all cancer cell lines, automatically loaded in scPharm.}
}
\value{
a seurat object added new meta data
}
\description{
Identify pharmacological cell subpopulations
}
\examples{
\dontrun{
result <- scPharmIdentify(seurat.object, type = "tissue", cancer = "LUAD")
}
}
