#' Generate the null distribution and calculate the thresholds for labeling sensitive and resistant cells
#'
#' @param object Seurat object
#' @param cancer the TCGA cancer type associated with your cells for constructing null distribution. A character or vector.(eg: BRCA or c('LUAD', 'LUSC')). cancer='pan' means calculating in the context of pan-cancer.
#' @param nmcs number of components to compute and store for MCA. Default:50
#' @param nfeatures number of genes used to make cell ID
#' @param cores number of CPU cores to use. This parameter can only be set to 1 on windows. Default:1
#' @param features character vector of feature names for MCA. If not specified all features will be taken.
#' @param slot slot of seurat object used to run MCA
#' @param assay assay of seurat object used to run MCA
#' @param bulkdata file for all cancer cell lines, automatically loaded in scPharm.
#' @param gdscdata pharmacological file for all cancer cell lines, automatically loaded in scPharm.
#'
#' @import Seurat
#' @importFrom sparseMatrixStats rowVars
#' @importFrom stringr str_length
#' @importFrom irlba irlba
#' @importFrom SeuratObject CreateDimReducObject
#' @importFrom CelliD GetCellGeneSet
#' @importFrom stats cor.test
#' @importFrom fgsea fgseaMultilevel
#' @import dplyr
#' @importFrom utils data
#' @importFrom tidyr gather
#' @importFrom mixtools normalmixEM
#'
#' @return a list contains the null distribution(NullDist), thresholds for labeling sensitive and resistant cells(threshold_s and threshold_r).
#' @export
#'
#' @useDynLib scPharm
#' @encoding UTF-8
#'
#' @examples
#' \dontrun{
#' result <- scPharmIdentify(seurat.object, type = "tissue", cancer = "LUAD")
#' }
scPharmGenNullDist <- function(object,
                               cancer,
                               nmcs = 50,
                               nfeatures = 200,
                               cores = 1,
                               features = NULL,
                               slot = "data",
                               assay = "RNA",
                               bulkdata = bulkdata,
                               gdscdata = gdscdata) {
  # classification class of object
  if (class(object)[1] == "Seurat") {
    InitAssay <- DefaultAssay(object)
    DefaultAssay(object) <- assay
    mat <- GetAssayData(object, slot = slot)
  } else {
    message("ERROR:: not a Seurat object!")
    return(NULL)
  }


  # MCA
  reduction.name <- "mca"
  if (!is.null(features)) {
    mat <- mat[features, ]
  }
  mat <- mat[rowVars(mat) != 0, ]
  mat <- mat[str_length(rownames(mat)) > 0, ]
  mat <- mat[!duplicated(rownames(mat)), ]
  cellsN <- colnames(mat)
  featuresN <- rownames(mat)
  # tic()
  message("Computing Fuzzy Matrix")
  MCAPrepRes <- SparseMCAStep1(mat)
  # toc()
  message("Computing SVD")
  # tic()
  SVD <- irlba(A = MCAPrepRes$Z, nv = nmcs + 1, nu = 1)[seq(3)]
  # toc()
  message("Computing Coordinates")
  # tic()
  MCA <- MCAStep2(Z = MCAPrepRes$Z, V = SVD$v[, -1], Dc = MCAPrepRes$Dc)
  component <- paste0(reduction.name, "_", seq(ncol(MCA$cellsCoordinates)))
  colnames(MCA$cellsCoordinates) <- component
  colnames(MCA$featuresCoordinates) <- component
  rownames(MCA$cellsCoordinates) <- cellsN
  rownames(MCA$featuresCoordinates) <- featuresN
  MCA$stdev <- SVD$d[-1]
  class(MCA) <- "MCA"
  # toc()

  # add MCA to Seurat
  geneEmb <- MCA$featuresCoordinates
  cellEmb <- MCA$cellsCoordinates
  stdev <- MCA$stdev

  colnames(cellEmb) <- paste0(reduction.name, "_", seq(ncol(cellEmb)))
  colnames(geneEmb) <- paste0(reduction.name, "_", seq(ncol(geneEmb)))
  DimReducObject <- CreateDimReducObject(embeddings = cellEmb, loadings = geneEmb, key = paste0(reduction.name, "_"), assay = assay)
  object@reductions[[reduction.name]] <- DimReducObject
  if (!is.null(stdev)) {
    object@reductions[[reduction.name]]@stdev <- sqrt(stdev)
  }
  object@reductions[[reduction.name]]@misc[["mca.flag"]] <- TRUE
  DefaultAssay(object) <- InitAssay

  # calculate ID signatures for cell
  message("Computing cell ID geneset")
  cell.geneset <- suppressMessages(GetCellGeneSet(object,
    reduction = "mca",
    dims = seq(1, nmcs),
    n.features = nfeatures
  ))
  bulk_data <- bulkdata
  GDSC <- gdscdata

  # decipher response
  if (is.character(cancer)) {
    if (cancer != "pan") {
      if (cancer %in% unique(GDSC$TCGA_DESC)) {
        GDSC <- GDSC[GDSC$TCGA_DESC == cancer, ]
      } else {
        message("ERROR: unrecognized cancer type")
      }
    }
  } else if (is.vector(cancer)) {
    if (cancer %in% unique(GDSC$TCGA_DESC)) {
      GDSC <- GDSC[which(GDSC$TCGA_DESC %in% cancer), ]
    } else {
      message("ERROR: contains the wrong cancer type")
    }
  } else {
    message("ERROR: unrecognized parameter")
  }

  drug_id <- GDSC[, c(8, 9, 10, 11)]
  drug_id <- drug_id[!duplicated(drug_id$DRUG_ID), ]
  drug_id <- dplyr::filter(drug_id, !dplyr::if_all(.fns = is.na))

  meta.data <- object@meta.data

  # scPharm
  message("GSEA")
  i <- 1
  while (i <= nrow(drug_id)) {
    message(paste(" ", drug_id$DRUG_ID[i], drug_id$DRUG_NAME[i], " ", sep = ""))
    drug_data <- GDSC[GDSC$DRUG_ID == drug_id$DRUG_ID[i], c(5, 17)]
    drug_data <- drug_data[which(drug_data$CELL_LINE_NAME %in% colnames(bulk_data)), ]
    drug_data <- drug_data[order(drug_data$AUC, decreasing = TRUE), ]
    if (nrow(drug_data) < 3) {
      message(paste0("Insufficient amount of data to ", drug_id$DRUG_ID[i], " ", drug_id$DRUG_NAME[i]))
      i <- i + 1
      next
    }
    sub_bulk_data <- bulk_data[, drug_data$CELL_LINE_NAME]
    sub_bulk_data <- sub_bulk_data[rowSums(sub_bulk_data) > 0, ]
    # compute drug response-determined ranking gene list
    geneList <- apply(sub_bulk_data, 1, function(expr) {
      return(cor.test(expr, drug_data$AUC, method = "pearson", alternative = "two.sided")$estimate)
    })
    geneList <- sort(geneList, decreasing = TRUE)
    enrich_drug <- tryCatch(
      {
        suppressWarnings(fgseaMultilevel(
          pathways = cell.geneset,
          stats = geneList,
          eps = 0,
          minSize = 10,
          maxSize = 500,
          nproc = cores,
          nPermSimple = 10000
        ))
      },
      error = function(e) {
        2
      }
    )
    if (is.numeric(enrich_drug)) {
      message("Re-do")
    } else {
      message("Right running")
      enrich_drug[is.na(enrich_drug)] <- 0
      # saveRDS(enrich_drug, file = "./result/fgsea.rds")
      if (nrow(enrich_drug) == 0) {
        label <- data.frame(cell = rownames(meta.data), nes = 0)
        rownames(label) <- label$cell
      } else {
        label <- data.frame(cell = rownames(meta.data), nes = 0)
        rownames(label) <- label$cell
        label[enrich_drug$pathway, 2] <- enrich_drug[, 6]
      }
      col_name <- paste("scPharm_nes", drug_id$DRUG_ID[i], drug_id$DRUG_NAME[i], sep = "_")
      meta.data[, col_name] <- label[, "nes"]
      i <- i + 1
    }
  }
  meta.data <- meta.data[, grep("scPharm_nes", colnames(meta.data), value = T)]
  meta.data.long <- gather(meta.data, key = "Drug", value = "NES")
  out.mix <- normalmixEM(meta.data.long$NES)
  threshold.s <- out.mix$mu[1] - out.mix$sigma[1]
  threshold.r <- out.mix$mu[2] + out.mix$sigma[2]

  return(list(NullDist = meta.data.long$NES, threshold_r = threshold.r, threshold_s = threshold.s))
}
