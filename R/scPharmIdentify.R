


#' Identify pharmacological cell subpopulations
#'
#' @param object Seurat object
#' @param type the source of cell, cell line or tumor tissue. Can be set to 'cellline' or 'tissue'. Default:tissue
#' @param cancer the cancer type of cell. cancer='pan' means calculating in the context of pan-cancer.
#' @param drug the drug name for identifying. If not specified all drugs from GDSC2 project will be taken.
#' @param nmcs number of components to compute and store for MCA, default set to 50
#' @param nfeatures number of genes used to make cell ID
#' @param cores number of CPU cores to use
#' @param features character vector of feature names. If not specified all features will be taken.
#' @param slot slot of seurat object used to run MCA
#' @param assay assay of seurat object used to run MCA
#'
#' @import Seurat
#' @importFrom copykat copykat
#' @importFrom sparseMatrixStats rowVars
#' @importFrom stringr str_length
#' @importFrom irlba irlba
#' @importFrom SeuratObject CreateDimReducObject
#' @importFrom CelliD GetCellGeneSet
#' @importFrom stats cor.test
#' @importFrom fgsea fgseaMultilevel
#' @import dplyr
#' @importFrom utils data
#'
#' @return a seurat object added new meta data
#' @export
#'
#' @useDynLib scPharm
#' @encoding UTF-8
#'
#' @examples
#' \dontrun{
#' result <- scPharmIdentify(seurat.object, type = "tissue", cancer = "LUAD")
#' }
scPharmIdentify <- function(object,
                            type,
                            cancer,
                            drug = NULL,
                            nmcs = 50,
                            nfeatures = 200,
                            cores = 4,
                            features = NULL,
                            slot = "data",
                            assay = "RNA") {

  # classification class of object
  if(class(object)[1]=="Seurat") {
    InitAssay <- DefaultAssay(object)
    DefaultAssay(object) <- assay
    mat <- GetAssayData(object, slot = slot)
  }else {
    message("ERROR:: Not a Seurat object!")
    return(NULL)
  }

  # identify tumor or adjacent cell
  if (type == "tissue") {
    message("Identify tumor or adjacent cell")
    # tic()
    exp.rawdata <- GetAssayData(object, slot = "counts")
    out <- file("scPharm_copykat_out.txt", open = "w")
    sink(out, type = "output")
    copykat.result <- copykat(exp.rawdata,
                              id.type = "S",
                              ngene.chr = 3,
                              win.size = 25,
                              KS.cut = 0.1,
                              sam.name = "scPharm",
                              distance = "euclidean",
                              n.cores = cores,
                              norm.cell.names = "",
                              output.seg = FALSE)
    sink()
    close(out)
    # toc()
    object@meta.data[,"cell.label"] = "adjacent"
    pred.label <- data.frame(copykat.result$prediction)
    tumor.cells <- pred.label$cell.names[which(pred.label$copykat.pred=="aneuploid")]
    object@meta.data[tumor.cells, "cell.label"] = "tumor"
    system("rm -f *_copykat_*")
  }else {
    object@meta.data[,"cell.label"] = "tumor"
  }

  # MCA
  reduction.name = "mca"
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
  cell.geneset = suppressMessages(GetCellGeneSet(object,
                                                 reduction = "mca",
                                                 dims = seq(1,nmcs),
                                                 n.features = nfeatures))
  bulk_data = bulkdata
  GDSC = gdscdata

  # decipher response
  if (cancer != 'pan') {
    if (cancer %in% unique(GDSC$TCGA_DESC)) {
      GDSC <- GDSC[GDSC$TCGA_DESC == cancer,]
    }else {
      message("ERROR: cancer type")
    }
  }
  drug_id <- GDSC[,c(8,9,10,11)]
  drug_id <- drug_id[!duplicated(drug_id$DRUG_ID),]
  drug_id <- dplyr::filter(drug_id, !dplyr::if_all(.fns = is.na))
  if (!is.null(drug)) {
    drug_id = drug_id[drug_id$DRUG_NAME == drug,]
    if (nrow(drug_id) == 0) {
      return(message("ERROR:: No drug data about such cancer!"))
    }
  }
  meta.data = object@meta.data

  #scPharm
  message("GSEA")
  i = 1
  while (i <= nrow(drug_id)) {
    message(paste(" ", drug_id$DRUG_ID[i], drug_id$DRUG_NAME[i], " ", sep = ""))
    drug_data = GDSC[GDSC$DRUG_ID == drug_id$DRUG_ID[i], c(5,17)]
    drug_data = drug_data[which(drug_data$CELL_LINE_NAME %in% colnames(bulk_data)),]
    drug_data = drug_data[order(drug_data$AUC, decreasing = TRUE),]
    if (nrow(drug_data) < 3) {
      message(paste0("Insufficient amount of data to ", drug_id$DRUG_ID[i], " ", drug_id$DRUG_NAME[i]))
      i = i+1
      next
    }
    sub_bulk_data = bulk_data[,drug_data$CELL_LINE_NAME]
    sub_bulk_data = sub_bulk_data[rowSums(sub_bulk_data) > 0, ]
    # compute drug response-determined ranking gene list
    geneList = apply(sub_bulk_data, 1, function(expr) {
      return(cor.test(expr, drug_data$AUC, method = "pearson", alternative = 'two.sided')$estimate)
    })
    geneList = sort(geneList, decreasing = TRUE)
    enrich_drug = tryCatch({
      suppressWarnings(fgseaMultilevel(pathways = cell.geneset,
                                       stats    = geneList,
                                       eps = 0,
                                       minSize  = 10,
                                       maxSize  = 500,
                                       nproc = cores,
                                       nPermSimple = 10000))
    }, error=function(e) {
      2
    })
    if (is.numeric(enrich_drug)) {
      message("Re-do")
    }else {
      message("Right running")
      enrich_drug[is.na(enrich_drug)] = 0
      # saveRDS(enrich_drug, file = "./result/fgsea.rds")
      if (nrow(enrich_drug) == 0) {
        label = data.frame(cell = rownames(meta.data), nes = 0)
        rownames(label) = label$cell
        label$label = "other"
      }else {
        label = data.frame(cell = rownames(meta.data), nes = 0)
        rownames(label) = label$cell
        label[enrich_drug$pathway, 2] = enrich_drug[,6]
        label$label = "other"
        label[label$nes > 1.518551, 3] = "resistant"
        label[label$nes < -1.751302, 3] = "sensitive"
      }
      col_name = c(paste("scPharm_label", drug_id$DRUG_ID[i], drug_id$DRUG_NAME[i], sep = "_"),
                   paste("scPharm_nes", drug_id$DRUG_ID[i], drug_id$DRUG_NAME[i], sep = "_"))
      meta.data[,col_name] = label[, c('label','nes')]
      i = i+1
    }
  }
  object@meta.data = meta.data
  return(object)
}
