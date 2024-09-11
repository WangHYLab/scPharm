#' Predict drug side effects
#'
#' @param object Seurat object after running scPharmIdentify function
#'
#' @return a data.frame consist of drug information and drug Dse
#' @export
#' @encoding UTF-8
#' @examples
#' \dontrun{
#' Dse <- scPharmDse(scPharmIdentify.result)
#' }
scPharmDse <- function(object) {
  # predict drug side effects
  message("computing Dse")
  meta.data <- object@meta.data
  result <- meta.data[, c("cell.label", grep("scPharm_label", colnames(meta.data), value = T))]
  result <- result[result$cell.label == "adjacent", !names(result) %in% "cell.label"]
  sensi.ratio <- colMeans(result == "sensitive")
  # resis.ratio = colMeans(result == "resistant")
  score <- data.frame(
    DRUG_ID = sapply(strsplit(names(sensi.ratio), split = "_"), function(x) {
      x[3]
    }),
    DRUG_NAME = sapply(strsplit(names(sensi.ratio), split = "_"), function(x) {
      if (length(x) == 4) {
        return(x[4])
      } else {
        return(paste(x[4:length(x)], collapse = "_"))
      }
    }),
    Dse = sensi.ratio
  )
  score <- score[order(-score[, 3]), ]
  # score$Rank = seq(1,nrow(score),1)
  return(score)
}
