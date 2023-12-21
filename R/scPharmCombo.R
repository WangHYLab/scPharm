#' Identify potential drug combinations
#'
#' @param object Seurat object after running scPharmIdentify function
#' @param score output of running scPhrmDr
#' @param drug name of drug as drug 1. Default:NULL
#' @param topN number of drugs at top of Dr list as drug 1. Default:1
#'
#' @importFrom mixtools normalmixEM
#' @importFrom utils data
#'
#' @return a list consist of data.frame consist of drug combinations and  its effects value
#' @export
#' @encoding UTF-8
#'
#' @examples
#' \dontrun{
#' combo <- scPharmCombo(scPharmIdentify.result, scPharmDr.result)
#' }
scPharmCombo <- function(object,
                         score,
                         drug = NULL,
                         topN = 1) {
  # identify potential drug combinations
  drug_id = drug_info
  if (is.null(drug)) {
    score = score[seq(1, topN),]
  }else {
    score = score[score$DRUG_NAME == drug,]
    if (nrow(score) < 1) {
      message(paste("ERROR: ", drug))
      return(NULL)
    }
  }
  combo.result = list()
  meta.data = object@meta.data
  nes.data = meta.data[, grep("scPharm_nes", colnames(meta.data), value = T)]
  meta.data = meta.data[, grep("scPharm_label", colnames(meta.data), value = T)]
  for (i in seq(1, nrow(score))) {
    col.name = c(paste("scPharm_label", score[i, "DRUG_ID"], score[i, "DRUG_NAME"], sep = "_"),
                 paste("scPharm_nes", score[i, "DRUG_ID"], score[i, "DRUG_NAME"], sep = "_"))
    # Compensation effects
    out.mix = normalmixEM(nes.data[,col.name[2]])
    sub.meta1 = meta.data[meta.data[,col.name[1]] == "resistant",]
    if (max(out.mix$mu) > 0 & nrow(sub.meta1) != 0) {
      combo.overlap = colMeans(sub.meta1 == "sensitive")
      if (combo.overlap[col.name[1]] != 0) {
        warning("calculate error")
      }
      combo.overlap = sort(combo.overlap, decreasing = T)
      combo.overlap = combo.overlap[1:5]
      combo.overlap = data.frame(DRUG_FIRST = score[i, "DRUG_NAME"],
                                 DRUG_ID = sapply(strsplit(names(combo.overlap), split = "_"), function(x) {x[3]}),
                                 DRUG_NAME = sapply(strsplit(names(combo.overlap), split = "_"), function(x) {
                                   if (length(x) ==4) {
                                     return(x[4])
                                   }else {
                                     return(paste(x[4:length(x)], collapse = "_"))
                                   }
                                 }),
                                 Effect = combo.overlap,
                                 Strategy = "compensation effects")
    }else {
      message("No resistant nes peak")
      combo.overlap = NULL
    }

    # Booster effects
    target.path = drug_id[drug_id$DRUG_ID == score[i, "DRUG_ID"],]$PATHWAY_NAME
    # print(target.path)
    if (!grepl("Other|Unclassified", target.path)) {
      sub.meta2 = meta.data[meta.data[,col.name[1]] == "sensitive",]
      combo.enhance = colMeans(sub.meta2 == "sensitive")
      if (combo.enhance[col.name[1]] != 1) {
        warning("calculate error")
      }
      same.path.drug = drug_id[which(drug_id$PATHWAY_NAME == target.path),c("DRUG_ID","DRUG_NAME")]
      # print(same.path.drug)
      to.drop = paste("scPharm_label", same.path.drug$DRUG_ID, same.path.drug$DRUG_NAME, sep = "_")
      # print(to.drop)
      combo.enhance = combo.enhance[-which(names(combo.enhance) %in% to.drop)]
      # print(length(combo.enhance))
      combo.enhance = sort(combo.enhance, decreasing = T)
      combo.enhance = combo.enhance[1:5]
      combo.enhance = data.frame(DRUG_FIRST = score[i, "DRUG_NAME"],
                                 DRUG_ID = sapply(strsplit(names(combo.enhance), split = "_"), function(x) {x[3]}),
                                 DRUG_NAME = sapply(strsplit(names(combo.enhance), split = "_"), function(x) {
                                   if (length(x) ==4) {
                                     return(x[4])
                                   }else {
                                     return(paste(x[4:length(x)], collapse = "_"))
                                   }
                                 }),
                                 Effect = combo.enhance,
                                 Strategy = "booster effects")
      if (!is.null(combo.overlap)) {
        combo = rbind(combo.enhance, combo.overlap)
      }else {
        combo = combo.enhance
      }
      combo = combo[order(combo$Effect, decreasing = T),]
      # combo = combo[1:5,]
    }else {
      combo = combo.overlap
    }
    combo.result[[paste(score$DRUG_ID[i], score$DRUG_NAME[i], sep = "_")]] = combo
  }
  return(combo.result)
}
