getAbundanceThreshold <- function(sites, numGenes){
  orderedAbundances <- sites[order(-sites$estAbundProp)]
  #'unique' functionally removes anything that's duplicated, thus preserving order
  abundCutoff <- orderedAbundances$nearest_refSeq_gene==tail(head(unique(orderedAbundances$nearest_refSeq_gene), numGenes),1)
  abundCutoff <- orderedAbundances[which(abundCutoff)[1]]$estAbundProp  
  abundCutoff
}

#' list of most abundant genes and abundance proportion cutoff
#'
#' @param sites df with columns: estAbundProp, nearest_refSeq_gene
#' @param numGenes number of genes with highest abundance proportions
#' @return list with 2 elements: abundanceCutoff( cut off value), 
#' topGenes(vector of gene names)
getMostAbundantGenes <- function(sites, numGenes) {
    cutoff <- getAbundanceThreshold(sites, numGenes)
    orderedAbundances <- sites[order(-sites$estAbundProp)]
    topGenes <- head(unique(orderedAbundances$nearest_refSeq_gene), numGenes)
    list(abundanceCutoff=cutoff, topGenes=topGenes)
}

filterLowAbund <- function(sites, abundCutoff){
  sites$maskedRefGeneName <- ifelse(sites$estAbundProp >= abundCutoff,
                                    sites$nearest_refSeq_gene,
                                    "LowAbund")
  sites
}

#' all the genes that are not from white list masked as 'LowAbund'
#'
#' @param df with column: nearest_refSeq_gene
#' @param genes_to_keep genes that do NOT need to be masked
#' @return dataframe the same as input with additional column 'maskedRefGeneName'
maskGenes <- function(sites, genes_to_keep) {
    sites$maskedRefGeneName <- ifelse(sites$nearest_refSeq_gene %in% genes_to_keep, 
        sites$nearest_refSeq_gene, 'LowAbund')
    sites
}

getAbundanceSums <- function(sites, splitBy){
  splitBy <- mcols(sites)[,splitBy]
  splitSites <- split(sites, apply(as.data.frame(splitBy), 1, paste, collapse=""))
  
  do.call(rbind, lapply(splitSites, function(sites){
    res <- aggregate(estAbundProp~maskedRefGeneName, mcols(sites), sum)
    res$Timepoint <- sites[1]$Timepoint
    res$CellType <- sites[1]$CellType
    
    res
  }))
}
