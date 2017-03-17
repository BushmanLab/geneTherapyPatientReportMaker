#    This source code file is a component of the larger INSPIIRED genomic analysis software package.
#    Copyright (C) 2016 Frederic Bushman
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

#' list of the most abundant genes and absolute abundance cutoff
#'
#' @param sites df with columns: estAbund, nearest_refSeq_gene
#' @param numGenes number of most abundant genes to keep
#' @return list with 2 elements: abundanceCutoff( cut off value - absolute abundance), 
#' topGenes(character vector of gene names)
getMostAbundantGenes <- function(sites, numGenes) {
    orderedSites <- sites[order(-sites$estAbund)]
    #'unique' functionally removes anything that's duplicated, thus preserving order
    topGenes <- head(unique(orderedSites$nearest_refSeq_gene), numGenes)
    cutoff <- .getAbsoluteAbundanceThreshold(orderedSites, topGenes)
    list(abundanceCutoff=cutoff, topGenes=topGenes)
}

.getAbsoluteAbundanceThreshold <- function(orderedSites, topGenes){
    abundCutoff <- orderedSites$nearest_refSeq_gene==tail(topGenes, 1)
    orderedSites[which(abundCutoff)[1]]$estAbund  
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
  
  ###browser()
  
  do.call(rbind, lapply(splitSites, function(sites){
    res <- aggregate(estAbundProp ~ maskedRefGeneName, mcols(sites), sum)
    est <- aggregate(estAbund ~ maskedRefGeneName, mcols(sites), sum)
    res <- join(res, est, by = "maskedRefGeneName")
    res$Timepoint <- sites[1]$Timepoint
    res$CellType <- sites[1]$CellType
    
    res
  }))
}
