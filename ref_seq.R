#' read Ref Seq gene names and locations from file or download it
#'
#' @param ref_seq_file filename(if not null use the file, if null download it)
#' @return GRanges object with UCSC ref seq genes
#' @note side effect: if download is done file 'refSeq.rds' is created
read_ref_seq <- function(ref_seq_file) {
    if ( ! is.null(ref_seq_file)) {
        .has_ref_seq_file(ref_seq_file)
        return(readRDS(ref_seq_file))
    }
    refSeq_genes <- makeGRanges(
      getUCSCtable("refGene", "RefSeq Genes", freeze=freeze),
      freeze=freeze
    )
    saveRDS(refSeq_genes, file="refSeq.rds")
    refSeq_genes
}

.has_ref_seq_file <- function(ref_seq_file) {
    if ( ! file.exists(ref_seq_file)) {
        stop(paste("COULD NOT FIND REF SEQ FILE: ", ref_seq_file))
    }
}
