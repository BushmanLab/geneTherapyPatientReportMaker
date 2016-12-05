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
