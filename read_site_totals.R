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

#' for a given df with 2 cols: sampleName and GTSP
#' gets counts from intSite DB and collapse replicates.
#' @param df with cols: sampleName, GTSP
get_read_site_totals <- function(sampleName_GTSP, connection) {
    reads <- get_count_per_GTSP(sampleName_GTSP, getUniqueSiteReadCounts, "TotalReads", connection)
    reads
    #sites <- get_count_per_GTSP(sampleName_GTSP, getUniqueSiteCounts, "UniqueSites", connection)
    #merge(reads, sites)
}


#' call database_function to get counts then collapse all samples with the same GTSP
#' numbers.
#' @return df with 2 cols: GTSP, column_name
get_count_per_GTSP <- function(sampleName_GTSP, database_function, column_name, connection) {
    # workaround suppressWarnings: type in DB of count int lead to warning when casting to num in R
    sample_count <- suppressWarnings(database_function(sampleName_GTSP, connection))
    if (nrow(sample_count) == 0) {
        message("None of the samples have integration sites.")
        message("Report will NOT be generated.")
        stop(0)
    }
    names(sample_count) <- c("sampleName", "refGenome", column_name)
    sample_count_GTSP <- merge(sample_count, sampleName_GTSP)
    aggregate(sample_count_GTSP[column_name], sample_count_GTSP['GTSP'], FUN=sum)
}

get_data_totals <- function(sampleName_GTSP, uniqueSites, multihitSites){
    uniqueReads <- sapply(uniqueSites, function(sites){sum(sites$counts)})
    multihitReads <- sapply(multihitSites, length)
    readTotals <- data.frame( 
	    "sampleName" = sampleName_GTSP$sampleName,
	    "refGenome" = ref_genome,
	    "TotalReads" = uniqueReads + multihitReads)
    
    sample_count_GTSP <- merge(readTotals, sampleName_GTSP)
    aggregate(sample_count_GTSP["TotalReads"], sample_count_GTSP["GTSP"], FUN=sum)
}
