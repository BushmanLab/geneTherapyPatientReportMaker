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
    uniqueReads <- sapply(uniqueSites, function(sites){length(unique(names(sites)))})
    multihitReads <- sapply(multihitSites, function(multihits){
	reads <- sum(sapply(multihits, function(hits){length(unique(names(hits)))}))
	reads
    })
    readTotals <- data.frame( 
	"sampleName" = sampleName_GTSP$sampleName,
	"refGenome" = refGenome,
	"TotalReads" = uniqueReads + multihitReads
    )
    sample_count_GTSP <- merge(readTotals, sampleName_GTSP)
    aggregate(sample_count_GTSP["TotalReads"], sample_count_GTSP["GTSP"])
}
