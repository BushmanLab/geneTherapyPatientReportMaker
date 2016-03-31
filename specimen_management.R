# access patient metadata for GTSP

#' @param GTSP vector of GTSP 
#' @return df with columns: Trial, GTSP, Patient, Timepoint, CellType, FragMethod, VCN.
#' @note all GTSP in th vector should be unique
get_metadata_for_GTSP <- function(GTSP, db_group) {
    stopifnot(length(GTSP) == length(unique(GTSP)))
    
    HIVSP <- unique(sapply(strsplit(GTSP, split = "-"), "[[", 1))
    query_selection <- "SELECT trial,parentAlias,patient,timepoint,cellType,vcn FROM hivsp"
    string <- paste(HIVSP, collapse = "' OR parentAlias = '")
    query_condition <- paste0("WHERE parentAlias = '", string, "'")  
    query_sets <- paste(query_selection, query_condition)

    query_selection <- "SELECT parentAlias,childAlias,prepMethod FROM hivsam"
    string <- paste(GTSP, collapse = "%' OR childAlias LIKE '")
    query_condition <- paste0("WHERE childAlias LIKE '", string, "%'")
    query_samples <- paste(query_selection, query_condition)
        
    #GTSP = paste0("^", GTSP, "$", collapse="|")
    
    GTSPDBconn <- dbConnect(MySQL(), group=db_group)
    
    #query = paste0("SELECT Trial, SpecimenAccNum, Patient, Timepoint, CellType,
    #                       SamplePrepMethod, VCN
    #               FROM gtsp
    #               WHERE SpecimenAccNum
    #               REGEXP ", dbQuoteString(GTSPDBconn, GTSP), ";")
    
    sets <- dbGetQuery(GTSPDBconn, query_sets)
    samples <- dbGetQuery(GTSPDBconn, query_samples)
    
    dbDisconnect(GTSPDBconn)
    
    sets <- sets[,c("trial", "parentAlias", "patient", "timepoint", "cellType", "vcn")]
    
    samples$GTSP <- mapply(paste0, sapply(strsplit(samples$childAlias, split = "-"), "[[", 1),
                           rep("-", nrow(samples)),
                           sapply(strsplit(samples$childAlias, split = "-"), "[[", 2))
    samples <- unique(samples[,c("parentAlias", "GTSP", "prepMethod")])
    
    sets <- left_join(samples, sets, by = "parentAlias")
    sets <- sets[, c("trial", "GTSP", "patient", "timepoint", "cellType", "prepMethod", "vcn")]
    names(sets) <- c("Trial", "GTSP", "Patient", "Timepoint", "CellType", "FragMethod", "VCN")
    sets    
}
