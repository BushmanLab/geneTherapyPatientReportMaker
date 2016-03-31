library("methods", quietly=TRUE)
library("RMySQL", quietly = TRUE) #also loads DBI

specimen_group <- "hiv_specimen.database"
intsites_group <- "hiv_intsites.database"

junk <- sapply(dbListConnections(MySQL()), dbDisconnect)
dbConn1 <- dbConnect(MySQL(), group = specimen_group)
dbConn2 <- dbConnect(MySQL(), group = intsites_group)

sql <- "select trial, patient, cellType, timepoint, parentAlias from hivsp"
specimens <- dbGetQuery(dbConn1, sql)
#names(specimens) <- tolower(names(specimens))

sql <- "select parentAlias, childAlias, uniqRegion, primerType from hivsam"
samples <- dbGetQuery(dbConn1, sql)
#names(samples) <- tolower(names(samples))

sql <- "select sampleName, refGenome, gender from samples"
intsites <- dbGetQuery(dbConn2,sql)
#names(intsites) <- tolower(names(intsites))
#set_ref_gen$gtsp <- sub("-\\d+$", "", set_ref_gen$samplename)

junk <- sapply(dbListConnections(MySQL()), dbDisconnect)

sp_sam <- merge(specimens, samples, by = "parentAlias")
sp_sam_int <- merge(sp_sam, intsites, by.x = "childAlias", by.y = "sampleName")

sp_sam_int$GTSP <- paste0(sp_sam_int$parentAlias, "-", sp_sam_int$primerType, sp_sam_int$uniqRegion)
sp_sam_int$cellType <- paste0(sp_sam_int$cellType, ":", sp_sam_int$primerType, sp_sam_int$uniqRegion)

#merged.tab <- merge(trial_pat_gtsp, set_ref_gen,
#                    by.x="specimenaccnum", by.y="gtsp")

#merged.tab <- plyr:::arrange(merged.tab, trial, patient, specimenaccnum, samplename, refgenome, gender)

#merged.tab <- subset(merged.tab, select=c("trial", "patient", "celltype", "timepoint", "specimenaccnum", "samplename", "refgenome", "gender"))

##message()

pat <- commandArgs(trailingOnly=TRUE)[1]
if( is.na(pat) | !(pat %in% sp_sam_int$patient) ) {
    write.table(sp_sam_int, "", sep = "\t", row.names=FALSE, quote=FALSE)
    if( is.na(pat) ) q()
    if( !(pat %in% sp_sam_int$patient) ) stop(pat, " patient not found in the above table")
}

df <- sp_sam_int[sp_sam_int$patient == pat, c("childAlias", "GTSP", "patient", "cellType", "timepoint")]
names(df) <- c("sampleName", "GTSP", "patient", "celltype", "timepoint")

write.csv(df, file="", row.names=FALSE, quote=FALSE)
q()

