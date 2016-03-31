library("methods", quietly=TRUE)
library("RMySQL", quietly = TRUE) #also loads DBI

group <- "hiv_specimen.database"

junk <- sapply(dbListConnections(MySQL()), dbDisconnect)
dbConn <- dbConnect(MySQL(), group=group)

sql <- "select trial, patient, cellType, timepoint, parentAlias from hivsp"
##message(sql)
trial_pat_gtsp <- dbGetQuery(dbConn,sql)
names(trial_pat_alias) <- tolower(names(trial_pat_alias))

sql <- "select sampleName, refGenome, gender from samples"
##message(sql)
set_ref_gen <- dbGetQuery(dbConn,sql)
names(set_ref_gen) <- tolower(names(set_ref_gen))
set_ref_gen$gtsp <- sub("-\\d+$", "", set_ref_gen$samplename)

merged.tab <- merge(trial_pat_gtsp, set_ref_gen,
                    by.x="specimenaccnum", by.y="gtsp")

merged.tab <- plyr:::arrange(merged.tab, trial, patient, specimenaccnum, samplename, refgenome, gender)

merged.tab <- subset(merged.tab, select=c("trial", "patient", "celltype", "timepoint", "specimenaccnum", "samplename", "refgenome", "gender"))

##message()

pat <- commandArgs(trailingOnly=TRUE)[1]
if( is.na(pat) | !(pat %in% merged.tab$patient) ) {
    write.table(merged.tab, "", sep = "\t", row.names=FALSE, quote=FALSE)
    if( is.na(pat) ) q()
    if( !(pat %in% merged.tab$patient) ) stop(pat, " patient not found in the above table")
}

df <- subset(merged.tab, patient==pat, select=c("samplename", "specimenaccnum", "patient", "celltype", "timepoint"))
colnames(df) <- c("sampleName", "GTSP", "patient", "celltype", "timepoint")

write.csv(df, file="", row.names=FALSE, quote=FALSE)
q()

