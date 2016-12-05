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

library("methods", quietly=TRUE)
library("RMySQL", quietly = TRUE) #also loads DBI

SPECIMEN_MANAGEMENT <- "specimen_management_2"

junk <- sapply(dbListConnections(MySQL()), dbDisconnect)
dbConn <- dbConnect(MySQL(), group=SPECIMEN_MANAGEMENT)

sql <- "select trial, patient, CellType, Timepoint, SpecimenAccNum from specimen_management.gtsp"
##message(sql)
trial_pat_gtsp <- dbGetQuery(dbConn,sql)
names(trial_pat_gtsp) <- tolower(names(trial_pat_gtsp))

INTSITES_GROUP <- "gt_intsites.database"

junk <- sapply(dbListConnections(MySQL()), dbDisconnect)
dbConn <- dbConnect(MySQL(), group=INTSITES_GROUP)

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

