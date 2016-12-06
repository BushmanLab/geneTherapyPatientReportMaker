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

group <- "intsites_miseq.read"

junk <- sapply(dbListConnections(MySQL()), dbDisconnect)
dbConn <- dbConnect(MySQL(), group=group)

sql <- "select trial, patient, CellType, Timepoint, SpecimenAccNum from specimen_management.gtsp"
##message(sql)
trial_pat_gtsp <- dbGetQuery(dbConn,sql)
names(trial_pat_gtsp) <- tolower(names(trial_pat_gtsp))

workDir <- commandArgs(trailingOnly=TRUE)[1]
if( is.na(workDir) ) workDir <- "."
workDir <- normalizePath(workDir, mustWork=TRUE)
gtspDir <- grep("^GTSP", list.dirs(path=workDir, full.names=FALSE), ignore.case=TRUE, value=TRUE)

gtsp <- sub("-\\d+$", "", gtspDir)

df <- data.frame(GTSP=gtsp, sampleName=gtspDir)

merged.tab <- merge(df, trial_pat_gtsp,
                    by.x="GTSP", by.y="specimenaccnum")
write.table(merged.tab, "", sep = "\t", row.names=FALSE, quote=FALSE)

q()

