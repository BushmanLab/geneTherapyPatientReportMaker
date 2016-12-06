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

warnlvl <- getOption("warn")
options(warn = -1)

dependancies <- c("methods", "DBI", "RMySQL", "yaml")
null <- sapply(dependancies, function(x){
  suppressPackageStartupMessages(
    try(library(x, character.only = TRUE), silent = TRUE))
})

options(warn = warnlvl)

args <- commandArgs(trailingOnly = TRUE)

pat <- args[1]
config <- yaml.load_file(args[grep("-c", args) + 1]) #Path to INSPIIRED config file

junk <- sapply(dbListConnections(MySQL()), dbDisconnect)

if( config$dataBase == 'mysql' ){
  stopifnot(file.exists("~/.my.cnf"))
  dbConn1 <- dbConnect(MySQL(), group = config$mysqlSpecimenManagementGroup)
#  info1 <- dbGetInfo(dbConn1)
#  dbConn1 <- src_sql("mysql", dbConn1, info = info1)
  dbConn2 <- dbConnect(MySQL(), group = config$mysqlConnectionGroup)
#  info2 <- dbGetInfo(dbConn2)
#  dbConn2 <- src_sql("mysql", dbConn2, info = info2)
}else if( config$dataBase == 'sqlite' ){
  dbConn1 <- dbConnect(RSQLite::SQLite(), dbname = config$sqliteSampleManagement)
#  info1 <- dbGetInfo(dbConn1)
#  dbConn1 <- src_sql("sqlite", dbConn1, info = info1)
  dbConn2 <- dbConnect(RSQLite::SQLite(), dbname = config$sqliteIntSitesDB)
#  info2 <- dbGetInfo(dbConn2)
#  dbConn2 <- src_sql("sqlite", dbConn2, info = info2)
}else{ 
  stop('Can not establish a connection to the database') 
}

sql <- "select trial, patient, cellType, timepoint, parentAlias from hivsp"
specimens <- dbGetQuery(dbConn1, sql)

sql <- "select parentAlias, childAlias, uniqRegion, primerType from hivsam"
samples <- dbGetQuery(dbConn1, sql)

sql <- "select sampleName, refGenome, gender from samples"
intsites <- dbGetQuery(dbConn2, sql)

junk <- sapply(dbListConnections(MySQL()), dbDisconnect)

sp_sam <- merge(specimens, samples, by = "parentAlias")
sp_sam_int <- merge(sp_sam, intsites, by.x = "childAlias", by.y = "sampleName")

sp_sam_int$GTSP <- paste0(sp_sam_int$parentAlias, "-", sp_sam_int$primerType, sp_sam_int$uniqRegion)
sp_sam_int$cellType <- paste0(sp_sam_int$cellType, ":", sp_sam_int$primerType, sp_sam_int$uniqRegion)

if( is.na(pat) | !(pat %in% sp_sam_int$patient) ) {
    write.table(sp_sam_int, "", sep = "\t", row.names=FALSE, quote=FALSE)
    if( is.na(pat) ) q()
    if( !(pat %in% sp_sam_int$patient) ) stop(pat, " patient not found in the above table")
}

df <- sp_sam_int[sp_sam_int$patient == pat, c("childAlias", "GTSP", "patient", "cellType", "timepoint")]
names(df) <- c("sampleName", "GTSP", "patient", "celltype", "timepoint")

write.csv(df, file="", row.names=FALSE, quote=FALSE)
q()

