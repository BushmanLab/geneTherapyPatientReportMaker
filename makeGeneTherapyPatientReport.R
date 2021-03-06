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

library(DBI, quietly=TRUE, verbose=FALSE)
library(yaml, quietly=TRUE, verbose=FALSE)
options(stringsAsFactors = FALSE, useFancyQuotes=FALSE)

#' set all argumentgs for the script
#' @return list of argumentgs
#' @example set_args()
#'          set_args(c("--ref_genome", "mm9"))
#' Rscript ~/geneTherapyPatientReportMaker/makeGeneTherapyPatientReport.R --ref_genome hg19 hs.csv
#' Rscript ~/geneTherapyPatientReportMaker/makeGeneTherapyPatientReport.R --ref_genome mm9 mm.csv
set_args <- function(...) {
    ## arguments from command line
    suppressMessages(library(argparse))
    parser <- ArgumentParser(description="Gene Therapy Patient Report for Single Patient")
    parser$add_argument("sample_gtsp", nargs='?', default='sampleName_GTSP.csv')
    parser$add_argument("-c", default="./INSPIIRED.yml", help="path to INSPIIRED configuration file.")
    parser$add_argument("-s", action='store_true', help="abundance by sonicLength package (Berry, C. 2012)")
    parser$add_argument("-r", "--ref_genome", default="hg18", help="reference genome used for all samples")
#    parser$add_argument("--sites_group", default="intsites_miseq.read", help="group to use for integration sites db from ~/.my.cnf")
#    parser$add_argument("--gtsp_group", default="specimen_management", help="group to use for specimen management GTSP db from ~/.my.cnf")
    parser$add_argument("--ref_seq", help="read Ref Seq genes from file")
    parser$add_argument("-o", "--output", help='HTML and MD file names instead of Trial.Patient.Date name')

    arguments <- parser$parse_args(...)

    ## gene files
    arguments$oncoGeneFile <- ""
    if(grepl("^hg", arguments$ref_genome)) arguments$oncoGeneFile <- "allonco_no_pipes.csv"
    if(grepl("^mm", arguments$ref_genome)) arguments$oncoGeneFile <- "allonco_no_pipes.mm.csv"
    stopifnot( arguments$ref_genome!="" )
    
    ## note this file was obtained by
    ## wget http://www.bushmanlab.org/assets/doc/humanLymph.tsv
    ## it contains 38 gene names associated with Lymphoma
    arguments$adverseGeneFile <- "humanLymph.tsv"
    
    ## code dir past to Rscript
    codeDir <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly=FALSE), value=T)))
    if( length(codeDir)!=1 ) codeDir <- list.files(path="~", pattern="geneTherapyPatientReportMaker$", recursive=TRUE, include.dirs=TRUE, full.names=TRUE)
    stopifnot(file.exists(file.path(codeDir, "GTSPreport.css")))
    stopifnot(file.exists(file.path(codeDir, "GTSPreport.Rmd")))
    stopifnot(file.exists(file.path(codeDir, arguments$oncoGeneFile)))
    arguments$codeDir <- codeDir
    
    return(arguments)
}
arguments <- set_args()
print(arguments)

arguments$gtsp_group
# Load configuration file
if (!file.exists(arguments$c)) stop("the configuration file can not be found.")
config <<- yaml.load_file(arguments$c)


## defaults:
use.sonicLength <-  ! arguments$s
#db_group_sites <- arguments$sites_group
#db_group_gtsp <- arguments$gtsp_group
ref_genome <- arguments$ref_genome
codeDir <- arguments$codeDir
ref_seq_filename <- arguments$ref_seq

#### INPUTS: csv file/table GTSP to sampleName ####
csvfile <- arguments$sample_gtsp

if( !file.exists(csvfile) ) stop(csvfile, "not found")

#### load up require packages + objects #### 
libs <- c("RMySQL", "plyr", "dplyr", "stringr", "reshape2",
          "scales", "ggplot2", "devtools", "reldist",
          "hiAnnotator", "sonicLength", "intSiteRetriever",
          "BiocParallel", "PubMedWordcloud", "markdown",
          "RColorBrewer", "magrittr", "knitr")
null <- suppressMessages(sapply(libs, library, character.only=TRUE))

R_source_files <- c("utilities.R", "estimatedAbundance.R", "read_site_totals.R", "ref_seq.R",
                    "populationInfo.R", "abundanceFilteringUtils.R")

null <- sapply(R_source_files, function(x) source(file.path(codeDir, x)))

url <- "https://raw.githubusercontent.com/BushmanLab/intSiteCaller/master/"
source_url(paste0(url, "hiReadsProcessor.R"))
source_url(paste0(url, "standardization_based_on_clustering.R"))

#### load datasets and process them before knit #### 
message("\nReading csv from ", csvfile)
sampleName_GTSP <- read.csv(csvfile)
stopifnot(all(c("sampleName", "GTSP") %in% colnames(sampleName_GTSP)))
sampleName_GTSP$refGenome <- ref_genome
message("\nGenerating report from the following sets")
print(sampleName_GTSP)

# Connect to my database
if (config$dataBase == 'mysql'){
   stopifnot(file.exists("~/.my.cnf"))
   dbConn <- dbConnect(MySQL(), group=config$mysqlConnectionGroup)
   info <- dbGetInfo(dbConn)
   dbConn <- src_sql("mysql", dbConn, info = info)
   dbConnSampleManagemnt <-  dbConnect(MySQL(), group=config$mysqlSpecimenManagementGroup)
}else if (config$dataBase == 'sqlite') {
   dbConn <- dbConnect(RSQLite::SQLite(), dbname=config$sqliteIntSitesDB)
   info <- dbGetInfo(dbConn)
   dbConn <- src_sql("sqlite", dbConn, info = info)
   dbConnSampleManagemnt <- dbConnect(RSQLite::SQLite(), dbname=config$sqliteSampleManagement)
} else { stop('Can not establish a connection to the database') }

if( !all(setNameExists(sampleName_GTSP, dbConn)) ) {
    sampleNameIn <- paste(sprintf("'%s'", sampleName_GTSP$sampleName), collapse=",")

    q <- sprintf("SELECT * FROM samples WHERE sampleName IN (%s)", sampleNameIn)
    message("\nChecking database:\n",q,"\n")

    t <- dbSendQuery(dbConn, q)
    write.table(t)

    stop("Was --ref_genome specified correctly or did query return all entries")
   } else {
    message("All samples are in DB.")
}

read_sites_sample_GTSP <- get_read_site_totals(sampleName_GTSP, dbConn)

get_metadata_for_GTSP <- function(GTSP, GTSPDBconn) {
    stopifnot(length(GTSP) == length(unique(GTSP)))

    GTSP = paste(sQuote(GTSP), collapse=',')

    query = paste0("SELECT Trial, SpecimenAccNum, Patient, Timepoint, CellType, SamplePrepMethod, VCN
                   FROM gtsp
                   WHERE SpecimenAccNum in (", GTSP, ");");
    
    sets <- dbGetQuery(GTSPDBconn, query)

    names(sets) <- c("Trial", "GTSP", "Patient", "Timepoint", "CellType", "FragMethod", "VCN")
    sets
}

sets <- get_metadata_for_GTSP(unique(sampleName_GTSP$GTSP), dbConnSampleManagemnt)

## some clean up for typos, dates, spaces etc
sets[sets$Timepoint=="NULL", "Timepoint"] <- "d0"
sets[sets$Timepoint=="", "Timepoint"] <- "d0"
sets$Timepoint <- gsub('_', '.', sets$Timepoint, fixed=TRUE)
for(col in which(!sapply(sets, class) %in% c("numeric", "integer"))) {
    sets[[col]] <-  gsub("\\s", '', sets[[col]])
}

# reports are for a single patient
stopifnot(length(unique(sets$Patient)) == 1)
patient <- sets$Patient[1]
# and for a single trial
stopifnot(length(unique(sets$Trial)) == 1)
trial <- sets$Trial[1]


RDataFile <- paste(trial, patient, format(Sys.Date(), format="%Y%m%d"), "RData", sep=".")

# all GTSP in the database
stopifnot(nrow(sets) == length(unique(sampleName_GTSP$GTSP)))


##end INPUTS

sets <- merge(sets, read_sites_sample_GTSP)
sets$Timepoint <- sortFactorTimepoints(sets$Timepoint)

# at present the whole report is done for one genome
stopifnot(length(unique(sampleName_GTSP$refGenome))==1)
freeze <- sampleName_GTSP[1, "refGenome"]

##==========GET AND PERFORM BASIC DEREPLICATION/SONICABUND ON SITES=============
message("Fetching unique sites and estimating abundance")


sites <- merge(getUniquePCRbreaks(sampleName_GTSP, dbConn), sampleName_GTSP)
names(sites)[names(sites)=="position"] <- "integration"

#we really don't care about seqinfo - we just want a GRange object for easy manipulation
uniqueSites.gr <- GRanges(seqnames=Rle(sites$chr),
                          ranges=IRanges(start=pmin(sites$integration, sites$breakpoint),
                                         end=pmax(sites$integration, sites$breakpoint)),
                          strand=Rle(sites$strand))
mcols(uniqueSites.gr) <- sites[,c("sampleName", "GTSP")]

#standardize sites across all GTSPs
isthere <- which("dplyr" == loadedNamespaces()) # temp work around  of
#Known conflict with package:dplyr::count(), need to unload package if present
# unloading and reloading the package
if(length(isthere) > 0){detach("package:dplyr", unload = TRUE)}
standardizedReplicatedSites <- standardizeSites(uniqueSites.gr)
if(length(isthere) > 0){suppressMessages(library("dplyr"))}

standardizedReplicatedSites$posid <- paste0(seqnames(standardizedReplicatedSites),
                                            strand(standardizedReplicatedSites),
                                            start(flank(standardizedReplicatedSites, -1, start=T)))
standardizedReplicatedSites <- split(standardizedReplicatedSites,
                                     standardizedReplicatedSites$GTSP)
standardizedReplicatedSites <- lapply(standardizedReplicatedSites, function(x){
    x$replicate <- as.integer(as.factor(x$sampleName))
    x$sampleName <- NULL
    x
})

#this is slow (~1.5min/sample), but would be easy to parallelize - just be
#careful about memory consumption!  sonic abundance could get 20GB+ per thread
standardizedReplicatedSites <- standardizedReplicatedSites[sapply(standardizedReplicatedSites, length)>0]

standardizedDereplicatedSites <- lapply(standardizedReplicatedSites, function(sites){
    res <- getEstimatedAbundance(sites, use.sonicLength = use.sonicLength)
    res$GTSP <- sites[1]$GTSP
    res$posid <- paste0(seqnames(res), strand(res), start(flank(res, -1, start=T)))
    res
})

standardizedReplicatedSites <- prepSiteList(standardizedReplicatedSites)
standardizedDereplicatedSites <- prepSiteList(standardizedDereplicatedSites)
standardizedDereplicatedSites <- flank(standardizedDereplicatedSites, -1, start=TRUE)

unique_sites_per_GTSP <- sapply(split(standardizedDereplicatedSites,
                                      standardizedDereplicatedSites$GTSP),
                                function(x){length(unique(x$posid))})
unique_sites_per_GTSP <- data.frame("GTSP" = names(unique_sites_per_GTSP),
                                    "UniqueSites" = unique_sites_per_GTSP)
sets <- merge(sets, unique_sites_per_GTSP, by = "GTSP")

cells_recovered <- (standardizedReplicatedSites %>%
                              as.data.frame %>%
                              select(GTSP, replicate, posid, width) %>%
                              distinct %>%
                              group_by(GTSP) %>%
                              count(GTSP))

unique_cells_per_GTSP <- data.frame("GTSP" = cells_recovered$GTSP,
                                      "InferredCells" = cells_recovered$n)                             

sets <- merge(sets, unique_cells_per_GTSP, by = "GTSP")
#============CALCULATE POPULATION SIZE/DIVERSITY INFORMATION=================
populationInfo <- getPopulationInfo(standardizedReplicatedSites,
                                    standardizedDereplicatedSites,
                                    "GTSP")
populationInfo$Replicates <- sapply(split(standardizedReplicatedSites$replicate,
                                          standardizedReplicatedSites$GTSP),
                                    max)

#========CALCULATE POPULATION SIZE/DIVERSITY INFORMATION BY TIMEPOINT==========
timepointPopulationInfo <- getPopulationInfo(standardizedReplicatedSites,
                                             standardizedDereplicatedSites,
                                             "Timepoint")

timepointPopulationInfo$UniqueSites <- sapply(split(standardizedDereplicatedSites, 
                                                    standardizedDereplicatedSites$Timepoint),
                                              function(x){length(unique(x$posid))})


#=======================ANNOTATE DEREPLICATED SITES==========================
#standard refSeq genes
message("Annotating unique hit sites")
refSeq_genes <- read_ref_seq(ref_seq_filename)

save.image(RDataFile)

standardizedDereplicatedSites <- getNearestFeature(standardizedDereplicatedSites,
                                                   refSeq_genes,
                                                   colnam="nearest_refSeq_gene",
                                                   feature.colnam="name2")

standardizedDereplicatedSites <- getNearestFeature(standardizedDereplicatedSites,
                                                   refSeq_genes,
                                                   colnam="nearest_refSeq_gene",
                                                   side="5p",
                                                   feature.colnam="name2")

standardizedDereplicatedSites <- getSitesInFeature(standardizedDereplicatedSites,
                                                   refSeq_genes,
                                                   colnam="inGene",
                                                   feature.colnam="name2")

#oncogenes
oncogenes <- scan(file= file.path(codeDir, arguments$oncoGeneFile), what='character')
oncogenes <- oncogenes[!grepl("geneName", oncogenes, ignore.case=TRUE)]

refSeq_oncogene <- refSeq_genes[toupper(refSeq_genes$name2) %in% toupper(oncogenes)]

standardizedDereplicatedSites <- getNearestFeature(standardizedDereplicatedSites,
                                                   refSeq_oncogene,
                                                   colnam="NrstOnco",
                                                   side="5p",
                                                   feature.colnam="name2")


wantedgenes <- as.character(
    read.csv(file.path(codeDir, arguments$adverseGeneFile),
             sep="\t",
             header=TRUE)$symbol)
stopifnot(length(wantedgenes)>=1)

## * in transcription units
## ~ within 50kb of a onco gene 
## ! nearest is a known bad gene 
standardizedDereplicatedSites$geneMark <- ""

## ~ nearest is a known bad gene 
isNearWanted <- standardizedDereplicatedSites$nearest_refSeq_gene %in% wantedgenes 
isInWanted <- sapply( standardizedDereplicatedSites$inGene,
                     function(txt) any(unlist(strsplit(txt, ',')) %in% wantedgenes) )
standardizedDereplicatedSites$geneMark <- ifelse(
    isNearWanted | isInWanted,
    paste0(standardizedDereplicatedSites$geneMark, "!"),
    standardizedDereplicatedSites$geneMark )

## ~ within 50kb of a onco gene 
isNearOnco <- (standardizedDereplicatedSites$nearest_refSeq_gene %in% oncogenes &
               abs(standardizedDereplicatedSites$nearest_refSeq_geneDist) < 50000 )
isInOnco <- sapply( standardizedDereplicatedSites$inGene,
                   function(txt) any(unlist(strsplit(txt, ',')) %in% oncogenes) )
standardizedDereplicatedSites$geneMark <- ifelse(
    isNearOnco | isInOnco,
    paste0(standardizedDereplicatedSites$geneMark, "~"),
    standardizedDereplicatedSites$geneMark )

## * in transcription units
standardizedDereplicatedSites$geneMark <- ifelse(
    toupper(standardizedDereplicatedSites$inGene)!="FALSE",
    paste0(standardizedDereplicatedSites$geneMark, "*"),
    standardizedDereplicatedSites$geneMark)

## attach gene marks
standardizedDereplicatedSites$nearest_refSeq_gene <- paste0(
    standardizedDereplicatedSites$nearest_refSeq_gene,
    standardizedDereplicatedSites$geneMark)

#===================GENERATE EXPANDED CLONE DATAFRAMES======================
standardizedDereplicatedSamples <- split(
  standardizedDereplicatedSites, 
  standardizedDereplicatedSites$CellType
)

cutoff_genes_barplot <- lapply(
  standardizedDereplicatedSamples,
  getMostAbundantGenes, 
  numGenes = 10
)

abundCutoff.barplots <- sapply(cutoff_genes_barplot, "[[", 1)
frequent_genes_barplot_by_sample <- lapply(cutoff_genes_barplot, "[[", 2)
frequent_genes_barplot <- unique(unlist(frequent_genes_barplot_by_sample))

barplotAbunds <- lapply(1:length(standardizedDereplicatedSamples), function(i){
  sites <- standardizedDereplicatedSamples[[i]]
  genes <- frequent_genes_barplot
  getAbundanceSums(maskGenes(sites, genes), c("CellType", "Timepoint"))
})
barplotAbundsBySample <- lapply(1:length(standardizedDereplicatedSamples), function(i){
  sites <- standardizedDereplicatedSamples[[i]]
  genes <- frequent_genes_barplot_by_sample[[i]]
  getAbundanceSums(maskGenes(sites, genes), c("CellType", "Timepoint"))
})

barplotAbunds <- bind_rows(lapply(barplotAbunds, order_barplot))
barplotAbundsBySample <- bind_rows(lapply(barplotAbundsBySample, order_barplot))
CellType_order <- unique(standardizedDereplicatedSites$CellType)
barplotAbunds$CellType <- factor(barplotAbunds$CellType, levels=CellType_order)
barplotAbundsBySample$CellType <- factor(barplotAbundsBySample$CellType, levels=CellType_order)

#detailed abundance plot
cutoff_genes <- getMostAbundantGenes(standardizedDereplicatedSites, 50)
abundCutoff.detailed <- cutoff_genes[[1]]
frequent_genes <- cutoff_genes[[2]]

detailedAbunds <- getAbundanceSums(maskGenes(
    standardizedDereplicatedSites, frequent_genes), c("CellType", "Timepoint"))

categorySums <- sapply(split(detailedAbunds$estAbundProp,
                             detailedAbunds$maskedRefGeneName),sum)

detailedAbunds$maskedRefGeneName <- factor(detailedAbunds$maskedRefGeneName,
                                           levels=names(sort(categorySums)))
detailedAbunds$CellType <- factor(detailedAbunds$CellType,
                                  levels=CellType_order)

#================Longitudinal Behaviour===============================
longitudinal <- as.data.frame(standardizedDereplicatedSites)[,c("Timepoint",
                                                                "CellType",
                                                                "estAbundProp",
                                                                "posid")]
longitudinal$CellType <- factor(longitudinal$CellType,
                                levels=CellType_order)

has_longitudinal_data <- length(unique(longitudinal$Timepoint)) > 1


#==================DETAILED REPORTS FOR BAD ACTORS=====================
badActors <- c("LMO2", "IKZF1", "CCND2", "HMGA2", "MECOM")

badActorData <- sapply(badActors, function(badActor){
  hasBadActor <- grepl(badActor, standardizedDereplicatedSites$X5pNrstOnco)
  badActorWithin100K <- abs(standardizedDereplicatedSites$X5pNrstOncoDist) <= 100000
  standardizedDereplicatedSites[hasBadActor & badActorWithin100K]
})

badActorData <- lapply(badActorData, function(x){
  x$CellType <- factor(x$CellType, levels=CellType_order)
  x
})

#==================SET VARIABLES FOR MARKDOWN REPORT=====================
timepoint <- levels(sets$Timepoint)

popSummaryTable <- merge(sets,  populationInfo, by.x="GTSP", by.y="group")
popSummaryTable <- arrange(popSummaryTable,Timepoint,CellType)

cols <- c("Trial", "GTSP", "Replicates", "Patient", "Timepoint", "CellType", 
          "TotalReads", "InferredCells", "UniqueSites", "FragMethod", "VCN", "S.chao1", "Gini", "Shannon", "UC50")
summaryTable <- popSummaryTable[,cols]

summaryTable$VCN <- ifelse(summaryTable$VCN == 0, NA, summaryTable$VCN)
    
timepointPopulationInfo <- melt(timepointPopulationInfo, "group")

#==================Get abundance for multihit events=====================
message("Fetching multihit sites and estimating abundance")

options(useFancyQuotes=FALSE)
sql <- paste0("select samples.sampleName, samples.refGenome, multihitpositions.multihitID, ",
              "multihitlengths.length from multihitlengths left join multihitpositions on ",
              "multihitpositions.multihitID = multihitlengths.multihitID left join samples on ",
              "samples.sampleID = multihitpositions.sampleID where sampleName in (",
               paste0(sQuote(unique(sampleName_GTSP$sampleName)), collapse=","),  ")")

# Create a new DB connection.
# This connection is not a dplyr src_sql() connection like the previous dbConns.


# Connect to my database
if (config$dataBase == 'mysql'){
   stopifnot(file.exists("~/.my.cnf"))
   dbConn <- dbConnect(MySQL(), group=config$mysqlConnectionGroup)
}else if (config$dataBase == 'sqlite') {
   dbConn <- dbConnect(RSQLite::SQLite(), dbname=config$sqliteIntSitesDB)
} else { stop('Can not establish a connection to the database') }


multiHitLengths <- unique(dbGetQuery(dbConn, sql))
sites.multi <- merge(multiHitLengths, sampleName_GTSP)

if( nrow(sites.multi) > 0 ) {
    sites.multi <- (sites.multi %>%
                    group_by(multihitID) %>%
                    mutate(replicate=as.integer(as.factor(sampleName))) )
    
    sites.multi <- (sites.multi %>%
                    group_by(sampleName, multihitID) %>%
                    mutate( estAbund = length(unique(length)) ) )
    
    sites.multi <- merge(sites.multi, sets, by="GTSP")
    
    sites.multi <- (sites.multi %>%
                    group_by(Patient, Timepoint, CellType) %>%
                    mutate(Rank=rank(-estAbund, ties.method="max")))
}


save.image(RDataFile)
##end setting variables for markdown report

fig.path <- paste(unique(trial), unique(patient),
                  format(Sys.Date(), format="%Y%m%d"), "Figures",
                  sep=".")

#### begin generating markdown ####
unlink(fig.path, force=TRUE, recursive=TRUE)
mdfile <- paste(unique(trial), unique(patient),
                format(Sys.Date(), format="%Y%m%d"), "md",
                sep=".")

if ( ! is.null(arguments$output)) {
    mdfile <- paste(arguments$output, 'md', sep='.')
}

htmlfile <- gsub("\\.md$",".html",mdfile)
options(knitr.table.format='html')
theme_set(theme_bw()) #for ggplot2
knit(file.path(codeDir, "GTSPreport.Rmd"), output=mdfile)
markdownToHTML(mdfile, htmlfile, extensions=c('tables'),
               options=c(markdownHTMLOptions(defaults=T), "toc"),
               stylesheet=file.path(codeDir, "GTSPreport.css") )


message("\nReport ", htmlfile, " is generated from ", csvfile)
save.image(RDataFile)

