{
# Set this flag to FALSE if you are developing with an IDE such as Rstudio.
# Command line arguments will be defined in the code below.
commandLine = FALSE




# Set R parameters and loaded required libraries.
options(stringsAsFactors = FALSE, useFancyQuotes=FALSE)

loadLibraries <- function(libs){
  cat('Loading libraries... ')
  null <- invisible(lapply(libs, function(x) suppressPackageStartupMessages(library(x, character.only = TRUE))))
  cat('done.\n')
}

libs <- c( "DBI", "yaml", "RMySQL", "plyr", "dplyr", "stringr", "reshape2","scales", "ggplot2",
           "devtools", "reldist", "hiAnnotator", "sonicLength", "intSiteRetriever", "BiocParallel",
           "PubMedWordcloud", "markdown", "RColorBrewer", "magrittr", "knitr")

null <- suppressMessages(sapply(libs, library, character.only=TRUE))





# Parse command line arguments OR manually supply required parameters.
if (commandLine){
   set_args <- function(...) {
      library(argparse)
      parser <- ArgumentParser(description="Gene Therapy Patient Report for Single Patient")
      parser$add_argument("sample_gtsp", nargs='?', default='sampleName_GTSP.csv')
      parser$add_argument("-c", default="./INSPIIRED.yml", help="path to INSPIIRED configuration file.")
      parser$add_argument("-s", action='store_true', help="abundance by sonicLength package (Berry, C. 2012)")
      parser$add_argument("-r", "--ref_genome", default="hg18", help="reference genome used for all samples")
      parser$add_argument("--ref_seq", help="read Ref Seq genes from file")
      parser$add_argument("-o", "--output", help='HTML and MD file names instead of Trial.Patient.Date name')

      arguments <- parser$parse_args(...)

      codeDir <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly=FALSE), value=T)))
      if( length(codeDir)!=1 ) codeDir <- list.files(path="~", pattern="geneTherapyPatientReportMaker$", recursive=TRUE, include.dirs=TRUE, full.names=TRUE)
      stopifnot(file.exists(file.path(codeDir, "GTSPreport.css")))
      stopifnot(file.exists(file.path(codeDir, "GTSPreport.Rmd")))
      stopifnot(file.exists(file.path(codeDir, arguments$oncoGeneFile)))
      arguments$codeDir <- codeDir
    
      return(arguments)
   }

   arguments <- set_args()
   
   if (!file.exists(arguments$c)) stop("The configuration file can not be found.")
   config <- yaml.load_file(arguments$c)
   
} else {
  
  arguments <- list()
  arguments$ref_genome      <- 'hg18'
  arguments$s               <- FALSE
  arguments$codeDir         <- '/home/everett/projects/geneTherapyPatientReportMaker'
  arguments$c               <- paste0(arguments$codeDir, '/INSPIIRED.yml')
  arguments$sample_gtsp     <- paste0(arguments$codeDir, '/sampleName_GTSP.csv')
  arguments$ref_seq         <- paste0(arguments$codeDir, '/refSeq.rds')
  
  config <- list()
  config$dataBase <-  'mysql'
  config$mysqlConnectionGroup <- 'intsites_miseq'
  config$mysqlSpecimenManagementGroup <- 'specimen_management'
}



# Prepare oncogene lists.
# The human onco gene list contains gene symbols, previous symbols and synonyms.
# All three fields will be split on commas and combined into a single character vector.
# The mouse gene symbols are a single list of ids that can be scanned into a vector.

if(grepl("^hg", arguments$ref_genome)){
  t <- suppressWarnings(read.table(paste0(arguments$codeDir, '/', 'allOnco.human.tsv'), 
                                   header=TRUE, 
                                   sep='\t', 
                                   as.is=T, 
                                   strip.white=T))
  
  oncoGenes <- toupper(unlist(apply(t, 1, function(x){
    synonyms    <- gsub('\\s', '', unname(unlist(strsplit( x['synonyms'],    ',', fixed=TRUE))))
    prevSymbols <- gsub('\\s', '', unname(unlist(strsplit( x['prevSymbols'], ',', fixed=TRUE))))
    c( unname(x['symbol']), prevSymbols, synonyms )
  })))
}else if(grepl("^mm", arguments$ref_genome)){
  oncoGenes <- toupper(scan(paste0(arguments$codeDir, '/', 'allOnco.mouse.list'), what='character', sep='\n', strip.white=TRUE))
} else {
  stop("The software only supports human and mouse ref genomes.")
}
  
humanLymphGenes <- toupper(scan(paste0(arguments$codeDir, '/', 'humanLymph.list'), what='character', sep='\n', strip.white=TRUE))




# Read in and check sample file.
if( !file.exists(arguments$sample_gtsp) ) stop(paste0(arguments$sample_gtsp, " not found"))
message("\nReading csv from ", arguments$sample_gtsp)
sampleName_GTSP <- read.csv(arguments$sample_gtsp)
stopifnot(all(c("sampleName", "GTSP") %in% colnames(sampleName_GTSP)))
sampleName_GTSP$refGenome <- arguments$ref_genome
message("\nGenerating report from the following sets")
print(sampleName_GTSP)


# Throw an error if the expexted columns are not in the provided CSV file.
# Throw an error if more than one GTSP has the same patient / cell type / timepoint combination.

stopifnot(all(c('sampleName', 'GTSP', 'patient', 'celltype', 'timepoint') %in% colnames(sampleName_GTSP)))
stopifnot(length(unique(sampleName_GTSP$GTSP)) == nrow(unique(sampleName_GTSP[, c('patient', 'celltype', 'timepoint')])))



# Source required files.
R_source_files <- c("utilities.R", "estimatedAbundance.R", "read_site_totals.R", "ref_seq.R",
                    "populationInfo.R", "abundanceFilteringUtils.R")

null <- sapply(R_source_files, function(x) source(file.path(arguments$codeDir, x)))

url <- "https://raw.githubusercontent.com/BushmanLab/intSiteCaller/master/"
source_url(paste0(url, "hiReadsProcessor.R"))
source_url(paste0(url, "standardization_based_on_clustering.R"))




# Connect to my database
if (config$dataBase == 'mysql'){
   stopifnot(file.exists("~/.my.cnf"))
   sapply(dbListConnections(MySQL()), dbDisconnect)
   dbConn0 <- dbConnect(MySQL(), group=config$mysqlConnectionGroup)
   info <- dbGetInfo(dbConn0)
   dbConn <- src_sql("mysql", dbConn0, info = info)
   dbConnSampleManagemnt <-  dbConnect(MySQL(), group=config$mysqlSpecimenManagementGroup)
}else if (config$dataBase == 'sqlite') {
   dbConn0 <- dbConnect(RSQLite::SQLite(), dbname=config$sqliteIntSitesDB)
   info <- dbGetInfo(dbConn0)
   dbConn <- src_sql("sqlite", dbConn0, info = info)
   dbConnSampleManagemnt <- dbConnect(RSQLite::SQLite(), dbname=config$sqliteSampleManagement)
} else { stop('Can not establish a connection to the database') }




# Check samples against the datbase.
if( !all(setNameExists(sampleName_GTSP, dbConn)) ) {
    sampleNameIn <- paste(sprintf("'%s'", sampleName_GTSP$sampleName), collapse=",")

    q <- sprintf("SELECT * FROM samples WHERE sampleName IN (%s)", sampleNameIn)
    message("\nChecking database:\n",q,"\n")

    t <- dbSendQuery0(con, q)
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
    res <- getEstimatedAbundance(sites, use.sonicLength = arguments$s)
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
refSeq_genes <- read_ref_seq(arguments$ref_seq)

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


refSeq_oncogene <- refSeq_genes[toupper(refSeq_genes$name2) %in% toupper(oncoGenes)]

standardizedDereplicatedSites <- getNearestFeature(standardizedDereplicatedSites,
                                                   refSeq_oncogene,
                                                   colnam="NrstOnco",
                                                   side="5p",
                                                   feature.colnam="name2")

## * in transcription units
## ~ within 50kb of a onco gene 
## ! nearest is a known bad gene 
standardizedDereplicatedSites$geneMark <- ""

## ! nearest is a known bad gene 

isNearWanted <- unlist(lapply(standardizedDereplicatedSites$nearest_refSeq_gene,
                              function(txt){ any(toupper(unlist(strsplit(txt, ','))) %in% toupper(humanLymphGenes)) }))

isInWanted <- unlist(lapply(standardizedDereplicatedSites$inGene,
                            function(txt){ any(toupper(unlist(strsplit(txt, ','))) %in% toupper(humanLymphGenes)) }))

standardizedDereplicatedSites$geneMark <- ifelse(
    isNearWanted | isInWanted,
    paste0(standardizedDereplicatedSites$geneMark, "!"),
    standardizedDereplicatedSites$geneMark )





# ~ within 50kb of a onco gene 
isNearOnco <- unname(unlist(mapply(function(gene, dist){
  any(toupper(unlist(strsplit(gene, ','))) %in% toupper(oncoGenes)) & abs(dist) < 50000
}, standardizedDereplicatedSites$nearest_refSeq_gene, standardizedDereplicatedSites$nearest_refSeq_geneDist)))

isInOnco <- unlist(lapply(standardizedDereplicatedSites$inGene,
                          function(txt) any(toupper(unlist(strsplit(txt, ','))) %in% toupper(oncoGenes))))

standardizedDereplicatedSites$geneMark <- ifelse(
    isNearOnco | isInOnco,
    paste0(standardizedDereplicatedSites$geneMark, "~"),
    standardizedDereplicatedSites$geneMark )


# * in transcription units
standardizedDereplicatedSites$geneMark <- ifelse(
    toupper(standardizedDereplicatedSites$inGene)!="FALSE",
    paste0(standardizedDereplicatedSites$geneMark, "*"),
    standardizedDereplicatedSites$geneMark)


# There will be instances where a comma delimited list of gene names are used by
# nearest_refSeq_gene. Go through each gene name, if 1 or more in the humanLymphGenes list
# then use those names. If not , check the oncoGenes. If none of the genes are in either of 
# those two lists, just use the first name in the list.

geneLabels <- unlist(lapply(standardizedDereplicatedSites$nearest_refSeq_gene, function(x){
   g <- toupper(unlist(strsplit(x, ',', fixed = TRUE)))
  
   if (any(g %in% toupper(humanLymphGenes))){
      r <- paste0(g[g %in% toupper(humanLymphGenes)], collapse=',')
   } else if (any(g %in% toupper(oncoGenes))){
      r <- paste0(g[g %in% toupper(oncoGenes)], collapse=',')
   } else {
      r <- g[1]
   }

   r
}))

standardizedDereplicatedSites$nearest_refSeq_gene <- geneLabels


# attach gene marks
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



order_barplot <- function(df){
  df$timepointCelltype <- paste0(df$Timepoint, '/', df$CellType)
  df <- do.call("rbind", by(df, df$timepointCelltype, function(x){
    a <- x[grepl('LowAbund', x$maskedRefGeneName),]
    b <- x[!grepl('LowAbund', x$maskedRefGeneName),]
    b <- b[order(b$estAbund, decreasing=TRUE),]
    rbind(a,b)
  }))
  df$timepointCelltype <- NULL
  df
} 


barplotAbunds <- bind_rows(lapply(barplotAbunds, order_barplot))
barplotAbundsBySample <- bind_rows(lapply(barplotAbundsBySample, order_barplot))

CellType_order <- unique(standardizedDereplicatedSites$CellType)
barplotAbunds$CellType <- factor(barplotAbunds$CellType, levels=CellType_order)
barplotAbundsBySample$CellType <- factor(barplotAbundsBySample$CellType, levels=CellType_order)

barplotAbunds$maskedRefGeneName <- factor(barplotAbunds$maskedRefGeneName,
                                          levels=unique(barplotAbunds$maskedRefGeneName))

if ('LowAbund' %in% levels(barplotAbunds$maskedRefGeneName)){
  barplotAbunds$maskedRefGeneName <- relevel(barplotAbunds$maskedRefGeneName, 'LowAbund')
}

barplotAbundsBySample$maskedRefGeneName <- factor(barplotAbundsBySample$maskedRefGeneName,
                                                  levels=unique(barplotAbundsBySample$maskedRefGeneName))

if ('LowAbund' %in% levels(barplotAbundsBySample$maskedRefGeneName)){
  barplotAbundsBySample$maskedRefGeneName <- relevel(barplotAbundsBySample$maskedRefGeneName, 'LowAbund')
}


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
knit(file.path(arguments$codeDir, "GTSPreport.Rmd"), output=mdfile)
markdownToHTML(mdfile, htmlfile, extensions=c('tables'),
               options=c(markdownHTMLOptions(defaults=T), "toc"),
               stylesheet=file.path(arguments$codeDir, "GTSPreport.css") )


message("\nReport ", htmlfile, " is generated from ", arguments$sample_gtsp)
save.image(RDataFile)
}
