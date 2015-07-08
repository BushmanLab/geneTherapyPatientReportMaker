options(stringsAsFactors = FALSE)

#### INPUTS: csv file/table GTSP to sampleName ####
csvfile <- "sampleName_GTSP.csv"
args <- commandArgs(trailingOnly=TRUE)

use.sonicLength <- FALSE

if( length(args)==1 ) {
  csvfile <- args[1]
}
if( length(args)==2 & args[2] == "-s"){
  csvfile <- args[1]
  use.sonicLength <- TRUE
}else if(length(args)==2 & args[2] != "-s"){
  message("Incorrect flags.")
  stop()
}
if( !file.exists(csvfile) ) stop(csvfile, "not found")

codeDir <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly=FALSE), value=T)))
if( length(codeDir)!=1 ) codeDir <- list.files(path="~", pattern="geneTherapyPatientReportMaker", recursive=TRUE, include.dirs=TRUE, full.names=TRUE)
if( length(codeDir)!=1 ) codeDir <- "~/geneTherapyPatientReportMaker"
stopifnot(file.exists(file.path(codeDir, "GTSPreport.css")))
stopifnot(file.exists(file.path(codeDir, "GTSPreport.Rmd")))

#### load up require packages + objects #### 
library("RMySQL") #also loads DBI
library("plyr")
library("dplyr")
library("markdown")
library("knitr")
library("PubMedWordcloud")
library("hiAnnotator")
library("ggplot2")
library("reldist")
library("sonicLength")
library("reshape2")
library("scales")
library("intSiteRetriever") 
## intSiteRetriever package was installed from github as follows
## git clone https://github.com/BushmanLab/intSiteRetriever.git
## cd intSiteRetriever
## git checkout remove_multihitClusterID
## R
## devtools::document()
## devtools::install()

source(file.path(codeDir, "utilities.R"))
source(file.path(codeDir, "specimen_management.R"))
source(file.path(codeDir, "estimatedAbundance.R"))
source(file.path(codeDir, "dereplicateSites.R"))
source(file.path(codeDir, "standardizeSites.R"))
source(file.path(codeDir, "read_site_totals.R"))
source(file.path(codeDir, "populationInfo.R"))
source(file.path(codeDir, "abundanceFilteringUtils.R"))

unlink("CancerGeneList", force=TRUE, recursive=TRUE)
cmd <- "git clone https://github.com/BushmanLab/CancerGeneList.git"
message(cmd)
stopifnot( system(cmd)==0 )
source("CancerGeneList/onco_genes.R")


#### load datasets and process them before knit #### 
message("\nReading csv from ", csvfile)
sampleName_GTSP <- read.csv(csvfile)
stopifnot(all(c("sampleName", "GTSP") %in% colnames(sampleName_GTSP)))
message("\nGenerating report from the following sets")
print(sampleName_GTSP)

junk <- sapply(dbListConnections(MySQL()), dbDisconnect)
dbConn <- dbConnect(MySQL(), group="intSitesDev237")
stopifnot(all(setNameExists(sampleName_GTSP$sampleName, dbConn)))

read_sites_sample_GTSP <- get_read_site_totals(sampleName_GTSP, dbConn)

sets <- get_metadata_for_GTSP(unique(sampleName_GTSP$GTSP))
# reports are for a single patient
stopifnot(length(unique(sets$Patient)) == 1)
patient <- sets$Patient[1]
# and for a single trial
stopifnot(length(unique(sets$Trial)) == 1)
trial <- sets$Trial[1]

# all GTSP in the database
stopifnot(nrow(sets) == length(unique(sampleName_GTSP$GTSP)))

##end INPUTS
sets[sets$Timepoint=="NULL", "Timepoint"] = "d0"

sets <- merge(sets, read_sites_sample_GTSP)
sets$Timepoint <- sortFactorTimepoints(sets$Timepoint)

junk <- sapply(dbListConnections(MySQL()), dbDisconnect)
dbConn <- dbConnect(MySQL(), group="intSitesDev237")
refGenomes <- getRefGenome(sampleName_GTSP$sampleName, dbConn)
# at present the whole report is done for one genome
stopifnot(length(unique(refGenomes$refGenome))==1)
freeze <- refGenomes[1, "refGenome"]

##==========GET AND PERFORM BASIC DEREPLICATION/SONICABUND ON SITES=============
message("Fetching unique sites and estimating abundance")
junk <- sapply(dbListConnections(MySQL()), dbDisconnect)
dbConn <- dbConnect(MySQL(), group="intSitesDev237")
sites <- merge(getUniquePCRbreaks(sampleName_GTSP$sampleName, dbConn), sampleName_GTSP)

#we really don't care about seqinfo - we just want a GRange object for easy manipulation
uniqueSites.gr <- GRanges(seqnames=Rle(sites$chr),
                          ranges=IRanges(start=pmin(sites$integration, sites$breakpoint),
                                         end=pmax(sites$integration, sites$breakpoint)),
                          strand=Rle(sites$strand))
mcols(uniqueSites.gr) <- sites[,c("sampleName", "GTSP")]

#standardize sites across all GTSPs
standardizedReplicatedSites <- standardizeSites(uniqueSites.gr)
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
                                              length)


#=======================ANNOTATE DEREPLICATED SITES==========================
#standard refSeq genes
refSeq_genes <- makeGRanges(
  getUCSCtable("refGene", "RefSeq Genes", freeze=freeze),
  freeze=freeze
)

standardizedDereplicatedSites <- getNearestFeature(standardizedDereplicatedSites,
                                                   refSeq_genes,
                                                   colnam="nearest_refSeq_gene",
                                                   feature.colnam="name2")

#oncogenes
oncogene_file <- "CancerGeneList/allonco_no_pipes.csv"
oncogenes <- get_oncogene_from_file(oncogene_file)

refSeq_oncogene <- refSeq_genes[is_onco_gene(refSeq_genes$name2, oncogenes)]

standardizedDereplicatedSites <- getNearestFeature(standardizedDereplicatedSites,
                                                   refSeq_oncogene,
                                                   colnam="NrstOnco",
                                                   side="5p",
                                                   feature.colnam="name2")


#===================GENERATE EXPANDED CLONE DATAFRAMES======================
#barplots
abundCutoff.barplots <- getAbundanceThreshold(standardizedDereplicatedSites, 10)

barplotAbunds <- getAbundanceSums(filterLowAbund(standardizedDereplicatedSites,
                                                 abundCutoff.barplots),
                                  c("CellType", "Timepoint"))

barplotAbunds <- arrange(barplotAbunds, estAbundProp)

#detailed abundance plot
abundCutoff.detailed <- getAbundanceThreshold(standardizedDereplicatedSites, 50)

detailedAbunds <- getAbundanceSums(filterLowAbund(standardizedDereplicatedSites,
                                                 abundCutoff.detailed),
                                  c("CellType", "Timepoint"))

categorySums <- sapply(split(detailedAbunds$estAbundProp,
                             detailedAbunds$maskedRefGeneName),sum)

detailedAbunds$maskedRefGeneName <- factor(detailedAbunds$maskedRefGeneName,
                                           levels=names(sort(categorySums)))


#================Longitudinal Behaviour===============================
longitudinal <- as.data.frame(standardizedDereplicatedSites)[,c("Timepoint",
                                                                "CellType",
                                                                "estAbundProp",
                                                                "posid")]
has_longitudinal_data <- length(unique(longitudinal$Timepoint)) > 1


#==================DETAILED REPORTS FOR BAD ACTORS=====================
badActors <- c("LMO2", "IKZF1", "CCND2", "HMGA2", "MECOM")

badActorData <- sapply(badActors, function(badActor){
  hasBadActor <- grepl(badActor, standardizedDereplicatedSites$X5pNrstOnco)
  badActorWithin100K <- abs(standardizedDereplicatedSites$X5pNrstOncoDist) <= 100000
  standardizedDereplicatedSites[hasBadActor & badActorWithin100K]
})


#==================SET VARIABLES FOR MARKDOWN REPORT=====================
timepoint <- levels(sets$Timepoint)

cols <- c("Trial", "GTSP", "Patient", "Timepoint", "CellType", 
          "TotalReads", "UniqueSites", "FragMethod", "VCN")
summaryTable <- arrange(sets,Timepoint,CellType)
summaryTable <- summaryTable[,cols]

##cols <- c("Patient", "Timepoint", "CellType", "UniqueSites",
##          "Replicates", "FragMethod", "VCN", "S.chao1", "Gini", "Shannon")
cols <- c("Patient", "Timepoint", "CellType", "UniqueSites",
          "Replicates", "FragMethod", "VCN", "Gini", "Shannon")
popSummaryTable <- merge(sets,  populationInfo, by.x="GTSP", by.y="group")
popSummaryTable <- arrange(popSummaryTable,Timepoint,CellType)
##popSummaryTable <- popSummaryTable[,cols]

cols <- c("Trial", "GTSP", "Replicates", "Patient", "Timepoint", "CellType", 
          "TotalReads", "UniqueSites", "FragMethod", "VCN", "Gini", "Shannon")
summaryTable <- popSummaryTable[,cols]
    
timepointPopulationInfo <- melt(timepointPopulationInfo, "group")

#==================Get abundance for multihit events=====================
message("Fetching multihit sites and estimating abundance")
junk <- sapply(dbListConnections(MySQL()), dbDisconnect)
dbConn <- dbConnect(MySQL(), group="intSitesDev237")
sites.multi <- merge( suppressWarnings(getMultihitLengths(sampleName_GTSP$sampleName, dbConn)), sampleName_GTSP)
if( nrow(sites.multi) > 0 ) {
    sites.multi <- sites.multi %>%
    group_by(multihitID) %>%
    mutate(replicate=as.integer(as.factor(sampleName)))
    
    dfr <- data.frame(ID=sites.multi$multihitID,
                      fragLength=sites.multi$length,
                      replicate=sites.multi$replicate)
    
    if(length(unique(dfr$replicate))==1){
        estimatedAbundances <- estAbund(dfr$ID, dfr$fragLength)
    }else{
        estimatedAbundances <- estAbund(dfr$ID, dfr$fragLength, dfr$replicate)
    }
    
    sites.multi <- subset(sites.multi, !duplicated(multihitID))
    sites.multi$estAbund <- round(estimatedAbundances$theta[as.character(sites.multi$multihitID)])
    
    sites.multi <- merge(sites.multi, sets, by="GTSP")
    sites.multi <- sites.multi %>% group_by(Patient, Timepoint, CellType) %>%
    mutate(Rank=rank(-estAbund, ties.method="max"))

}


write.csv(as.data.frame(standardizedDereplicatedSites),
          file=paste(trial, patient, "uniquehit.csv", sep="."))

##save.image("debug.RData")
##end setting variables for markdown report

fig.path <- paste(unique(trial), unique(patient),
                  format(Sys.Date(), format="%Y%m%d"), "Figures",
                  sep=".")

#### begin generating markdown ####
unlink(fig.path, force=TRUE, recursive=TRUE)
mdfile <- paste(unique(trial), unique(patient),
                format(Sys.Date(), format="%Y%m%d"), "md",
                sep=".")

htmlfile <- gsub("\\.md$",".html",mdfile)
options(knitr.table.format='html')
theme_set(theme_bw()) #for ggplot2
knit(file.path(codeDir, "GTSPreport.Rmd"), output=mdfile)
markdownToHTML(mdfile, htmlfile, extensions=c('tables'),
               options=c(markdownHTMLOptions(defaults=T), "toc"),
               stylesheet=file.path(codeDir, "GTSPreport.css") )

#### clean up ####
#unlink("figureByPatient", force=TRUE, recursive=TRUE)
unlink("CancerGeneList", force=TRUE, recursive=TRUE)
unlink(mdfile, force=TRUE, recursive=TRUE)

message("\nReport ", htmlfile, " is generated from ", csvfile)

