#### set variables for the analysis ####
freeze <- "hg18"
col.keys <- c("patient"="Pat", "timepoint"="Tp", "celltype"="Cell", "posID"="posID", 
              "estAbundance1"="BreakSites", "estAbundance1Prop"="RelativeAbundance",
              "estAbundance1Rank"="SiteRank", "ccPropPropRank"="SiteRank", 
              "nrst5pOncoGeneDist"="DistTo5pCancerGene", "geneType"="SiteType",
              "isMultiHit2"="isMultiHit", "siteType2"="SiteType", 
              "ccProp"="RelativeClonecounts", "clonecount"="Clonecount")
valueVars <- c("estAbundance1", "clonecount", "primeridcounts")
sampleVars <- c("patient", "timepoint", "celltype", "Trial")
positionVars <- c("posID")
positionVars.multi <- c("otuID")

# windows to test for the overlap analysis
windows <- c(5)

# top ranking sites to show #
abundCutoff <- 0.03

#### functions to streamline analysis calculations ####
source("helper_functions.R")
require(gridExtra)

#### find patient data objects ####
message("Loading previously saved data.")
dataDir <- "PatientData"
oldData <- list.files(path=dataDir, pattern=".*.RData", full.names=T, recursive=T)
if(length(oldData)==0) {  
  stop("No Patient Data files found in ", dataDir)
}

##### make reports by patients #####
for(patData in oldData) {
  # make a playing copy of the data objects #
  load(patData)
  sites.qc <- droplevels(allpatdata$sites.qc)
  sites.qc$GTSP <- sub(".+-(GTSP\\d+)-.+","\\1",sites.qc$setName)
  sites.qc$GTSP[sites.qc$enzyme!="FRAG"] <- ""
  sites.qc$geneType[sites.qc$geneType=="HMGA2,AK128707,AY387666*~!"] <- "HMGA2*~!"  
  
  if(!is.null(allpatdata$sites.multi)) {
    sites.multi <- droplevels(allpatdata$sites.multi)
    sites.multi$GTSP <- sub(".+-(GTSP\\d+)-.+","\\1",sites.multi$setName)
    sites.multi$GTSP[sites.multi$enzyme!="FRAG"] <- ""
  }

  sites.all <- droplevels(allpatdata$sites.all)
  rm(allpatdata)
  
  sites.qc$Alias <- with(sites.qc, paste(timepoint,celltype,sep=":"))  
  
  # make global variables for the patient #
  for(f in sampleVars) {
    assign(f, sort(unique(sites.qc[,f])))
  }
  
  timepoint <- sortTimePoints(timepoint)

  # make the reports #
  if(!file.exists(paste0("Reports/",Trial))) {
  	system(paste0("mkdir Reports/",Trial))
  }
  
  filename <- paste("Report",Trial,patient,"md",sep=".")
  outFilename <- gsub("\\.md",".html",filename)
  options(knitr.table.format = 'html')
  knit("GTSPReportByPatient.Rmd", output = filename)
  markdownToHTML(filename, outFilename, extensions=c('tables'),
                 options=c(markdownHTMLOptions(defaults=T),"toc"),
                 stylesheet="markdown_custom.css")  
  cmd <- sprintf("cp %s %s", outFilename, file.path("Reports",Trial,""))
  system(cmd)
  plotfiles <- list.files(path="figureByPatient", full.names=TRUE, pattern="*.png")
  file.remove(filename,outFilename, plotfiles)  
}