#### load up require packages + objects #### 
source("helper_functions.R")
source("../intSiteRetriever/intSiteRetriever.R")
source("../genomicHeatmapMaker/CancerGeneList/onco_genes.R")

#INPUTS: either GTSP numbers or patient name
#        csv file/table GTSP to sampleName


#### Update Mysql table for each sample within a trial for newly processed samples ####
# updateSetstable <- function(dbConn) {
# 	# extract patient, timepoint, celltype, and GTSP identifiers from the setname #
# 	sqls <- c("UPDATE intsites.genetherapy_samples SET
#             celltype=SUBSTRING_INDEX(SUBSTRING_INDEX(setname,'-',-5),'-',1),
#             time=SUBSTRING_INDEX(SUBSTRING_INDEX(setname,'-',-4),'-',1),
#             patient=replace(SUBSTRING_INDEX(SUBSTRING_INDEX(setname,'-',-6),'-',1),
#                                             'XSCIDp',''),
#             SpecimenAccNum=SUBSTRING_INDEX(SUBSTRING_INDEX(setname,'-',-3),'-',1),
#             enzyme='FRAG', gender='M' WHERE celltype = '' 
#             AND setname like '%GTSP%' AND trial='FirstSCID'",
#             
#             "UPDATE intsites.genetherapy_samples SET
#             celltype=SUBSTRING_INDEX(SUBSTRING_INDEX(setname,'-',-5),'-',1),
#             time=SUBSTRING_INDEX(SUBSTRING_INDEX(setname, '-',-4),'-',1),
#             patient=replace(SUBSTRING_INDEX(SUBSTRING_INDEX(setname,'-',-6),'-',1),
#                                             'XSCID',''),
#             SpecimenAccNum=SUBSTRING_INDEX(SUBSTRING_INDEX(setname,'-',-3),'-',1),
#             enzyme='FRAG', gender='M' WHERE celltype = '' 
#             AND setname like '%GTSP%' AND trial='sinSCID'",
#             
#             "UPDATE intsites.genetherapy_samples SET
#             celltype=SUBSTRING_INDEX(SUBSTRING_INDEX(setname,'-',-5),'-',1),
#             time=SUBSTRING_INDEX(SUBSTRING_INDEX(setname,'-',-4),'-',1),
#             patient=SUBSTRING_INDEX(SUBSTRING_INDEX(setname,'-',-6), '-',1),
#             SpecimenAccNum=SUBSTRING_INDEX(SUBSTRING_INDEX(setname,'-',-3),'-',1),
#             enzyme='FRAG', gender='M' WHERE celltype = '' 
#             AND setname like '%GTSP%' AND trial='WasLenti'",			  
#             
#             "UPDATE intsites.genetherapy_samples SET
#             celltype=SUBSTRING_INDEX(SUBSTRING_INDEX(setname,'-',-5),'-',1),
#             time=SUBSTRING_INDEX(SUBSTRING_INDEX(setname,'-',-4),'-',1),
#             patient=replace(SUBSTRING_INDEX(SUBSTRING_INDEX(setname,'-',-6),'-',1),
#                                             'bThal',''),
#             SpecimenAccNum=SUBSTRING_INDEX(SUBSTRING_INDEX(setname,'-',-3),'-',1), 
#             enzyme='FRAG' WHERE celltype = '' 
#             AND setname like '%GTSP%' AND trial='betaThal'",
#             
#             "UPDATE intsites.genetherapy_samples SET gender = if(patient='pPLB','M','F'),
#             time=if(time like '%preinfusio%','d0',time)",
#             
#             "UPDATE intsites.genetherapy_samples SET time = REPLACE(time,'_','.')",
#             
#             "UPDATE intsites.geneTherapy_samples set time='d0' WHERE time like '%hit%' or   
#             time like '%d0.%'",
#             
#             "UPDATE intsites.genetherapy_samples JOIN specimen_management.GTSP 
#             USING (SpecimenAccNum) 
#             SET intsites.genetherapy_samples.location=specimen_management.GTSP.location",
#             
# 	          "UPDATE intsites.genetherapy_samples JOIN specimen_management.GTSP 
# 	          USING (SpecimenAccNum) 
# 	          SET intsites.genetherapy_samples.VCN=specimen_management.GTSP.VCN",
#             
# 	          "UPDATE intsites.genetherapy_samples JOIN intsites.sets on name=setname 
#             SET intsites.genetherapy_samples.passingsetsize=intsites.sets.passingsetsize 
#             WHERE intsites.genetherapy_samples.passingsetsize=0")
# 	
# 	sapply(sqls, function(x) dbSendQuery(dbConn,x))
# 	
# 	# fix the new patient identifier to match the old patient naming schema for FIRST SCID TRIAL
# 	sqls <- c("UPDATE intsites.genetherapy_samples SET patient='1' 
#             WHERE patient='ML' and setname like '%GTSP%' AND trial='FirstSCID'", 
#             "UPDATE intsites.genetherapy_samples SET patient='2' 
#             WHERE patient='RN' and setname like '%GTSP%' AND trial='FirstSCID'", 
#             "UPDATE intsites.genetherapy_samples SET patient='4' 
#             WHERE patient='CW' and setname like '%GTSP%' AND trial='FirstSCID'", 
#             "UPDATE intsites.genetherapy_samples SET patient='5' 
#             WHERE patient='RD' and setname like '%GTSP%' AND trial='FirstSCID'", 
#             "UPDATE intsites.genetherapy_samples SET patient='6' 
#             WHERE patient='BA' and setname like '%GTSP%' AND trial='FirstSCID'", 
#             "UPDATE intsites.genetherapy_samples SET patient='7' 
#             WHERE patient='FM' and setname like '%GTSP%' AND trial='FirstSCID'", 
#             "UPDATE intsites.genetherapy_samples SET patient='8' 
#             WHERE patient='HR' and setname like '%GTSP%' AND trial='FirstSCID'", 
#             "UPDATE intsites.genetherapy_samples SET patient='9' 
#             WHERE patient='RE' and setname like '%GTSP%' AND trial='FirstSCID'", 
#             "UPDATE intsites.genetherapy_samples SET patient='10' 
#             WHERE patient='KA' and setname like '%GTSP%' AND trial='FirstSCID'")
# 	sapply(sqls, function(x) dbSendQuery(dbConn,x))
# }
# updateSetstable(dbConn)

####  get set names with other metadata #### 
# sets <- dbGetQuery(dbConn,"SELECT setName, celltype, patient, time as timepoint, 
#                    Trial, upper(enzyme) as enzyme, VCN FROM intsites.genetherapy_samples")

# sets$Trial[grepl("betaThal",sets$Trial,ignore.case=T)] <- "betaThal"
# fix month/day mislabelling of old betaThal samples from pPLB #
# rows <- with(sets,Trial=="betaThal" & patient=="pPLB" & grepl("m",timepoint))
# sets$timepoint[rows] <- sub("(\\d+)m","d\\1",sets$timepoint[rows])

# let's make sure there are no patients with the same name in different trials! #
# stopifnot(all(count(unique(sets[,c("patient","Trial")]))$freq==1))

####  add replicate information for abundance estimation by soniclength #### 
## fix d0 timepoints ##
sets$timepoint[sets$timepoint %in% c("d0","0")] <- "d0" #can be done in db as REGEXP in future
sets$timepointDay <- timepointDays(sets$timepoint)

#This information will be provided by the GTSP_sampleName table
#GTSP duplications mean there were >1 replicate for that GTSP

# sets$Replicate <- 1
# rows <- grepl("GTSP",sets$setName)
# sets$Replicate[rows] <- as.numeric(sub(".+-(\\d+)-.+","\\1",sets$setName[rows]))
# replicates <- sapply(with(sets,
#                           tapply(Replicate,
#                                  paste(celltype,patient,timepointDay,Trial),
#                                  unique)),
#                      function(x) 1:length(x), simplify=F)

## get qc passed sites for all sets ## 
#getUniqueSites from intSiteRetriever
sites.qc <- getSitesFromDB(dbConn, setName=sets$setName[rows], 
                           freeze="hg18", nrstOnco=T, nrst5pOnco=T, nrstGene=T,
                           nrst5pGene=T, inGene=T)
sites.qc$Position <- as.integer(sites.qc$Position)
sites.qc$setposid <- with(sites.qc, paste0(setName,Chr,strand,Position)) #new PK
sites.qc <- merge(sites.qc,sets,all.x=T,by="setName") #adding metadata

## get all sequences per site for all sets ##
#getUniqueBreakpoints from intSiteRetriever
sites <- getSitesFromDB(dbConn,setName=sets$setName[rows], freeze="hg18", allsites=T)
sites$Position <- as.integer(sites$Position)
sites$setposid <- with(sites,paste0(setName,Chr,Ort,Position))
sites$isqc98 <- sites$setposid %in% sites.qc$setposid
sites <- merge(sites,sets,all.x=T,by="setName")

sites$setposid <- NULL; sites.qc$setposid <- NULL;

## get the multihits for dangerously abundant hits in repeat elements ##
#will be handled by intSiteRetriever, eventually...
sites.multi <- getSitesFromDB(dbConn, setName=sets$setName[rows], 
                              freeze="hg18", multihit=T)
sites.multi <- subset(sites.multi, isMultiHit)
sites.multi$Position <- as.integer(sites.multi$Position)
sites.multi <- merge(sites.multi,sets,all.x=T,by="setName")
sites.rd <- makeGRanges(sites.multi, soloStart=TRUE, chromCol="Chr", strandCol="Ort")

## remove know vector contaminants ##
badass <- data.frame(Chr=c("chr11"), Ort=c("+"), Position=c(61776245),
                     stringsAsFactors=FALSE)
badass$posID <- with(badass, paste0(Chr,Ort,Position))

sites.qc <- droplevels(subset(sites.qc, 
                              !paste0(Chr,strand,Position) %in% badass$posID))
sites <- droplevels(subset(sites, 
                           !paste0(Chr,Ort,Position) %in% badass$posID))

##  get genes of interest ##
wantedgenes <- toupper(c('HIVEP3', 'VAV3', 'NOTCH2', 'ITPR1', 'FOXP1', 'MDS1', 
                         'C3ORF21', 'MAML3', 'TRIO', 'LOC441108', 'CXXC5', 'JARID2',
                         'IKZF1', 'ANGPT1', 'MLLT3', 'CDH23', 'C10ORF54', 'SWAP70', 
                         'ADRBK1', 'NINJ2', 'ETV6', 'HMGA2', 'CRADD', 'C14ORF43', 
                         'PRKCB1', 'GAS7', 'RAI1', 'FMNL1', 'MSI2', 'SEPT9', 'PREX1', 
                         'RUNX1', 'LMO2', 'CCND2', 'BMI1', 'EVI1'))

##  get all oncogenes ##
#this will be from the CancerGeneList R project
oncos <- dbGetQuery(dbConn,"select * from oncogenelists.allonco")$geneName
oncos <- unlist(strsplit(gsub(" ","",toupper(oncos)),"\\|"))
load("oncogenes.rl.Rdata")
oncos <- union(oncos,oncogenes.rl$Symbol)
oncogenes.rl <- as(oncogenes.rl,"GRanges")

allgenes.rd <- makeGRanges(
  dbGetQuery(dbConn,"(select distinct geneName, Chrom, strand, txStart, txEnd from hg18.refflat) UNION (select distinct kgXref.geneSymbol as geneName, knownGene.chrom as Chrom, knownGene.strand as strand, txStart, txEnd FROM hg18.knownGene JOIN hg18.kgXref ON knownGene.name=kgXref.kgID)"))

#ACTUALLY START DOING ANNOTATIONS HERE
sites.rd <- getNearestFeature(sites.rd, allgenes.rd, colnam="nrstGene")  
sites.rd <- getNearestFeature(sites.rd, allgenes.rd, colnam="nrstGene", side="5p") 
sites.rd <- getSitesInFeature(sites.rd, allgenes.rd, colnam="inGene")
rm(allgenes.rd)

sites.rd <- getNearestFeature(sites.rd, oncogenes.rl, 
                              colnam="nrstOncoGene", feature.colnam="Symbol")  
sites.rd <- getNearestFeature(sites.rd, oncogenes.rl, side="5p", 
                              colnam="nrstOncoGene", feature.colnam="Symbol") 
rm("oncogenes.rl")

allIntSites <- as.data.frame(sites.rd)
allIntSites$seqnames <- NULL;
allIntSites$start <- NULL;
allIntSites$end <- NULL;
allIntSites$width <- NULL;
sites.multi <- allIntSites
rm("sites.rd","allIntSites")
#has single hits and multihits at this point
sites.multi <- merge(sites.multi, count(sites.multi,"Sequence"))

####  setup gene columns for later analysis #### 
sites.qc$nrstOncoGeneName <- toupper(as.character(sites.qc$nrstOncoGeneName))
sites.qc$nrstGene <- toupper(as.character(sites.qc$nrstGene))
sites.qc$X5pnrstGene <- toupper(as.character(sites.qc$X5pnrstGene))
sites.qc$inGene <- toupper(as.character(sites.qc$inGene))

#think about alternate schema for this, using colors, fonts (bold, etc.)
sites.qc$geneType <- with(sites.qc,ifelse(inGene=="FALSE",
                                          nrstGene, 
                                          paste0(inGene,"*")))
sites.qc$geneType <- with(sites.qc,
                          ifelse(unlist(lapply(lapply(strsplit(inGene,","),
                                                      "%in%",oncos),any)) | 
                                   unlist(lapply(lapply(strsplit(nrstGene,","),
                                                        "%in%",oncos),any)), 
                                 paste0(geneType,"~"), geneType))
sites.qc$geneType <- with(sites.qc,
                          ifelse(unlist(lapply(lapply(strsplit(inGene,","),
                                                      "%in%",wantedgenes),any)) | 
                                   unlist(lapply(lapply(strsplit(nrstGene,","),
                                                        "%in%",wantedgenes),any)), 
                                 paste0(geneType,"!"), geneType))

sites.qc$inTheGene <- with(sites.qc,ifelse(inGene=="FALSE",FALSE,TRUE))
sites.qc$oncoWithin50kb <- abs(sites.qc$nrstOncoGeneDist)<=50000 | 
  abs(sites.qc$nrst5pOncoGeneDist)<=50000
sites.qc$wantedWithin50kb <- unlist(lapply(lapply(strsplit(sites.qc$inGene,","),
                                                  "%in%",wantedgenes),any)) | 
  with(sites.qc,
       (nrstOncoGeneName %in% wantedgenes & abs(nrstOncoGeneDist)<=50000) | 
         (nrst5pOncoGeneName %in% wantedgenes & abs(nrst5pOncoGeneDist)<=50000) |
         (nrstGene %in% wantedgenes & abs(nrstGeneDist)<=50000) | 
         (X5pnrstGene %in% wantedgenes & abs(X5pnrstGeneDist)<=50000))

sites.multi$nrstGene <- toupper(as.character(sites.multi$nrstGene))
sites.multi$X5pnrstGene <- toupper(as.character(sites.multi$X5pnrstGene))
sites.multi$nrstOncoGene <- toupper(as.character(sites.multi$nrstOncoGene))
sites.multi$X5pnrstOncoGene <- toupper(as.character(sites.multi$X5pnrstOncoGene))
sites.multi$inGene <- toupper(as.character(sites.multi$inGene))
sites.multi$geneType <- with(sites.multi,ifelse(inGene=="FALSE",
                                                nrstGene,
                                                paste0(inGene,"*")))
sites.multi$geneType <- with(sites.multi,
                             ifelse(unlist(lapply(lapply(strsplit(inGene,","),
                                                         "%in%",oncos),any)) | 
                                      unlist(lapply(lapply(strsplit(nrstGene,","),
                                                           "%in%",oncos),any)), 
                                    paste0(geneType,"~"), geneType))
sites.multi$geneType <- with(sites.multi,
                             ifelse(unlist(lapply(lapply(strsplit(inGene,","),
                                                         "%in%",wantedgenes),any)) |
                                      unlist(lapply(lapply(strsplit(nrstGene,","),
                                                           "%in%",wantedgenes),any)),
                                    paste0(geneType,"!"), geneType))
sites.multi$inTheGene <- with(sites.multi, ifelse(inGene=="FALSE", FALSE, TRUE))
sites.multi$wantedWithin50kb <- unlist(lapply(lapply(strsplit(sites.multi$inGene,","),
                                                     "%in%",wantedgenes),any)) |
  with(sites.multi, 
       (nrstGene %in% wantedgenes & abs(nrstGeneDist)<=50000) | 
         (X5pnrstGene %in% wantedgenes & abs(X5pnrstGeneDist)<=50000))
sites.multi$oncoWithin50kb <- abs(sites.multi$nrstOncoGeneDist)<=50000 | 
  abs(sites.multi$X5pnrstOncoGeneDist)<=50000

#### create an alias + alias-posid column to pool sites and perform subsequent analysis #### 
sites$posID <- with(sites,paste0(Chr,Ort,Position))
sites.qc$posID <- with(sites.qc,paste0(Chr,strand,Position))
sites.multi$posID <- with(sites.multi,paste0(Chr,strand,Position))

sites$Alias <- with(sites,paste(timepointDay,celltype,sep=":"))
sites.qc$Alias <- with(sites.qc,paste(timepointDay,celltype,sep=":"))
sites.multi$Alias <- with(sites.multi,paste(timepointDay,celltype,sep=":"))

sites$Aliasposid <- with(sites,paste0(Alias,posID))
sites.qc$Aliasposid <- with(sites.qc,paste0(Alias,posID))

pats.to.do <- as.character(unique(sites.qc$patient))

sites <- arrange(sites, Trial,patient,timepointDay,celltype,Chr,Position) 
sites.qc <- arrange(sites.qc, Trial,patient,timepointDay,celltype,Chr,Position)
sites.multi <- arrange(sites.multi, Trial,patient,timepointDay,celltype,Chr,Position) 

message("*** done updating...doing analysis ***")

##### some analysis: do this after updating/adding new data since abundance relies on # of sites per pool. #####
## do it by patients to put less stress on computing ##

cleanit <- gc()
cols <- c('posID','Chr','strand','Position','Sequence','Alias',
          'patient','celltype','timepoint','timepointDay','isMultiHit')
sites.all <- unique(rbind(sites.multi[,cols], 
                          cbind(sites.qc,isMultiHit=FALSE)[,cols]))

if(length(pats.to.do)>0) {
  sites <- split(sites, sites$patient)[pats.to.do]
  sites.qc <- split(sites.qc, sites.qc$patient)[pats.to.do]
  sites.multi <- split(sites.multi, sites.multi$patient)[pats.to.do]
  sites.all <- split(sites.all, sites.all$patient)[pats.to.do]
  
  # write out the data by patient to save on memory! #    
  message("*** Saving new patient data objects ***")
  for(pat in pats.to.do) {
    trial <- unique(sites.qc[[pat]]$Trial)    
    if(!file.exists(file.path(dataDir, trial))) {
      dir.create(file.path(dataDir, trial))
    }
    filename <- file.path(dataDir, trial, paste0(paste(trial,pat,sep="_"),".RData"))
    allpatdata <- list("sites"=sites[[pat]], "sites.qc"=sites.qc[[pat]], 
                       "sites.multi"=sites.multi[[pat]], "sites.all"=sites.all[[pat]])
    save(allpatdata, file=filename)
  }
  
  rm("sites","sites.qc","sites.multi","sites.all","allpatdata")
  cleanit <- gc()

  ## find all the new patient data files ##
  oldData <- list.files(path=dataDir, pattern=".*.RData", full.names=T, recursive=T)
	
  for(pat in pats.to.do) {
    sec.num.messg <- 1
    message(pat)
    
    message(sec.num.messg,") Loading newly saved data.")
    filename <- grep(paste0("_",pat,".RData"), oldData, value=T,fixed=T)
    load(filename)
      
      # re-create alias-posid column to pool sites and perform subsequent analysis #
      allpatdata$sites$posID <- with(allpatdata$sites, 
                                     paste0(Chr,Ort,Position))
      allpatdata$sites.qc$posID <- with(allpatdata$sites.qc, 
                                        paste0(Chr,strand,Position))
      
      allpatdata$sites$Aliasposid <- with(allpatdata$sites, 
                                          paste0(Alias,posID))
      allpatdata$sites.qc$Aliasposid <- with(allpatdata$sites.qc, 
                                             paste0(Alias,posID))
    }
    
    ## combine unique sites with multihits & get OTUs ##	
#     sec.num.messg <- sec.num.messg + 1
#     message(sec.num.messg,") Making OTUs by Alias.")
#     cl <- makeCluster(3)
#     registerDoParallel(cl)
#     
#     # lets isolate really big Aliases since they can occupy crap load of memory #
#     aliasCounts <- table(allpatdata$sites.all$Alias) > 200000
#     if(any(aliasCounts)) {
#       message("Skipping OTU step for: ",
#               paste(names(which(aliasCounts)),collapse=","))
#       #       sites.all.otus <- lapply(names(which(aliasCounts)), function(x)
#       #         with(droplevels(subset(allpatdata$sites.all,Alias==x)), 
#       #              otuSites2(posID=paste0(Chr,strand), value=Position, 
#       #                        readID=Sequence, grouping=Alias, parallel=F))
#       #       )
#       #       sites.all.otus <- do.call(rbind, sites.all.otus)
#       sites.all.otus <- with(droplevels(subset(allpatdata$sites.all,
#                                                !Alias %in% names(which(aliasCounts)))), 
#                              otuSites2(posID=paste0(Chr,strand), value=Position, 
#                                        readID=Sequence, grouping=Alias))    
#     } else {
#       sites.all.otus <- with(allpatdata$sites.all, 
#                              otuSites2(posID=paste0(Chr,strand), value=Position, 
#                                        readID=Sequence, grouping=Alias))
#     }
#     stopCluster(cl)
#     
#     allpatdata$sites.all <- merge(arrange(allpatdata$sites.all, Sequence), 
#                                   arrange(unique(sites.all.otus[,c("readID","otuID")]),
#                                           readID), 
#                                   by.x="Sequence", by.y="readID")
#     rm(sites.all.otus)
#     
#     # need to remove this since an OTU should uniquely define a row!
#     #allpatdata$sites.all$Sequence <- NULL     
#     #allpatdata$sites.all <- unique(allpatdata$sites.all)
#     
#     # remove multihit indicator for sites that merged with singleton!
#     bore <- subset(melt(with(allpatdata$sites.all, 
#                              tapply(isMultiHit, list(Alias,otuID), 
#                                     function(x) !any(!x)))),
#                    !is.na(value))
#     names(bore) <- c('Alias','otuID','isMultiHit2')
#     allpatdata$sites.all <- merge(allpatdata$sites.all,bore,all.x=T)
#     allpatdata$sites.all$AliasOTUid <- with(allpatdata$sites.all, paste0(Alias,otuID))
#     rm(bore)
#     
#     # trickle OTUs back to original multi and all sites frames for estAbund calcs #
#     if(!is.null(allpatdata$sites.multi)) {
#       allpatdata$sites.multi <- arrange(allpatdata$sites.multi, Alias, posID)
#       bore <- arrange(unique(allpatdata$sites.all[,c('Alias','posID','otuID')]), 
#                       Alias, posID)
#       
#       # lets make sure there is only one otuID for each posID per Alias
#       #test <- with(bore, tapply(otuID, paste0(Alias,posID), unique))
#       #stopifnot(!any(sapply(test,length)>1))
#       
#       allpatdata$sites.multi <- merge(allpatdata$sites.multi, bore, all.x=TRUE)      
#       allpatdata$sites.multi$AliasOTUid <- with(allpatdata$sites.multi, 
#                                                 paste0(Alias,otuID))
#       rm(bore)
#     }
#     
#     allpatdata$sites <- 
#       merge(arrange(allpatdata$sites,Alias,posID),
#             arrange(unique(allpatdata$sites.all[,c('Alias','posID','otuID')]),
#                     Alias,posID), 
#             all.x=TRUE)
#     allpatdata$sites$AliasOTUid <- with(allpatdata$sites,paste0(Alias,otuID))
#         
    # get estimated abundance by break points per site #
    sec.num.messg <- sec.num.messg + 1
    message(sec.num.messg,") Doing SonicLength abundance")	
    allpatdata$sites.qc$estAbundance1 <- 0
    allpatdata$sites.all$estAbundance1 <- 0
    if(!is.null(allpatdata$sites.multi)) {
      allpatdata$sites.multi$estAbundance <- 0
    }
    
  #talk to Eric when re-working sonic abundance code
# it's better to sapply through a list of data frames rather than do a single call
# to estAbund using Alias as sampleName
    dfr <- droplevels(unique(subset(allpatdata$sites,isqc98,
                                    select=c(Aliasposid,qEnd,Replicate,Alias), 
                                    drop=T)))
    replicates <- if(length(unique(dfr$Replicate))>1) {dfr$Replicate} else {NULL}
    siteAbund <- with(dfr, getEstAbund(posID=Aliasposid, fragLen=qEnd, group=Alias, 
                                       replicate=replicates, parallel=TRUE, 
                                       clusterfragLen=FALSE))
    stopifnot(all(table(siteAbund$posID)==1))
    rows <- match(allpatdata$sites.qc$Aliasposid, siteAbund$posID)
    stopifnot(!any(is.na(rows)))
    allpatdata$sites.qc$estAbundance1 <- siteAbund$estAbund[rows]
    cleanit <- gc()  
    
    # get estimated abundance by break points per multihit site/OTU #		
#     dfr <- droplevels(unique(
#       rbind(allpatdata$sites[,c("AliasOTUid","qEnd","Replicate","Alias")],
#             allpatdata$sites.multi[,c("AliasOTUid","qEnd","Replicate","Alias")])))
#     replicates <- if(length(unique(dfr$Replicate))>1) {dfr$Replicate} else {NULL}
#     siteAbund <- with(dfr, getEstAbund(posID=AliasOTUid, fragLen=qEnd, group=Alias, 
#                                        replicate=replicates, parallel=TRUE, 
#                                        clusterfragLen=FALSE))      
#     stopifnot(all(table(siteAbund$posID)==1))
#     rows <- match(allpatdata$sites.all$AliasOTUid, siteAbund$posID)
#     allpatdata$sites.all$estAbundance1 <- siteAbund$estAbund[rows]    
#     if(!is.null(allpatdata$sites.multi)) {
#       rows <- match(allpatdata$sites.multi$AliasOTUid, siteAbund$posID)
#       allpatdata$sites.multi$estAbundance <- siteAbund$estAbund[rows]
#     }
#     rm(rows)
#     cleanit <- gc()  
#     
    # get estAbundance1 proportions & ranks#
    # since estAbundance1 is same for a posID in different replicate..take unique! #
    res <- with(allpatdata$sites.qc, 
                getPropsAndRanks(estAbundance1, Aliasposid, Alias, "estAbundance1"))
    stopifnot(all(table(res$posID)==1))
    rows <- match(allpatdata$sites.qc$Aliasposid, res$posID)
    allpatdata$sites.qc$estAbundance1Prop <- res$estAbundance1Prop[rows]
    allpatdata$sites.qc$estAbundance1PropRank <- res$estAbundance1PropRank[rows]
    
    # get estAbundance proportions for multihits #
#     if(!is.null(allpatdata$sites.multi)) {
#       res <- with(allpatdata$sites.multi, 
#                   getPropsAndRanks(estAbundance, AliasOTUid, Alias, "estAbundance"))
#       stopifnot(all(table(res$posID)==1))
#       rows <- match(allpatdata$sites.multi$AliasOTUid, res$posID)
#       allpatdata$sites.multi$estAbundanceProp <- res$estAbundanceProp[rows]
#     }
    
    # get estAbundance1 proportions for all sites: multi + unique #
    res <- with(allpatdata$sites.all, 
                getPropsAndRanks(estAbundance1, AliasOTUid, Alias, "estAbundance1"))
    stopifnot(all(table(res$posID)==1))
    rows <- match(allpatdata$sites.all$AliasOTUid, res$posID)
    allpatdata$sites.all$estAbundance1Prop <- res$estAbundance1Prop[rows]

    # add site ranks for plotting and reports 
    sec.num.messg <- sec.num.messg + 1
    message(sec.num.messg,") Adding Ranks.")
    test <- droplevels(allpatdata$sites.qc)
     
    test$estAbundance1Rank <- NULL
    test <- merge(test,
                  with(unique(test[,c("estAbundance1","posID","Alias")]),
                       getRanks(estAbundance1, posID, Alias, 
                                "estAbundance1Rank")),
                  by.x=c("posID","Alias"), by.y=c("posID","grouping"))
    
    allpatdata$sites.qc <- test
    rm(test)
    cleanit <- gc()
    
#     if(!is.null(allpatdata$sites.multi)) {
#       allpatdata$sites.multi$estAbundance1Rank <- NULL
#       allpatdata$sites.multi <- 
#         merge(allpatdata$sites.multi, 
#               with(unique(allpatdata$sites.multi[,c("estAbundance","otuID","Alias")]), 
#                    getRanks(estAbundance, otuID, Alias, "estAbundance1Rank")),
#               by.x=c("otuID","Alias"), by.y=c("posID","grouping"))
#     }
    
    allpatdata$sites.all$estAbundance1Rank <- NULL
    allpatdata$sites.all <- 
      merge(allpatdata$sites.all, 
            with(unique(allpatdata$sites.all[,c("estAbundance1","otuID","Alias")]), 
                 getRanks(estAbundance1, otuID, Alias, "estAbundance1Rank")),
            by.x=c("otuID","Alias"), by.y=c("posID","grouping"))
    
    sec.num.messg <- sec.num.messg + 1
    message(sec.num.messg,") Saving objects")
    save(allpatdata, file=filename)
    rm("allpatdata")
    cleanit <- gc()   
  }
}

message("*** Done analyzing & updating Trial objects! ***")