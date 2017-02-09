#### load all the libraries for plotting and other calculations ####
libs <- c('hiAnnotator','knitr','ggplot2','scales','plyr','pheatmap','xtable',
          'reldist','vegan','reshape','reshape2','hiReadsProcessor','sonicLength',
          'doParallel','RColorBrewer','Vennerable','directlabels','markdown',
          'PubMedWordcloud')
sapply(libs, require, character.only=T)
theme_set(theme_bw())

#### functions to streamline analysis calculations ####
gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}

sanitize <- function(string) {
  result <- gsub("&", "\\&", string, fixed = TRUE)
  result <- gsub("_", "\\_", result, fixed = TRUE)
  result
}

timepointDays <- function(y) {
  suppressWarnings(ifelse(grepl("^m|m$",y),
                          as.numeric(gsub("m","",y))*30.5,
                          as.numeric(gsub("d(\\d+).*","\\1",y,ignore.case=T))))
}

sortTimePoints <- function(x) {
  sort(sapply(as.character(unique(x)), 
              function(y) timepointDays(y)
  ))
}

getPopEstimates <- function(counts) {
  if(all(dim(counts)>2) & is.matrix(counts)) {
    # do jackknife correction if more than 1 samples/replicates are available
    pseudo <- ncol(counts)*
      estimateR(rowSums(counts))-(ncol(counts)-1)*
      sapply(1:ncol(counts), function(y) estimateR(rowSums(counts[,-y])))
    # rowMeans(pseudo)[c(1,2,4)] #drop invalid se's
    round(rowMeans(pseudo)["S.chao1"])
  } else {
    pseudo <- estimateR(rowSums(counts))
    # pseudo[c(1,2,4)]
    round(pseudo["S.chao1"])
  }
}

getEstAbund <- function(posID, fragLen, group, replicate=NULL, parallel=TRUE, 
                        clusterfragLen=FALSE) {
  if(parallel) {
    cl <- makeCluster(3)
    registerDoParallel(cl)
  }
  
  dfr <- data.frame(posID, fragLen, group, row.names=NULL, stringsAsFactors=FALSE)
  if(!is.null(replicate)) {
    dfr$replicate <- replicate
  } else {
    dfr$replicate <- 1
  }
 
  if(clusterfragLen) {
    corrected <- with(dfr, 
                      clusterSites(paste(posID,replicate,sep=":"), 
                                   fragLen, group, parallel=parallel))
    newValues <- with(corrected,
                      split(clusteredValue, paste(posID,value)))
    dfr$fragLen <- as.numeric(newValues[with(dfr, 
                                             paste(paste(posID,replicate,sep=":"),
                                                   fragLen))])
  }
  
  counts.fragLen <- count(count(dfr, c("group","posID","fragLen"))[,-4],
                          c("group","posID"))
  names(counts.fragLen)[3] <- "fragLenCounts"
  
  dfr <- unique(dfr)
  dfr <- split(dfr, dfr$group)
  
  res <- mclapply(dfr, function(x) {
    dummy.theta <- structure(rep(1,length(x$posID)), names=x$posID)
    if(length(unique(x$replicate))>1) {
      siteAbund <- tryCatch(with(x, estAbund(factor(posID), fragLen, 
                                             factor(replicate))),  
                            error = function(z) list("theta"=dummy.theta))     
    } else {
      siteAbund <- tryCatch(with(x, estAbund(factor(posID), fragLen)),  
                            error = function(z) list("theta"=dummy.theta))
    }
    x$estAbund <- round(siteAbund$theta)[x$posID]
    x
  }, mc.cores = ifelse(parallel,getDoParWorkers(),1))
  
  if(parallel) {
    stopCluster(cl)
  }
  res <- do.call(rbind, res)
  rownames(res) <- NULL
  merge(unique(res[,c("group","posID","estAbund")]), counts.fragLen, all.x=TRUE)
}

getRanks <- function(value, posID, grouping, colname="valueRanks") {
  sums <- ddply(data.frame(value,posID,grouping,stringsAsFactors=FALSE),
                .(posID,grouping), summarize, value.sum=sum(value))
  sums <- arrange(sums, grouping, plyr::desc(value.sum))
  ranks <- with(sums, tapply(value.sum, grouping, rank))
  stopifnot(class(ranks)=="array")
  maxranks <- lapply(ranks, max)
  ranks <- mapply(function(x,y) (y+1)-x, ranks, maxranks)
  ranks <- lapply(ranks, function(x) as.numeric(as.factor(x)))
  sums[,colname] <- as.numeric(unlist(ranks, use.names=FALSE))
  return(sums[,c(1,2,4)])
}

getPropsAndRanks <- function(value, posID, grouping, coreColname) {
  
  stopifnot(!any(is.na(posID)))
  stopifnot(!any(is.na(grouping)))
  if(any(is.na(value))) {
    message("Found NAs in value variable...using 1 as replacement!")
    value[is.na(value)] <- 1
  }
  
  # IMPORTANT: work with unique sites only...
  # sum of Proportions for each combination of grouping should be 1! #
  dat <- data.frame(value, posID, grouping, stringsAsFactors=FALSE)
  if(any(duplicated(dat))) {
    message("taking unique values to perform calculations")
    dat <- unique(dat)
  }
  
  # get proportions #
  valueSums <- with(dat, tapply(value, paste0(grouping,posID), sum))
  dat$valueSums <- valueSums[with(dat, paste0(grouping,posID))]
  
  groupSums <- with(dat, tapply(value, grouping, sum))
  dat$groupSums <- groupSums[as.character(dat$grouping)]
  
  propColname <- paste0(coreColname,"Prop")
  dat[, propColname] <- as.numeric(with(dat, valueSums/groupSums))
  dat$valueSums <- NULL
  dat$groupSums <- NULL
  
  # take unique proportions when doing ranking since data has already been summed!
  dat <- unique(dat[,c("posID", "grouping", propColname)])
  
  stopifnot(all(round(tapply(dat[,propColname], dat$grouping, sum), digits=3)==1))
  
  # get ranks #     
  dat <- merge(dat, with(dat,
                         getRanks(get(propColname), posID, grouping, 
                                  paste0(propColname,"Rank"))))
  
  dat
}

aggregateSiteTypes <- function(dat, groupingVar="Alias", posIDvar="posID",
                               siteTypeVar="siteType",
                               valueVar="estAbundance1") {
  stopifnot(groupingVar %in% names(dat))        
  stopifnot(posIDvar %in% names(dat))
  stopifnot(siteTypeVar %in% names(dat))
  
  # IMPORTANT: work with unique sites only...    
  uniques <- unique(dat[,c(posIDvar, valueVar, siteTypeVar, groupingVar)])
  props <- getPropsAndRanks(value=uniques[,valueVar], posID=uniques[,posIDvar],
                            grouping=uniques[,groupingVar], coreColname="")
  
  props <- merge(props, unique(uniques[,c(posIDvar, groupingVar, siteTypeVar)]), 
                 by.y=c(posIDvar, groupingVar), 
                 by.x=c("posID","grouping"))
  names(props)[names(props)=="grouping"] <- groupingVar
  
  sums <- ddply(.data=props, c(groupingVar, siteTypeVar), summarize,
                Props=sum(Prop), Ranks=sum(PropRank))
  sums$siteType <- sapply(strsplit(sums[,siteTypeVar],","), tail, 1)
#   sums$siteType[nchar(sums$siteType)>10] <- 
#     with(sums,
#          abbreviate(siteType[nchar(siteType)>10], 
#                     minlength = mean(nchar(siteType))))
#   

  # sum of proportions for each combination of alias should be 1! #
  test <- ddply(.data=sums, groupingVar, summarize, Props=sum(Props))  
  stopifnot(all(round(test$Props)==1))    
  
  sums
}
