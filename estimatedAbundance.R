

getEstimatedAbundance <- function(sites, use.sonicLength=FALSE){
  
  if (use.sonicLength){
     sites$posid = paste0(seqnames(sites), strand(sites), start(flank(sites, -1, start=T)))
     dfr <- data.frame("ID"=sites$posid,
                       "fragLength"=width(sites),
                       "replicate"=sites$replicate)
  
     # Sonic abundance will crash if given a vector of replicates with one value
     if(length(unique(dfr$replicate)) == 1){
       estimatedAbundances <- estAbund(dfr$ID, dfr$fragLength)
     }else{
       estimatedAbundances <- estAbund(dfr$ID, dfr$fragLength, dfr$replicate)
     }
  
     dereplicatedSites = granges(dereplicateSites(sites)) #instead of mcols()=NULL

     # Regenerate posid which is guaranteed to match since the input data was already run through the dereplicator
     dereplicatedSites$posid = paste0(seqnames(dereplicatedSites),
                                      strand(dereplicatedSites),
                                      start(flank(dereplicatedSites, -1, start=T)))
     
     dereplicatedSites$estAbund <- round(estimatedAbundances$theta) #estAbund preserves order
  } else {
    
    sites$posid = paste0(seqnames(sites), strand(sites), start(flank(sites, -1, start=T)))
    sites$w <- paste0(sites$replicate, '/', width(sites))
    
    o <- rbind.fill(lapply(split(sites, sites$posid), function(g){
      data.frame(posid=g$posid[1], reads=length(g), frags=length(unique(g$w)))
    }))
 
    sites.reduced <- flank(sites, -1, start=TRUE)
    dereplicatedSites <- unlist(reduce(sites.reduced, min.gapwidth = 0L))
    
    dereplicatedSites$posid <- paste0(seqnames(dereplicatedSites), strand(dereplicatedSites), start(dereplicatedSites))
    dereplicatedSites <- dereplicatedSites[match(o$posid, dereplicatedSites$posid)]
    
    dereplicatedSites$reads <- o$reads
    dereplicatedSites$frags <- o$frags
    dereplicatedSites$estAbund <- dereplicatedSites$frags
  }

  dereplicatedSites$estAbundProp <- dereplicatedSites$estAbund/sum(dereplicatedSites$estAbund)
  dereplicatedSites$estAbundRank <- rank(-1*dereplicatedSites$estAbundProp, ties.method="max")
  dereplicatedSites
}
