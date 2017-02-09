

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
    
    # Convert the GRanges to a data frame.
    # The strandedness of the range determines which position is the integration site and which is the breakpoint.
    df <- rbind.fill(
      lapply(split(sites, as.vector(strand(sites))), function(x){
        if (as.vector(strand(x[1])) == '+'){
          d <- data.frame(seqnames=as.vector(seqnames(x)), strand=as.vector(strand(x)), intSite=start(x), breakPoints=end(x), GTSP=x$GTSP)
        } else {
          d <- data.frame(seqnames=as.vector(seqnames(x)), strand=as.vector(strand(x)), intSite=end(x), breakPoints=start(x), GTSP=x$GTSP)
        }
        d
      }))
    
    df$posid <- paste0(df$seqnames, df$strand, df$intSite)
    df$samplePosid <- paste0(df$GTSP, '/', df$posid)
    
    # Aggregate the data frame on intSites while counting both the number of reads and unique fragments.
    df2 <- rbind.fill(by(df, df$samplePosid, function(x){
      x$reads <- nrow(x)
      aggregate(breakPoints ~ intSite+GTSP+seqnames+strand+reads, FUN=function(x){ length(unique(x)) }, data=x)
    }))
    
    # Convert data frame back to a GRange object
    dereplicatedSites <- GRanges(seqnames=df2$seqnames, strand=df2$strand, ranges= IRanges(start=df2$intSite, end=df2$intSite), GTSP=df2$GTSP, reads=df2$reads, breakPoints=df2$breakPoints)
    
    # Set the default estAbund to the number of reads
    dereplicatedSites$estAbund <- dereplicatedSites$breakPoints
  }

  dereplicatedSites$estAbundProp <- dereplicatedSites$estAbund/sum(dereplicatedSites$estAbund)
  dereplicatedSites$estAbundRank <- rank(-1*dereplicatedSites$estAbundProp, ties.method="max")
  
  dereplicatedSites
}
