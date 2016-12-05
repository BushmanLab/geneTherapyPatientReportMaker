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

getEstimatedAbundance <- function(sites, use.sonicLength=FALSE){
  #inputs at this point have been standardized, so it's ok to use posid as PK
  if(use.sonicLength){
    estAbund.uniqueFragLen <- function(location, fragLen, replicate=NULL){
      if(is.null(replicate)){replicate <- 1}  #Need for downstream workflow
      dfr <- data.frame(location = location, fragLen = fragLen, 
                      replicate = replicate)
      dfr_dist <- distinct(dfr)
      site_list <- split(dfr_dist, dfr_dist$location)
      theta <- sapply(site_list, function(x){nrow(x)})
      theta <- theta[unique(dfr$location)]
      list(theta=theta)
    }
    estAbund <- estAbund.uniqueFragLen
  }
  
  sites$posid = paste0(seqnames(sites), strand(sites), start(flank(sites, -1, start=T)))
  dfr <- data.frame("ID"=sites$posid,
                    "fragLength"=width(sites),
                    "replicate"=sites$replicate)
  
  #sonic abundance will crash if given a vector of replicates with one value
  if(length(unique(dfr$replicate)) == 1){
    estimatedAbundances <- estAbund(dfr$ID, dfr$fragLength)
  }else{
    estimatedAbundances <- estAbund(dfr$ID, dfr$fragLength, dfr$replicate)
  }
  
  dereplicatedSites = granges(dereplicateSites(sites)) #instead of mcols()=NULL

  #regenerate posid which is guaranteed to match since the input data was
  #already run through the dereplicator
  dereplicatedSites$posid = paste0(seqnames(dereplicatedSites),
                                   strand(dereplicatedSites),
                                   start(flank(dereplicatedSites, -1, start=T)))
  
  dereplicatedSites$estAbund <- round(estimatedAbundances$theta) #estAbund preserves order
  dereplicatedSites$estAbundProp <- dereplicatedSites$estAbund/sum(dereplicatedSites$estAbund)
  dereplicatedSites$estAbundRank <- rank(-1*dereplicatedSites$estAbundProp, ties.method="max")
  
  dereplicatedSites
}
