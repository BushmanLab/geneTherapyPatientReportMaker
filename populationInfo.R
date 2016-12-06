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

calculateShannon <- function(dereplicated){
  shannonInput <- as.data.frame(mcols(dereplicated)[,c("posid", "Timepoint", "estAbund")])
  
  diversity(acast(shannonInput, Timepoint~posid, fill=0, value.var="estAbund",
                  fun.aggregate=sum))
}

jackIID <- function(specie,jrep=NULL,nrep=10L){
    if (is.null(jrep)) jrep <-
        sample(rep(1:nrep,length=length(specie)))
    est0 <- estimateR(table(specie))
    jackrep <- jrep
    urepl <- unique(jrep)
    jackmat <-
        sapply(urepl,
               function(x) estimateR(table(specie[jackrep!=x])))
    pseudo <- length(urepl)*est0 - (length(urepl)-1)*jackmat
    rowMeans(pseudo)
}

#' jackknife biased or unbiased Chao estimator
#'
#' @param replicatedSites df with column posid
#' @return number population size estimate
calculateChao <- function(replicatedSites, biased=TRUE){
    if ( ! biased) { #regular Chao
        cluster.tab <- table(replicatedSites$posid)
        return(round(estimateR(cluster.tab)["S.chao1"]))
    }
    round(jackIID(replicatedSites$posid)["S.chao1"])
}

calculateGini <- function(dereplicated){
  stopifnot(require(reldist))
  gini(dereplicated$estAbundProp)
}

getPopulationInfo <- function(replicated, dereplicated, splitBy){
  stopifnot(require(vegan))
  stopifnot((splitBy %in% names(mcols(replicated))) &
              (splitBy %in% names(mcols(dereplicated))))

  replicated <- split(replicated, mcols(replicated)[,splitBy])
  dereplicated <- split(dereplicated, mcols(dereplicated)[,splitBy])

  stopifnot(length(replicated) == length(dereplicated))

  populationInfo <- lapply(names(replicated), function(name){
    #can iterate through standardizedReplicatedSites and standardizedDereplicatedSites using GTSP#
    replicatedSites <- replicated[[name]]
    dereplicatedSites <- dereplicated[[name]]
    
    data.frame("group"=name,
               "S.chao1"=calculateChao(replicatedSites),
               "Gini"=calculateGini(dereplicatedSites),
               "Shannon"=calculateShannon(dereplicatedSites),
               "UC50"=calculateUC50(dereplicatedSites$estAbund))
  })

  populationInfo <- do.call(rbind, populationInfo)
  rownames(populationInfo) <- NULL
  populationInfo
}

calculateUC50 <- function(abund){
  stopifnot(is.vector(abund) & is.numeric(abund))
  abund <- abund[order(abund)]
  accum <- sapply(1:length(abund), function(i){sum(abund[1:i])})
  length(accum[accum >= sum(abund)/2])
}
