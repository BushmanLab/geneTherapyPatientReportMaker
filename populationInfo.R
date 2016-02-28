calculateShannon <- function(dereplicated){
  stopifnot(require(vegan))
  shannonInput <- as.data.frame(mcols(dereplicated)[,c("posid", "Timepoint", "estAbund")])
  
  diversity(acast(shannonInput, Timepoint~posid, fill=0, value.var="estAbund",
                  fun.aggregate=sum))
}

calculateChao <- function(replicatedSites, dereplicatedSites){
  stopifnot(require(vegan))
  chaoInput <- as.data.frame(merge(mcols(replicatedSites[,c("posid", "replicate")]),
                                          mcols(dereplicatedSites[,c("posid", "estAbund")])))
  
  counts <- acast(chaoInput, posid~replicate, fun.aggregate=sum,
                  value.var="estAbund", fill=0)

  if(all(dim(counts)>2) & is.matrix(counts)) {
    # do jackknife correction if more than 1 samples/replicates are available
    pseudo <- ncol(counts)*
      estimateR(rowSums(counts))-(ncol(counts)-1)*
      sapply(1:ncol(counts), function(y) estimateR(rowSums(counts[,-y])))
    round(rowMeans(pseudo)["S.chao1"])
  } else {
    pseudo <- estimateR(rowSums(counts))
    round(pseudo["S.chao1"])
  }  
}

calculateGini <- function(dereplicated){
  stopifnot(require(reldist))
  gini(dereplicated$estAbundProp)
}

getPopulationInfo <- function(replicated, dereplicated, splitBy){
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
               "S.chao1"=calculateChao(replicatedSites, dereplicatedSites), 
               "Gini"=calculateGini(dereplicatedSites),
               "Shannon"=calculateShannon(dereplicatedSites))
  })

  populationInfo <- do.call(rbind, populationInfo)
  rownames(populationInfo) <- NULL
  populationInfo
}
