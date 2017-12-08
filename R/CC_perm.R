#' Generates a matrix of permutation indeces. Can be used across simulations with the same number of CC lines
#'
#' This function takes the output from sim.CC.data() to produce a matrix of permutation indeces, allowing
#' the same permutations to be performed for different phenotypes but the same CC lines.
#'
#' @param sim.CC.object Simulated CC data output from sim.CC.data().
#' @param num.perm The number of permutations.
#' @param seed DEFAULT: 1. A seed is necessary to produce the same results over multiple runs and different machines.
#' @export
#' @examples generate.permutation.index.matrix()
generate.permutation.index.matrix <- function(sim.CC.object, 
                                              num.perm, 
                                              seed=1){
  
  n <- nrow(sim.CC.object$data)
  
  set.seed(seed)
  perm.index.matrix <- replicate(n=num.perm, sample(1:n, replace=FALSE))
  colnames(perm.index.matrix) <- paste0("perm.", 1:num.perm)
  return(perm.index.matrix)
}

#' Runs permutation scans from a permutation index matrix, simulated CC data, and simulated CC scans.
#'
#' This function takes the outputs from generate.permutation.index.matrix(), sim.CC.data(), and run.sim.scans() to perform
#' permutation scans and determine a significance threshold. 
#'
#' @param perm.index.matrix Permutation index matrix that is output from generate.permutation.index.matrix().
#' @param sim.CC.scans Genome scans from simulated CC data output from run.sim.scans().
#' @param sim.CC.object Simulated CC data output from sim.CC.data().
#' @param phenotype.index The phenotype index of simulated phenotype that correspond to those found in sim.CC.scans
#' and sim.CC.objects.
#' @param all.sim.qr DEFAULT: NULL. Allows qr decompositions to only be saved once. If NULL, it expects that they are
#' stored in sim.CC.scans$all.sim.qr, which can be specified in run.sim.scans().
#' @param keep.full.scans DEFAULT: TRUE. If TRUE, the full scans for the permutations of a phenotype are kept. If FALSE,
#' only the minimum p-values are kept, saving space.
#' @param scan.index DEFAULT: NULL. If NULL, it performs scans for all permuted data sets based on permutation index.
#' If given a vector of integers, representing the permutation index, it will only run scans for those permutations.
#' @param chr DEFAULT: "all". Specifies which chromosomes to scan.
#' @param just.these.loci DEFAULT: NULL. If NULL, all loci in genome cache are scanned.
#' @param use.progress.bar DEFAULT: FALSE. Specifies whether to use a progress bar for qr decompositions 
#' and genome scans. 
#' @export
#' @examples run.permutation.threshold.scans()
run.permutation.threshold.scans <- function(perm.index.matrix, 
                                            sim.CC.scans, 
                                            sim.CC.object,
                                            phenotype.index,
                                            all.sim.qr=NULL,
                                            keep.full.scans=TRUE, 
                                            scan.index=NULL,
                                            chr="all", 
                                            just.these.loci=NULL, 
                                            use.progress.bar=FALSE,
                                            ...){
  
  if(is.null(scan.index)){ scan.index <- 1:ncol(perm.index.matrix) }
  
  if(is.null(all.sim.qr)){
    all.sim.qr <- sim.CC.scans$all.sim.qr
  }
  if(sim.CC.scans$properties$vary.lines){
    loci <- names(all.sim.qr[[1]]$qr.list)
    rh.formula <- all.sim.qr[[1]]$formula
    loci.chr <- all.sim.qr[[1]]$chr
  }
  else{
    loci <- names(all.sim.qr$qr.list)
    rh.formula <- all.sim.qr$formula
    loci.chr <- all.sim.qr$chr
  }
  
  if(chr[1] != "all"){
    loci <- loci[loci.chr %in% chr]
  }
  if(!is.null(just.these.loci)){
    loci <- loci[loci %in% just.these.loci]
    loci.chr <- loci.chr[loci %in% just.these.loci]
  }
  data <- sim.CC.object$data[,c(paste0("sim.y.", phenotype.index), 
                                paste0("SUBJECT.NAME.", ifelse(sim.CC.scans$properties$vary.lines, phenotype.index, 1)))]
  
  full.p <- these.pos <- NULL
  if(keep.full.scans){
    full.p <- matrix(NA, nrow=length(scan.index), ncol=length(loci))
    colnames(full.p) <- loci
    if(sim.CC.scans$properties$vary.lines){
      these.pos <- list(Mb=all.sim.qr[[1]]$pos$Mb[loci],
                        cM=all.sim.qr[[1]]$pos$cM[loci])
    }
    else{
      these.pos <- list(Mb=all.sim.qr$pos$Mb[loci],
                        cM=all.sim.qr$pos$cM[loci])
    }
  }
  min.p <- rep(NA, length(scan.index))
  
  if(sim.CC.scans$properties$vary.lines){
    this.qr <- all.sim.qr[[phenotype.index]]
  }
  else{
    this.qr <- all.sim.qr
  }
  for(i in 1:length(scan.index)){
    perm.data <- data
    perm.data[,1] <- data[perm.index.matrix[,scan.index[i]], 1]
    names(perm.data) <- c("perm_y", "SUBJECT.NAME")
    
    this.scan <- miqtl::scan.qr(qr.object=this.qr, data=perm.data, 
                                phenotype="perm_y", chr=chr,
                                return.allele.effects=FALSE, use.progress.bar=use.progress.bar,
                                ...)
    if(keep.full.scans){
      full.p[i,] <- this.scan$p.value
    }
    min.p[i] <-  min(this.scan$p.value)
    cat("\n", "Threshold scan: index", scan.index[i], "complete ---------- final index of this run:", scan.index[length(scan.index)], "\n")
  }
  return(list(full.results=list(LOD=NULL,
                                p.value=full.p,
                                chr=loci.chr, 
                                pos=these.pos), 
              max.statistics=list(LOD=NULL,
                                  p.value=min.p)))
}