#' Runs the genome scans of the simulated data output by sim.CC.data()
#'
#' This function takes the output from sim.CC.data() and performs the genome scans. Internally it runs
#' the require QR decompositions, which it can save for later if specified.
#'
#' @param sim.data Output simulated data from sim.CC.data()
#' @param scan.index DEFAULT: NULL. If NULL, it performs scans for all simulated data sets stored in sim.data.
#' If given a vector of integers, representing the simulation index, it will only run scans for those phenotypes.
#' @param chr DEFAULT: "all". Specifies which chromosomes to scan.
#' @param just.these.loci DEFAULT: NULL. If NULL, all loci in genome cache are scanned.
#' @param use.progress.bar DEFAULT: FALSE. Specifies whether to use a progress bar for qr decompositions 
#' and genome scans. 
#' @param all.sim.qr DEFAULT: NULL. If NULL, necessary qr decompositions are performed, which will slow down
#' the function significantly. If non-NULL, it expects either output from miqtl::extract.qr for a single sample
#' of CC lines, or a list of output for multiple samples.
#' @param return.all.sim.qr DEFAULT: TRUE. Saves the qr decompositions for later use, though will make output
#' substantially larger. Could then be used for later scans with the same set of CC lines.
#' @export
#' @examples run.sim.scans()
run.sim.scans <- function(sim.data,
                          scan.index=NULL, 
                          chr="all", 
                          just.these.loci=NULL,
                          use.progress.bar=FALSE,
                          all.sim.qr=NULL,
                          return.all.sim.qr=TRUE,
                          ...){
  
  if(is.null(scan.index)){ 
    scan.index <- 1:length(sim.data$data)
  }
  h <- miqtl::DiploprobReader$new(sim.data$genomecache)
  loci <- h$getLoci()
  loci.chr <- h$getChromOfLocus(loci)
  if(chr[1] != "all"){
    loci <- loci[loci.chr %in% chr]
  }
  if(!is.null(just.these.loci)){
    loci <- loci[loci %in% just.these.loci]
    loci.chr <- loci.chr[loci %in% just.these.loci]
  }
  
  full.p <- matrix(NA, nrow=length(scan.index), ncol=length(loci))
  colnames(full.p) <- loci
  these.pos <- list(Mb=h$getLocusStart(loci, scale="Mb"),
                    cM=h$getLocusStart(loci, scale="cM"))
  
  if(sim.data$properties$vary.lines){
    all.sim.qr <- NULL
  }
  else{
    if(is.null(all.sim.qr)){
      this.qr <- miqtl::extract.qr(genomecache=sim.data$genomecache, model="additive",
                                   formula=~1, data=sim.data$data, 
                                   id="SUBJECT.NAME.1",
                                   just.these.loci=just.these.loci, chr=chr,
                                   use.progress.bar=use.progress.bar)
    }
    else{
      this.qr <- all.sim.qr
    }
  }
  
  if(return.all.sim.qr & is.null(all.sim.qr)){ all.sim.qr <- list() }
  for(i in 1:length(scan.index)){
    if(sim.data$properties$vary.lines){
      this.qr <- miqtl::extract.qr(genomecache=sim.data$genomecache, model="additive",
                                   formula=~1, data=sim.data$data,
                                   id=paste0("SUBJECT.NAME.", i),
                                   just.these.loci=just.these.loci, chr=chr,
                                   use.progress.bar=use.progress.bar)
      this.id <- paste0("SUBJECT.NAME.", i)
    }
    else{ this.id <- "SUBJECT.NAME.1" }
    this.scan <- miqtl::scan.qr(qr.object=this.qr, data=sim.data$data, 
                                phenotype=paste0("sim.y.", i), id=this.id, chr=chr, 
                                return.allele.effects=FALSE, use.progress.bar=use.progress.bar,
                                ...)
    full.p[i,] <- this.scan$p.value
    if(return.all.sim.qr & sim.data$properties$vary.lines){ all.sim.qr[[i]] <- this.qr }
    cat("\n", "Simulation scan: index", scan.index[i], "complete ---------- final index of this run:", scan.index[length(scan.index)], "\n")
  }
  output <- list(LOD=NULL,
                 p.value=full.p,
                 chr=loci.chr, 
                 pos=these.pos,
                 locus=sim.data$locus,
                 locs.pos=sim.data$locus.pos,
                 properties=sim.data$properties)
  if(return.all.sim.qr){
    output$all.sim.qr <- all.sim.qr
  }
  return(output)
}

#' Simulate Collaborative Cross (CC) phenotype data from the realized CC lines founder haplotype probabilities stored in a genome cache directory
#'
#' This function takes various input parameters to simulate CC data to be used in power calculations or as input data for 
#' other tools that analyze CC data.
#'
#' @param genomecache The path to the genome cache directory. The genome cache is a particularly structured
#' directory that stores the haplotype probabilities/dosages at each locus. It has an additive model
#' subdirectory and a full model subdirectory. Each contains subdirectories for each chromosome, which then
#' store .RData files for the probabilities/dosages of each locus.
#' @param CC.lines DEFAULT: NULL. If NULL is specified, sim.CC.data() will randomly draw samples of the 
#' available CC lines of the size specified in num.lines.
#' @param num.lines DEFAULT: NULL. If NULL, sim.CC.data() expects that CC.lines is non-NULL. If CC.lines
#' is NULL, num.lines determines the number of CC lines that are sampled.
#' @param vary.lines DEFAULT: TRUE. If CC.lines is NULL and vary.lines is TRUE, then sim.CC.data() will
#' sample multiple sets of CC lines of the size specified in num.lines. If CC.lines is NULL and vary.lines is
#' FALSE, then just one set of CC lines will be sampled.
#' @param locus DEFAULT: NULL. If NULL is specified, sim.CC.data() will randomly draw a locus stored in
#' the genome cache. 
#' @param vary.locus DEFAULT: TRUE. If locus is NULL and vary.locus is TRUE, then sim.CC.data() will sample
#' as many loci as specified in num.sim. If locus is NULL and vary.locus is FALSE, then sim.CC.data() will
#' only sample one locus.
#' @param num.replicates The number of replicates of each CC line that will be simulated. Mapping for power calculations
#' will use strain means. Currently requires that all lines have the same number of replicates.
#' @param num.sim The number of phenotypes to be simulated for a given parameter setting.
#' @param M.ID DEFAULT: NULL. M is a matrix that maps from counts of founder haplotypes to counts of functional 
#' alleles. Mapping will be based on haplotype association, but potentially there are two to the number of
#' founders alleles at the QTL. M.ID is a string that codifies this mapping. One potential balanced two allele
#' M.ID would be "c(0,0,0,0,1,1,1,1)". With 8 functional alleles, on per founder, the only M.ID is "c(0,1,2,3,4,5,6,7)".
#' If M.ID is NULL, M.ID will be sampled.
#' @param num.alleles DEFAULT: 8. The number of functional alleles. Must be less than or equal to the number of 
#' founders.
#' @param num.founders DEFAULT: 8. The number of founders, which must correspond to the genome cache. The CC has
#' eight.
#' @param qtl.effect.size The size of the simulated QTL effect. The scale of the input is in proportion of
#' the phenotypic variance due to the QTL, thus should be greater than or equal to zero, and less than one.
#' @param beta DEFAULT: NULL. Allows for the manual specification of QTL effect. Is expected to be a vector
#' the length of the number of alleles. It will be scaled based on qtl.effect.size.
#' @param strain.effect.size The size of the simulated strain effect, which represents something akin to a polygenic 
#' effect. Other variants specific to CC lines will result in overall strain-specific effects. The scale of the input 
#' is in proportion of the phenotypic variance due to the strain, thus should be greater than or equal to zero, and less than one.
#' @param impute DEFAULT: TRUE. If TRUE, the QTL portion of the design matrix in the simulation is a realized sampling
#' of haplotypes from the probabilities. If FALSE, the simulations are based on the probabilities, which is flawed in
#' terms of biological reality.
#' @param scale.by.effect DEFAULT: FALSE. If TRUE, the effects are scaled by var(QTL effect) so that the effect matches
#' the stated proportion of variance in a balanced population.
#' @param scale.by.varp DEFAULT: TRUE. If TRUE, the effects are scaled by var(X %*% QTL effect) so that the effect matches
#' the stated proportion of variance in the observed population.
#' @export
#' @examples sim.CC.data()
sim.CC.data <- function(genomecache, 
                        CC.lines=NULL, 
                        num.lines=NULL, 
                        vary.lines=TRUE,
                        locus=NULL, 
                        vary.locus=TRUE,
                        num.replicates, 
                        num.sim,
                        M.ID=NULL,
                        num.alleles=8, 
                        num.founders=8,
                        qtl.effect.size,
                        beta=NULL,
                        strain.effect.size,
                        impute=TRUE,
                        scale.by.effect=FALSE,
                        scale.by.varp=TRUE){
  
  h <- miqtl::DiploprobReader$new(genomecache)
  
  ## Sampling lines
  if(is.null(CC.lines)){
    CC.lines <- sapply(1:ifelse(vary.lines, num.sim, 1), function(i) sample(x=h$getSubjects(), size=num.lines, replace=FALSE))
    if(vary.lines){
      cc.index <- 1:num.sim
    }
    else{
      cc.index <- rep(1, num.sim)
    }
  }
  else{
    vary.lines <- FALSE
    num.lines <- length(CC.lines)
    CC.lines <- matrix(CC.lines, ncol=1)
    cc.index <- rep(1, num.sim)
  }
  
  ## Setting number of alleles to match pre-specified QTL effect - convenience
  if(!is.null(beta)){
    num.alleles <- length(beta)
  }
  
  ## Sampling loci
  if(is.null(locus)){
    locus <- sample(h$getLoci(), size=ifelse(vary.locus, num.sim, 1), replace=TRUE)
    if(vary.locus){
      locus.index <- 1:num.sim
    }
    else{
      locus.index <- rep(1, num.sim)
    }
  }
  else{
    vary.locus <- FALSE
    locus.index <- rep(1, num.sim)
  }
  
  M <- NULL
  if(!is.null(M.ID)){
    M <- model.matrix.from.ID(M.ID)
    num.alleles <- length(unique(unlist(strsplit(x=M.ID, split=","))))
  }
  sim.matrix <- matrix(NA, nrow=nrow(CC.lines), ncol=num.sim)
  id.matrix <- matrix(NA, nrow=nrow(CC.lines), ncol=ifelse(vary.lines, num.sim, 1))
  for(i in 1:num.sim){
    this.sim <- simulate.CC.qtl(CC.lines=CC.lines[,cc.index[i]], 
                                num.replicates=num.replicates,
                                M=M,
                                num.alleles=num.alleles, 
                                num.founders=num.founders,
                                qtl.effect.size=qtl.effect.size, 
                                strain.effect.size=strain.effect.size,
                                impute=impute, 
                                scale.by.varp=scale.by.varp,
                                locus.matrix=h$getLocusMatrix(locus=locus[locus.index[i]], model="full"), 
                                num.sim=1)$data
    sim.matrix[,i] <- this.sim[,1]
    if(vary.lines | i == 1){
      id.matrix[,i] <- as.character(this.sim[,2])
    }
  }
  colnames(sim.matrix) <- paste0("sim.y.", 1:num.sim)
  colnames(id.matrix) <- paste0("SUBJECT.NAME.", 1:ncol(id.matrix))
  
  outcome <- data.frame(sim.matrix, id.matrix)
  
  return(list(data=outcome,
              locus=locus,
              locus.pos=list(Mb=h$getLocusStart(locus, scale="Mb"),
                             cM=h$getLocusStart(locus, scale="cM")),
              genomecache=genomecache,
              properties=list(num.alleles=num.alleles,
                              num.replicates=num.replicates,
                              num.founders=num.founders,
                              qtl.effect.size=qtl.effect.size, 
                              strain.effect.size=strain.effect.size,
                              num.lines=num.lines,
                              impute=impute,
                              scale.by.varp=scale.by.varp,
                              M.ID=M.ID,
                              vary.lines=vary.lines,
                              vary.locus=vary.locus)))

simulate.CC.qtl <- function(CC.lines, 
                            num.replicates,
                            M=NULL,
                            qtl.effect.size, 
                            beta=NULL,
                            strain.effect.size,
                            num.alleles=8, 
                            num.founders=8,
                            locus.matrix, 
                            num.sim,
                            impute=TRUE, 
                            scale.by.varp=TRUE,
                            scale.by.effect=FALSE,
                            ...){
  
  this.locus.matrix <- locus.matrix[CC.lines,]
  this.locus.matrix <- this.locus.matrix[rep(1:nrow(this.locus.matrix), each=num.replicates),]
  
  full.to.add.matrix <- t(miqtl::straineff.mapping.matrix(M=num.founders))
  
  ## QTL
  if(qtl.effect.size != 0){
    QTL.effect <- simulate.QTL.model.and.effects(num.alleles=num.alleles, 
                                                 num.founders=num.founders, 
                                                 M=M,
                                                 effect.var=qtl.effect.size, 
                                                 beta=beta,
                                                 scale.by.effect=scale.by.effect,
                                                 ...)
    if(impute){
      this.locus.matrix <- t(apply(this.locus.matrix, 1, function(x) rmultinom(1, 1, x)))
    }
    this.locus.matrix <- tcrossprod(this.locus.matrix, full.to.add.matrix)
    QTL.predictor <- tcrossprod(tcrossprod(this.locus.matrix, t(QTL.effect$M)), matrix(QTL.effect$beta, nrow=1))
    if(scale.by.varp){
      QTL.predictor <- as.vector(scale(QTL.predictor))
      QTL.predictor <- as.vector(QTL.predictor*sqrt(qtl.effect.size))
    }
  }
  else{ QTL.predictor <- rep(0, length(CC.lines))  }
  
  ## Strain
  if(strain.effect.size != 0){
    this.strain.matrix <- incidence.matrix(factor(CC.lines, levels=CC.lines))
    this.strain.matrix <- this.strain.matrix[rep(1:nrow(this.strain.matrix), each=num.replicates),]
    strain.effect <- rnorm(mean=0, sd=sqrt(strain.effect.size), n=length(CC.lines))
    strain.effect <- as.vector(scale(strain.effect))
    if(scale.by.effect){
      strain.effect <- as.vector(strain.effect*sqrt(strain.effect.size)) 
    }
    strain.predictor <- this.strain.matrix %*% matrix(strain.effect, ncol=1)

    if(scale.by.varp){
      strain.predictor <- as.vector(scale(strain.predictor))
      strain.predictor <- as.vector(strain.predictor*sqrt(strain.effect.size))
    }
  }
  else{ strain.predictor <- rep(0, length(CC.lines)) }

  
  sim.outcome <- sapply(1:num.sim, 
                        function(i) QTL.predictor + strain.predictor + calc.scaled.residual(qtl.effect.size=qtl.effect.size, 
                                                                                            strain.effect.size=strain.effect.size, 
                                                                                            n=nrow(this.locus.matrix)))
  mean.sim.outcome <- apply(sim.outcome, 2, function(x) 
    tapply(x, INDEX=as.factor(rep(CC.lines, each=num.replicates)), FUN=mean))
  colnames(mean.sim.outcome) <- paste0("sim.y.", 1:ncol(mean.sim.outcome))
  mean.sim.outcome <- data.frame(mean.sim.outcome, "SUBJECT.NAME"=rownames(mean.sim.outcome))
  
  return(list(data=mean.sim.outcome,
              properties=list(qtl.effect.size=qtl.effect.size,
                              strain.effect.size=strain.effect.size,
                              num.alleles=num.alleles,
                              num.replicates=num.replicates,
                              num.lines=length(CC.lines),
                              impute=impute)))
}

## Draws and scales residuals in single function
calc.scaled.residual <- function(qtl.effect.size, 
                                 strain.effect.size, 
                                 n){

  residual <- rnorm(mean=0, sd=sqrt(1 - qtl.effect.size - strain.effect.size), n=n)
  residual <- as.vector(scale(residual))
  residual <- as.vector(residual*sqrt(1 - qtl.effect.size - strain.effect.size))
  return(residual)
}

## From Wes. Returns SDP matrix and QTL effects
simulate.QTL.model.and.effects <- function(num.alleles=8, 
                                           num.founders=8, 
                                           M=NULL,
                                           model.method="crp", 
                                           effect.type="random", 
                                           effect.var, 
                                           beta=NULL,
                                           scale.by.effect=TRUE,
                                           ...){
  
  if(is.null(M)){
    M <- matrix(0, num.founders, num.alleles)
    M[cbind(sample(1:num.founders, num.alleles), 1:num.alleles)] <- 1
    
    if(model.method == "crp"){
      for(i in which(rowSums(M)==0)){
        M[i,] <- rmultinom(1, 1, colSums(M)/sum(M))
      }
    }
    else if(model.method == "uniform"){
      M[which(rowSums(M)==0),] <- t(rmultinom(num.founders - num.alleles, 1, rep(1/num.alleles, num.alleles)))
    }
  }
  
  if(is.null(beta)){
    if(effect.type == "random"){
      beta <- rnorm(num.alleles)
    } 
    else if(effect.type == "even"){
      beta <- 1:num.alleles
    }  
  }
  
   beta <- as.vector(scale(beta))
  if(scale.by.effect){
    beta <- as.vector(0.5*beta*sqrt(effect.var)) 
  }

  effect <- list(M=M, 
                 beta=beta)
  return(effect)
}

#' Plots a single simulated scan. Built on miqtl::genome.plotter.whole()
#'
#' This function takes the output from run.sim.scans(), and plots a single genome scane of a simulated phenotype.
#'
#' @param sim.scans Genome scans of simulated data, output by run.sim.scans().
#' @param phenotype.index DEFAULT: 1. The index of the simulated phenotype to plot. Cannot exceed the number of 
#' phenotypes contained in sim.scans.
#' @param scale DEFAULT: "Mb". The scale in which to plot genomic position. "cM" is the other option.
#' @param thresh DEFAULT: NULL. If NULL, no significance is threshold is plotted. Expects a numeric value, commonly
#' output from get.gev.thresholds() based on permutation scans.
#' @export
#' @examples single.sim.plot()
single.sim.plot <- function(sim.scans, 
                            phenotype.index=1, 
                            scale=c("Mb", "cM"), 
                            thresh=NULL, 
                            ...){
  
  scale <- scale[1]
  dummy.scan <- sim.scans
  if(is.matrix(dummy.scan$p.value)){
    dummy.scan$p.value <- sim.scans$p.value[phenotype.index,]
  }
  else{
    dummy.scan$p.value <- sim.scans$p.value
  }
  if(length(sim.scans$locus) > 1){
    locus <- sim.scans$locus[phenotype.index]
  }
  else{
    locus <- sim.scans$locus
  }
  
  this.title <- c(paste0("QTL Effect=",  sim.scans$properties$qtl.effect.size*100, "%"), 
                  paste0("Strain Effect=", sim.scans$properties$strain.effect.size*100, "%"),
                  paste(ifelse(sim.scans$properties$impute, "Impute", "ROP"), "sim,",
                        sim.scans$properties$num.alleles, "functional alleles,",
                        sim.scans$properties$num.lines, "lines,",
                        sim.scans$properties$num.replicates, ifelse(sim.scans$properties$num.replicates == 1, "replicate", "replicates")))
  genome.plotter.whole(scan.list=list(dummy.scan), override.title=this.title, distinguish.chr.type="color",
                       scale=scale, mark.locus=locus, hard.thresholds=thresh, ...)
}

## Very simply checks what proportion of simulations correctly called the QTL. Needs thresholds
#' @export
pull.power <- function(sim.scans, locus, thresholds){
  p.values <- sim.scans$p.value[,colnames(sim.scans$p.value) == locus]
  power <- mean(-log10(p.values) > thresholds)
  return(power)
}

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

## Produce SDP mapping matrix from SDP string
model.matrix.from.ID <- function(M.ID){
  m <- as.numeric(unlist(strsplit(M.ID, ","))) + 1
  J <- length(m)
  M <- matrix(0, J, max(m))
  M[cbind(1:J, m)] <- 1
  return(M)
}

## Produces an incidence matrix from a factor
incidence.matrix <- function(fact){
  m <- diag(nlevels(fact))[fact,]
  colnames(m) <- levels(fact)
  return(m)
}
  
  

