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
#' @param print.scans.progress DEFAULT: FALSE. Specifies whether to output a message that indicates the number of
#' scans completed.
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
                          print.scans.progress=FALSE,
                          all.sim.qr=NULL,
                          return.all.sim.qr=TRUE,
                          ...){
  
  if(is.null(scan.index)){ 
    scan.index <- 1:sum(grepl("sim.y", colnames(sim.data$data)))
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
  
  if(!sim.data$properties$vary.lines){
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
    this.id <- "SUBJECT.NAME.1"
  }
  compute.qr.boolean <- is.null(all.sim.qr)
  if(return.all.sim.qr & is.null(all.sim.qr)){ all.sim.qr <- list() }
  for(i in 1:length(scan.index)){
    if(sim.data$properties$vary.lines){
      if(compute.qr.boolean){
        this.qr <- miqtl::extract.qr(genomecache=sim.data$genomecache, model="additive",
                                     formula=~1, data=sim.data$data,
                                     id=paste0("SUBJECT.NAME.", i),
                                     just.these.loci=just.these.loci, chr=chr,
                                     use.progress.bar=use.progress.bar)
      }
      else{
        this.qr <- pull.qr.from.compact(compact.qr.list=all.sim.qr, 
                                        qr.index=scan.index[i])
      }
      this.id <- paste0("SUBJECT.NAME.", i)
    }
    this.scan <- miqtl::scan.qr(qr.object=this.qr, data=sim.data$data, 
                                phenotype=paste0("sim.y.", i), id=this.id, chr=chr, 
                                return.allele.effects=FALSE, use.progress.bar=use.progress.bar,
                                ...)
    full.p[i,] <- this.scan$p.value
    if(return.all.sim.qr & sim.data$properties$vary.lines){ all.sim.qr[[i]] <- this.qr }
    if(print.scans.progress){
      cat("\n", "Simulation scan: index", scan.index[i], "complete ---------- final index of this run:", scan.index[length(scan.index)], "\n")
    }
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
#' @param sample.sadp.method DEFAULT: "uniform". The procedure used for sampling the strain allele distribution pattern. If
#' every strain has its own allele, this option does not matter. Alternatively, a Chinese restaurant process ("crp") can be used,
#' which is possibly more biologically accurate, and will favor strain allele distribution patterns that are less balanced (1 vs 7).
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
#' @param strain.reduction.fraction DEFAULT: NULL. This argument allows the strain effect to be reduced by a specified fraction.
#' The DEFAULT of NULL for return.value="raw" and return.means=TRUE sets strain.reduction.fraction to (num.replicates - 1)/num.replicates, which
#' corresponds to the expected reduction of the strain effect if we could estimate the BLUE of the QTL. A value of 1 reflects the strain effect being 
#' perfectly controlled, and thus removed. The parameterization is such that strain.effect.size = strain.effect.size * (1 - strain.reduction.fraction). 
#' Setting strain.reduction.fraction to 0 reflects strain not being modeled at all, and simply being absorbed into the residual. Argument is forced to NULL
#' for residuals or return.means=FALSE.
#' @param impute DEFAULT: TRUE. If TRUE, the QTL portion of the design matrix in the simulation is a realized sampling
#' of haplotypes from the probabilities. If FALSE, the simulations are based on the probabilities, which is flawed in
#' terms of biological reality.
#' @param scale.by.varp DEFAULT: TRUE. If TRUE, the effects are scaled by var(X %*% QTL effect) so that the effect matches
#' the stated proportion of variance in the observed population.
#' @param return.value DEFAULT: "raw". If "raw", residuals are not taken. If "fixef.resid", then the data
#' are residuals after regressing phenotype on strain. If "ranef.resid", then the data have had the strain BLUP
#' effect subtracted.
#' @param return.means DEFAULT: TRUE. If TRUE, strain means are returned. If FALSE, the full data
#' with replicate observations of strains are returned.
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
                        sample.sadp.method=c("uniform", "crp"),
                        num.alleles=8, 
                        num.founders=8,
                        qtl.effect.size,
                        beta=NULL,
                        strain.effect.size,
                        strain.reduction.fraction=NULL,
                        impute=TRUE,
                        scale.by.varp=FALSE,
                        return.value=c("raw", "fixef.resid", "ranef.resid"),
                        return.means=TRUE){
  
  h <- miqtl::DiploprobReader$new(genomecache)
  return.value <- return.value[1]
  sample.sdp.method <- sample.sdp.method[1]
  
  original.effects <- list(qtl.effect.size=qtl.effect.size,
                           strain.effect.size=strain.effect.size)
  
  if (return.value == "raw" & return.means) {
    if (is.null(strain.reduction.fraction)) { strain.reduction.fraction <- (num.replicates - 1)/num.replicates }
    qtl.effect.size <- qtl.effect.size/(1 - strain.effect.size*strain.reduction.fraction)
    strain.effect.size <- strain.effect.size*(1 - strain.reduction.fraction)
  }
  else { strain.reduction.fraction <- NULL }

  ## Sampling lines
  if (is.null(CC.lines)) {
    CC.lines <- sapply(1:ifelse(vary.lines, num.sim, 1), function(i) sample(x=h$getSubjects(), size=num.lines, replace=FALSE))
    if (vary.lines) {
      cc.index <- 1:num.sim
    }
    else {
      cc.index <- rep(1, num.sim)
    }
  }
  else{
    if(is.vector(CC.lines)){
      CC.lines <- matrix(CC.lines, ncol=1)
      vary.lines <- FALSE
      num.lines <- length(CC.lines)
      cc.index <- rep(1, num.sim)
    }
    else{
      vary.lines <- TRUE
      num.lines <- nrow(CC.lines)
      num.sim <- ncol(CC.lines)
      cc.index <- 1:num.sim
    }
  }
  num.ind <- ifelse(return.means, num.lines, num.lines*num.replicates)
  
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
  sim.matrix <- matrix(NA, nrow=num.ind, ncol=num.sim)
  id.matrix <- matrix(NA, nrow=num.ind, ncol=ifelse(vary.lines, num.sim, 1))
  for(i in 1:num.sim){
    this.sim <- simulate.CC.qtl(CC.lines=CC.lines[,cc.index[i]], 
                                num.replicates=num.replicates,
                                M=M,
                                sample.sadp.method=sample.sadp.method,
                                num.alleles=num.alleles, 
                                num.founders=num.founders,
                                qtl.effect.size=qtl.effect.size, 
                                strain.effect.size=strain.effect.size,
                                impute=impute,
                                scale.by.varp=scale.by.varp,
                                locus.matrix=h$getLocusMatrix(locus=locus[locus.index[i]], model="full"),
                                return.value=return.value,
                                return.means=return.means,
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
                              sample.sadp.method=sample.sadp.method,
                              num.replicates=num.replicates,
                              num.founders=num.founders,
                              qtl.effect.size=original.effects$qtl.effect.size, 
                              strain.effect.size=original.effects$strain.effect.size,
                              strain.reduction.fraction=strain.reduction.fraction,
                              num.lines=num.lines,
                              impute=impute,
                              scale.by.varp=scale.by.varp,
                              M.ID=M.ID,
                              vary.lines=vary.lines,
                              vary.locus=vary.locus,
                              return.value=return.value,
                              return.means=return.means)))
}

simulate.CC.qtl <- function(CC.lines, 
                            num.replicates,
                            M=NULL,
                            sample.sadp.method=c("uniform", "crp"),
                            qtl.effect.size, 
                            beta=NULL,
                            strain.effect.size,
                            num.alleles=8, 
                            num.founders=8,
                            locus.matrix, 
                            num.sim,
                            impute=TRUE, 
                            scale.by.varp=FALSE,
                            return.value=c("raw", "fixef.resid", "ranef.resid"),
                            return.means=TRUE,
                            ...){
  
  return.value <- return.value[1]
  sample.sadp.method <- sample.sadp.method[1]
  this.locus.matrix <- locus.matrix[CC.lines,]
  this.locus.matrix <- this.locus.matrix[rep(1:nrow(this.locus.matrix), each=num.replicates),]
  
  full.to.add.matrix <- t(miqtl::straineff.mapping.matrix(M=num.founders))
  
  ## QTL
  if(qtl.effect.size != 0){
    QTL.effect <- simulate.QTL.model.and.effects(num.alleles=num.alleles, 
                                                 num.founders=num.founders, 
                                                 M=M,
                                                 sample.sadp.method=sample.sadp.method,
                                                 effect.var=qtl.effect.size, 
                                                 beta=beta,
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
    
    strain.effect <- rnorm(n=length(CC.lines))
    strain.effect <- as.vector(scale(strain.effect))
    strain.effect <- as.vector(strain.effect*sqrt(strain.effect.size)) 

    strain.predictor <- this.strain.matrix %*% matrix(strain.effect, ncol=1)

    if(scale.by.varp){
      strain.predictor <- as.vector(scale(strain.predictor))
      strain.predictor <- as.vector(strain.predictor*sqrt(strain.effect.size))
    }
  }
  else{ strain.predictor <- rep(0, length(CC.lines)) }
  
  sim.data <- sapply(1:num.sim, function(i) QTL.predictor + strain.predictor + calc.scaled.residual(qtl.effect.size=qtl.effect.size, 
                                                                                                    strain.effect.size=strain.effect.size,
                                                                                                    n=nrow(this.locus.matrix)))
  colnames(sim.data) <- paste0("sim.y.", 1:ncol(sim.data))
  sim.data <- data.frame(sim.data, "SUBJECT.NAME"=rep(CC.lines, each=num.replicates))
  
  ## Taking residuals
  if (return.value == "fixef.resid") {
    sim.data$sim.y.1 <- take.fixef.residuals(data=sim.data, 
                                             locus.matrix=this.locus.matrix, 
                                             strain.matrix=this.strain.matrix)
  }
  else if (return.value == "ranef.resid") {
    sim.data$sim.y.1 <- take.ranef.residuals(data=sim.data,
                                             locus.matrix=this.locus.matrix)
  }
  if(return.means){
    sim.data <- apply(sim.data[,-ncol(sim.data), drop=FALSE], 2, function(x) tapply(x, INDEX=as.factor(rep(CC.lines, each=num.replicates)), FUN=mean))
    colnames(sim.data) <- paste0("sim.y.", 1:ncol(sim.data))
    sim.data <- data.frame(sim.data, "SUBJECT.NAME"=rownames(sim.data))
  }
  return(list(data=sim.data,
              properties=list(qtl.effect.size=qtl.effect.size,
                              strain.effect.size=strain.effect.size,
                              num.alleles=num.alleles,
                              sample.sadp.method=sample.sadp.method,
                              num.replicates=num.replicates,
                              num.lines=length(CC.lines),
                              impute=impute,
                              return.value=return.value,
                              return.means=return.means)))
}

take.fixef.residuals <- function(data, locus.matrix, strain.matrix){
  locus.matrix <- locus.matrix[,-which.max(colSums(locus.matrix))]
  rownames(locus.matrix) <- NULL
  colnames(locus.matrix) <- paste("allele", 1:ncol(locus.matrix), sep=".")
  
  strain.matrix <- strain.matrix[,-1]
  
  temp.data <- data.frame(data, locus.matrix, strain.matrix)
  this.formula <- formula(paste("sim.y.1~1", 
                               paste(paste("allele", 1:ncol(locus.matrix), sep="."), collapse="+"),
                               paste(colnames(strain.matrix), collapse="+"),
                               sep="+"))
  this.fit <- lm(this.formula, data=temp.data)
  
  drop.na <- is.na(this.fit$coefficients[colnames(strain.matrix)])
  
  these.residuals <- data$sim.y.1 - strain.matrix[,!drop.na] %*% this.fit$coefficients[colnames(strain.matrix)][!drop.na]
  return(these.residuals)
}

take.ranef.residuals <- function(data, locus.matrix){
  locus.matrix <- locus.matrix[,-which.max(colSums(locus.matrix))]
  rownames(locus.matrix) <- NULL
  colnames(locus.matrix) <- paste("allele", 1:ncol(locus.matrix), sep=".")
  
  temp.data <- data.frame(data, locus.matrix)
  this.formula <- formula(paste0("sim.y.1~", 
                                paste(paste("allele", 1:ncol(locus.matrix), sep="."), collapse="+"),
                                "+(1|SUBJECT.NAME)"))

  this.fit <- lme4::lmer(this.formula, data=temp.data)
  these.residuals <- data$sim.y.1 - lme4::ranef(this.fit)$SUBJECT.NAME[data$SUBJECT.NAME,]
  return(these.residuals)
}

## Draws and scales residuals in single function
calc.scaled.residual <- function(qtl.effect.size, 
                                 strain.effect.size, 
                                 n){

  residual <- rnorm(n=n)
  residual <- as.vector(scale(residual))
  residual <- as.vector(residual*sqrt(1 - qtl.effect.size - strain.effect.size))
  return(residual)
}

## From Wes. Returns SDP matrix and QTL effects
simulate.QTL.model.and.effects <- function(num.alleles=8, 
                                           num.founders=8, 
                                           M=NULL,
                                           sample.sadp.method=c("uniform", "crp"), 
                                           effect.var, 
                                           beta=NULL,
                                           ...){
  
  model.method <- model.method[1]
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
    beta <- 1:num.alleles
  }
  beta <- as.vector(scale(beta))
  beta <- as.vector(0.5*beta*sqrt(effect.var)) 
  
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
  dummy.scan$p.value[dummy.scan$p.value == 0] <- .Machine$double.xmin
  
  if(length(sim.scans$locus) > 1){
    locus <- sim.scans$locus[phenotype.index]
  }
  else{
    locus <- sim.scans$locus
  }
  
  outcome.type <- NULL
  if (sim.scans$properties$return.value != "raw") { outcome.type <- sim.scans$properties$return.value }
  
  if (!is.null(sim.scans$properties$strain.reduction.fraction)) { 
    effect.title <- paste(paste0("QTL=",  sim.scans$properties$qtl.effect.size*100, "%"),
                          paste0("Strain=", sim.scans$properties$strain.effect.size*100, "%"),
                          paste0("Strain reduction=", round(sim.scans$properties$strain.reduction.fraction*100, 2), "%"),
                          sep=", ")
  } 
  else {
    effect.title <- paste(paste0("QTL=",  sim.scans$properties$qtl.effect.size*100, "%"),
                          paste0("Strain=", sim.scans$properties$strain.effect.size*100, "%"),
                          sep=", ")
  }
  
  this.title <- c(paste(outcome.type, 
                        ifelse(sim.scans$properties$return.means, "means:", "replicates:"),
                        effect.title),
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

#' @export
extract.compact.qr <- function(genomecache, 
                               CC.lines.matrix,
                               use.progress.bar=TRUE){
  
  results <- list(qr.alt.list=list(), CC.lines.matrix=CC.lines.matrix)
  for(i in 1:ncol(CC.lines.matrix)){
    this.qr <- extract.qr(genomecache=genomecache, 
                          model="additive",
                          formula=~1, 
                          id="SUBJECT.NAME.1", 
                          data=data.frame(SUBJECT.NAME.1=CC.lines.matrix[,i]),
                          use.progress.bar=use.progress.bar)
    results$qr.alt.list[[i]] <- this.qr$qr.list
    if(i == 1){
      results$shared <- list(intercept.allele=this.qr$intercept.allele,
                             condition.loci=this.qr$condition.loci,
                             qr.0=this.qr$qr.0,
                             chr=this.qr$chr,
                             pos=this.qr$pos,
                             model=this.qr$model,
                             founders=this.qr$founders,
                             formula=this.qr$formula)
    }
  }
  return(results)
}

#' @export
pull.qr.from.compact <- function(compact.qr.list, 
                                 qr.index){
  
  output.qr <- list(qr.list=compact.qr.list$qr.alt.list[[qr.index]],
                    intercept.allele=compact.qr.list$shared$intercept.allele,
                    condition.loci=compact.qr.list$shared$condition.loci,
                    qr.0=compact.qr.list$shared$qr.0,
                    chr=compact.qr.list$shared$chr,
                    pos=compact.qr.list$shared$pos,
                    model=compact.qr.list$shared$model,
                    founders=compact.qr.list$shared$founders,
                    subjects=compact.qr.list$CC.lines.matrix[,qr.index],
                    formula=compact.qr.list$shared$formula)
  return(output.qr)
}

  