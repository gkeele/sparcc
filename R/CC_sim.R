#' Runs the genome scans of the simulated data output by sim.CC.data()
#'
#' This function takes the output from sim.CC.data() and performs the genome scans. Internally it runs
#' the require QR decompositions, which it can save for later if specified.
#'
#' @param sim.data Output simulated data from sim.CC.data().
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
  if (return.all.sim.qr & is.null(all.sim.qr)){ all.sim.qr <- list() }
  for (i in 1:length(scan.index)) {
    if (sim.data$properties$vary.lines) {
      if (compute.qr.boolean) {
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
    if (return.all.sim.qr & sim.data$properties$vary.lines) { all.sim.qr[[i]] <- this.qr }
    else if (return.all.sim.qr & !sim.data$properties$vary.lines) { all.sim.qr <- this.qr }
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
#' @param sample.as.method DEFAULT: "uniform". The procedure used for sampling the allelic series. If
#' every strain has its own allele, this option does not matter. Alternatively, a Chinese restaurant process ("crp") can be used,
#' which is possibly more biologically accurate, and will favor allelic series that are less balanced (1 vs 7).
#' @param num.alleles DEFAULT: 8. The number of functional alleles. Must be less than or equal to the number of 
#' founders.
#' @param num.founders DEFAULT: 8. The number of founders, which must correspond to the genome cache. The CC has
#' eight.
#' @param qtl.effect.size The size of the simulated QTL effect. The scale of the input is in proportion of
#' the phenotypic variance due to the QTL, thus should be greater than or equal to zero, and less than one.
#' @param beta DEFAULT: NULL. Allows for the manual specification of QTL effect. Is expected to be a vector
#' the length of the number of alleles. It will be scaled based on qtl.effect.size.
#' @param strain.effect.size DEFAULT: 0. The size of the simulated strain effect, which represents something akin to a polygenic 
#' effect. Other variants specific to CC lines will result in overall strain-specific effects. The scale of the input 
#' is in proportion of the phenotypic variance due to the strain, thus should be greater than or equal to zero, and less than one.
#' @param impute DEFAULT: TRUE. If TRUE, the QTL portion of the design matrix in the simulation is a realized sampling
#' of haplotypes from the probabilities. If FALSE, the simulations are based on the probabilities, which is flawed in
#' terms of biological reality.
#' @param scale.qtl.mode DEFAULT: "B". Specifies how the QTL effect is scaled. If "B", then the variance of the qtl effect
#' vector beta is scaled to the effect size specified in qtl.effect.size, which would be the variance explained in a population
#' perfectly balanced in terms of the functional alleles. If "MB", the the variance of M %*% beta is scaled, which would be
#' the variance explained in a population balanced in terms of founder strains with the allelic series that is specified in M. 
#' If "DAMB", the variance of D %*% A %*% M %*% beta is scaled, which would be the variance explained in a population of a 
#' specific set of CC strains (specified in D). If "ZAMB", the variance of Z %*% D %*% A %*% M %*% beta is scaled, which would
#' be the variance explained in the specific population of CC strains the number of replicates (specified in Z).
#' @param scale.strain.mode DEFAULT: "d". Specifies how the strain effect is scaled. If "d", then the variance of
#' delta is scaled to the effect size specified in strain.effect.size, which would be the varianced explained in the
#' population without replicates. If "Zd", then Z %*% delta is scaled, which would be the variance explained in a population with
#' the number of replicates specified in Z. This option more strongly controls the variation from strain effect.
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
                        sample.as.method=c("uniform", "crp"),
                        num.alleles=8, 
                        num.founders=8,
                        qtl.effect.size,
                        beta=NULL,
                        strain.effect.size=0,
                        impute=TRUE,
                        scale.qtl.mode=c("B", "MB", "DAMB", "ZDAMB", "none"),
                        scale.strain.mode=c("d", "Zd"),
                        return.value=c("raw", "fixef.resid", "ranef.resid"),
                        return.means=TRUE){
  
  h <- miqtl::DiploprobReader$new(genomecache)
  scale.qtl.mode <- scale.qtl.mode[1]
  scale.strain.mode <- scale.strain.mode[1]
  return.value <- return.value[1]
  sample.as.method <- sample.as.method[1]
  
  original.effects <- list(qtl.effect.size=qtl.effect.size,
                           strain.effect.size=strain.effect.size)
  
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
  var.table <- matrix(NA, nrow=num.sim, ncol=8)
  sim.matrix <- matrix(NA, nrow=num.ind, ncol=num.sim)
  id.matrix <- matrix(NA, nrow=num.ind, ncol=ifelse(vary.lines, num.sim, 1))
  for(i in 1:num.sim){
    this.sim <- simulate.CC.qtl(CC.lines=CC.lines[,cc.index[i]], 
                                num.replicates=num.replicates,
                                M=M,
                                sample.as.method=sample.as.method,
                                num.alleles=num.alleles, 
                                num.founders=num.founders,
                                qtl.effect.size=qtl.effect.size, 
                                strain.effect.size=strain.effect.size,
                                impute=impute,
                                scale.qtl.mode=scale.qtl.mode,
                                scale.strain.mode=scale.strain.mode,
                                locus.matrix=h$getLocusMatrix(locus=locus[locus.index[i]], model="full"),
                                return.value=return.value,
                                return.means=return.means,
                                num.sim=1)
    this.var.table <- this.sim$properties$var.table
    this.sim <- this.sim$data
    sim.matrix[,i] <- this.sim[,1]
    var.table[i,] <- this.var.table[1,]
    if(vary.lines | i == 1){
      id.matrix[,i] <- as.character(this.sim[,2])
    }
  }
  colnames(sim.matrix) <- paste0("sim.y.", 1:num.sim)
  colnames(id.matrix) <- paste0("SUBJECT.NAME.", 1:ncol(id.matrix))
  colnames(var.table) <- colnames(this.var.table)
  
  outcome <- data.frame(sim.matrix, id.matrix)
  return(list(data=outcome,
              locus=locus,
              locus.pos=list(Mb=h$getLocusStart(locus, scale="Mb"),
                             cM=h$getLocusStart(locus, scale="cM")),
              genomecache=genomecache,
              properties=list(num.alleles=num.alleles,
                              sample.as.method=sample.as.method,
                              num.replicates=num.replicates,
                              num.founders=num.founders,
                              qtl.effect.size=original.effects$qtl.effect.size, 
                              strain.effect.size=original.effects$strain.effect.size,
                              num.lines=num.lines,
                              impute=impute,
                              scale.qtl.mode=scale.qtl.mode,
                              scale.strain.mode=scale.strain.mode,
                              M.ID=M.ID,
                              vary.lines=vary.lines,
                              vary.locus=vary.locus,
                              return.value=return.value,
                              return.means=return.means,
                              var.table=var.table)))
}

simulate.CC.qtl <- function(CC.lines, 
                            num.replicates,
                            M=NULL,
                            sample.as.method=c("uniform", "crp"),
                            qtl.effect.size, 
                            beta=NULL,
                            strain.effect.size,
                            num.alleles=8, 
                            num.founders=8,
                            locus.matrix, 
                            num.sim,
                            impute=TRUE, 
                            scale.qtl.mode=c("B", "MB", "DAMB", "ZDAMB", "none"),
                            scale.strain.mode=c("d", "Zd"),
                            return.value=c("raw", "fixef.resid", "ranef.resid"),
                            return.means=TRUE,
                            ...){
  scale.qtl.mode <- scale.qtl.mode[1]
  scale.strain.mode <- scale.strain.mode[1]
  return.value <- return.value[1]
  sample.as.method <- sample.as.method[1]

  this.locus.matrix <- locus.matrix[CC.lines,]
  ## Imputation for simulation
  if(impute){
    this.locus.matrix <- t(apply(this.locus.matrix, 1, function(x) rmultinom(1, 1, x)))
    rownames(this.locus.matrix) <- CC.lines
  }
  D <- this.locus.matrix ## Saved for variance calculations
  this.locus.matrix <- this.locus.matrix[rep(1:length(CC.lines), each=num.replicates),]
  
  full.to.add.matrix <- t(miqtl::straineff.mapping.matrix(M=num.founders))
  
  ## QTL
  this.locus.matrix <- tcrossprod(this.locus.matrix, full.to.add.matrix)
  if (qtl.effect.size != 0) {
    QTL.effect <- simulate.QTL.model.and.effects(num.alleles=num.alleles, 
                                                 num.founders=num.founders, 
                                                 M=M,
                                                 sample.as.method=sample.as.method,
                                                 beta=beta,
                                                 ...)
  }
  else {
    QTL.effect <- list(M=diag(8),
                       beta=rep(0, 8))
  }
  if (scale.qtl.mode != "none") {
    beta <- as.vector(scale(QTL.effect$beta))
  }
  if (scale.qtl.mode == "B") {
    beta <- as.vector(0.5*beta*sqrt(qtl.effect.size)) 
    
    QTL.predictor <- tcrossprod(tcrossprod(this.locus.matrix, t(QTL.effect$M)), matrix(beta, nrow=1))
  }
  else if (scale.qtl.mode == "MB") {
    var.ratio <- var(2*QTL.effect$M %*% beta)/var(2*beta)

    beta <- as.vector(0.5*beta*sqrt(qtl.effect.size))*sqrt(1/var.ratio)
  }
  else if (scale.qtl.mode == "DAMB") {
    var.ratio <- var(D %*% t(full.to.add.matrix) %*% QTL.effect$M %*% beta)/var(2*beta)
    
    if (var.ratio != 0) { # Case when more than one allele is observed
      beta <- as.vector(0.5*beta*sqrt(qtl.effect.size))*sqrt(1/var.ratio)
    }
    else { # Case when only one allele is observed, should result in a null scan
      beta <- as.vector(0.5*beta*sqrt(qtl.effect.size))
    }
  }
  else if (scale.qtl.mode == "ZDAMB") {
    var.ratio <- var(this.locus.matrix %*% QTL.effect$M %*% beta)/var(2*beta)
    
    if (var.ratio != 0) { # Case when more than one allele is observed
      beta <- as.vector(0.5*beta*sqrt(qtl.effect.size))*sqrt(1/var.ratio)
    }
    else { # Case when only one allele is observed, should result in a null scan
      beta <- as.vector(0.5*beta*sqrt(qtl.effect.size))
    }
  }
  QTL.predictor <- tcrossprod(tcrossprod(this.locus.matrix, t(QTL.effect$M)), matrix(beta, nrow=1))

  ## Strain
  if (strain.effect.size != 0) {
    this.strain.matrix <- incidence.matrix(factor(CC.lines, levels=CC.lines))
    this.strain.matrix <- this.strain.matrix[rep(1:nrow(this.strain.matrix), each=num.replicates),]
    
    strain.effect <- rnorm(n=length(CC.lines))
    strain.effect <- as.vector(scale(strain.effect))

    strain.predictor <- this.strain.matrix %*% matrix(strain.effect, ncol=1)

    if (scale.strain.mode == "Zd") {
      strain.predictor <- as.vector(scale(strain.predictor))
      strain.predictor <- as.vector(strain.predictor*sqrt(strain.effect.size))
    }
  }
  else{ strain.predictor <- rep(0, length(CC.lines)) }
  
  var.table <- matrix(NA, nrow=num.sim, ncol=8)
  sim.data <- matrix(NA, nrow=nrow(this.locus.matrix), ncol=num.sim)
  for (i in 1:num.sim) {
    scaled.resid <- calc.scaled.residual(qtl.effect.size=qtl.effect.size, 
                                         strain.effect.size=strain.effect.size,
                                         n=nrow(this.locus.matrix))
    this.denominator <- 
    sim.data[,i] <- QTL.predictor + strain.predictor + scaled.resid
    var.table[i,] <- c(var(2*beta), 
                       var(2*QTL.effect$M %*% beta),
                       var(D %*% t(full.to.add.matrix) %*% QTL.effect$M %*% beta),
                       var(this.locus.matrix %*% QTL.effect$M %*% beta),
                       var(2*beta)/(var(2*beta) + var(strain.predictor) + var(scaled.resid)),
                       var(2*QTL.effect$M %*% beta)/(var(2*QTL.effect$M %*% beta) + var(strain.predictor) + var(scaled.resid)),
                       var(D %*% t(full.to.add.matrix) %*% QTL.effect$M %*% beta)/(var(D %*% t(full.to.add.matrix) %*% QTL.effect$M %*% beta) + var(strain.predictor) + var(scaled.resid)),
                       var(this.locus.matrix %*% QTL.effect$M %*% beta)/(var(this.locus.matrix %*% QTL.effect$M %*% beta) + var(strain.predictor) + var(scaled.resid)))
  }
  
  colnames(sim.data) <- paste0("sim.y.", 1:ncol(sim.data))
  sim.data <- data.frame(sim.data, "SUBJECT.NAME"=rep(CC.lines, each=num.replicates))
  colnames(var.table) <- c("B.effect", "MB.effect", "DAMB.effect", "ZDAMB.effect",
                           "B.ve", "MB.ve", "DAMB.ve", "ZDAMB.ve")
  
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
                              sample.as.method=sample.as.method,
                              num.replicates=num.replicates,
                              num.lines=length(CC.lines),
                              impute=impute,
                              return.value=return.value,
                              return.means=return.means,
                              var.table=var.table)))
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
                                           sample.as.method=c("uniform", "crp"), 
                                           beta=NULL,
                                           ...){
  
  sample.as.method <- sample.as.method[1]
  if(is.null(M)){
    M <- matrix(0, num.founders, num.alleles)
    M[cbind(sample(1:num.founders, num.alleles), 1:num.alleles)] <- 1
    
    if(sample.as.method == "uniform"){
      M[which(rowSums(M)==0),] <- t(rmultinom(num.founders - num.alleles, 1, rep(1/num.alleles, num.alleles)))
    }
    else if(sample.as.method == "crp"){
      for(i in which(rowSums(M)==0)){
        M[i,] <- rmultinom(1, 1, colSums(M)/sum(M))
      }
    }
  }
  
  if(is.null(beta)){
    beta <- rnorm(num.alleles)
  }
  #beta <- as.vector(scale(beta))
  #beta <- as.vector(0.5*beta*sqrt(effect.var)) 
  
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
                       scale=scale, mark.locus=locus, hard.thresholds=thresh[phenotype.index], ...)
}

#' @export
sim.data.effect.plot <- function(sim.data,
                             x.var=c("B.effect", "MB.effect", "DAMB.effect", "ZDAMB.effect"),
                             y.var=c("MB.effect", "DAMB.effect", "ZDAMB.effect", "B.effect"),
                             qtl.col="red", title="") {
  x.var <- x.var[1]
  y.var <- y.var[1]
  nominal.effect.size <- sim.data$properties$qtl.effect.size
  var.table <- sim.data$properties$var.table
  
  x.val <- var.table[,x.var]
  y.val <- var.table[,y.var]
  
  x.max <- max(x.val, nominal.effect.size)
  y.max <- max(y.val, nominal.effect.size)
  
  plot(x.val, y.val,
       xlab=x.var, ylab=y.var, xlim=c(0, x.max), ylim=c(0, y.max),
       frame.plot=FALSE, pch=20, las=1, main=title)
  abline(0, 1, lty=3)
  abline(v=nominal.effect.size, lty=2, col=qtl.col)
  abline(h=nominal.effect.size, lty=2, col=qtl.col)
}

#' @export
sim.data.ve.plot <- function(sim.data,
                             x.var=c("B.ve", "MB.ve", "DAMB.ve", "ZDAMB.ve"),
                             y.var=c("MB.ve", "DAMB.ve", "ZDAMB.ve", "B.ve"),
                             qtl.col="red", title="") {
  x.var <- x.var[1]
  y.var <- y.var[1]
  nominal.effect.size <- sim.data$properties$qtl.effect.size
  var.table <- sim.data$properties$var.table
  
  x.val <- var.table[,x.var]
  y.val <- var.table[,y.var]
  
  plot(x.val, y.val,
       xlab=x.var, ylab=y.var, xlim=c(0, 1), ylim=c(0, 1),
       frame.plot=FALSE, pch=20, las=1, main=title)
  abline(0, 1, lty=3)
  abline(v=nominal.effect.size, lty=2, col=qtl.col)
  abline(h=nominal.effect.size, lty=2, col=qtl.col)
}

#' @export
sim.data.all.ve.plot <- function(sim.data,
                                 x.var=c("B.ve", "MB.ve", "DAMB.ve", "ZDAMB.ve"),
                                 include.types=c("B.ve", "MB.ve", "DAMB.ve"),
                                 ve.col=c("#A6CEE3", "#B2DF8A", "#FB9A99", "#FDBF6F"), # alts #1F78B4 #33A02C #E31A1C #FF7F00
                                 ve.pch=c(15, 17, 18, 19),
                                 qtl.col="red",
                                 transparency=0.5,
                                 legend.pos="topleft",
                                 title="") {
  x.var <- x.var[1]
  nominal.effect.size <- sim.data$properties$qtl.effect.size
  original.types <- colnames(sim.data$properties$var.table)[grep(x=colnames(sim.data$properties$var.table), pattern="ve")]

  var.table <- sim.data$properties$var.table[,unique(c(x.var, include.types))]
  
  ve.col <- c(scales::alpha(ve.col[1], transparency),
              scales::alpha(ve.col[2], transparency),
              scales::alpha(ve.col[3], transparency),
              scales::alpha(ve.col[4], transparency))
  ve.col <- ve.col[original.types %in% include.types]
  ve.pch <- ve.pch[original.types %in% include.types]
  
  x.val <- var.table[,x.var]
  keep <- !(colnames(var.table) %in% x.var)
  y.tab <- var.table[,keep]
  ve.col <- ve.col[keep]
  ve.pch <- ve.pch[keep]
  
  plot(x.val, y.tab[,ncol(y.tab)],
       xlab=x.var, ylab="Proportion of variance", xlim=c(0, 1), ylim=c(0, 1),
       frame.plot=FALSE, pch=ve.pch[ncol(y.tab)], las=1, col=ve.col[ncol(y.tab)], cex=1.5, main=title)
  if (ncol(y.tab) > 1) {
    for (i in rev(2:ncol(y.tab)-1)) {
      points(x.val, y.tab[,i], col=ve.col[i], pch=ve.pch[i], cex=1.5)
    }
  }
  abline(0, 1, lty=3)
  abline(v=nominal.effect.size, lty=2, col=qtl.col)
  abline(h=nominal.effect.size, lty=2, col=qtl.col)
  
  legend(legend.pos, legend=colnames(y.tab), col=ve.col, pch=ve.pch, cex=1.5, bty="n")
}

#' @export
sim.data.all.effect.plot <- function(sim.data,
                                 x.var=c("B.effect", "MB.effect", "DAMB.effect", "ZDAMB.effect"),
                                 include.types=c("B.effect", "MB.effect", "DAMB.effect"),
                                 ve.col=c("#A6CEE3", "#B2DF8A", "#FB9A99", "#FDBF6F"), # alts #1F78B4 #33A02C #E31A1C #FF7F00
                                 ve.pch=c(15, 17, 18, 19),
                                 qtl.col="red",
                                 transparency=0.5,
                                 legend.pos="topleft",
                                 title="") {
  x.var <- x.var[1]
  nominal.effect.size <- sim.data$properties$qtl.effect.size
  original.types <- colnames(sim.data$properties$var.table)[grep(x=colnames(sim.data$properties$var.table), pattern="effect")]
  
  var.table <- sim.data$properties$var.table[,unique(c(x.var, include.types))]
  
  ve.col <- c(scales::alpha(ve.col[1], transparency),
              scales::alpha(ve.col[2], transparency),
              scales::alpha(ve.col[3], transparency),
              scales::alpha(ve.col[4], transparency))
  ve.col <- ve.col[original.types %in% include.types]
  ve.pch <- ve.pch[original.types %in% include.types]
  
  x.val <- var.table[,x.var]
  keep <- !(colnames(var.table) %in% x.var)
  y.tab <- var.table[,keep]
  ve.col <- ve.col[keep]
  ve.pch <- ve.pch[keep]
  
  x.max <- max(x.val, nominal.effect.size)
  y.max <- max(y.tab, nominal.effect.size)
  
  plot(x.val, y.tab[,ncol(y.tab)],
       xlab=x.var, ylab="Proportion of variance", xlim=c(0, x.max), ylim=c(0, y.max),
       frame.plot=FALSE, pch=ve.pch[ncol(y.tab)], las=1, col=ve.col[ncol(y.tab)], cex=1.5, main=title)
  if (ncol(y.tab) > 1) {
    for (i in rev(2:ncol(y.tab)-1)) {
      points(x.val, y.tab[,i], col=ve.col[i], pch=ve.pch[i], cex=1.5)
    }
  }
  abline(0, 1, lty=3)
  abline(v=nominal.effect.size, lty=2, col=qtl.col)
  abline(h=nominal.effect.size, lty=2, col=qtl.col)
  
  legend(legend.pos, legend=colnames(y.tab), col=ve.col, pch=ve.pch, cex=1.5, bty="n")
}

#' Calculate the QTL mapping power from SPARCC genome scans output from run.sim.scans()
#'
#' This function takes output genome scans from run.sim.scans() and calculates the power to map the simulated
#' QTL. A window around the QTL can be specified to allow peaks not immediately at the simulated locus but likely
#' tagging the QTL to also be included.
#'
#' @param sim.scans Output simulated genome scans from run.sim.scans().
#' @param thresh A list of threshold, calculated from get.gev.thresh(), that correspond to the scans in sim.scans. 
#' @param window.mb DEFAULT: 5. Loci upstream and downstream the specified window.mb in Mb will also be checked 
#' for statistically significant signals. Sometimes the statistical score will not pass at the simulated QTL, but
#' does at nearby loci.
#' @export
#' @examples pull.power()
pull.power <- function(sim.scans, thresh, window.mb=5){
  map <- rep(NA, nrow(sim.scans$p.value))
  total <- length(sim.scans$p.value)
  for (i in 1:nrow(sim.scans$p.value)) {
    this.index <- which(colnames(sim.scans$p.value) == sim.scans$locus[i])
    this.chr <- sim.scans$chr[this.index]
    this.pos <- sim.scans$pos$Mb[this.index]
    p.value <- sim.scans$p.value[i, sim.scans$pos$Mb >= this.pos - window.mb & sim.scans$pos$Mb <= this.pos + window.mb & sim.scans$chr == this.chr]
    map[i] <- any(-log10(p.value) > thresh[i])
  }
  return(mean(map))
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