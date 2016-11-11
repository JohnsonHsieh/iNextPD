check_datatype <- function(datatype){
  TYPE <- c("abundance", "incidence", "incidence_freq", "incidence_raw")
  if((is.na(pmatch(datatype, TYPE))) | (pmatch(datatype, TYPE) == -1))
    stop("invalid datatype")
  datatype <- match.arg(datatype, TYPE)
  if(datatype=="incidence_freq") datatype <- "incidence"
  return(datatype)
}

iNextPD.Ind <- function(x, labels, phy, q, size=NULL, endpoint=2*sum(x), knots=40, se=FALSE, conf=0.95, nboot=50){
  m <- size
  n <- sum(x)
  if(is.null(m)) {
    if(endpoint <= n) {
      m <- floor(seq(1, endpoint, length=floor(knots)))
    } else {
      m <- c(floor(seq(1, sum(x)-1, length.out=floor(knots/2)-1)), sum(x), floor(seq(sum(x)+1, to=endpoint, length.out=floor(knots/2))))
    }
    m <- c(1, m[-1])
  }else if(!is.null(m)) {	
    if(max(m)>n & length(m[m==n])==0){
      m <- c(m, n-1, n, n+1)
    }
    m <- c(1, m)
    m <- sort(unique(m))
  }
  
  if(se==TRUE & !is.null(conf)){
    if (!is.numeric(conf) || conf > 1 || conf < 0) {
      warning("\"conf\"(confidence level) must be a numerical value between 0 and 1, We use \"conf\" = 0.95 to calculate!")
      conf <- 0.95
    }
  }
  
  tmp <- ExpandData_(x, labels, phy, datatype="abundance")
  U <- tmp$branch_abun
  L <- tmp$branch_length
  names(x) <- labels
  x <- x[names(phy$leaves)]
  phd.hat <- PhD.m(n, x, U, L, q, m)
  C.hat <- Coverage_(x, "abundance", m)
  
  if(se==TRUE & nboot > 0 & length(x) > 1) {
    n <- sum(x)
    
    Prob.hat <- SPBoot_(x, "abundance")
    Abun.Mat <- rmultinom(nboot, n, Prob.hat)
    
    BootComm <- PDBoot_(x, labels, phy, "abundance")
    Prob.hat <- BootComm$branch_abun
    L.hat <- BootComm$branch_length
    U.Mat <- t(sapply(Prob.hat, function(p) rbinom(nboot, n, p)))
    
    error <-  qnorm(1-(1-conf)/2) * apply(apply(U.Mat, 2, function(U) PhD.m(n, x, U, L.hat, q, m)), 1, sd, na.rm=TRUE)
    left  <- phd.hat - error
    right <- phd.hat + error
    
    error.C <-  qnorm(1-(1-conf)/2) * apply(apply(Abun.Mat, 2, function(x) Coverage_(x, "abundance", m)), 1, sd, na.rm=TRUE)
    left.C  <- C.hat - error.C
    right.C <- C.hat + error.C
    
    right.C[right.C>=1] <- 1
    left.C[left.C<0] <- 0
    
    # left.C <- apply(apply(X.Mat, 2, function(x) Cvg.m(x, m)), 1, quantile, na.rm=TRUE, probs=0.025)
    # right.C <- apply(apply(X.Mat, 2, function(x) Cvg.m(x, m)), 1, quantile, na.rm=TRUE, probs=0.975)
    out <- data.frame("m"=m, "qPD"=phd.hat, "qPD.95.LCL"=left, "qPD.95.UCL"=right, "SC"=C.hat, "SC.95.LCL"=left.C, "SC.95.UCL"=right.C)
  }else{
    out <- data.frame("m"=m, "qPD"=phd.hat, "SC"=C.hat)
  }
  out <- data.frame(out)
  out$method <- ifelse(out$m<n, "interpolated", ifelse(out$m==n, "observed", "extrapolated"))
  out$order <- q
  id <- match(c("m", "method", "order", "qPD", "qPD.95.LCL", "qPD.95.UCL", "SC", "SC.95.LCL", "SC.95.UCL"), names(out), nomatch = 0)
  out <- out[, id]
  return(out)
}


iNextPD.Sam <- function(x, labels, phy, q, size=NULL, endpoint=2*ncol(x), knots=40, se=FALSE, conf=0.95, nboot=50){
  
  Spec <- iNEXT::as.incfreq(x)
  #if(which.max(Spec)!=1) 
  #  stop("invalid data structure!, first element should be number of sampling units")
  t <- size
  nT <- Spec[1]
  endpoint <- ifelse(is.null(endpoint), 2*ncol(x), endpoint)
  if(is.null(t)) {
    if(endpoint <= nT) {
      t <- floor(seq(1, endpoint, length.out=floor(knots)))
    } else {
      t <- c(floor(seq(1, nT-1, length.out=floor(knots/2)-1)), nT, floor(seq(nT+1, to=endpoint, length.out=floor(knots/2))))
    }
    t <- c(1, t[-1])
    t <- sort(unique(t))
  } else if(is.null(t)==FALSE) {	
    if(max(t)>nT & length(t[t==nT])==0)  t <- c(t, nT-1, nT, nT+1)
    t <- c(1, t)
    t <- sort(unique(t))
  }
  
  if(se==TRUE & !is.null(conf)){
    if (!is.numeric(conf) || conf > 1 || conf < 0) {
      warning("\"conf\"(confidence level) must be a numerical value between 0 and 1, We use \"conf\" = 0.95 to calculate!")
      conf <- 0.95
    }
  }
  
  tmp <- ExpandData_(x, labels, phy, datatype="incidence_raw")
  U <- tmp$branch_abun
  L <- tmp$branch_length
  xx <- rowSums(x)
  names(xx) <- labels
  xx <- xx[names(phy$leaves)]
  phd.hat <- PhD.m(nT, xx, U, L, q, t)
  C.hat <- Coverage_(Spec, "incidence",t)
  
  if(se==TRUE & nboot > 0 & length(x) > 1) {
    Prob.hat <- SPBoot_(Spec, "incidence")
    Abun.Mat <- t(sapply(Prob.hat, function(p) rbinom(nboot, nT, p)))
    
    BootComm <- PDBoot_(x, labels, phy, "incidence_raw")
    Prob.hat <- BootComm$branch_abun
    L.hat <- BootComm$branch_length
    U.Mat <- t(sapply(Prob.hat, function(p) rbinom(nboot, nT, p)))
    
    error <-  qnorm(1-(1-conf)/2) * apply(apply(U.Mat, 2, function(U) PhD.m(nT, x, U, L.hat, q, t)), 1, sd, na.rm=TRUE)
    left  <- phd.hat - error
    right <- phd.hat + error
    
    error.C <-  qnorm(1-(1-conf)/2) * apply(apply(Abun.Mat, 2, function(x) Coverage_(x, "incidence", t)), 1, sd, na.rm=TRUE)
    left.C  <- C.hat - error.C
    right.C <- C.hat + error.C
    
    right.C[right.C>=1] <- 1
    left.C[left.C<0] <- 0
    
    # left.C <- apply(apply(X.Mat, 2, function(x) Cvg.m(x, m)), 1, quantile, na.rm=TRUE, probs=0.025)
    # right.C <- apply(apply(X.Mat, 2, function(x) Cvg.m(x, m)), 1, quantile, na.rm=TRUE, probs=0.975)
    out <- data.frame("t"=t, "qPD"=phd.hat, "qPD.95.LCL"=left, "qPD.95.UCL"=right, "SC"=C.hat, "SC.95.LCL"=left.C, "SC.95.UCL"=right.C)
  }else{
    out <- data.frame("t"=t, "qPD"=phd.hat, "SC"=C.hat)
  }
  out <- data.frame(out)
  out$method <- ifelse(out$t<nT, "interpolated", ifelse(out$t==nT, "observed", "extrapolated"))
  out$order <- q
  id <- match(c("t", "method", "order", "qPD", "qPD.95.LCL", "qPD.95.UCL", "SC", "SC.95.LCL", "SC.95.UCL"), names(out), nomatch = 0)
  out <- out[, id]
  return(out)
}


#' iNterpolation and EXTrapolation of Hill number
#' 
#' \code{iNextPD}: Interpolation and extrapolation of Hill number with order q.
#' 
#' @param x a matrix, data.frame (species by sites), or list of species abundances or incidence data.
#' @param labels species names for object x.
#' @param phy a phylog objcet for input phylo-tree.
#' @param q a numeric value specifying the diversity order of Hill number .
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}),  
#' or species by sampling-units incidence matrix (\code{datatype = "incidence_raw"}).
#' @param size an integer vector of sample sizes (number of individuals or sampling units) for which diversity estimates will be computed. 
#' If NULL, then diversity estimates will be computed for those sample sizes determined by the specified/default \code{endpoint} and \code{knots} .
#' @param endpoint an integer specifying the sample size that is the \code{endpoint} for rarefaction/extrapolation. 
#' If NULL, then \code{endpoint} \code{=} double reference sample size.
#' @param knots an integer specifying the number of equally-spaced \code{knots} (say K, default is 40) between size 1 and the \code{endpoint};
#' each knot represents a particular sample size for which diversity estimate will be calculated.  
#' If the \code{endpoint} is smaller than the reference sample size, then \code{iNextPD()} computes only the rarefaction esimates for approximately K evenly spaced \code{knots}. 
#' If the \code{endpoint} is larger than the reference sample size, then \code{iNextPD()} computes rarefaction estimates for approximately K/2 evenly spaced \code{knots} between sample size 1 and the reference sample size, and computes extrapolation estimates for approximately K/2 evenly spaced \code{knots} between the reference sample size and the \code{endpoint}.
#' @param se a logical variable to calculate the bootstrap standard error and \code{conf} confidence interval.
#' @param conf a positive number < 1 specifying the level of confidence interval, default is 0.95.
#' @param nboot an integer specifying the number of replications.
#' @return a list of three objects: \code{$DataInfo} for summarizing data information; 
#' \code{$iNextPDEst} for showing diversity estimates for rarefied and extrapolated samples along with related statistics;
#' \code{$AsyPDEst} for showing asymptotic diversity estimates along with related statistics, and \code{$ExpandData} (xi, Li, i=1,2,...,B).
#' @examples
#' data(bird)
#' bird.abu <- bird$abun
#' bird.lab <- rownames(bird$abun)
#' bird.phy <- ade4::newick2phylog(bird$tre)
#' iNextPD(bird.abu, labels=bird.lab, phy=bird.phy, q=0, datatype="abundance")
#' @export
#' 
iNextPD <- function(x, labels, phy, q=0, datatype="abundance", size=NULL, endpoint=NULL, knots=40, se=FALSE, conf=0.95, nboot=50){
  
  datatype <- check_datatype(datatype)
  
  if(datatype=="incidence_freq" | datatype=="incidence") 
    stop('only support datatype="incidence_raw"')
  
  Fun <- function(x, q){
    if(datatype == "abundance"){
      if(sum(x)==0) stop("Zero abundance counts in one or more sample sites")
      x <- as.numeric(unlist(x))
      out <- iNextPD.Ind(x, labels, phy, q, size, endpoint=ifelse(is.null(endpoint), 2*sum(x), endpoint), knots, se, conf, nboot)
    }
    if(datatype == "incidence_raw"){
      y <- iNEXT::as.incfreq(x)
      t <- y[1]
      y <- y[-1]
      if(t>sum(y)){
        warning("Insufficient data to provide reliable estimators and associated s.e.") 
      }
      if(sum(y)==0) stop("Zero incidence frequencies in one or more sample sites")
      
      out <- iNextPD.Sam(x, labels, phy, q, size, endpoint=ifelse(is.null(endpoint), 2*max(x), endpoint), knots, se, conf, nboot)  
    }
    out
  }
  
  if(class(q) != "numeric")
    stop("invlid class of order q, q should be a postive value/vector of numeric object")
  if(min(q) < 0){
    warning("ambigous of order q, we only compute postive q")
    q <- q[q >= 0]
  }
  
  if(datatype=="abundance"){
    if(class(x)=="matrix" | class(x)=="data.frame"){
      name_x <- colnames(x)
      x <- lapply(seq_len(ncol(x)), function(i) x[,i])
      names(x) <- name_x
    }
  }
  
  if(datatype=="abundance"){
    
    if(class(x)=="numeric" | class(x)=="integer" | class(x)=="double"){
      out <- do.call("rbind", lapply(q, function(q) Fun(x, q)))
      out[,-(1:3)] <- round(out[,-(1:3)],3)
      index <- rbind(as.matrix(estPD(x, labels, phy, q=0, datatype, se=TRUE, conf=conf)), 
                     as.matrix(estPD(x, labels, phy, q=1, datatype, se=TRUE, conf=conf)),
                     as.matrix(estPD(x, labels, phy, q=2, datatype, se=TRUE, conf=conf)))
      rownames(index) <- c("q = 0", "q = 1", "q = 2")
      
    }else if(class(x)=="list"){
      out <- lapply(x, function(x) {
        tmp <- do.call("rbind", lapply(q, function(q) Fun(x,q)))
        tmp[,-(1:3)] <- round(tmp[,-(1:3)],3)
        tmp
      })
      arr <- array(0, dim = c(3, 5, length(x)))
      arr[1,,] <- t(as.matrix(estPD(x, labels, phy, q=0, datatype, se=TRUE, conf=conf)))
      arr[2,,] <- t(as.matrix(estPD(x, labels, phy, q=1, datatype, se=TRUE, conf=conf)))
      arr[3,,] <- t(as.matrix(estPD(x, labels, phy, q=2, datatype, se=TRUE, conf=conf)))  
      dimnames(arr)[[3]] <- names(x)
      dimnames(arr)[[1]] <- c("q = 0", "q = 1", "q = 2")
      dimnames(arr)[[2]] <- c("Observed", "Estimator", "Est_s.e.", "95% Lower", "95% Upper")
      index <- ftable(arr, row.vars = c(3,1))
    }else{
      stop("invalid class of x, x should be a object of numeric, matrix, data.frame, or list")
    }
  }else if(datatype=="incidence_raw"){
    if(class(x)=="matrix" | class(x)=="data.frame"){
      out <- do.call("rbind", lapply(q, function(q) Fun(x, q)))
      out[,-(1:3)] <- round(out[,-(1:3)],3)
      index <- rbind(as.matrix(estPD(x, labels, phy, q=0, datatype, se=TRUE, conf=conf)), 
                     as.matrix(estPD(x, labels, phy, q=1, datatype, se=TRUE, conf=conf)),
                     as.matrix(estPD(x, labels, phy, q=2, datatype, se=TRUE, conf=conf)))
      rownames(index) <- c("q = 0", "q = 1", "q = 2")
    }else if(class(x)=="list"){
      out <- lapply(x, function(x) {
        tmp <- do.call("rbind", lapply(q, function(q) Fun(x,q)))
        tmp[,-(1:3)] <- round(tmp[,-(1:3)],3)
        tmp
      })
      
      arr <- array(0, dim = c(3, 5, length(x)))
      arr[1,,] <- t(as.matrix(estPD(x, labels, phy, q=0, datatype, se=TRUE, conf=conf)))
      arr[2,,] <- t(as.matrix(estPD(x, labels, phy, q=1, datatype, se=TRUE, conf=conf)))
      arr[3,,] <- t(as.matrix(estPD(x, labels, phy, q=2, datatype, se=TRUE, conf=conf)))  
      dimnames(arr)[[3]] <- names(x)
      dimnames(arr)[[1]] <- c("q = 0", "q = 1", "q = 2")
      dimnames(arr)[[2]] <- c("Observed", "Estimator", "Est_s.e.", "95% Lower", "95% Upper")
      index <- ftable(arr, row.vars = c(3,1))
    }
  }
  
  info <- DataInfo(x, datatype)
  tree <- ExpandData(x, labels, phy, datatype)
  
  z <- list("DataInfo"=info, "iNextPDEst"=out, "AsyPDEst"=index, "ExpandData"=tree)
  class(z) <- c("iNextPD")
  return(z)
}

