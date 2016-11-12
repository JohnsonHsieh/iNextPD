invChat_ <- function(x, datatype="abundance", C=NULL, digits=4){
  datatype <- check_datatype(datatype)
  
  if(datatype=="abundance"){
    n <- sum(x)
  }else if(datatype=="incidence"){
    n <- x[1]
  }else if(datatype=="incidence_raw"){
    n <- ncol(x)    
  }

  refC <- Coverage_(x, datatype, n)
  
  if(is.null(C)){ C <- refC }
  
  invChatSub_ <- function(C){
    f <- function(m, C) abs(Coverage_(x, datatype, m)-C)
    if(refC > C){
      opt <- optimize(f, C=C, lower=0, upper=sum(x))
      mm <- opt$minimum
      #mm <- round(mm)
    }else if(refC <= C){
      f1 <- sum(x==1)
      f2 <- sum(x==2)
      if(f1>0 & f2>0){A <- (n-1)*f1/((n-1)*f1+2*f2)}
      if(f1>1 & f2==0){A <- (n-1)*(f1-1)/((n-1)*(f1-1)+2)}
      if(f1==1 & f2==0){A <- 1}
      if(f1==0 & f2==0){A <- 1}
      mm <- (log(n/f1)+log(1-C))/log(A)-1
      mm <- n+mm
      # mm <- round(mm)
    }
    if(mm > 2*n) 
      warning("The maximum size of the extrapolation exceeds double reference sample size, the results for q = 0 may be subject to large prediction bias.")
    mm
  }
  out <- sapply(C, invChatSub_)
  return(data.frame(m=out, SC=C))
}

# Compute sample size with fixed sample coverage
invChat <- function(x, datatype="abundance", C=NULL, digits=4){
  if(is.null(C)){
    C <- min(DataInfo(x, datatype)[,4])
  }
  if(class(x)=="numeric" | class(x)=="integer"){
    data.frame(site="Site1", invChat_(x, datatype, C, digits))
  }else if(class(x)=="list"){
    if(is.null(names(x))){
      names(x) <- paste0("Site", 1:length(x))
    }
    out <- do.call(rbind, lapply(x, function(x) invChat_(x, datatype, C, digits)))
    out <- data.frame(site=rep(names(x), each=length(C)), out)
    rownames(out) <- NULL
    out
  }else if(class(x)=="data.frame" | class(x)=="matrix"){
    if(is.null(colnames(x))){
      colnames(x) <- paste0("Site", 1:ncol(x))
    }
    out <- do.call(rbind, apply(x, 2, function(x) invChat_(x, datatype, C, digits)))
    out <- data.frame(site=rep(colnames(x), each=length(C)), out)
    rownames(out) <- NULL
    out
  }
}

# Compute sample coverage with fixed sample size
invSize <- function(x, datatype="abundance", size=NULL, digits=4){
  if(is.null(size)){
    size <- min(DataInfo(x, datatype)[,2])
  }
  if(class(x)=="numeric" | class(x)=="integer"){
    data.frame(site="Site1", m=size, SC=round(Coverage_(x, datatype, size), digits))
  }else if(class(x)=="list"){
    if(is.null(names(x))){
      names(x) <- paste0("Site", 1:length(x))
    }
    out <- do.call(cbind, lapply(x, function(x) {
      data.frame(m=size, SC=round(Coverage_(x, datatype, size), digits))
      }))
    out <- data.frame(site=rep(names(x), each=length(size)), out)
    rownames(out) <- NULL
    out
  }else if(class(x)=="data.frame" | class(x)=="matrix"){
    if(is.null(colnames(x))){
      colnames(x) <- paste0("Site", 1:ncol(x))
    }
    out <- do.call(rbind, apply(x, 2, function(x) {
      data.frame(m=size, SC=round(Coverage_(x, datatype, size), digits))
      }))
    out <- data.frame(site=rep(colnames(x), each=length(size)), out)
    rownames(out) <- NULL
    out
  }
}





###############################################
#' Compute phylogenetic diversity with a particular of sample size/coverage 
#' 
#' \code{estimatePD}: computes phylogenetic diversity (Hill numbers with q = 0, 1 and 2) with a particular user-specified level of sample size or sample coverage.
#' @param x a \code{data.frame} or \code{list} of species abundances or incidence frequencies.\cr 
#' If \code{datatype = "incidence"}, then the first entry of the input data must be total number of sampling units, followed 
#' by species incidence frequencies in each column or list.
#' @param labels species names for object x
#' @param phy a phylog objcet for input phylo-tree
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}),  
#' or species by sampling-units incidence matrix (\code{datatype = "incidence_raw"}).
#' @param base comparison base: sample-size-based (\code{base="size"}) or coverage-based \cr (\code{base="coverage"}).
#' @param level an integer specifying a particular sample size or a number (between 0 and 1) specifying a particular value of sample coverage. 
#' If \code{base="size"} and \code{level=NULL}, then this function computes the diversity estimates for the minimum sample size among all sites. 
#' If \code{base="coverage"} and \code{level=NULL}, then this function computes the diversity estimates for the minimum sample coverage among all sites. 
#' @param conf a positive number < 1 specifying the level of confidence interval, default is 0.95. Remove C.I. by setting conf=NULL.
#' @param digits integer indicating the number of decimal places \code{round} to be used.
#' @return a \code{data.frame} of phylogenetic diversity table including the sample size, sample coverage,
#' method (rarefaction or extrapolation), and diversity estimates with q = 0, 1, and 2 for the user-specified sample size or sample coverage.
#' @examples
#' data(bird)
#' bird.abu <- bird$abun
#' bird.lab <- rownames(bird$abun)
#' bird.phy <- ade4::newick2phylog(bird$tre)
#' estimatePD(bird.abu, bird.lab, bird.phy, "abundance", base="size", level=NULL, conf=NULL)
#' \dontrun{
#' estimatePD(bird.abu, bird.lab, bird.phy, "abundance", base="size", level=NULL)
#' }
#' @export
estimatePD <- function(x, labels, phy, datatype="abundance", base="size", level=NULL, conf=0.95, digits=4){
  datatype <- check_datatype(datatype)
  site <- m <- NULL
  
  if (!inherits(phy, "phylog")) 
    stop("Non convenient data")
  
  datatype <- check_datatype(datatype)
  
  if(datatype=="incidence_freq" | datatype=="incidence") 
    stop('only support datatype="incidence_raw"')
  
  BASE <- c("size", "coverage")
  if(is.na(pmatch(base, BASE)))
    stop("invalid datatype")
  if(pmatch(base, BASE) == -1)
    stop("invalid datatype")
  base <- match.arg(base, BASE)
  
  if(base=="size"){
    out1 <- invSize(x, datatype, size=level, digits=digits)
  }else if(base=="coverage"){
    out1 <- invChat(x, datatype, C=level, digits=digits)
  }
  size <- sort(unique(out1$m))
  
  q <- c(0, 1, 2)
  Fun <- function(x, q){
    m <- NULL
    se <- ifelse(is.null(conf), FALSE, TRUE)
    if(datatype == "abundance"){
      if(sum(x)==0) stop("Zero abundance counts in one or more sample sites")
      x <- as.numeric(unlist(x))
      tmp <- iNextPD.Ind(x, labels, phy, q, size=size, se=se, conf=conf, nboot=50)
      out <- subset(tmp, m%in%size)
      
    }else if(datatype == "incidence_raw"){
      y <- as.incfreq(x)
      t <- y[1]
      y <- y[-1]
      if(t>sum(y)){
        warning("Insufficient data to provide reliable estimators and associated s.e.") 
      }
      if(sum(y)==0) stop("Zero incidence frequencies in one or more sample sites")
      
      tmp <- iNextPD.Sam(x, labels, phy, q, size=size, se=se, conf=conf, nboot=50)  
      out <- subset(tmp, t%in%size)
    }
    out
  }
  
  if(class(x)=="numeric" | class(x)=="integer" | class(x)=="double"){
    out <- do.call("rbind", lapply(q, function(q) Fun(x, q)))
    out[,-(1:3)] <- round(out[,-(1:3)], digits)
    out <- data.frame("site"="Site1", out)
  }else if(class(x)=="list"){
    if(is.null(names(x))){
      names(x) <- paste0("Site", 1:length(x))
    }
    out <- lapply(x, function(x) {
      tmp <- do.call("rbind", lapply(q, function(q) Fun(x,q)))
      tmp[,-(1:3)] <- round(tmp[,-(1:3)], digits)
      tmp
    })
    out <- do.call(rbind, out)
    out <- data.frame(site=rep(names(x), each=length(size)*length(q)), out)
    if(base=="coverage"){
      tmp2 <- sapply(1:length(size), function(i){
        subset(out, site==unique(out$site)[i], m==size[i])
      })
      out <- do.call(rbind, tmp2)
    }
    rownames(out) <- NULL
  }else if(class(x)=="data.frame" | class(x)=="matrix"){
    if(is.null(colnames(x))){
      colnames(x) <- paste0("Site", 1:ncol(x))
    }
    out <- apply(x, 2, function(x) {
      tmp <- do.call("rbind", lapply(q, function(q) Fun(x,q)))
      tmp[,-(1:3)] <- round(tmp[,-(1:3)], digits)
      tmp
    })
    out <- do.call(rbind, out)
    out <- data.frame(site=rep(names(x), each=length(size)*length(q)), out)
    rownames(out) <- NULL
    
    if(base=="coverage"){
      size <- unique(out$m)
      tmp2 <- lapply(1:length(size), function(i){
        subset(out, out$site%in%unique(out$site)[i] & out$m%in%size[i])
      })
      out <- do.call(rbind, tmp2)
    }
    
   
  }

  if(!is.null(conf)){
    out[,c(1,2,3,4,8,5,6,7)]
  }else{
    out[,c(1,2,3,4,6,5)]
  }
}

# out <- estimatePD(bird.abu,bird.lab,bird.phy, datatype="abundance", level = c(100,200,300.5), base="size")
#out <- estimatePD(bird.abu$North.site,bird.lab,bird.phy, datatype="abundance", level = c(200), base="size")



# -----------------
# 2015-12-27, add transformation function 
# from incidence raw data to incidence frequencies data (iNEXT input format)
# 

###############################################
# Transform incidence raw data to incidence frequencies (iNEXT input format) 
# 
# \code{as.incfreq}: transform incidence raw data (a species by sites presence-absence matrix) to incidence frequencies data (iNEXT input format, a row-sum frequencies vector contains total number of sampling units).
# @param x a \code{data.frame} or \code{matirx} of species by sites presence-absence matrix.
# @return a \code{vector} of species incidence frequencies, the first entry of the input data must be total number of sampling units.
# @examples
# data(plant)
# lapply(plant, as.incfreq)
# @export 
# 
as.incfreq <- function(x){
  if(class(x) == "data.frame" | class(x) == "matrix"){
    a <- sort(unique(c(unlist(x))))
    if(!identical(a, c(0,1))){
      warning("Invalid data type, the element of species by sites presence-absence matrix should be 0 or 1. Set nonzero elements as 1.")
      x <- (x > 0)
    }
    nT <- ncol(x)
    y <- rowSums(x)
    y <- c(nT, y)
    # names(y) <- c("nT", rownames(x))
    y
  }else if(class(x)=="numeric" | class(x)=="integer" | class(x)=="double"){
    warnings("Ambiguous data type, the input object is a vector. Set total number of sampling units as 1.")
    c(1, x) 
  }else{
    stop("Invalid data type, it should be a data.frame or matrix.")
  }
}

###############################################
# Transform abundance raw data to abundance row-sum counts (iNEXT input format) 
# 
# \code{as.abucount}: transform abundance raw data (a species by sites matrix) to abundance rwo-sum counts data (iNEXT input format).
# @param x a \code{data.frame} or \code{matirx} of species by sites matrix.
# @return a \code{vector} of species abundance row-sum counts.
# @examples
# data(plant)
# lapply(plant, as.abucount)
# @export
# 
as.abucount <- function(x){
  if(class(x) == "data.frame" | class(x) == "matrix"){
    y <- rowSums(x)
    y
  }else if(class(x)=="numeric" | class(x)=="integer" | class(x)=="double"){
    warnings("Ambiguous data type, the input object is a vector. Set total number of sampling units as 1.")
    x 
  }else{
    stop("invalid data type, it should be a data.frame or matrix.")
  }
}