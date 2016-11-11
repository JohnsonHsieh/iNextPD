Coverage_ <- function(data, datatype, m){
  # add ckeck_datatype()
  datatype <- check_datatype(datatype)
  ifelse(class(data)=="matrix" || class(data)=="data.frame", type <- "raw", type <- "numeric") 
  ifelse(type == "raw", x <- rowSums(data), x <- data )
  #if(type=="raw" || datatype=='incidence_raw') u <- sum(data)
  
  # checking ambiguous case
  if(datatype=="abundance"){
    # if(sum(x)!=n) warning("argument n is not equal to total sample size.")
    n <- sum(x)
  }else if(datatype=="incidence"){
    # if(max(x)!=n) warning("argument n is not equal to total sampling units.")
    n <- x[1]
    x <- x[-1]
    u <- sum(x)
  }else if(datatype=="incidence_raw"){
    # if(ncol(data)!=n) warning("argument n is not equal to total sampling units.")
    n <- ncol(data)
    u <- sum(data)
  }
  
  f1 <- sum(x == 1)
  f2 <- sum(x == 2)
  f0.hat <- ifelse(f2 == 0, (n - 1) / n * f1 * (f1 - 1) / 2, (n - 1) / n * f1 ^ 2/ 2 / f2)  #estimation of unseen species via Chao1
  A <- ifelse(f1>0, n*f0.hat/(n*f0.hat+f1), 1)
  
  Coverage.Ind <- function(m){
    if(m <= (n-1)) {
      xx <- x[(n-x)>=m]
      out <- 1-sum(xx / n * exp(lgamma(n-xx+1)-lgamma(n-xx-m+1)-lgamma(n)+lgamma(n-m)))
    }else if(m > (n-1) & m <= n){
      xx <- x[(n-x)>=(n-1)]
      a <- 1-sum(xx / n * exp(lgamma(n-xx+1)-lgamma(n-xx-n+2)-lgamma(n)+lgamma(1)))
      b <- 1-f1/n*A
      dat <- data.frame(x=c(n-1, n), y=c(a,b))
      fit <- stats::lm(log(y)~x, dat)
      a[1] <- exp(stats::predict(fit, newdata=data.frame(x=m)))
      out <- a
    }else if(m == n){
      out <- 1-f1/n*A
    }else if(m > n) {
      out <- 1-f1/n*A^(m-n+1)
    }
    out
  }
  
  Coverage.Sam <- function(m){
    if(m <= (n-1)) {
      xx <- x[(n-x)>=m]
      out <- 1-sum(xx / u * exp(lgamma(n-xx+1)-lgamma(n-xx-m+1)-lgamma(n)+lgamma(n-m)))
    }else if(m > (n-1) & m < n){
      xx <- x[(n-x)>=(n-1)]
      a <- 1-sum(xx / u * exp(lgamma(n-xx+1)-lgamma(n-xx-n+2)-lgamma(n)+lgamma(1)))
      b <- 1-f1/u*A
      dat <- data.frame(x=c(n-1, n), y=c(a,b))
      fit <- stats::lm(log(y)~x, dat)
      a[1] <- exp(stats::predict(fit, newdata=data.frame(x=m)))
      out <- a
    }else if(m == n){
      out <- 1-f1/u*A
    }else if(m > n){
      out <- 1-f1/u*A^(m-n+1)
    }
    out
  }
  sapply(m, FUN = function(i){
    if(datatype=="abundance") Coverage.Ind(i)
    else if(datatype=="incidence") Coverage.Sam(i)
    else if(datatype=="incidence_raw") Coverage.Sam(i)
  })		
}


###########################################
#' Compute sample coverage
#' 
#' \code{Coverage}: Compute sample coverage
#' 
#' @param x a vector/matrix/list of species abundances/incidence frequencies or a matrix of incidence table.\cr 
#' If \code{datatype = "incidence_freq"}, then the first entry of the input data must be total number of sampling units, followed by species incidence frequencies.
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}),  
#' sampling-unit-based incidence frequencies data (\code{datatype = "incidence_freq"}) or species by sampling-units incidence matrix (\code{datatype = "incidence_raw"}).
#' @param m an integer vector of sample sizes (number of individuals or sampling units) for which diversity estimates will be computed.
#' @return a data.frame with sample size and sample coverage.
#' @examples 
#' data(bird)
#' Coverage(bird$abun, datatype="abundance", m=c(10,50,100,150,200))
#' @export
Coverage <- function(x, datatype, m){
  # add ckeck_datatype()
  datatype <- check_datatype(datatype)
  
  if(class(x)=="list"){
    tmp <- lapply(x, function(x) Coverage_(x, datatype, m))
    data.frame(m=m, do.call(cbind, tmp))
  }else if(class(x) %in% c("matrix","data.frame") & datatype=="abundance"){
    tmp <- apply(x, 2, function(x) Coverage_(x, datatype, m))
    data.frame(m=m, tmp)
  }else{
    data.frame(m=m, coverage=Coverage_(x, datatype, m))
  }
}