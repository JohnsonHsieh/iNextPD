# ###########################################
#' Estimate asymptotic phylogenetic diversity
#' 
#' \code{estPD}: Estimate asymptotic phylogenetic diversity with order q = 0, 1, 2.
#' 
#' @param x a vector/matrix/list of species abundances or a matrix of raw incidence table.\cr 
#' @param labels a vector of species name for input data.\cr 
#' @param phy a phylogenetic tree with \code{"phylog"} class.\cr 
#' @param q a numeric value specifying the diversity order of Hill number.\cr
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}), species by sampling-units incidence matrix (\code{datatype = "incidence_raw"}).
#' @param se a logical variable to calculate the bootstrap standard error and conf confidence interval.
#' @param conf a positive number < 1 specifying the level of confidence interval, default is 0.95.
#' @return a data.frame with sample size and sample coverage.
#' @examples 
#' # abundance-based example
#' data(bird)
#' bird.abu <- bird$abun
#' bird.lab <- rownames(bird$abun)
#' bird.phy <- ade4::newick2phylog(bird$tre)
#' estPD(bird.abu, labels=bird.lab, phy=bird.phy, q=0, datatype="abundance")
#' 
#' # incidence_based example
#' bird.inc <- bird$inci
# bird.lab <- rownames(bird$inci[[1]])
# bird.phy <- ade4::newick2phylog(bird$tre)
#' estPD(bird.inc, labels=bird.lab, phy=bird.phy, q=0, datatype="incidence_raw", se=TRUE)
#' @export
estPD <- function(x, labels, phy, q=0, datatype="abundance", se=FALSE, conf=0.95){
  
  datatype <- check_datatype(datatype)
  
  if(datatype=="incidence_freq" | datatype=="incidence") 
    stop('only support datatype="incidence_raw"')
  
  if (!is.numeric(conf) || conf > 1 || conf < 0) {
    warning("\"conf\"(confidence level) must be a numerical value between 0 and 1, We use \"conf\" = 0.95 to calculate!")
    conf <- 0.95
  }
  
  z <- qnorm(1 - (1 - conf)/2)
  
  if (!inherits(phy, "phylog")) 
    stop("Non convenient data")
  
  if(!is.numeric(x) & !is.matrix(x) & !is.data.frame(x) & is.integer(x))
    stop("invalid data structure")
  if(is.matrix(x) | is.data.frame(x)){
    if (ncol(x) == 1 & nrow(x) == 1)
      stop("invalid data structure")
  }
  
  
  myFun0 <- function(x,U,L,datatype){
    if(datatype=="abundance"){
      n <- sum(x)
    }else if(datatype=="incidence_raw"){
      y <- x
      x <- as.incfreq(x)
      n <- x[1]
      x <- x[-1]
    }
    L.obs <- L[U>0]
    U.obs <- U[U>0]
    T <- sum(L*U/n)
    f1 <- sum(U.obs==1)
    f2 <- sum(U.obs==2)
    f0.hat <- ceiling((n-1)/n * ifelse(f2 > 0, f1^2/(2*f2), f1*(f1-1)/2))
    g1 <- sum(L.obs[U.obs==1])
    g2 <- sum(L.obs[U.obs==2])
    g0.hat <- ifelse((2*g2*sum(x==1))>g1*sum(x==2),
                     (n-1)/n*g1^2/2/g2,
                     (n-1)/n*g1*(sum(x==1)-1)/2/(sum(x==2)+1)) 
    est <- sum(L.obs) + g0.hat
    obs <- sum(L.obs)
    out <- data.frame(round(t(c(obs, est)),3))
    colnames(out) <- c("Observed", "Estimator")
    if(se==TRUE){
      if(datatype=="abundance"){
        bse <-  PdBootstrapFun(x, labels, phy, q=0, datatype, B=200)
      }else if(datatype=="incidence_raw"){
        bse <-  PdBootstrapFun(y, labels, phy, q=0, datatype, B=200)
      }
      
      CI <- c(max(est - z * bse, obs), est + z * bse)
      out <- round(t(c(obs, est, bse, CI)),3)
      out <- data.frame(out)
      colnames(out) <- c("Observed", "Estimator", "Est_s.e",
                         paste(conf*100, "% Lower", sep=""), paste(conf*100, "% Upper", sep=""))
    }
    out
  }
  
  myFun1 <- function(x,U,L,datatype){
    
    if(datatype=="abundance"){
      n <- sum(x)
    }else if(datatype=="incidence_raw"){
      y <- x
      x <- as.incfreq(x)
      n <- x[1]
      x <- x[-1]
    }
    
    L.obs <- L[U>0]
    U.obs <- U[U>0]
    T <- sum(L*U/n)
    UE <- sum(L.obs/T*U.obs/n*(digamma(n)-digamma(U.obs)), na.rm=1)
    Hn <- sum(-L.obs/T*U.obs/n*log(U.obs/n))
    
    f1 <- sum(U.obs==1)
    g1 <- sum(L.obs[U.obs==1])
    g2 <- sum(L.obs[U.obs==2])
    g0.hat <- ifelse((2*g2*sum(x==1))>g1*sum(x==2),
                     (n-1)/n*g1^2/2/g2,
                     (n-1)/n*g1*(sum(x==1)-1)/2/(sum(x==2)+1)) 
    A <- g1/(n*g0.hat+g1)
    B <- ifelse(A<1,
                g1/T/n*(1-A)^(1-n)*(-log(A)-sum(sapply(1:(n-1), function(r){(1-A)^r/r}))),
                0)
    D.hat <- exp(UE + B)
    est <- D.hat*T
    obs <- T*exp(Hn)
    out <- data.frame(round(t(c(obs, est)),3))
    colnames(out) <- c("Observed", "Estimator")
    if(se==TRUE){
      if(datatype=="abundance"){
        bse <-  PdBootstrapFun(x, labels, phy, q=1, datatype, B=200)
      }else if(datatype=="incidence_raw"){
        bse <-  PdBootstrapFun(y, labels, phy, q=1, datatype, B=200)
      }
      CI <- c(max(est - z * bse, obs), est + z * bse)
      out <- round(t(c(obs, est, bse, CI)),3)
      out <- data.frame(out)
      colnames(out) <- c("Observed", "Estimator", "Est_s.e",
                         paste(conf*100, "% Lower", sep=""), paste(conf*100, "% Upper", sep=""))
    }
    out
  }
  
  myFun2 <- function(x,U,L,datatype){
    
    if(datatype=="abundance"){
      n <- sum(x)
    }else if(datatype=="incidence_raw"){
      y <- x
      x <- as.incfreq(x)
      n <- x[1]
      x <- x[-1]
    }
    
    L.obs <- L[U>0]
    U.obs <- U[U>0]
    T <- sum(L*U/n)
    lam <- sum(L.obs*(U.obs*(U.obs-1)/n/(n-1)))
    est <- T/(lam/T)
    obs <- T/sum(L.obs*(U.obs^2/n^2)/T)
    out <- data.frame(round(t(c(obs, est)),3))
    colnames(out) <- c("Observed", "Estimator")
    if(se==TRUE){
      if(datatype=="abundance"){
        bse <-  PdBootstrapFun(x, labels, phy, q=2, datatype, B=200)
      }else if(datatype=="incidence_raw"){
        bse <-  PdBootstrapFun(y, labels, phy, q=2, datatype, B=200)
      }
      
      CI <- c(max(est - z * bse, obs), est + z * bse)
      out <- round(t(c(obs, est, bse, CI)),3)
      out <- data.frame(out)
      colnames(out) <- c("Observed", "Estimator", "Est_s.e",
                         paste(conf*100, "% Lower", sep=""), paste(conf*100, "% Upper", sep=""))
    }
    out
    
  }
  
  myFun <- function(x,labels, phy, datatype){
    
    tmp <- ExpandData_(x, labels, phy, datatype)  
    U <- tmp$branch_abun
    L <- tmp$branch_length
    
    if(q==0) myFun0(x, U, L, datatype)
    else if(q==1) myFun1(x, U, L, datatype)
    else if(q==2) myFun2(x, U, L, datatype)
  }
  
  
  if(class(x) == "numeric" | class(x) == "integer"){
    out <- myFun(x,labels, phy, datatype)
  }else if(class(x) == "list"){
    out <- do.call("rbind", lapply(x, myFun, labels, phy, datatype))
  } else if(class(x) == "matrix" | class(x) == "data.frame"){
    if(datatype=="abundance"){
      out <- do.call("rbind", apply(as.matrix(x), 2, myFun,labels, phy, datatype))  
    }else if(datatype=="incidence_raw"){
      out <- myFun(x,labels, phy, datatype)
    }
  }
  return(out)
  
}

# sub function to build estPD bootstrap s.e.
estPDsub <- function(n, x, U, L, q=0, datatype="abundance"){
  datatype <- check_datatype(datatype)
  
  if(datatype=="incidence_freq" | datatype=="incidence") 
    stop('only support datatype="incidence_raw"')
  
  if(datatype=="abundance"){
    L.obs <- L[U>0]
    U.obs <- U[U>0]
    T <- sum(L.obs*U.obs/n)
    
    if(q==0){
      g1 <- sum(L.obs[U.obs==1])
      g2 <- sum(L.obs[U.obs==2])
      g0.hat <- ifelse((2*g2*sum(x==1))>(g1*sum(x==2)),
                       (n-1)/n*g1^2/2/g2,
                       (n-1)/n*g1*(sum(x==1)-1)/2/(sum(x==2)+1))  
      PD.hat <- sum(L.obs) + g0.hat
      PD.hat
    }else if(q==1){
      UE <- sum(L.obs/T*U.obs/n*(digamma(n)-digamma(U.obs)), na.rm=1)
      Hn <- sum(-L.obs/T*U.obs/n*log(U.obs/n))
      f1 <- sum(U.obs==1)
      g1 <- sum(L.obs[U.obs==1])
      g2 <- sum(L.obs[U.obs==2])
      g0.hat <- ifelse((2*g2*sum(x==1))>(g1*sum(x==2)),
                       (n-1)/n*g1^2/2/g2,
                       (n-1)/n*g1*(sum(x==1)-1)/2/(sum(x==2)+1)) 
      A <- g1/(n*g0.hat+g1)
      B <- ifelse(A<1,
                  g1/T/n*(1-A)^(1-n)*(-log(A)-sum(sapply(1:(n-1), function(r){(1-A)^r/r}))),
                  0)
      D.hat <- exp(UE + B)
      c(D.hat*T)
    }else if(q==2){
      lam <- sum(L.obs*(U.obs*(U.obs-1)/n/(n-1)))
      T/(lam/T)
    }else{
      warning("Not support")
    }
  }else if(datatype=="incidence_raw"){
    x <- as.incfreq(x)
    n <- x[1]
    x <- x[-1]
    
    L.obs <- L[U>0]
    U.obs <- U[U>0]
    T <- sum(L.obs*U.obs/n)
    
    if(q==0){
      g1 <- sum(L.obs[U.obs==1])
      g2 <- sum(L.obs[U.obs==2])
      g0.hat <- ifelse((2*g2*sum(x==1))>(g1*sum(x==2)),
                       (n-1)/n*g1^2/2/g2,
                       (n-1)/n*g1*(sum(x==1)-1)/2/(sum(x==2)+1))  
      PD.hat <- sum(L.obs) + g0.hat
      PD.hat
    }else if(q==1){
      UE <- sum(L.obs/T*U.obs/n*(digamma(n)-digamma(U.obs)), na.rm=1)
      Hn <- sum(-L.obs/T*U.obs/n*log(U.obs/n))
      f1 <- sum(U.obs==1)
      g1 <- sum(L.obs[U.obs==1])
      g2 <- sum(L.obs[U.obs==2])
      g0.hat <- ifelse((2*g2*sum(x==1))>(g1*sum(x==2)),
                       (n-1)/n*g1^2/2/g2,
                       (n-1)/n*g1*(sum(x==1)-1)/2/(sum(x==2)+1)) 
      A <- g1/(n*g0.hat+g1)
      # B <- ifelse(A<1,
      #             f1/n*(1-A)^(1-n)*(-log(A)-sum(sapply(1:(n-1), function(r){(1-A)^r/r}))),
      #             0)
      # 
      B <- ifelse(A<1,
                  g1/T/n*(1-A)^(1-n)*(-log(A)-sum(sapply(1:(n-1), function(r){(1-A)^r/r}))),
                  0)
      D.hat <- exp(UE + B)
      c(D.hat*T)
    }else if(q==2){
      lam <- sum(L.obs*(U.obs*(U.obs-1)/n/(n-1)))
      T/(lam/T)
    }else{
      warning("Not support")
    }
  }
}

# using bootstrap assemblage to compute se of qPD
PdBootstrapFun <-  function(x, labels, phy, q, datatype="abundance", B){
  if (!inherits(phy, "phylog")) 
    stop("Non convenient data")
  
  datatype <- check_datatype(datatype)
  
  if(datatype=="incidence_freq" | datatype=="incidence") 
    stop('only support datatype="incidence_raw"')
  
  if(datatype=="abundance"){
    n <- sum(x)
    BootComm <- PDBoot_(x, labels, phy, "abundance")
    Prob.hat <- BootComm$branch_abun
    L.hat <- BootComm$branch_length
    U.Mat <- t(sapply(Prob.hat, function(p) rbinom(B, n, p)))
    se <- sd(apply(U.Mat, 2, function(U) estPDsub(n, x, U, L.hat, q, datatype)), na.rm = TRUE)
    se
  }else if(datatype=="incidence_raw"){
    t <- as.incfreq(x)[1]
    BootComm <- PDBoot_(x, labels, phy, "incidence_raw")
    Prob.hat <- BootComm$branch_abun
    L.hat <- BootComm$branch_length
    U.Mat <- t(sapply(Prob.hat, function(p) rbinom(B, t, p)))
    se <- sd(apply(U.Mat, 2, function(U) estPDsub(n, x, U, L.hat, q, datatype)), na.rm = TRUE)
    se
  }
}