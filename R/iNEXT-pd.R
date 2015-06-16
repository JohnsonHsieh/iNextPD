DataInfo <- function(x, datatype="abundance"){
  
  Fun.abun <- function(x) {
    n <- sum(x)
    fk <- sapply(1:10, function(k) sum(x == k))
    f1 <- fk[1]
    f2 <- fk[2]
    Sobs <- sum(x > 0)
    f0.hat <- ifelse(f2 == 0, (n - 1)/n * f1 * (f1 - 1)/2, 
                     (n - 1)/n * f1^2/2/f2)
    A <- ifelse(f1 > 0, n * f0.hat/(n * f0.hat + f1), 1)
    Chat <- round(1 - f1/n * A, 4)
    c(n, Sobs, Chat, fk)
  }
  
  if(class(x) == "numeric" | class(x) == "integer"){
    out <- matrix(Fun.abun(x), nrow=1)
  }else if(class(x) == "list"){
    out <- do.call("rbind", lapply(x, Fun.abun))
  } else if(class(x) == "matrix" | class(x) == "data.frame"){
    out <- t(apply(as.matrix(x), 2, Fun.abun))  
  }
  colnames(out) <-  c("n", "S.obs", "C.hat", paste("f",1:10, sep=""))
  as.data.frame(out)
}

ViewTree <- function(abun, labels, phy){
<<<<<<< HEAD
  require(ade4)
=======
  library(ade4)
>>>>>>> 7958e4ab2fd7dba6eda75cecc2110aafda24ae6d
  if (!inherits(phy, "phylog")) 
    stop("Non convenient data")
  
  n <- sum(abun)
  names(abun) <- labels
  a <- abun[names(phy$leaves)]
  
  
  for(i in 1:length(phy$parts)){
    a[1+length(a)] <- sum(a[phy$parts[[i]]])
    names(a)[length(a)] <- names(phy$parts)[i]
  }
  out <- data.frame("branch_abun"=a, "branch_length"=c(phy$leaves, phy$nodes))
  out
}



Stirling2 <- function(n,m)
{
  ## Purpose:  Stirling Numbers of the 2-nd kind
  ##     S^{(m)}_n = number of ways of partitioning a set of
  ##                      $n$ elements into $m$ non-empty subsets
  ## Author: Martin Maechler, Date:  May 28 1992, 23:42
  ## ----------------------------------------------------------------
  ## Abramowitz/Stegun: 24,1,4 (p. 824-5 ; Table 24.4, p.835)
  ## Closed Form : p.824 "C."
  ## ----------------------------------------------------------------
  
  if (0 > m || m > n) stop("'m' must be in 0..n !")
  k <- 0:m
  sig <- rep(c(1,-1)*(-1)^m, length= m+1)# 1 for m=0; -1 1 (m=1)
  ## The following gives rounding errors for (25,5) :
  ## r <- sum( sig * k^n /(gamma(k+1)*gamma(m+1-k)) )
  ga <- gamma(k+1)
  round(sum( sig * k^n /(ga * rev(ga))))
}


PdBootstrapFun <- function(x, labels, phy, FunName, q, datatype="abundance", B){
  if(!is.numeric(x) & !is.matrix(x) & !is.data.frame(x))
    stop("invalid data structure")
  if(is.matrix(x) | is.data.frame(x)){
    if (ncol(x) != 1 & nrow(x) != 1)
      stop("invalid data structure")
  }
  n <- sum(x)
  BootComm <- SpBoot(x, labels, phy)
  Prob.hat <- BootComm$branch_abun
  L.hat <- BootComm$branch_length
  U.Mat <- t(sapply(Prob.hat, function(p) rbinom(B, n, p)))
  se <- sd(apply(U.Mat, 2, function(U) FunName(n=n, x=x, U=U, L=L.hat, q=q)))
  se
}

estPDsub <- function(n, x, U, L, q=0){
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
                f1/n*(1-A)^(1-n)*(-log(A)-sum(sapply(1:(n-1), function(r){(1-A)^r/r}))),
                0)
    D.hat <- exp(UE + B)
    D.hat*T
  }else if(q==2){
    lam <- sum(L.obs*(U.obs*(U.obs-1)/n/(n-1)))
    T/(lam/T)
  }else{
    warning("Not support")
  }
}


estPD <- function(x, labels, phy, q=0, datatype="abundance", conf=0.95){
  TYPE <- c("abundance", "incidence")
  if(is.na(pmatch(datatype, TYPE)))
    stop("invalid data type")
  if(pmatch(datatype, TYPE) == -1)
    stop("ambiguous data type")
  datatype <- match.arg(datatype, TYPE)
  
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
  
  
  myFun0 <- function(x,U,L){
<<<<<<< HEAD
    n <- sum(x)
=======
>>>>>>> 7958e4ab2fd7dba6eda75cecc2110aafda24ae6d
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
    se <-  PdBootstrapFun(x, labels, phy, estPDsub, q=0, B=200)
    CI <- c(max(est - z * se, obs), est + z * se)
    out <- round(t(c(obs, est, se, CI)),3)
    out <- data.frame(out)
    colnames(out) <- c("Observed", "Estimator", "Est_s.e",
                       paste(conf*100, "% Lower", sep=""), paste(conf*100, "% Upper", sep=""))
    out
  }
  
  myFun1 <- function(x,U,L){
<<<<<<< HEAD
    n <- sum(x)
=======
>>>>>>> 7958e4ab2fd7dba6eda75cecc2110aafda24ae6d
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
                f1/n*(1-A)^(1-n)*(-log(A)-sum(sapply(1:(n-1), function(r){(1-A)^r/r}))),
                0)
    D.hat <- exp(UE + B)
    est <- D.hat*T
    obs <- T*exp(Hn)
    se <-  PdBootstrapFun(x, labels, phy, estPDsub, q=0, B=200)
    CI <- c(max(est - z * se, obs), est + z * se)
    out <- round(t(c(obs, est, se, CI)),3)
    out <- data.frame(out)
    colnames(out) <- c("Observed", "Estimator", "Est_s.e",
                       paste(conf*100, "% Lower", sep=""), paste(conf*100, "% Upper", sep=""))
    out
  }
  
  myFun2 <- function(x,U,L){
<<<<<<< HEAD
    n <- sum(x)
=======
>>>>>>> 7958e4ab2fd7dba6eda75cecc2110aafda24ae6d
    L.obs <- L[U>0]
    U.obs <- U[U>0]
    T <- sum(L*U/n)
    lam <- sum(L.obs*(U.obs*(U.obs-1)/n/(n-1)))
    est <- T/(lam/T)
    obs <- T/sum(L.obs*(U.obs^2/n^2)/T)
    se <-  PdBootstrapFun(x, labels, phy, estPDsub, q=0, B=200)
    CI <- c(max(est - z * se, obs), est + z * se)
    out <- round(t(c(obs, est, se, CI)),3)
    out <- data.frame(out)
    colnames(out) <- c("Observed", "Estimator", "Est_s.e",
                       paste(conf*100, "% Lower", sep=""), paste(conf*100, "% Upper", sep=""))
    out
    
  }
  
  myFun <- function(x){
    
    n <- sum(x)        #sample size  
    names(x) <- labels
    a <- x[names(phy$leaves)]  
    for(i in 1:length(phy$parts)){
      a[1+length(a)] <- sum(a[phy$parts[[i]]])
      names(a)[length(a)] <- names(phy$parts)[i]
    }
    
    tmp <- data.frame("branch_abun"=a, "branch_length"=c(phy$leaves, phy$nodes))
    U <- tmp$branch_abun
    L <- tmp$branch_length
    T <- sum(L*U/n)
    
    if(q==0) myFun0(x, U, L)
    else if(q==1) myFun1(x, U, L)
    else if(q==2) myFun2(x, U, L)
  }
  
  
  if(class(x) == "numeric"){
    out <- myFun(x)
  }else if(class(x) == "integer"){
    out <- myFun(x)
  }else if(class(x) == "list"){
    out <- do.call("rbind", lapply(x, myFun))
  } else if(class(x) == "matrix" | class(x) == "data.frame"){
    out <- do.call("rbind", apply(as.matrix(x), 2, myFun))  
  }
  return(out)
  
}


PhD.m.mle <- function(p, phy, q, m){
  if (!inherits(phy, "phylog")) 
    stop("Non convenient data")
  names(p) <- names(phy$leaves)
  a <- p[names(phy$leaves)]
  
  for(i in 1:length(phy$parts)){
    a <- c(a, sum(a[phy$parts[[i]]]))
    names(a)[length(a)] <- names(phy$parts)[i]
  }
  L <- c(phy$leaves, phy$nodes)
  T <- sum(L*a)
  sub <- function(m){
    if(q==0) sum(L*(1-(1-a)^m)) 
    else if(q==1) {
      a <- a[a>0]
      sub <- function(k) {
        sum(L*choose(m,k)*a^k*(1-a)^(m-k))
      }
      k <- 1:m
      gk <- sapply(k, sub)
      Hm <- sum(-k/m/T*log(k/m/T)*gk)
      exp(Hm)
    } else if(q==2) {
      lam <- sum(L*((m-1)/m*a^2+1/m*a))
      T/(lam/T)
    } else{
      sub <- function(x, j){
        a1 <- Stirling2(q, j) * exp(lgamma(m+1) - lgamma(m-j+1) - q*log(m))
        b1 <- x^j
        a1*b1
      }
      tmp <- rep(0, length(a))
      for(i in 1:length(a)){
        tmp[i] <- L[i] * sum(sapply(1:q, function(j) sub(a[i], j)))
      }
      lam <- sum(tmp)
      T*(lam/T)^(1/(1-q))
    }
  }
  sapply(m, sub)
}

PhD.m <- function(n, x, U, L, q, m){
  L.obs <- L[U>0]
  U.obs <- U[U>0]
  T <- sum(L.obs*U.obs/n)
  
  if(q == 0){
    Fun0 <- function(m){
      PD.obs <- sum(L.obs)
      if(m < n){
        ga <- sum(L.obs * exp(lchoose(n - U.obs, m) - lchoose(n, m)))
        PD.obs - ga
      } else {
        #         f1 <- sum(U.obs==1)
        #         f2 <- sum(U.obs==2)
        #         g1 <- sum(L.obs[U.obs==1])
        #         g2 <- sum(L.obs[U.obs==2])
        #         g0 <- (n-1)/n * ifelse((2*g2*f1) > (g1*f2), g1^2/(2*g2), g1*(f1-1)/(2*(f2+1)))
        g1 <- sum(L.obs[U.obs==1])
        g2 <- sum(L.obs[U.obs==2])
        g0.hat <- ifelse((2*g2*sum(x==1))>g1*sum(x==2),
                         (n-1)/n*g1^2/2/g2,
                         (n-1)/n*g1*(sum(x==1)-1)/2/(sum(x==2)+1)) 
        
        A <- g1/(n*g0.hat+g1)
        PD.obs + g0.hat*(1 - (1-A)^(m-n))
      }
    }
    sapply(m, Fun0)
  } else if(q == 1){
    Fun1 <- function(m){
      if(m < n){
        sub <- function(k) sum(L.obs*exp(lchoose(U.obs,k)+lchoose(n-U.obs,m-k)-lchoose(n,m)))
        k <- 1:m
        gk <- sapply(k, sub)
        Hm <- sum(-k/m/T*log(k/m/T)*gk)
        exp(Hm)
      }else if(m==n){
        exp(sum(-L.obs*U.obs/n/T*log(U.obs/T/n)))
      } else {
        UE <- sum(L.obs/T*U.obs/n*(digamma(n)-digamma(U.obs)), na.rm=1)
        Hn <- sum(-L.obs/T*U.obs/n*log(U.obs/n))
        g1 <- sum(L.obs[U.obs==1])
        g2 <- sum(L.obs[U.obs==2])
        g0.hat <- ifelse((2*g2*sum(x==1))>g1*sum(x==2),
                         (n-1)/n*g1^2/2/g2,
                         (n-1)/n*g1*(sum(x==1)-1)/2/(sum(x==2)+1)) 
        
        f1 <- sum(U.obs==1)
        A <- g1/(n*g0.hat+g1)        
        B <- ifelse(A<1,
                    f1/n*(1-A)^(1-n)*(-log(A)-sum(sapply(1:(n-1), function(r){(1-A)^r/r}))),
                    0)
        D.hat <- exp(UE + B)
        est <- T*D.hat
        obs <- T*exp(Hn)
        
        sub <- function(k) sum(L.obs*exp(lchoose(U.obs,k)+lchoose(n-U.obs,n-1-k)-lchoose(n,n-1)))
        k <- 1:(n-1)
        gk <- sapply(k, sub)
        Dm <- T*exp(sum(-k/(n-1)/T*log(k/(n-1))*gk))
        A <- (obs-Dm)/(est-obs)
        
        ifelse(A!=0, obs + (est-obs)*(1-(1-A)^(m-n)), obs)
        
        #w <- (m-n)/m
        #T*exp(w*(UE+B)+(1-w)*Hn)
      }
    }
    sapply(m, Fun1)
  } else if(q == 2){
    Fun2 <- function(m){
      lam <- sum(L.obs*((m-1)/m*U.obs*(U.obs-1)/n/(n-1)+1/m*U.obs/n))
      T/(lam/T)
    }
    sapply(m, Fun2)
  } else{
    Funq <- function(m){
      
      sub <- function(x, j){
        a <- Stirling2(q, j) * exp(lgamma(m+1) - lgamma(m-j+1) - q*log(m))
        if(x > j) {
          b <- exp(lgamma(x+1) - lgamma(x-j+1) - lgamma(n+1) + lgamma(n-j+1))
        } else {b <- 0}
        a*b
      }
      
      tmp <- rep(0, length(U.obs))
      for(i in 1:length(U.obs)){
        tmp[i] <- L.obs[i] * sum(sapply(1:q, function(j) sub(U.obs[i], j)))
      }
      lam <- sum(tmp)
      T*(lam/T)^(1/(1-q))
    }
    sapply(m, Funq)
  }
}

Cvg.m <- function(x, m){
  x <- x[x>0]
  n <- sum(x)
  f1 <- sum(x == 1)
  f2 <- sum(x == 2)
  f0.hat <- ifelse(f2 == 0, (n - 1) / n * f1 * (f1 - 1) / 2, (n - 1) / n * f1 ^ 2/ 2 / f2)  #estimation of unseen species via Chao1
  A <- ifelse(f1>0, n*f0.hat/(n*f0.hat+f1), 1)
  Sub <- function(m){
    if(m < n) out <- 1-sum(x / n * exp(lchoose(n - x, m)-lchoose(n - 1, m)))
    if(m == n) out <- 1-f1/n*A
    if(m > n) out <- 1-f1/n*A^(m-n+1)
    out
  }
  sapply(m, Sub)
}

EstiBootComm.Ind <- function(Spec)
{
  Sobs <- sum(Spec > 0)   #observed species
  n <- sum(Spec)  	  	#sample size
  f1 <- sum(Spec == 1) 	#singleton 
  f2 <- sum(Spec == 2) 	#doubleton
  a <- ifelse(f1 == 0, 0, (n - 1) * f1 / ((n - 1) * f1 + 2 * f2) * f1 / n)
  b <- sum(Spec / n * (1 - Spec / n) ^ n)
  w <- a / b  			#adjusted factor for rare species in the sample
  f0.hat <- ceiling(ifelse(f2 == 0, (n - 1) / n * f1 * (f1 - 1) / 2, (n - 1) / n * f1 ^ 2/ 2 / f2))	#estimation of unseen species via Chao1
  Prob.hat <- Spec / n * (1 - w * (1 - Spec / n) ^ n)					#estimation of relative abundance of observed species in the sample
  Prob.hat.Unse <- rep(2 * f2/((n - 1) * f1 + 2 * f2), f0.hat)		#estimation of relative abundance of unseen species in the sample
  return(c(Prob.hat, Prob.hat.Unse))									#Output: a vector of estimated relative abundance
}

SpBoot <- function(abun, labels, phy){
  n <- sum(abun)        #sample size
  if (!inherits(phy, "phylog")) 
    stop("Non convenient data")
  names(abun) <- labels
  a <- abun[names(phy$leaves)]
  
  f1 <- sum(a==1)
  f2 <- sum(a==2)
  f0.hat <- (n-1)/n * ifelse(f2 > 0, f1^2/(2*f2), f1*(f1-1)/2)
  A <- n*f0.hat/(n*f0.hat+f1)
  a1 <- f1/n*A
  b1 <- sum(a / n * (1 - a / n) ^ n)
  w <- a1 / b1    
  aa <- a/n*(1-w*(1-a/n)^n)
  
  for(i in 1:length(phy$parts)){
    a[1+length(a)] <- sum(a[phy$parts[[i]]])
    names(a)[length(a)] <- names(phy$parts)[i]
    
    aa[1+length(aa)] <- sum(aa[phy$parts[[i]]])
    names(aa)[length(aa)] <- names(phy$parts)[i]
  }
  tmp <- data.frame("branch_abun"=a, "branch_length"=c(phy$leaves, phy$nodes))
  U <- tmp$branch_abun
  L <- tmp$branch_length
  U.obs <- U[U>0]
  L.obs <- L[U>0]
  a.obs <- aa[U>0]
  
  a.und <- rep(a1/ceiling(f0.hat), times = ceiling(f0.hat))
  
  g1 <- sum(L.obs[U.obs==1])
  g2 <- sum(L.obs[U.obs==2])
  g0 <- (n-1)/n * ifelse((2*g2*f1) > (g1*f2), g1^2/(2*g2), g1*(f1-1)/(2*(f2+1)))
  L0 <-  rep(g0/ceiling(f0.hat), ceiling(f0.hat))
  
  ai <- c(a.obs, a.und)
  Li <- c(L.obs, L0)
  data.frame("branch_abun"=ai, "branch_length"=Li)
}


phy2abun <- function(phy, abun){
  a <- abun[names(phy$leaves)]  
  for(i in 1:length(phy$parts)){
    a[1+length(a)] <- sum(a[phy$parts[[i]]])
    names(a)[length(a)] <- names(phy$parts)[i]
  }
  tmp <- data.frame("branch_abun"=a, "branch_length"=c(phy$leaves, phy$nodes))
  tmp
}


PhBoot <- function(phy, abun){
  n <- sum(abun)        #sample size
  if (!inherits(phy, "phylog")) 
    stop("Non convenient data")
  a <- abun[names(phy$leaves)]  
  for(i in 1:length(phy$parts)){
    a[1+length(a)] <- sum(a[phy$parts[[i]]])
    names(a)[length(a)] <- names(phy$parts)[i]
  }
  tmp <- data.frame("branch_abun"=a, "branch_length"=c(phy$leaves, phy$nodes))
  U <- tmp$branch_abun
  L <- tmp$branch_length
  T <- sum(L*U/n)
  
  L.obs <- L[U>0]
  U.obs <- U[U>0]
  f1 <- sum(U.obs==1)
  f2 <- sum(U.obs==2)
  f0.hat <- ceiling((n-1)/n * ifelse(f2 > 0, f1^2/(2*f2), f1*(f1-1)/2))
  g1 <- sum(L.obs[U.obs==1]) / T
  g2 <- sum(L.obs[U.obs==2]) / T
  g2 <- ifelse(g2>0, g2, g1/f1)
  a1.hat <- 2*g2/((n-1)*g1+2*g2)
  A <- g1/n*(1-a1.hat)
  w <- A / sum(L.obs / T * U.obs / n * (1 - U.obs / n) ^ n)
  B <- f1/n*(n-1)*f1/((n-1)*f1+2*max(f2,1))
  L0 <- T * A / B
  a0 <- ifelse(f0.hat>0, B / f0.hat, 0)
  Prob.hat <- U.obs / n * (1 - w * (1 - U.obs / n) ^ n)      		#estimation of relative abundance of observed species in the sample
  Prob.hat.Unse <- rep(a0, f0.hat)		#estimation of relative abundance of unseen species in the sample
  ai <- c(Prob.hat, Prob.hat.Unse)
  Li <- c(L.obs, rep(L0, f0.hat))
  return(data.frame("branch_abun"=ai, "branch_length"=Li))	
}


iNextPD <- function(abun, labels, phy, size=NULL, endpoint=2*sum(abun), knots=40, se=FALSE, nboot=50){
  
  require(ade4)
  if (!inherits(phy, "phylog")) 
    stop("Non convenient data")
  datatype <- "abundance"  
  q <- 0
  
  Fun <- function(x, q){
    n <- sum(x)      	#sample size
    
    if(is.null(size)) {
      if(endpoint <= n) {
        m <- floor(seq(1, endpoint, length=floor(knots)))
      } else {
        m <- c(floor(seq(1, sum(x)-1, length.out=floor(knots/2)-1)), sum(x), floor(seq(sum(x)+1, to=endpoint, length.out=floor(knots/2))))
      }
      m <- c(1, m[-1])
    } else if(is.null(size)==FALSE) {  
      if(max(m)>n & length(m[m==n])==0)  m <- c(m, n, n+1)
      m <- sort(m)
    }
    
    names(x) <- labels
    a <- x[names(phy$leaves)]
    
    for(i in 1:length(phy$parts)){
      a[1+length(a)] <- sum(a[phy$parts[[i]]])
      names(a)[length(a)] <- names(phy$parts)[i]
    }
    
    tmp <- data.frame("branch_abun"=a, "branch_length"=c(phy$leaves, phy$nodes))
    U <- tmp$branch_abun
    L <- tmp$branch_length
    T <- sum(L*U/n)
    
    phd.hat <- PhD.m(n, x, U, L, q, m)
    cvg.hat <- Cvg.m(x, m)
    
    if(se==TRUE & nboot > 0 & length(x) > 1){
      BootComm <- SpBoot(abun = x, labels = labels, phy=phy)
      Prob.hat <- BootComm$branch_abun
      L.hat <- BootComm$branch_length
      U.Mat <- t(sapply(Prob.hat, function(p) rbinom(nboot, n, p)))
      
      p.hat <- EstiBootComm.Ind(x)
      X.Mat <- rmultinom(nboot, n, p.hat)
      
      error <-  qnorm(0.975) * apply(apply(U.Mat, 2, function(U) PhD.m(n, x, U, L.hat, q, m)), 1, sd, na.rm=TRUE)
      left  <- phd.hat - error
      right <- phd.hat + error
      
      error.C <-  qnorm(0.975) * apply(apply(X.Mat, 2, function(x) Cvg.m(x, m)), 1, sd, na.rm=TRUE)
      left.C  <- cvg.hat - error.C
      right.C <- cvg.hat + error.C
      
      right.C[right.C>=1] <- 1
      left.C[left.C<0] <- 0
      
      # left.C <- apply(apply(X.Mat, 2, function(x) Cvg.m(x, m)), 1, quantile, na.rm=TRUE, probs=0.025)
      # right.C <- apply(apply(X.Mat, 2, function(x) Cvg.m(x, m)), 1, quantile, na.rm=TRUE, probs=0.975)
      out <- data.frame("m"=m, "qPD"=phd.hat, "qPD.95.LCL"=left, "qPD.95.UCL"=right, "SC"=cvg.hat, "SC.95.LCL"=left.C, "SC.95.UCL"=right.C)
    }else{
      out <- data.frame("m"=m, "qPD"=phd.hat, "SC"=cvg.hat)
    }
    out$method <- ifelse(out$m<n, "interpolated", ifelse(out$m==n, "observed", "extrapolated"))
    out$order <- q
    id <- match(c("m", "method", "order", "qPD", "qPD.95.LCL", "qPD.95.UCL", "SC", "SC.95.LCL", "SC.95.UCL"), names(out), nomatch = 0)
    out <- out[, id]
    out
  }
  
  if(class(q) != "numeric")
    stop("invlid class of order q, q should be a postive value/vector of numeric object")
  if(min(q) < 0){
    warning("ambigous of order q, we only compute postive q")
    q <- q[q >= 0]
  }
  
  
  x <- abun
  if (class(x) == "numeric" | class(x) == "integer" | class(x) == 
        "double") {
    out <- do.call("rbind", lapply(q, function(q) Fun(x, q)))
    out[, -(1:3)] <- round(out[, -(1:3)], 3)
    index <- rbind(as.matrix(estPD(x, labels, phy, q=0)), 
                   as.matrix(estPD(x, labels, phy, q=1)), 
                   as.matrix(estPD(x, labels, phy, q=2)))
    rownames(index) <- c("q = 0", "q = 1", "q = 2")
    tree <- ViewTree(x, labels, phy)
  }
  else if (class(x) == "matrix" | class(x) == "data.frame") {
    out <- apply(as.matrix(x), 2, function(x) {
      tmp <- do.call("rbind", lapply(q, function(q) Fun(x, q)))
      tmp[, -(1:3)] <- round(tmp[, -(1:3)], 3)
      tmp
    })
    arr <- array(0, dim = c(3, 5, ncol(x)))
    arr[1, , ] <- t(as.matrix(estPD(x, labels, phy, q=0)))
    arr[2, , ] <- t(as.matrix(estPD(x, labels, phy, q=1)))
    arr[3, , ] <- t(as.matrix(estPD(x, labels, phy, q=2)))
    dimnames(arr)[[3]] <- colnames(x)
    dimnames(arr)[[1]] <- c("q = 0", "q = 1", "q = 2")
    dimnames(arr)[[2]] <- c("Observed", "Estimator", "Est_s.e.", 
                            "95% Lower", "95% Upper")
    index <- ftable(arr, row.vars = c(3, 1))
    tmp <- apply(abun, 2, ViewTree, labels, phy)
    tree <- data.frame(do.call(data.frame,lapply(tmp, "[",1)), tmp[[1]][,2])
    colnames(tree) <- c(paste("branch_abun_",names(tmp),sep=""), "branch_length")
  }
  else if (class(x) == "list") {
    out <- lapply(x, function(x) {
      tmp <- do.call("rbind", lapply(q, function(q) Fun(x, q)))
      tmp[, -(1:3)] <- round(tmp[, -(1:3)], 3)
      tmp
    })
    arr <- array(0, dim = c(3, 5, length(x)))
    arr[1, , ] <- t(as.matrix(estPD(x, labels, phy, q=0)))
    arr[2, , ] <- t(as.matrix(estPD(x, labels, phy, q=1)))
    arr[3, , ] <- t(as.matrix(estPD(x, labels, phy, q=2)))
    dimnames(arr)[[3]] <- names(x)
    dimnames(arr)[[1]] <- c("q = 0", "q = 1", "q = 2")
    dimnames(arr)[[2]] <- c("Observed", "Estimator", "Est_s.e.", 
                            "95% Lower", "95% Upper")
    index <- ftable(arr, row.vars = c(3, 1))
    
    tmp <- lapply(abun, ViewTree, labels, phy)
    tree <- data.frame(do.call(data.frame,lapply(tmp, "[",1)), tmp[[1]][,2])
    colnames(tree) <- c(paste("branch_abun_",names(tmp),sep=""), "branch_length")
    
  }
  else {
    stop("invlid class of x, x should be a object of numeric, matrix, data.frame, or list")
  }
  
  
  
  z <- list("DataInfo"=DataInfo(abun, datatype = "abundance"), 
            "ViewTree"=tree, 
            "BasicIndex"=index,
            "Accumulation"=out)
  
  class(z) <- c("iNEXT")      
  return(z)
}


ggiNEXT <- function (x, type = 1, se = TRUE, facet.var = "none", color.var = "order") {
  require("ggplot2")
  if (class(x) != "iNEXT") 
    stop("invalid object class")
  TYPE <- c(1, 2, 3)
  SPLIT <- c("none", "order", "site", "both")
  if (is.na(pmatch(type, TYPE)) | pmatch(type, TYPE) == -1) 
    stop("invalid plot type")
  if (is.na(pmatch(facet.var, SPLIT)) | pmatch(facet.var, SPLIT) == 
        -1) 
    stop("invalid facet variable")
  if (is.na(pmatch(color.var, SPLIT)) | pmatch(color.var, SPLIT) == 
        -1) 
    stop("invalid color variable")
  type <- pmatch(type, 1:3)
  facet.var <- match.arg(facet.var, SPLIT)
  color.var <- match.arg(color.var, SPLIT)
  if (facet.var == "order") 
    color.var <- "site"
  if (facet.var == "site") 
    color.var <- "order"
  y <- method <- site <- y.lwr <- y.upr <- NULL
  site <<- NULL
  z <- x$Accumulation
  if (class(z) == "list") {
    z <- data.frame(do.call("rbind", z), site = rep(names(z), sapply(z, nrow)))
    rownames(z) <- NULL
  }else{
    z$site <- "site"
  }
  
  if(names(z)[4]!="qD") is.phylo <-TRUE
  
  if(ncol(z)==10){
    names(z) <- c("m", "method", "order", "qD", "qD.95.LCL", "qD.95.UCL", "SC", "SC.95.LCL", "SC.95.UCL", "site")    
  }else{
    names(z) <- c("m", "method", "order", "qD", "SC", "site")    
  }
  
  if ("qD.95.LCL" %in% names(z) == FALSE & se) {
    warning("invalid se setting, the iNEXT object do not consist confidence interval")
    se <- FALSE
  }
  else if ("qD.95.LCL" %in% names(z) & se) {
    se <- TRUE
  }
  else {
    se <- FALSE
  }
  if (type == 1L) {
    z$x <- z[, 1]
    z$y <- z$qD
    if (se) {
      z$y.lwr <- z$qD.95.LCL
      z$y.upr <- z$qD.95.UCL
    }
  }
  else if (type == 2L) {
    if (length(unique(z$order)) > 1) {
      z <- subset(z, order == unique(z$order)[1])
    }
    z$x <- z[, 1]
    z$y <- z$SC
    if (se) {
      z$y.lwr <- z$SC.95.LCL
      z$y.upr <- z$SC.95.UCL
    }
  }
  else if (type == 3L) {
    z$x <- z$SC
    z$y <- z$qD
    if (se) {
      z$y.lwr <- z$qD.95.LCL
      z$y.upr <- z$qD.95.UCL
    }
  }
  if (color.var == "none") {
    if (levels(factor(z$order)) > 1 & "site" %in% names(z)) {
      warning("invalid color.var setting, the iNEXT object consists multiple sites and orders, change setting as both")
      color.var <- "both"
      z$col <- paste(z$site, z$order, sep = "-")
    }
    else if ("site" %in% names(z)) {
      warning("invalid color.var setting, the iNEXT object consists multiple orders, change setting as order")
      color.var <- "site"
      z$col <- z$site
    }
    else if (levels(factor(z$order)) > 1) {
      warning("invalid color.var setting, the iNEXT object consists multiple sites, change setting as site")
      color.var <- "order"
      z$col <- z$order
    }
    else {
      z$col <- rep(1, nrow(z))
    }
  }
  else if (color.var == "order") {
    z$col <- z$order
  }
  else if (color.var == "site") {
    if (!"site" %in% names(z)) {
      warning("invalid color.var setting, the iNEXT object do not consist multiple sites, change setting as order")
      z$col <- order
    }
    z$col <- z$site
  }
  else if (color.var == "both") {
    if (!"site" %in% names(z)) {
      warning("invalid color.var setting, the iNEXT object do not consist multiple sites, change setting as order")
      z$col <- order
    }
    z$col <- paste(z$site, z$order, sep = "-")
  }
  if ("site" %in% names(z)) {
    g <- ggplot(z, aes(x = x, y = y, colour = factor(col))) + 
      geom_point(size = 5, data = subset(z, method == "observed"))
  }
  else {
    g <- ggplot(z, aes(x = x, y = y, colour = factor(col), shape = site)) + 
      geom_point(size = 5, data = subset(z, method == "observed"))
  }
  g <- g + geom_line(aes(linetype = factor(method, c("interpolated", 
                                                     "extrapolated"), c("interpolation", "extrapolation"))), 
                     size = 1.5) + guides(linetype = guide_legend(title = "Method"), 
                                          colour = guide_legend(title = "Order"), fill = guide_legend(title = "Order"), 
                                          shape = guide_legend(title = "Site")) + theme(legend.position = "bottom", 
                                                                                        text = element_text(size = 18))
  if (type == 2L) 
    g <- g + labs(x = "Number of sampling units", y = "Sample coverage")
  else if (type == 3L) 
    g <- g + labs(x = "Sample coverage", y = ifelse(is.phylo, "Phylogenetic diversity", "Species diversity"))
  else g <- g + labs(x = "Number of sampling units", y = ifelse(is.phylo, "Phylogenetic diversity", "Species diversity"))
  if (se) 
    g <- g + geom_ribbon(aes(ymin = y.lwr, ymax = y.upr, 
                             fill = factor(col), colour = NULL), alpha = 0.2)
  if (facet.var == "order") {
    if (length(levels(factor(z$order))) == 1) {
      warning("invalid facet.var setting, the iNEXT object do not consist multiple orders.")
    }
    else {
      g <- g + facet_grid(. ~ order)
      if (color.var == "both") {
        g <- g + guides(colour = guide_legend(title = "Site-Order", 
                                              ncol = length(levels(factor(z$order))), byrow = TRUE), 
                        fill = guide_legend(title = "Site-Order"))
      }
    }
  }
  if (facet.var == "site") {
    if (!"site" %in% names(z)) {
      warning("invalid facet.var setting, the iNEXT object do not consist multiple sites.")
    }
    else {
      g <- g + facet_grid(. ~ site, )
      if (color.var == "both") {
        g <- g + guides(colour = guide_legend(title = "Site-Order", 
                                              nrow = length(levels(factor(z$order)))), fill = guide_legend(title = "Site-Order"))
      }
    }
  }
  if (facet.var == "both") {
    if (length(levels(factor(z$order))) == 1 | !"site" %in% 
          names(z)) {
      warning("invalid facet.var setting, the iNEXT object do not consist multiple sites or orders.")
    }
    else {
      g <- g + facet_wrap(site ~ order)
      if (color.var == "both") {
        g <- g + guides(colour = guide_legend(title = "Site-Order", 
                                              nrow = length(levels(factor(z$site))), byrow = TRUE), 
                        fill = guide_legend(title = "Site-Order"))
      }
    }
  }
  return(g)
}
