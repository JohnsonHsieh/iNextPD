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

# accumulation function for qPD(m)
PhD.m.mle <- function(p, labels=names(p), phy, q, m, datatype="abundance"){
  if (!inherits(phy, "phylog")) 
    stop("Non convenient data")
  tmp <- ExpandData_(p, labels, phy, datatype)
  a <- tmp$branch_abun
  L <- tmp$branch_length
  T <- sum(L*a)
  sub <- function(m){
    if(q==0) sum(L*(1-(1-a)^m)) 
    else if(q==1) {
      a1 <- a[a>0 & a<1]
      sub1 <- function(k) {
        #sum(L[a>0 & a<1]*choose(m,k)*a1^k*(1-a1)^(m-k))
        sum(L[a>0 & a<1]*exp(lchoose(m,k)+log(a1)*k+log(1-a1)*(m-k)))
      }
      k <- 1:m
      gk <- sapply(k, sub1)
      Hm <- sum(-k/m/T*log(k/m/T)*gk)
      exp(Hm)
    } else if(q==2) {
      lam <- sum(L*((m-1)/m*a^2+1/m*a))
      T/(lam/T)
    } else{
      sub2 <- function(x, j){
        a1 <- Stirling2(q, j) * exp(lgamma(m+1) - lgamma(m-j+1) - q*log(m))
        b1 <- x^j
        a1*b1
      }
      tmp <- rep(0, length(a))
      for(i in 1:length(a)){
        tmp[i] <- L[i] * sum(sapply(1:q, function(j) sub2(a[i], j)))
      }
      lam <- sum(tmp)
      T*(lam/T)^(1/(1-q))
    }
  }
  sapply(m, sub)
}



PhD.m.Rcpp <- function(n, x, U, L, q, m){
  #require("Rcpp")
  # no visible binding for global variable [variable name]
  RPD <- NULL
  tmp <- as.matrix(cbind(U,L))
  data <- x
  t_bar <- sum(tmp[,1] * tmp[,2] / n)
  tmp <- as.matrix(tmp)
  Rcpp::cppFunction('
              double RPD(NumericMatrix x , int n  , int m , int q) {
              int nrow = x.nrow();
              double tbar=0;
              NumericVector ghat(m);
              for (int i = 0; i < nrow; i++) {
              tbar += x(i, 0)*x(i, 1)/n;
              }
              for (int k = 0; k < m ; k++) {
              for (int i = 0; i < nrow; i++) {
              if ( x(i,0) >= k+1 && x(i,0) <= n-m+k+1 )
              {
              ghat[k] +=  x(i,1)*exp(Rf_lchoose(x(i,0), k+1)+Rf_lchoose(n-x(i,0), m-k-1)-Rf_lchoose(n, m)) ;
              //ghat[k] +=  x(i,1)*exp(-Rf_lbeta(x(i,0)-k, k+2)-log(x(i,0)+1)-Rf_lbeta(n-x(i,0)-m+k+2, m-k)-log(n-x(i,0)+1)+Rf_lbeta(n-m+1, m+1)+log(n+1)) ;
              }
              else
              {
              ghat[k] += 0 ;
              }
              }
              }
              double out=0;
              if(q == 0){
              for(int j = 0; j < m; j++) {
              out += ghat[j];
              }
              }
              if(q == 1){
              for (int j = 0; j < m; j++) {
              out += -( (j+1) / (m*tbar) ) * log ( (j+1) / (m*tbar) ) * ghat[j]  ;
              }           
              out = exp(out) ;
              }
              if(q == 2){
              for (int j = 0; j < m; j++) {
              out += pow( ( (j+1) / (m*tbar) ),2) * ghat[j] ;
              }
              out = 1 / out ;
              }
              return out ;
              }')
  
  RPD_m <- RPD(tmp, n, n-1, q)
  obs <- RPD(tmp, n, n, q)
  
  #asymptotic value
  Asy <- function(q){

    f1 <- sum(tmp[, 1] == 1) 
    f2 <- sum(tmp[, 1] == 2)
    g1 <- sum(tmp[tmp[, 1] == 1, 2])  
    g2 <- sum(tmp[tmp[, 1] == 2, 2])  
    if(q==0){
      g0_hat <- ifelse( g2>((g1*f2)/(2*f1)) , ((n-1)/n)*(g1^2/(2*g2)) , ((n-1)/n)*(g1*(f1-1)/(2*(f2+1))) )
      asy <- obs+g0_hat
    }
    if(q==1){
      # A <- 0    
      # if(f2>0){
      #   A <- 2 * f2 / ((n-1) * f1 + 2 * f2)
      # }else if(f2 == 0&&f1 > 0){
      #   A <- 2/((n-1)*(f1-1)+2)  
      # }else{
      #   A <- 0 
      # } 
      g0.hat <- ifelse((2*g2*sum(x==1))>g1*sum(x==2),
                       (n-1)/n*g1^2/2/g2,
                       (n-1)/n*g1*(sum(x==1)-1)/2/(sum(x==2)+1)) 
      A <- g1/(n*g0.hat+g1)
      
      i <- 1:(n-1)
      qq <- sum( ((1-A)^i)/i )
      h2 <- (g1/n)*((1-A)^(-n+1))*(-log(A)-qq)
      tmp2 <- tmp[which((1 <= tmp[, 1] ) & ( tmp[, 1] <= ( n-1 ) ) ),]
      h1 <- sum(apply(X = tmp2, MARGIN = 1, FUN = function(x){
        a <- sum( 1/(x[1]:(n-1)) )
        x[2]*x[1]*a/n
      }))
      h <- h1+h2
      asy <- t_bar*exp(h/t_bar)
    }
    if(q==2){
      tmp3 <- tmp[tmp[, 1] >= 2, ]
      asy <- (sum( tmp3[, 2]*tmp3[, 1]*(tmp3[, 1]-1) / ( ((t_bar)^2)*n*(n-1) ) ) )^(-1)
    }
    return(asy)
  }
  asy <- Asy(q)
  #beta
  beta <- (obs-RPD_m)/(asy-RPD_m)
  
  #Extrapolation 
  EPD <- function(m,q){
    m <- m-n
    if( q == 0 | q == 1 ){
      EPD <- obs+(asy-obs)*(1-(1-beta)^m)
    }
    if( q == 2 ){
      g <- sum( (tmp[,2]/(t_bar)^2)*((1/(n+m))*(tmp[,1]/n)+((n+m-1)/(n+m))*(tmp[,1]*(tmp[,1]-1)/(n*(n-1)))) )
      EPD <- 1/g
    }
    return(EPD)
  }
  PD <- sapply(m,FUN = function(m){
    if(m<n){
      if(identical(m, floor(m))) RPD(tmp,n,m,q)
      else{
        a <- RPD(tmp,n,floor(m),q)
        b <- RPD(tmp,n,floor(m)+1,q)
        dat <- data.frame(x=c(floor(m), floor(m)+1),y=c(a,b))
        fit <- glm(y~x, data=dat, family=quasipoisson(link = "log"))
        a[1] <- predict(fit, data.frame(x=m), type = "response")
        a
      }
    }else if(m==n){
      obs
    }else{
      EPD(m,q)
    }
  })
  return(PD)
} 

# iNEXTPD estimator for q = 0, 1, 2, ...
PhD.m <- function(n, x, U, L, q, m){
  if(q%in%c(0,1,2)){
    PhD.m.Rcpp(n, x, U, L, q, m)
  }else{
    L.obs <- L[U>0]
    U.obs <- U[U>0]
    tbar <- sum(L.obs*U.obs/n)
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
      tbar*(lam/tbar)^(1/(1-q))
    }
    sapply(m, Funq)
  }
}