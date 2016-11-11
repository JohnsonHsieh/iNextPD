# Estimation of species reletive abundance distribution
SPBoot_ <- function(x, datatype){
  datatype <- check_datatype(datatype)
  
  jade1 <- function(Spec, n){
    Sobs <- sum(Spec > 0)   #observed species
    f1 <- sum(Spec == 1)   #singleton 
    f2 <- sum(Spec == 2)   #doubleton
    f0.hat <- ifelse(f2 == 0, (n - 1) / n * f1 * (f1 - 1) / 2, (n - 1) / n * f1 ^ 2/ 2 / f2)  #estimation of unseen species via Chao1
    A <- ifelse(f1>0, n*f0.hat/(n*f0.hat+f1), 1)
    a <- f1/n*A
    b <- sum(Spec / n * (1 - Spec / n) ^ n)
    w <- a / b      	#adjusted factor for rare species in the sample
    Prob.hat <- Spec / n * (1 - w * (1 - Spec / n) ^ n)					#estimation of relative abundance of observed species in the sample
    Prob.hat.Unse <- rep(a/ceiling(f0.hat), ceiling(f0.hat))  	#estimation of relative abundance of unseen species in the sample
    return(sort(c(Prob.hat, Prob.hat.Unse), decreasing=TRUE))
  }
  
  if(datatype=="abundance"){
    jade1(x, sum(x))
  }else if(datatype=="incidence"){
    jade1(x[-1], x[1])
  }else if(datatype=="incidence_raw"){
    xx <- as.incfreq(x)
    jade1(xx[-1], xx[1])
  }
}

#' Estimation of species relative abundance or detection probability distribution
#' 
#' \code{SPBoot}: Expand bootstraping species relative abundance or detection probability
#' 
#' @param x a vector/matrix/list of species abundances or a matrix of raw incidence table.\cr 
#' @param datatype of input data: individual-based abundance data (datatype = "abundance"), sampling-unit-based incidence frequencies data (datatype = "incidence_freq") or species by sampling-units incidence matrix (datatype = "incidence_raw")..
#' @return a list of vector with species relative abundance or detection probability distribution.
#' @export
#' @examples 
#' data(bird)
#' bird.inc <- bird$inci
#' SPBoot(bird$abun, datatype="abundance")
#' SPBoot(bird$inci, datatype="incidence_raw")
SPBoot <- function(x, datatype="abundance"){

  datatype <- check_datatype(datatype)
  
  if(class(x)=="list"){
    lapply(x, function(x) SPBoot_(x, datatype))
  }else if(class(x) %in% c("matrix","data.frame") & datatype=="abundance"){
    apply(x, 2, function(x) SPBoot_(x, datatype))
  }else{
    SPBoot_(x, datatype)
  }
}



# Show {(ai_hat, Li_hat); i=1, 2,.., Bost, ..., Best}
PDBoot_ <- function(abun, labels, phy, datatype="abundance"){
  if (!inherits(phy, "phylog")) 
    stop("Non convenient data")
  
  datatype <- check_datatype(datatype)
  
  if(datatype=="incidence_freq" | datatype=="incidence") 
    stop('only support datatype="incidence_raw"')
  
  if(datatype=="abundance"){
    n <- sum(abun)        #sample size
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
  }else if(datatype=="incidence_raw"){
    y <- iNEXT::as.incfreq(abun)
    t <- y[1]
    y <- y[-1]
    names(y) <- labels
    y <- y[names(phy$leaves)]
    Ut <- sum(y)
    
    #abun <- data.frame(abun)
    rownames(abun) <- labels
    abun <- abun[names(phy$leaves),]
    
    Q1 <- sum(y==1)
    Q2 <- sum(y==2)
    Q0.hat <- ifelse(Q2==0, (t-1)/t*Q1*(Q1-1)/2, (t-1)/t*Q1^2/2/Q2)
    A <- ifelse(Q1>0, t*Q0.hat/(t*Q0.hat+Q1), 1)
    Chat <- 1 - Q1/Ut*A
    Q0 <- max(round(Q0.hat), 1)
    tau <- Ut / t * (1 - Chat) / sum(y / t * (1 - y / t)^t)
    aa <- y / t * (1 - tau * (1 - y / t)^t)
    
    for(i in 1:length(phy$parts)){
      abun <- rbind(abun, colSums(abun[phy$parts[[i]],])>0)
      rownames(abun)[nrow(abun)] <- names(phy$parts)[i]
      
      aa[1+length(aa)] <- 1-prod(1-aa[phy$parts[[i]]])
      names(aa)[length(aa)] <- names(phy$parts)[i]
    }
    yy <- rowSums(abun)
    tmp <- data.frame("branch_abun"=yy, "branch_length"=c(phy$leaves, phy$nodes))
    U <- tmp$branch_abun
    L <- tmp$branch_length
    U.obs <- U[U>0]
    L.obs <- L[U>0]
    a.obs <- aa[U>0]
    
    a.und <- rep(Ut / t * (1 - Chat) / ceiling(Q0.hat), times=ceiling(Q0.hat))
    
    g1 <- sum(L.obs[U.obs==1])
    g2 <- sum(L.obs[U.obs==2])
    g0 <- (t-1)/t * ifelse((2*g2*Q1) > (g1*Q2), g1^2/(2*g2), g1*(Q1-1)/(2*(Q2+1)))
    L0 <-  rep(g0/ceiling(Q0.hat), ceiling(Q0.hat))
    
    ai <- c(a.obs, a.und)
    Li <- c(L.obs, L0)
    data.frame("branch_abun"=ai, "branch_length"=Li)
  }
}

###########################################
#' Expand bootstraping branch abundance/incience and branch length
#' 
#' \code{PDBoot}: Expand bootstraping branch abundance/incience and branch length
#' 
#' @param x a vector/matrix/list of species abundances or a matrix of raw incidence table.\cr 
#' @param labels a vector of species name for input data.\cr 
#' @param phy a phylogenetic tree with \code{"phylog"} class.\cr 
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}), species by sampling-units incidence matrix (\code{datatype = "incidence_raw"}).
#' @return a list of data.frame with bootstraping branch abundance/incience and branch length.
#' @export
#' @examples 
#' data(bird)
#' bird.lab <- rownames(bird$abun)
#' bird.phy <- ade4::newick2phylog(bird$tre)
#' bird.inc <- bird$inci
#' PDBoot(bird.inc, labels=bird.lab, phy=bird.phy, datatype="incidence_raw")
#' 
PDBoot <- function(x, labels, phy, datatype="abundance"){
  if (!inherits(phy, "phylog")) 
    stop("Non convenient data")
  datatype <- check_datatype(datatype)
  
  # no visible binding for global variable [variable name]
  abun <- NULL
  
  if(datatype=="incidence_freq" | datatype=="incidence") 
    stop('only support datatype="incidence_raw"')
  
  if(class(x)=="list"){
    lapply(x, function(x) PDBoot_(x, labels, phy, datatype))
  }else if(class(abun) %in% c("matrix","data.frame") & datatype=="abundance"){
    apply(x, 2, function(x) PDBoot_(x, labels, phy, datatype))
  }else{
    PDBoot_(x, labels, phy, datatype)
  }
}
