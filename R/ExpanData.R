ExpandData_ <- function(x, labels, phy, datatype="abundance"){
  if (!inherits(phy, "phylog")) 
    stop("Non convenient data")
  
  datatype <- check_datatype(datatype)
  
  if(datatype=="incidence_freq" | datatype=="incidence") 
    stop('only support datatype="incidence_raw"')
  
  my_match <- match(labels, names(phy$leaves))
  if(sum(is.na(my_match)) > 0) stop("Argument labels and tree leaves not matach")
  
  if(datatype=="abundance"){
    n <- sum(x)  #sample size
    names(x) <- labels
    a <- x[names(phy$leaves)]
    
    for(i in 1:length(phy$parts)){
      a[1+length(a)] <- sum(a[phy$parts[[i]]])
      names(a)[length(a)] <- names(phy$parts)[i]
    }
    data.frame("branch_abun"=a, "branch_length"=c(phy$leaves, phy$nodes))
    
  }else if(datatype=="incidence_raw"){
    y <- iNEXT::as.incfreq(x)
    t <- y[1]
    y <- y[-1]
    names(y) <- labels
    y <- y[names(phy$leaves)]
    Ut <- sum(y)
    
    rownames(x) <- labels
    x <- x[names(phy$leaves),]
    
    for(i in 1:length(phy$parts)){
      x <- rbind(x, colSums(x[phy$parts[[i]],])>0)
      rownames(x)[nrow(x)] <- names(phy$parts)[i]
    }
    yy <- rowSums(x)
    data.frame("branch_abun"=yy, "branch_length"=c(phy$leaves, phy$nodes))
  }
}


###########################################
#' Expand branch abundance/incience and branch length
#' 
#' \code{ExpandData}: Expand branch abundance/incience and branch length
#' 
#' @param x a vector/matrix/list of species abundances or a matrix of raw incidence table.\cr 
#' @param labels a vector of species name for input data.\cr 
#' @param phy a phylogenetic tree with \code{"phylog"} class.\cr 
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}), species by sampling-units incidence matrix (\code{datatype = "incidence_raw"}).
#' @return a data.frame with sample size and sample coverage.
#' @examples 
#' data(bird)
#' bird.abu <- bird$abun
#' bird.inc <- bird$inci
#' bird.lab <- rownames(bird$abun)
#' bird.phy <- ade4::newick2phylog(bird$tre)
#' ExpandData(bird.abu, labels=bird.lab, phy=bird.phy, datatype="abundance")
#' ExpandData(bird.inc, labels=bird.lab, phy=bird.phy, datatype="incidence_raw")
#' @export
ExpandData <- function(x, labels, phy, datatype="abundance"){
  if (!inherits(phy, "phylog")) 
    stop("Non convenient data")
  
  datatype <- check_datatype(datatype)
  
  if(datatype=="incidence_freq" | datatype=="incidence") 
    stop('only support datatype="incidence_raw"')
  
  if(class(x)=="list"){
    lapply(x, function(x) ExpandData_(x, labels, phy, datatype))
  }else if(class(x) %in% c("matrix","data.frame") & datatype=="abundance"){
    apply(x, 2, function(x) ExpandData_(x, labels, phy, datatype))
  }else{
    ExpandData_(x, labels, phy, datatype)
  }
}