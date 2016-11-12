#' Printing iNextPD object
#' 
#' \code{print.iNextPD}: Print method for objects inheriting from class "iNextPD"
#' @param x an \code{iNextPD} object computed by \code{\link{iNextPD}}.
#' @param ... additional arguments.
#' @export
print.iNextPD <- function(x, ...){
  site.n <- nrow(x$DataInfo)
  order.n <- ifelse(site.n > 1, 
                    paste(unique(x$iNextPDEst[[1]]$order), collapse = ", "),
                    paste(unique(x$iNextPDEst$order), collapse = ", "))
  cat("Compare ", site.n, " assemblages with Hill number order q = ", order.n,".\n", sep="")
  cat("$class: iNEXT\n\n")
  cat("$DataInfo: basic data information\n")
  print(x$DataInfo)
  cat("\n")
  cat("$iNextPDEst: phylogenetic diversity estimates with rarefied and extrapolated samples.\n")
  if(class(x$iNextPDEst)=="data.frame"){
    y <- x$iNextPDEst
    m <- quantile(y[,1], type = 1)
    res <- y[y[,1]%in%m,]
  }else{
    res <- lapply((x$iNextPDEst), function(y){
      m <- quantile(y[,1], type = 1)
      y[y[,1]%in%m,]
    })
  }
  print(res)
  cat("\n")
  cat("$AsyPDEst: asymptotic phylogenetic diversity estimates along with related statistics.\n")
  print(x$AsyPDEst)
  cat("\n")
  cat("NOTE1: Only show five estimates, call iNextPD.object$iNextPDEst to show complete output.\n")
  cat("NOTE2: call iNextPD.object$ExpandData to show the completed branch abundance/incience and branch length (Ui, Li), i = 1, 2, ..., B.\n")
  return(invisible())
}


###############################################
#' Plotting iNextPD object
#' 
#' \code{plot.iNextPD}: Plotting method for objects inheriting from class "iNextPD"
#' @param x an \code{iNextPD} object computed by \code{\link{iNextPD}}.
#' @param type three types of plots: sample-size-based rarefaction/extrapolation curve (\code{type = 1}); 
#' sample completeness curve (\code{type = 2}); coverage-based rarefaction/extrapolation curve (\code{type = 3}).                 
#' @param se a logical variable to display confidence interval around the estimated sampling curve.
#' @param show.legend a logical variable to display legend.
#' @param show.main a logical variable to display title.
#' @param col a vector for plotting color.
#' @param xlab a title for the x axis.
#' @param ylab a title for the y axis.
#' @param ... arguments to be passed to methods, such as graphical parameters (\code{\link{par}}).
#' @examples
#' 
#' # single-assemblage abundance data
#' data(bird)
#' bird.phy <- ade4::newick2phylog(bird$tre)
#' bird.lab <- rownames(bird$abun)
#' out1 <- iNextPD(bird$abun$North.site, bird.lab, bird.phy, 
#'         q=1, datatype="abundance", endpoint=500)
#' plot(x=out1, type=1)
#' plot(x=out1, type=2)
#' plot(x=out1, type=3)
#' 

#' @export
plot.iNextPD <- function(x, type=1, se=TRUE, show.legend=TRUE, show.main=TRUE, col=NULL, xlab=NULL, ylab=NULL,...){
  
  if(class(x) != "iNextPD")
    stop("invalid object class")
  TYPE <-  c(1, 2, 3)
  SPLIT <- c("none", "order", "site", "both")
  if(is.na(pmatch(type, TYPE)) | pmatch(type, TYPE) == -1)
    stop("invalid plot type")
  
  type <- pmatch(type, 1:3)
  
  y <- method <- site <- shape <- y.lwr <- y.upr <- NULL
  site <<- NULL
  
  z <- x$iNextPDEst
  if(class(z) == "list"){
    z <- data.frame(do.call("rbind", z), site=rep(names(z), sapply(z, nrow)))
    rownames(z) <- NULL
  }else{
    z$site <- ""
    z$site <- factor(z$site)
  }
  
  if(ncol(z) == 6  & se) {
    warning("invalid se setting, the iNextPD object do not consist confidence interval")
    se <- FALSE
  }else if(ncol(z) > 6 & se) {
    se <- TRUE
  }else{
    se <- FALSE
  }
  
  if(type==1L) {
    z$x <- z[,1]
    z$y <- z$qPD
    if(is.null(xlab)) xlab <- ifelse(names(x$DataInfo)[2]=="n", "Number of individuals", "Number of sampling units")
    if(is.null(ylab)) ylab <- "Phylogenetic diversity"
    if(se){
      z$y.lwr <- z[,5]
      z$y.upr <- z[,6]
    }
  }else if(type==2L){
    if(length(unique(z$order))>1){
      z <- subset(z, order==unique(z$order)[1])
    }
    z$x <- z[,1]
    z$y <- z$SC
    if(is.null(xlab)) xlab <- ifelse(names(x$DataInfo)[2]=="n", "Number of individuals", "Number of sampling units")
    if(is.null(ylab)) ylab <- "Sample coverage"
    if(se){
      z$y.lwr <- z[,8]
      z$y.upr <- z[,9]
    }
  }else if(type==3L){
    z$x <- z$SC
    z$y <- z$qPD
    if(is.null(xlab)) xlab <- "Sample coverage"
    if(is.null(ylab)) ylab <- "Phylogenetic diversity"
    if(se){
      z$y.lwr <- z[,5]
      z$y.upr <- z[,6]
    }
  }
  
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length=n+1)
    hcl(h=hues, l=65, c=100)[1:n]
  }
  
  conf.reg=function(x,LCL,UCL,...) {
    x.sort <- order(x)
    x <- x[x.sort]
    LCL <- LCL[x.sort]
    UCL <- UCL[x.sort]
    polygon(c(x,rev(x)),c(LCL,rev(UCL)), ...)
  }
  
  SITE <- levels(z$site)
  ORDER <- unique(z$order)
  
  if(is.null(col)){
    col <- gg_color_hue(length(SITE))
  }else{
    col <- rep(col,length(SITE))[1:length(SITE)]
  }
  pch <- (16+1:length(SITE))%%25
  
  for(j in 1:length(ORDER)){
    if(se==TRUE){
      tmp.sub <- subset(z, order==ORDER[j])
      tmp.j <- data.frame(site=tmp.sub$site, order=tmp.sub$order,
                          method=tmp.sub$method, 
                          x=tmp.sub$x, y=tmp.sub$y,
                          y.lwr=tmp.sub$y.lwr, y.upr=tmp.sub$y.upr)
      
      plot(y.upr~x, data=tmp.j, type="n", xlab="", ylab="", ...)
    }else{
      tmp.sub <- subset(z, order==ORDER[j])
      
      tmp.j <- data.frame(site=tmp.sub$site, order=tmp.sub$order,
                          method=tmp.sub$method, 
                          x=tmp.sub$x, y=tmp.sub$y)
      
      plot(y~x, data=tmp.j, type="n", xlab="", ylab="", ...)
    }
    
    for(i in 1:length(SITE)){
      tmp <- subset(tmp.j, site==SITE[i])
      if(se==TRUE){
        conf.reg(x=tmp$x, LCL=tmp$y.lwr, UCL=tmp$y.upr, border=NA, col=adjustcolor(col[i], 0.25))
      }
      lines(y~x, data=subset(tmp, method=="interpolated"), lty=1, lwd=2, col=col[i])
      lines(y~x, data=subset(tmp, method=="extrapolated"), lty=2, lwd=2, col=col[i])
      points(y~x, data=subset(tmp, method=="observed"), pch=pch[i], cex=2, col=col[i])
      
    }
    if(show.legend==TRUE){
      if(type==3L){
        legend("topleft", legend=paste(SITE), col=col, lty=1, lwd=2, pch=pch, cex=1, bty="n")
      }else{
        legend("bottomright", legend=paste(SITE), col=col, lty=1, lwd=2, pch=pch, cex=1, bty="n")
      }
    }
    title(xlab=xlab, ylab=ylab)
    if(show.main==TRUE) title(main=paste("Order q =", ORDER[j]))
    par(ask=TRUE)
  }
  par(ask=FALSE)
}




#
#
###############################################
#' ggplot2 extension for an iNextPD object
#' 
#' \code{ggiNEXT}: the \code{ggplot} extension for \code{\link{iNextPD}} Object to plot sample-size- and coverage-based rarefaction/extrapolation curves along with a bridging sample completeness curve
#' @param x an \code{iNextPD} object computed by \code{\link{iNextPD}}.
#' @param type three types of plots: sample-size-based rarefaction/extrapolation curve (\code{type = 1}); 
#' sample completeness curve (\code{type = 2}); coverage-based rarefaction/extrapolation curve (\code{type = 3}).                 
#' @param se a logical variable to display confidence interval around the estimated sampling curve.
#' @param facet.var create a separate plot for each value of a specified variable: 
#'  no separation \cr (\code{facet.var="none"}); 
#'  a separate plot for each diversity order (\code{facet.var="order"}); 
#'  a separate plot for each site (\code{facet.var="site"}); 
#'  a separate plot for each combination of order x site (\code{facet.var="both"}).              
#' @param color.var create curves in different colors for values of a specified variable:
#'  all curves are in the same color (\code{color.var="none"}); 
#'  use different colors for diversity orders (\code{color.var="order"}); 
#'  use different colors for sites (\code{color.var="site"}); 
#'  use different colors for combinations of order x site (\code{color.var="both"}).  
#' @param grey a logical variable to display grey and white ggplot2 theme. 
#' @param ... other arguments passed on to methods. Not currently used.
#' @return a ggplot2 object
#' @examples
#' # single-assemblage abundance data
#' data(bird)
#' bird.phy <- ade4::newick2phylog(bird$tre)
#' bird.lab <- rownames(bird$abun)
#' out1 <- iNextPD(bird$abun$North.site, bird.lab, bird.phy, 
#'         q=1, datatype="abundance", endpoint=400, se=TRUE)
#' ggiNEXT(x=out1, type=1)
#' ggiNEXT(x=out1, type=2)
#' ggiNEXT(x=out1, type=3)
#' 
#'\dontrun{
#' # single-assemblage incidence data with three orders q
#' out2 <- iNextPD(bird$inci$North.site, bird.lab, bird.phy, 
#'         q=c(0,1,2), datatype="incidence_raw", endpoint=25)
#' ggiNEXT(out2, se=FALSE, color.var="order")
#' 
#' # multiple-assemblage abundance data with three orders q
#' out3 <-  iNextPD(bird$abun, bird.lab, bird.phy, 
#'         q=c(0,1,2), datatype="abundance", endpoint=400)
#' ggiNEXT(out3, facet.var="site", color.var="order")
#' ggiNEXT(out3, facet.var="both", color.var="both")
#'}
#' @export
#' 
ggiNEXT <- function(x, type=1, se=TRUE, facet.var="none", color.var="site", grey=FALSE){  
  UseMethod("ggiNEXT", x)
}

#' @export
#' @rdname ggiNEXT
ggiNEXT.iNextPD <- function(x, type=1, se=TRUE, facet.var="none", color.var="site", grey=FALSE){
  TYPE <-  c(1, 2, 3)
  SPLIT <- c("none", "order", "site", "both")
  if(is.na(pmatch(type, TYPE)) | pmatch(type, TYPE) == -1)
    stop("invalid plot type")
  if(is.na(pmatch(facet.var, SPLIT)) | pmatch(facet.var, SPLIT) == -1)
    stop("invalid facet variable")
  if(is.na(pmatch(color.var, SPLIT)) | pmatch(color.var, SPLIT) == -1)
    stop("invalid color variable")
  
  type <- pmatch(type, 1:3)
  facet.var <- match.arg(facet.var, SPLIT)
  color.var <- match.arg(color.var, SPLIT)
  if(facet.var=="order") color.var <- "site"
  if(facet.var=="site") color.var <- "order"
  
  options(warn = -1)
  z <- fortify(x, type=type)
  options(warn = 0)
  if(ncol(z) ==7) {se <- FALSE}
  datatype <- unique(z$datatype)
  if(color.var=="none"){
    if(levels(factor(z$order))>1 & "site"%in%names(z)){
      warning("invalid color.var setting, the iNextPD object consists multiple sites and orders, change setting as both")
      color.var <- "both"
      z$col <- z$shape <- paste(z$site, z$order, sep="-")
      
    }else if("site"%in%names(z)){
      warning("invalid color.var setting, the iNextPD object consists multiple orders, change setting as order")
      color.var <- "site"
      z$col <- z$shape <- z$site
    }else if(levels(factor(z$order))>1){
      warning("invalid color.var setting, the iNextPD object consists multiple sites, change setting as site")
      color.var <- "order"
      z$col <- z$shape <- factor(z$order)
    }else{
      z$col <- z$shape <- rep(1, nrow(z))
    }
  }else if(color.var=="order"){     
    z$col <- z$shape <- factor(z$order)
  }else if(color.var=="site"){
    if(!"site"%in%names(z)){
      warning("invalid color.var setting, the iNextPD object do not consist multiple sites, change setting as order")
      z$col <- z$shape <- factor(z$order)
    }
    z$col <- z$shape <- z$site
  }else if(color.var=="both"){
    if(!"site"%in%names(z)){
      warning("invalid color.var setting, the iNextPD object do not consist multiple sites, change setting as order")
      z$col <- z$shape <- factor(z$order)
    }
    z$col <- z$shape <- paste(z$site, z$order, sep="-")
  }
  
  z$lty <- factor(z$method, c("interpolated", "observed", "extrapolated"), c("interpolation", "interpolation", "extrapolation"))
  z$col <- factor(z$col)
  data.sub <- z[which(z$method=="observed"),]
  
  g <- ggplot(z, aes_string(x="x", y="y", colour="col")) + 
    geom_point(aes_string(shape="shape"), size=5, data=data.sub)
  
  
  g <- g + geom_line(aes_string(linetype="lty"), lwd=1.5) +
    guides(linetype=guide_legend(title="Method"),
           colour=guide_legend(title="Guides"), 
           fill=guide_legend(title="Guides"), 
           shape=guide_legend(title="Guides")) +
    theme(legend.position = "bottom", 
          legend.title=element_blank(),
          text=element_text(size=18)) 
  
  if(type==2L) {
    g <- g + labs(x="Number of sampling units", y="Sample coverage")
    if(datatype=="abundance") g <- g + labs(x="Number of individuals", y="Sample coverage")
  }
  else if(type==3L) {
    g <- g + labs(x="Sample coverage", y="Phylogenetic diversity")
  }
  else {
    g <- g + labs(x="Number of sampling units", y="Phylogenetic diversity")
    if(datatype=="abundance") g <- g + labs(x="Number of individuals", y="Phylogenetic diversity")
  }
  
  if(se)
    g <- g + geom_ribbon(aes_string(ymin="y.lwr", ymax="y.upr", fill="factor(col)", colour="NULL"), alpha=0.2)
  
  
  if(facet.var=="order"){
    if(length(levels(factor(z$order))) == 1 & type!=2){
      warning("invalid facet.var setting, the iNextPD object do not consist multiple orders.")      
    }else{
      g <- g + facet_wrap(~order, nrow=1)
      if(color.var=="both"){
        g <- g + guides(colour=guide_legend(title="Guides", ncol=length(levels(factor(z$order))), byrow=TRUE),
                        fill=guide_legend(title="Guides"))
      }
    }
  }
  
  if(facet.var=="site"){
    if(!"site"%in%names(z)) {
      warning("invalid facet.var setting, the iNextPD object do not consist multiple sites.")
    }else{
      g <- g + facet_wrap(~site, nrow=1)
      if(color.var=="both"){
        g <- g + guides(colour=guide_legend(title="Guides", nrow=length(levels(factor(z$order)))),
                        fill=guide_legend(title="Guides"))
      }
    }
  }
  
  if(facet.var=="both"){
    if(length(levels(factor(z$order))) == 1 | !"site"%in%names(z)){
      warning("invalid facet.var setting, the iNextPD object do not consist multiple sites or orders.")
    }else{
      g <- g + facet_wrap(site~order) 
      if(color.var=="both"){
        g <- g +  guides(colour=guide_legend(title="Guides", nrow=length(levels(factor(z$site))), byrow=TRUE),
                         fill=guide_legend(title="Guides"))
      }
    }
  }
  
  if(grey){
    g <- g + theme_bw(base_size = 18) +
      scale_fill_grey(start = 0, end = .4) +
      scale_colour_grey(start = .2, end = .2) +
      guides(linetype=guide_legend(title="Method"), 
             colour=guide_legend(title="Guides"), 
             fill=guide_legend(title="Guides"), 
             shape=guide_legend(title="Guides")) +
      theme(legend.position="bottom",
            legend.title=element_blank())
  }
  return(g)
  
}

#' @export
#' @rdname ggiNEXT
ggiNEXT.default <- function(x, ...){
  stop(
    "iNextPD doesn't know how to deal with data of class ",
    paste(class(x), collapse = "/"),
    call. = FALSE
  )
}


#' Fortify method for classes from the iNextPD package.
#'
#' @name fortify.iNextPD
#' @param model \code{iNextPD} to convert into a dataframe.
#' @param data not used by this method
#' @param type three types of plots: sample-size-based rarefaction/extrapolation curve (\code{type = 1}); 
#' sample completeness curve (\code{type = 2}); coverage-based rarefaction/extrapolation curve (\code{type = 3}).                 
#' @param se a logical variable to calculate the bootstrap standard error and \code{conf} confidence interval.
#' @param ... not used by this method
#' @export
#' @examples
#' # single-assemblage abundance data
#' data(bird)
#' bird.lab <- rownames(bird$abun)
#' bird.phy <- ade4::newick2phylog(bird$tre)
#' out1 <- iNextPD(bird$abun$North.site, bird.lab, bird.phy, 
#'         q=0, datatype="abundance", endpoint=400)
#' ggplot2::fortify(out1, type=1)

fortify.iNextPD <- function(model, data = model$iNextPDEst, type = 1, se=TRUE,...) {
  datatype <- ifelse(names(model$DataInfo)[2]=="n","abundance","incidence")
  z <- data
  if(class(z) == "list"){
    z <- data.frame(do.call("rbind", z), site=rep(names(z), sapply(z, nrow)))
    rownames(z) <- NULL
  }else{
    z$site <- "Site1"
  }
  
  if(ncol(z)==6 & se==TRUE) {
    warning("invalid se setting, the iNEXT object do not consist confidence interval")
    se <- FALSE
  }else if(ncol(z)>6) {
    se <- TRUE
  }
  
  if(type==1L) {
    z$x <- z[,1]
    z$y <- z$qPD
    if(se){
      z$y.lwr <- z[,5]
      z$y.upr <- z[,6]
    }
  }else if(type==2L){
    if(length(unique(z$order))>1){
      z <- subset(z, order==unique(z$order)[1])
    }
    z$x <- z[,1]
    z$y <- z$SC
    if(se){
      z$y.lwr <- z[,8]
      z$y.upr <- z[,9]
    }
  }else if(type==3L){
    z$x <- z$SC
    z$y <- z$qPD
    if(se){
      z$y.lwr <- z[,5]
      z$y.upr <- z[,6]
    }
  }
  z$datatype <- datatype
  z$plottype <- type
  if(se){
    data <- z[,c("datatype","plottype","site","method","order","x","y","y.lwr","y.upr")]
  }else{
    data <- z[,c("datatype","plottype","site","method","order","x","y")]
  }
  data
}

