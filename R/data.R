#' Bird phylogeny and survey data 
#' 
#' This data set describes the phylogeny of 41 birds as reported by Jetz et al. (2012). 
#' It also gives the two sites of species abundance and incidence data to these 41 species in November 2012 at Barrington Tops National Park, Australia.
#'
#' @format A data frame with 53940 rows and 10 variables:
#' \describe{
#'   \item{tre}{is a character string giving the phylogenetic tree in Newick format.} 
#'   \item{abun}{is a data frame with 41 species and two sites: North and South sites.}
#'   \item{inic}{is a list of two site (data.frame) for species by sampling-units incidence matrix.}
#' }
#' @source Jetz, W., Thomas, G.H., Joy, J.B., Hartmann, K. & Mooers A.O. (2012). The global diversity of birds in space and time. Nature, 491, 444-448.
#' @examples
#' data(bird)
#' bird.phy <- ade4::newick2phylog(bird$tre)
#' plot(bird.phy)
#' bird.abun <- bird$abun
#' bird.lab <- rownames(bird$abun)
#' ade4::table.phylog(bird.abun, bird.phy, csize=4, f.phylog=0.7)
"bird"