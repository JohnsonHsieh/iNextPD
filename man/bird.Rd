\name{bird}
\alias{bird}
\docType{data}
\title{Phylogeny and abundance of bird}
\description{
This data set describes the phylogeny of 41 birds as reported by Jetz et al. (2012). It also gives the two sites of species abundance to these 41 species in November 2012 at Barrington Tops National Park, Australia.
}
\usage{data(bird)}
\format{
\code{bird} is a list containing the 2 following objects: 
\describe{
   \item{tre}{is a character string giving the phylogenetic tree in Newick format.} 
   \item{tab}{is a data frame with 41 species and two sites: North and South sites.}    
}}
\source{
Jetz, W., Thomas, G.H., Joy, J.B., Hartmann, K. & Mooers A.O. (2012). The global diversity of birds in space and time. Nature, 491, 444–448.
}
\examples{
\dontrun{
data(bird)
bird.phy <- newick2phylog(bird$tre)
plot(bird.phy)

bird.abun <- bird$abun
bird.lab <- rownames(bird$abun)
table.phylog(bird.abun, bird.phy, csize=4, f.phylog=0.7)
}}
\keyword{datasets}
