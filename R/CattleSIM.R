#' Cattle dataset with simulated management units
#'
#' The cattle dataset is simulated with QMSim software (Sargolzaei and Schenkel, 2009).   
#' This dataset includes 2500 individuals across six generations (from founder to generation 5), 
#' each with 10,000 single nucleotide polymorphisms spread over 29 autosomes. Single phenotype with heritability of 0.6 was simulated.
#' In addition, 6 fixed effects including one Cluster effect and 5 management unit simulation scenarios are simulated. 
#' Fixed effect Cluster is simulated using K-medoid algorithm (Kaufman and Rousseeuw, 1990) to assign 2,500 individuals into 8 clusters. 
#' MUSC1 is simulated by randomly allocate eighter clusters into two sets, which is regarded as the least connected design. From MUSC2 to 5,  
#' randomly sampled individuals (i.e., 140, 210, 280 and 350)  are exchanged between the two sets in MUSC1, which aims to steadily increase 
#' the degree of relatedness. 
#'
#' @docType data
#' @name CattleSIM
#' @keywords datasets
#' @usage 
#' data(CattleSIM)
#' 
#' @author Haipeng Yu and Gota Morota 
#' 
#' Maintainer: Haipeng Yu \email{haipengyu@@vt.edu}
#'
#' @example man/examples/CattleSIM.R 
#' @references Sargolzaei, M., and F. S. Schenkel. 2009. Qmsim: a large-scale
#'   genome simulator for livestock. Bioinformatics 25:680–681. doi:10.1093/bioinformatics/btp045
#' @references Kaufman, L. and P. Rousseeuw. 1990. Finding groups in data: 
#' an introduction to cluster analysis. John Wiley and Sons, New York.
NULL

