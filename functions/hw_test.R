library(HardyWeinberg)

#' HWChisq Wrapper function for autosomal loci
#'
#' @param homoref Integer. The number of homozygous reference genotypes.
#' @param het Integer. The number of heterozygous genotypes.
#' @param homoalt Integer. The number of homozygous alternative genotypes.
#'
#' @return Float. P-value from a Chi-square test for Hardy-Weinberg Equilibrium.
#' @export
#'
#' @examples
#' hw_test_autosomes(298, 489, 213)
hw_test_autosomes <- function(homoref=0, het=0, homoalt=0) {
  
  x <- c(AA=homoref, AB=het, BB=homoalt)
  
  HW.test <- HWChisq(x, verbose = FALSE)
  
  HW.test$pval
}

#' HWChisq Wrapper function for X-linked loci
#'
#' @param ref_males Integer. The number of reference alleles in males.
#' @param alt_males Integer. The number of alternative alleles in males.
#' @param homoref_females Integer. The number of homozygous reference genotypes in females.
#' @param het_females Integer. The number of heterozygous genotypes in females.
#' @param homoalt_females Integer. The number of homozygous alternative genotypes in females.
#'
#' @return Float. P-value from a Chi-square test for Hardy-Weinberg Equilibrium.
#' @export
#'
#' @examples
#' hw_test_xlinked(392, 212, 275, 296, 80)
hw_test_xlinked <- function(ref_males=0, alt_males=0, homoref_females=0, het_females=0, homoalt_females=0) {
  
  x <- c(A=ref_males, B=alt_males, AA=homoref_females,  AB=het_females, BB=homoalt_females)
  
  HW.test <- HWChisq(x, x.linked = TRUE, verbose = FALSE)
  
  HW.test$pval
}

hw_test_triallelic_autosome <- function(AA=0, AB=0, BB=0, BC=0, AC=0, CC=0) {
  
  x <- c(AA=AA, AB=AB, AC=AC, BB=BB, BC=BC, CC=CC)
  
  x <- toTriangular(x)
  
  HW.test <- HWPerm.mult(x, verbose = TRUE, nperm = 100)
  
  HW.test$pval
  
}

hw_test_triallelic_xlink <- function(AA=0, AB=0, BB=0, BC=0, AC=0, CC=0, A=0, B=0, C=0){
  
  x <- c(AA, AB, AC, BB, BC, CC)
  
  y <- c(A, B, C)
  
  HW.test <- HWTriExact(x, y)
  
  HW.test$pval
  
}