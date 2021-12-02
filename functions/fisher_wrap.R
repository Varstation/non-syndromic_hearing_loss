library(glue)

fisher_wrap <- function(A, B, C, D) {
  
  A=A
  B=B
  C=C
  D=D
  
  dat <- matrix(c(A, B, C, D), 2)
  
  result <- fisher.test(dat)
  
  result$p.value
  
  glue("{result$p.value},{result$estimate},{result$conf.int[1]},{result$conf.int[2]}")
  
}