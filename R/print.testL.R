#' Printing method for `testL`
#'
#' Prettify the output of `testL` when printed.
#' @export
#' @keywords internal

print.testL <- function(x, ...){
  cat(x$test, "test\n")
  cat("\nEstimated difference L(algorithm1) - L(algorithm2):\n", x$estimate)
  cat("\nAlternative hypothesis:\n L(algorithm1) - L(algorithm2) ", x$Ha)
  cat("\nt-stat:\n", x$tval)
  cat("\np-val:\n", x$pval)
}
