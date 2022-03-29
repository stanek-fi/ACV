#' Printing method for `estimateL`
#'
#' Prettify the output of `estimateL` when printed.
#' @export
#' @keywords internal

print.estimateL <- function(x, ...){
  temp=data.frame(
    "Estimate" = x$estimate,
    "Std.Error"=sqrt(x$var),
    "CI.2.5"=stats::qnorm(0.025, mean=x$estimate, sd=sqrt(x$var)),
    "CI.97.5"=stats::qnorm(0.975, mean=x$estimate, sd=sqrt(x$var))
  )
  row.names(temp)="Loss"
  print(temp)
}


