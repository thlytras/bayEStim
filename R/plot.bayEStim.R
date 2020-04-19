#' @import graphics
#' 
#' @export
#' @noRd
plot.bayEStim <- function(x, ...) {
  summ <- summary(x)
  plot(summ, ...)
}


#' @import graphics
#' 
#' @export
#' @noRd
plot.summary.bayEStim <- function(x, 
    col.inc = c("limegreen", "skyblue"), col.R = "red", col.ci = "grey", 
    lwd.bar = 5, lwd.line = 3, 
    xlab = "Date", ylab.inc = "Number of cases", ylab.R = "R_effective", ...) {

  par(mfrow=c(2,1))
  with(x$I, plot(date, est, type="n", lwd=5, lend=1, 
    ylim=c(0, max(x$I$hi)), xlim=range(date),
    xlab="Date", ylab="Number of cases"
  ))
  with(x$I, polygon(x=c(date,rev(date)), y=c(lo,rev(hi)), border=NA, col=col.ci))
  with(x$I, points(date, est, type="h", lwd=lwd.bar, lend=1, col=col.inc[2]))
  with(x$I_local, points(date, est, type="h", lwd=lwd.bar, lend=1, col=col.inc[1]))
  
  with(x$I, plot(date, est, type="n", 
    ylim=c(0, max(x$R$hi)), xlim=range(date),
    xlab="Date", ylab="R_effective"
  ))
  with(x$R, polygon(x=c(date,rev(date)), y=c(lo,rev(hi)), border=NA, col=col.ci))
  abline(h=1)
  with(x$R, points(date, est, type="l", lwd=lwd.line, col="red"))
}

