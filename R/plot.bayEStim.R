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
    xlab = "Date", ylab.inc = "Number of cases", ylab.R = "R_effective", 
    no.mfrow = FALSE, ...) {

  if (!no.mfrow) par(mfrow=c(2,1))
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
  drawBand(x$R, "red", alpha=0.2)
#   with(x$R, polygon(x=c(date,rev(date)), y=c(lo,rev(hi)), border=NA, col=col.ci))
  abline(h=1)
#   with(x$R, points(date, est, type="l", lwd=lwd.line, col="red"))
}


addalpha <- function(colors, alpha=1.0) {
  r <- col2rgb(colors, alpha=T)
  # Apply alpha
  r[4,] <- alpha*255
  r <- r/255.0
  return(rgb(r[1,], r[2,], r[3,], r[4,]))
}

drawBand <- function(x, col, alpha=0.2) {
  with(x, {
    polygon(c(date,rev(date)), c(lo,rev(hi)), border=NA, col=addalpha(col, alpha))
    points(date, est, type="l", lwd=2, col=col)
  })
}

