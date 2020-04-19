#' @importFrom coda varnames
#' @importFrom stats quantile
#' 
#' @export
#' @noRd
summary.bayEStim <- function(object, ...) {
  res <- list(
    I_local = as.data.frame(summary(object$mcmcSamples[, grep("I_local", 
        varnames(object$mcmcSamples))])[[2]][,c(3,1,5)]),
    I_imported = as.data.frame(summary(object$mcmcSamples[, grep("I_imported", 
        varnames(object$mcmcSamples))])[[2]][,c(3,1,5)]),
    I = as.data.frame(t(apply(
        sapply(grep("I_local", varnames(object$mcmcSamples)), 
            function(i) unlist(object$mcmcSamples[,i])) + 
        sapply(grep("I_imported", varnames(object$mcmcSamples)), 
            function(i) unlist(object$mcmcSamples[,i])), 
        2, function(y) quantile(y, c(0.5, 0.025, 0.975))))),
    R = as.data.frame(summary(object$mcmcSamples[, grep("^R", 
        varnames(object$mcmcSamples))])[[2]][,c(3,1,5)])
  )
  rownames(res$I) <- sprintf("I[%s]", 1:nrow(res$I))
  for (i in 1:length(res)) names(res[[i]]) <- c("est", "lo", "hi")
  for (i in 1:3) {
    res[[i]]$date <- seq.Date(object$config$date0+1, object$config$maxDate, by="days")
    res[[i]] <- res[[i]][,c("date","est","lo","hi")]
  }
  if (is.null(object$config$t_breaks)) {
    res$R$date <- seq.Date(object$config$t_from + object$config$t_window - 1, object$config$maxDate, by="days")
    res$R <- res$R[,c("date","est","lo","hi")]
  } else {
    res$R_periods <- res$R
    res$R_periods$from <- c(object$config$t_from, object$config$t_breaks)
    res$R_periods$to <- c(object$config$t_breaks-1, object$config$maxDate)
    res$R_periods <- res$R_periods[,c("from","to","est","lo","hi")]
    res$R <- data.frame(date = seq.Date(object$config$t_from, object$config$maxDate, by="days"))
    res$R <- cbind(res$R, res$R_periods[sapply(res$R$date, function(d) which(d<=res$R_periods$to)[1]), c("est","lo","hi")])
  }
  
  res$data <- list(donset=
    as.data.frame.matrix(with(object$config$data, table(donset, local))))
  names(res$data$donset) <- c("I_imported", "I_local")
  res$data$donset$I <- rowSums(res$data$donset)
  res$data$donset$date <- as.Date(rownames(res$data$donset))
  rownames(res$data$donset) <- NULL
  res$data$donset <- res$data$donset[,c(4:1)]
  if (!is.null(object$config$data$ddiag)) {
    res$data$ddiag <- as.data.frame.matrix(with(object$config$data, table(ddiag, local)))
    names(res$data$ddiag) <- c("I_imported", "I_local")
    res$data$ddiag$I <- rowSums(res$data$ddiag)
    res$data$ddiag$date <- as.Date(rownames(res$data$ddiag))
    rownames(res$data$ddiag) <- NULL
    res$data$ddiag <- res$data$ddiag[,c(4:1)]
  }
  
  class(res) <- "summary.bayEStim"
  return(res)
}

