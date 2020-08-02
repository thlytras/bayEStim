#' Estimate effective reproduction number in a Bayesian framework
#'
#' \code{estimate_R_Bayes} estimates the effective reproduction number (R_e)
#' using the method of Cori et al, but in a Bayesian framework allowing 
#' (a) uncertainty in the serial interval distribution, 
#' (b) imputation of missing symptom onset dates, 
#' (c) extrapolation of infection dates, and (d) adjusting for 
#' delays between diagnosis, symptom onset and infection ("nowcasting")
#'
#' @param donset A vector of symptom onset dates (of class \code{Date}).
#' @param ddiag A vector of diagnosis dates (of class \code{Date}), 
#'     of equal length to \code{donset}; alternatively \code{NULL} 
#'     if working only with symptom onset dates, see 'Details' section.
#' @param local A logical vector, of equal length to \code{donset}
#'     indicating whether the case is local or imported.
#' @param t_from Date (of class \code{Date}) from which to start estimating
#'     effective reproduction numbers (instantaneous or period-specific);
#'    see 'Details' section.
#' @param SI_mean Mean for the serial interval distribution, either a single
#'     value or two values defining a range; see 'Details' section.
#' @param SI_std Standard Deviation for the serial interval distribution, 
#'     either a single value or two values defining a range; see 'Details' section.
#' @param t_window Window length for estimating a rolling instantaneous 
#'     effective reproduction number, see 'Details' section.
#' @param t_breaks Dates defining specific time windows for estimating the 
#'     effective reproduction number, or \code{NULL} for a rolling 
#'     instantaneous reproduction number, see 'Details' section.
#' @param maxDate Maximum date for which to estimate the effective 
#'     reproduction number; if \code{NULL}, the maximum diagnosis date or 
#'     symptom onset date is used.
#' @param weights An optional vector of weights for the cases. If \code{NULL} 
#'     (the default), all cases are given equal weights of 1.
#' @param withInfectTimes Extrapolate infection times and use those for 
#'     estimating the effective reproduction number. If \code{TRUE} 
#'     (the default), a Gamma distribution with mean \code{incub_mean} and 
#'     standard deviation \code{incub_std} will be used to impute the
#'     infection times.
#' @param delayAdjust Adjust for delays in case ascertainment, using the 
#'     distribution of times from symptom onset to diagnosis; 
#'     see 'Details' section
#' @param incub_mean Mean incubation period, used for extrapolating 
#'     infection times from symptom onset dates.
#' @param incub_std Standard deviation of the incubation period distribution,
#'     used for extrapolating infection times from symptom onset dates.
#' @param date0 "Date zero" from which to start counting when converting the
#'     provided Date vectors to integers for MCMC fitting. Also, imputation
#'     of symptom onset and infection dates will start after this date.
#' @param burn Number of burn-in MCMC iterations
#' @param iter Number of sampling MCMC iterations
#'
#' @details \code{estimate_R_Bayes()} allows flexibility in how to estimate
#'     the effective reproduction number R_e, using different function arguments.
#' 
#' At the most basic, a vector of symptom onset dates must be provided 
#' (argument \code{donset}) and on this basis the R_e will be estimated, if 
#' \code{withInfectTimes = FALSE}. If \code{withInfectTimes = TRUE} (the 
#' default), the R_e is estimated on the basis of the infection times 
#' imputed using a Gamma distribution with mean \code{incub_mean} and 
#' standard deviation \code{incub_std}. 
#' 
#' If a vector of diagnosis dates is also provided (argument \code{ddiag}), 
#' then \code{donset} can have missing values, and these will be imputed 
#' using the distribution of the time between symptom onset and diagnosis
#' (assumed to be Gamma, with parameters estimated from the data). 
#' In addition, if diagnosis dates are provided, it is possible to use this
#' distribution (of the time between symptom onset and diagnosis) to 
#' "nowcast" the epidemic, if the argument \code{delayAdjust} is set to 
#' \code{TRUE}: the number of cases with symptom onset at each time is 
#' divided by the probability of ascertainment by date \code{maxDate}.
#' 
#' If argument \code{t_breaks} is \code{NULL} (the default), the function
#' estimates the instantaneous effective reproduction number, using rolling
#' weekly windows between dates \code{t_start} and \code{maxDate}. The size
#' of the window is adjusted by the argument \code{t_window}. 
#' Alternatively, \code{t_breaks} can contain a vector of dates that define
#' (together with \code{t_breaks} and \code{maxDate}) specific periods for
#' which the effective reproduction number will be estimated.
#' 
#' Arguments \code{SI_mean} and \code{SI_std} can receive either one or two
#' values. If one value is provided, the parameter is assumed fixed. 
#' If two values are provided, these define a range from which the parameter
#' is randomly sampled (using a uniform distribution). This allows 
#' to incorporate uncertainty in estimating the serial interval 
#' distribution parameters.
#'
#' @return An object of class 'bayEStim'. This is a list containing the following elements:
#'   \describe{
#'     \item{$model}{The fitted model; an object of class \code{rjags}}
#'
#'     \item{$mcmcSamples}{An object of class \code{coda} containing the MCMC samples for
#'         the various model parameters.}
#'
#'     \item{$config}{A list of all configuration parameters used in the model.}
#'   }
#' 
#' Objects of class 'bayEStim' can be \link[summary.bayEStim]{summarized} and 
#' \link[plot.bayEStim]{plotted}, see the respective methods for details.
#'
#' @references \itemize{
#'  \item Cori A, Ferguson NM, Fraser C, Cauchemez S. A new framework and 
#'  software to estimate time-varying reproduction numbers during epidemics. 
#'  \href{https://academic.oup.com/aje/article-lookup/doi/10.1093/aje/kwt133}{Am J Epidemiol} 
#'  2013;178(9):1505–12 (\href{https://www.ncbi.nlm.nih.gov/pubmed/24043437}{PubMed})
#'  
#'  \item Thompson RN, Stockwin JE, van Gaalen RD, et al. Improved inference of 
#'  time-varying reproduction numbers during infectious disease outbreaks. 
#'  \href{https://linkinghub.elsevier.com/retrieve/pii/S1755-4365(19)30035-0}{Epidemics} 
#'  2019;29:100356 (\href{https://www.ncbi.nlm.nih.gov/pubmed/31624039}{PubMed})
#' }
#'
#' @importFrom rjags jags.model coda.samples
#'
#' @examples
#' # We'll eventually put some examples here
#'
#' @export
estimate_R_Bayes <- function(donset, ddiag, local, t_from, SI_mean, SI_std, 
    t_window = 7, t_breaks = NULL, maxDate = NULL, weights = NULL,
    withInfectTimes = TRUE, delayAdjust = TRUE, incub_mean = 5.1, incub_std = 3.0,
    date0 = as.Date("2020-1-31"), burn = 500, iter = 5000) {
  if (is.null(weights)) weights <- rep(1, length(donset))
  
  # Check data
  if ((!is.null(ddiag) && length(ddiag)!=length(donset)) || 
        length(donset)!=length(local) || length(donset)!=length(weights)) {
    stop("arguments 'ddiag', 'donset', 'local' and 'weights' must have the same length.")
  }
  if (!is.null(ddiag) && !inherits(ddiag, "Date")) stop("Argument 'ddiag' must be of class \"Date\".")
  if (!inherits(donset, "Date")) stop("Argument 'donset' must be of class \"Date\".")
  if (!is.null(ddiag) && sum(is.na(ddiag))>0) stop("Missing values found in argument 'ddiag'.")
  if (sum(is.na(local))>0) stop("Missing values found in argument 'local'.")
  if (sum(!is.na(donset))==0) stop("No non-missing values found in argument 'donset'.")
  if (is.null(ddiag) && sum(is.na(donset))>0) stop("Missing values found in argument 'donset'.")
  
  # Config 
  if (is.null(maxDate)) {
    maxDate <- if (is.null(ddiag)) max(donset, na.rm=TRUE) else max(ddiag, na.rm=TRUE)
  }
  config <- list(
    data = data.frame(donset, local),
    SI = list(mean = SI_mean, std = SI_std),
    t_from = t_from, t_window = t_window, t_breaks = t_breaks, 
    withInfectTimes = withInfectTimes, delayAdjust = delayAdjust,
    date0 = date0, maxDate = maxDate
  )
  if (!is.null(ddiag)) {
    config$data$ddiag <- ddiag
    if (withInfectTimes) {
      config$incub <- list(mean = incub_mean, std = incub_std)
    }
  }

  # Serial interval
  if (length(SI_mean)==1) SI_mean <- rep(SI_mean, 2)
  if (length(SI_std)==1) SI_std <- rep(SI_std, 2)
  
  # Reorder, so that missing onset dates come last
  if (!is.null(ddiag)) ddiag <- as.integer(ddiag[order(donset)] - date0)
  local <- local[order(donset)]
  weights <- weights[order(donset)]
  donset <- as.integer(donset[order(donset)] - date0)
  maxDate <- as.integer(config$maxDate - date0)
  
  # Set up data for JAGS
  datJ <- list(
    donset = donset, weights = weights,
    local = local, maxDate = maxDate,
    SI_mean_lim = SI_mean, SI_std_lim = SI_std,
    SI_k = c(0:(maxDate-1), 0),
    R_mean_prior = 5, R_std_prior = 5
  )
  if (!is.null(ddiag)) {
    datJ$ddiag <- ddiag
    datJ$onsetToDiag <- ddiag - donset
    datJ$N <- length(donset)
    datJ$firstMissing <- which(is.na(donset))[1]
  }
  if (is.null(t_breaks)) {
    datJ$t_start <- as.integer(t_from - date0):(maxDate + 1 - t_window)
    datJ$t_end <- (as.integer(t_from - date0) + t_window - 1):maxDate
  } else {
    datJ$t_start <- as.integer(c(t_from, t_breaks) - date0)
    datJ$t_end <- c(as.integer(t_breaks - date0)-1, maxDate)
  }
  datJ$Nper <- length(datJ$t_start)
  
  if (withInfectTimes) {
    datJ$incub_mean <- incub_mean
    datJ$incub_std <- incub_std
    datJ$N <- length(donset)
  }

  jagsModel <- if (is.null(ddiag)) {
    JAGSmodel(withInfectTimes = withInfectTimes, single = TRUE)
  } else {
    JAGSmodel(withInfectTimes, delayAdjust, withMissing = sum(is.na(donset))>0)
  }
  
  model <- jags.model(file = textConnection(jagsModel), 
    data = datJ, n.chains = 4, n.adapt = burn)
  
  parameters_to_sample <- c("R", "I_local", "I_imported")
  if (!is.null(ddiag) & delayAdjust) {
    parameters_to_sample <- c(parameters_to_sample, "shp_onsetToDiag", "rate_onsetToDiag")
    if (withInfectTimes) {
      parameters_to_sample <- c(parameters_to_sample, "shp_infectToDiag", "rate_infectToDiag")
    }
  }
  
  mcmcSamples <- coda.samples(model, variable.names=parameters_to_sample, n.iter = iter)

  res <- list(
    model = model, mcmcSamples = mcmcSamples, config = config
  )
  class(res) <- "bayEStim"
  
  return(res)
}

