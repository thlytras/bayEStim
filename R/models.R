#' Return the appropriate JAGS model to use
#'
#' This function formats and returns the appropriate JAGS model according
#' to the details of the estimation desired, i.e. whether the infection times 
#' or symptom times are used, whether a delay adjustment is desired, whether 
#' missing values exist in the symptom onset dates vector, and whether we 
#' work only with symptom onset dates and not diagnosis dates.
#'
#' @param withInfectTimes Use extrapolated infection times for estimating R.
#' @param delayAdjust Make a delay adjustment, based on the distribution of 
#'     times from symptom onset to diagnosis, or infection to diagnosis.
#' @param withMissing Does the symptom onset dates vector contain missing 
#'     values? If so, dates of diagnosis should also be supplied.
#' @param single Do we work only with symptom onset dates? (No diagnosis 
#'     dates). If \code{TRUE}, then arguments \code{delayAdjust} and 
#'     \code{withMissing} are ignored (no delay adjustment possible and no
#'     missing symptom onset dates allowed).
#'
#' @return An character vector of length 1, containing the JAGS model to use
#'
#' @export
JAGSmodel <- function(withInfectTimes=TRUE, delayAdjust=TRUE, withMissing=TRUE, single=FALSE) {

  chunkA_1 <- sprintf(
'data {
# Calculate R distribution priors
  R_shp_prior <- (R_mean_prior / R_std_prior)^2
  R_rate_prior <- R_mean_prior / R_std_prior^2
  
}
model {
  shp_onsetToDiag ~ dunif(0, 10)
  rate_onsetToDiag ~ dunif(0, 10)
  
  for (i in 1:N) {
    onsetToDiag[i] ~ dgamma(shp_onsetToDiag, rate_onsetToDiag)
  }
%s  

# Tabulate local and imported cases
  for (t in 1:maxDate) {
    I_local[t] <- round(sum(donset==t && local==1)%s)
    I_imported[t] <- round(sum(donset==t && local==0)%s)
  }
# Make a reversed vector of total cases 
# (used below, in the overall infectivity calculation)
  for(t in 1:maxDate) {
    revI[t] <- I_local[maxDate-t+1] + I_imported[maxDate-t+1]
  }',
  ifelse(withMissing, "  for (i in firstMissing:N) {
    donset[i] <- ddiag[i] - round(onsetToDiag[i])
  }\n", ""),
  ifelse(delayAdjust, " / pgamma((maxDate - t + 1), shp_onsetToDiag, rate_onsetToDiag)", ""),
  ifelse(delayAdjust, " / pgamma((maxDate - t + 1), shp_onsetToDiag, rate_onsetToDiag)", ""))
  

  chunkA_2 <- sprintf('data {
# Calculate R distribution priors
  R_shp_prior <- (R_mean_prior / R_std_prior)^2
  R_rate_prior <- R_mean_prior / R_std_prior^2
# Calculate incubation period shape/rate
  incub_shp <- (incub_mean / incub_std)^2
  incub_rate <- incub_mean / incub_std^2
}
model {
  shp_onsetToDiag ~ dunif(0, 10)
  rate_onsetToDiag ~ dunif(0, 10)
%s

  for (i in 1:N) {
    infectToOnset[i] ~ dgamma(incub_shp, incub_rate)
    dinfect[i] <- donset[i] - round(infectToOnset[i])
    onsetToDiag[i] ~ dgamma(shp_onsetToDiag, rate_onsetToDiag)
  }
%s

# Tabulate local and imported cases
  for (t in 1:maxDate) {
    I_local[t] <- round(sum(dinfect==t && local==1)%s)
    I_imported[t] <- round(sum(dinfect==t && local==0)%s)
  }
# Make a reversed vector of total cases 
# (used below, in the overall infectivity calculation)
  for(t in 1:maxDate) {
    revI[t] <- I_local[maxDate-t+1] + I_imported[maxDate-t+1]
  }',
  ifelse(delayAdjust, '
  shp_infectToDiag <- (mean(ddiag-dinfect) / sd(ddiag-dinfect))^2
  rate_infectToDiag <- mean(ddiag-dinfect) / sd(ddiag-dinfect)^2', ""),
  ifelse(withMissing, "  for (i in firstMissing:N) {
    donset[i] <- ddiag[i] - round(onsetToDiag[i])
  }\n", ""),
  ifelse(delayAdjust, " / pgamma((maxDate - t + 1), shp_infectToDiag, rate_infectToDiag)", ""),
  ifelse(delayAdjust, " / pgamma((maxDate - t + 1), shp_infectToDiag, rate_infectToDiag)", ""))
  
  
  chunkA_3 <- 'data {
# Calculate R distribution priors
  R_shp_prior <- (R_mean_prior / R_std_prior)^2
  R_rate_prior <- R_mean_prior / R_std_prior^2
  
# Tabulate local and imported cases
  for (t in 1:maxDate) {
    I_local[t] <- round(sum(donset==t && local==1))
    I_imported[t] <- round(sum(donset==t && local==0))
  }
  
# Make a reversed vector of total cases 
# (used below, in the overall infectivity calculation)
  for(t in 1:maxDate) {
    revI[t] <- I_local[maxDate-t+1] + I_imported[maxDate-t+1]
  }
}
model {'
  
  
  chunkA_4 <- 'data {
# Calculate R distribution priors
  R_shp_prior <- (R_mean_prior / R_std_prior)^2
  R_rate_prior <- R_mean_prior / R_std_prior^2
# Calculate incubation period shape/rate
  incub_shp <- (incub_mean / incub_std)^2
  incub_rate <- incub_mean / incub_std^2
}
model {
  shp_onsetToDiag ~ dunif(0, 10)
  rate_onsetToDiag ~ dunif(0, 10)

  for (i in 1:N) {
    infectToOnset[i] ~ dgamma(incub_shp, incub_rate)
    dinfect[i] <- donset[i] - round(infectToOnset[i])
  }

# Tabulate local and imported cases
  for (t in 1:maxDate) {
    I_local[t] <- round(sum(dinfect==t && local==1))
    I_imported[t] <- round(sum(dinfect==t && local==0))
  }
  
# Make a reversed vector of total cases 
# (used below, in the overall infectivity calculation)
  for(t in 1:maxDate) {
    revI[t] <- I_local[maxDate-t+1] + I_imported[maxDate-t+1]
  }
'
  
  
  chunkB <- '

# Distributions for SI_mean and SI_std
  SI_mean ~ dnorm((SI_mean_lim[1]+SI_mean_lim[2])/2, 0.00001)T(SI_mean_lim[1], SI_mean_lim[2])
  SI_std ~ dnorm((SI_std_lim[1]+SI_std_lim[2])/2, 0.00001)T(SI_std_lim[1],min(SI_std_lim[2],SI_mean))

# Calculate the discretized serial interval
  SIdGshp <- ((SI_mean - 1)/SI_std)^2  # Gamma distribution shape of the SI
  SIdGrate <- (SI_mean - 1)/SI_std^2  # Gamma distribution rate of the SI
  si_distr <-  SI_k * pgamma(SI_k, SIdGshp, SIdGrate) + 
    (SI_k - 2) * pgamma(SI_k - 2, SIdGshp, SIdGrate) - 
    2 * (SI_k - 1) * pgamma(SI_k - 1, SIdGshp, SIdGrate) +
    SIdGshp * (1/SIdGrate) * (2 * pgamma(SI_k - 1, SIdGshp + 1, SIdGrate) - 
        pgamma(SI_k - 2, SIdGshp + 1, SIdGrate) - 
        pgamma(SI_k, SIdGshp + 1, SIdGrate))

# Calculate the "overall infectivity" due to previously infected individuals
  for (t in 2:maxDate) { # NOTE: needs a NA in position 1
    lambda[t] <- sum(si_distr[1:t] * revI[(maxDate+1-t):maxDate])
  }
  
# Calculate Gamma shape and rate of R_t
  final_mean_si <- sum(si_distr * (0:maxDate))
  for (t in 1:Nper) {
    R_shp_posterior[t] <- ifelse(t_end[t]>final_mean_si,
        R_shp_prior + sum(I_local[t_start[t]:t_end[t]]), -100)
    R_rate_posterior[t] <- ifelse(t_end[t]>final_mean_si,
        (R_rate_prior + sum(lambda[t_start[t]:t_end[t]])), -100)
  }

# Calculate the distribution of R_t
  for (t in 1:Nper) {
    R[t] ~ dgamma(R_shp_posterior[t], R_rate_posterior[t])
  }
}'
  
  model <- if (single) {
    if (withInfectTimes) paste0(chunkA_4, chunkB) else paste0(chunkA_3, chunkB)
  } else {
    if (withInfectTimes) paste0(chunkA_2, chunkB) else paste0(chunkA_1, chunkB)
  }
  
  return(model)
}


