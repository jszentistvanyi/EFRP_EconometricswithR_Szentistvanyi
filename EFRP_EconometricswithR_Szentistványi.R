# EFRP
# First R homework
# Szentistványi János

library(forecast)

# Parameters for simulation:
nModels <- 4        # Number of different models to be simulated
nSim <- list(1000, 1000, 1000, 1000)      # Number of simulations for each model
tslength <- list(1000, 1000, 500, 750)        # Length of each simulated time series for each model
ARc <- list(c(0.4, -0.3), c(0.4, -0.3), c(0.6, -0.3), c(0.33, -0.12, 0.41, -0.11))      # Vectors of AR coefficients for each model
MAc <- list(c(-0.2, 0.4, 0.3), c(-0.2, 0.4, 0.3), c(0.5, 0.25), c(0.23, -0.1, -0.2, 0.17))       # Vectors of MA coefficients for each model
noise <- list(0.1, 0.5, 0.25, 0.2)            # Noise parameter for each model
order <- list()     # Generating order vector from the model specification (for comparison of the fitted and specified model
for (i in 1:nModels) {
  order[[i]] <- c(length(ARc[[i]]), 0, length(MAc[[i]]))
}  


# Fit parameters: maximum p (for AR) and q (for MA) values, for each model
arMax <- list(3, 3, 4, 5)
maMax <- list(3, 3, 4, 4)

# Defining the simulation function
# Generating a list containing lists: for example 3*1000 dimension, where each element is a time series object
simulateModels <- function(nModels, ARc, MAc, noise, tslength, nSim) {
  simTseries <- list()
  for (i in 1:nModels) {
    simTseries[[i]] <- list()
    for (j in 1:nSim[[i]]) {
      simTseries[[i]][[j]] <- arima.sim(n = tslength[[i]], 
                                        model = list(ar = ARc[[i]], ma = MAc[[i]]),
                                        rand.gen = rnorm,
                                        sd = noise[[i]])
    }
  }
  invisible(simTseries)
}

# Defining the model fitting function
# Loops through the simulations, finds the best model for each time series by AIC and BIC as well
# Counting how many times it gives the same results as the specified order of the time series

fitBestModel <- function(simTseries, arMax, maMax) {
  modelFit <- list()
  for (i in 1:nModels) {
    modelFit[[i]] <- list(0,0)
    for (j in 1:nSim[[i]]) {
      fit_aic <- auto.arima(simTseries[[i]][[j]], ic = "aic", max.p = arMax[[i]], max.q = maMax[[i]], seasonal=FALSE)
      if (isTRUE(all.equal(c(fit_aic$arma[1], fit_aic$arma[6], fit_aic$arma[2]), order[[i]]))) {
        modelFit[[i]][[1]] = modelFit[[i]][[1]] + 1
      }
      fit_bic <- auto.arima(simTseries[[i]][[j]], ic = "bic", max.p = arMax[[i]], max.q = maMax[[i]], seasonal=FALSE)
      if (isTRUE(all.equal(c(fit_bic$arma[1], fit_bic$arma[6], fit_bic$arma[2]), order[[i]]))) {
        modelFit[[i]][[2]] = modelFit[[i]][[2]] + 1
      }
    }  
  }
  invisible(modelFit)
}

# Defining a main function that runs both the sumulation and fit functions and generates a table with the results.
simulateAndFitResults <- function(nModels, ARc, MAc, noise, tslength, nSim, arMax, maMax) {
  # Calling the simulation function with the given parameters
  simTseries <- simulateModels(nModels=nModels, ARc = ARc, MAc = MAc, noise = noise, 
                               nSim = nSim, tslength = tslength)
  # Calling the fit function with the given parameters
  modelFit <- fitBestModel(simTseries= simTseries, arMax=arMax, maMax=maMax)
  
  # Generating table with the results
  results_list <- list()
  for (i in 1:nModels) {
    results_list[[i]] <- list()
    results_list[[i]][[1]] <- i
    results_list[[i]][[2]] <- order[[i]]
    results_list[[i]][[3]] <- tslength[[i]]
    results_list[[i]][[4]] <- noise[[i]]
    results_list[[i]][[5]] <- nSim[[i]]
    results_list[[i]][[6]] <- modelFit[[i]][[1]]
    results_list[[i]][[7]] <- round(modelFit[[i]][[1]]/nSim[[i]]*100,2)
    results_list[[i]][[8]] <- modelFit[[i]][[2]]
    results_list[[i]][[9]] <- round(modelFit[[i]][[2]]/nSim[[i]]*100,2)
    results_list[[i]][[10]] <- arMax[[i]]
    results_list[[i]][[11]] <- maMax[[i]]
  }
  results <- matrix(unlist(results_list), nrow = nModels, byrow= TRUE)
  results <- data.frame(results)
  colnames(results) <- c("Model number", "Specified AR order", "(Specified D order)=0",
                         "Specified MA order", "Lenght of each time series","Noise in the simulation",
                         "Number of simulations", "Number of accurate AIC fits",
                         "% of accurate AIC fits", "Number of accurate BIC fits", "% of accurate BIC fits",
                         "Max p when fitting", "Max q when fitting")
  invisible(results)
}

# Running the main functions that calls both functions froma above, and also creates a table with the results
# Using the parameters specified in the first part on the top
results<- simulateAndFitResults(nModels, ARc, MAc, noise, tslength, nSim, arMax, maMax)
View(results)