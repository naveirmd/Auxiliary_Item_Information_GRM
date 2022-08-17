# The following code uses classical item discrimination and thresholds to obtain the linear regression
# RSD estimates necessary to utilize the method selection guideline provided in Figure F1.

########################################
###### RSD estimation function #########
########################################

RSD_function <- function(RSD_function_input){
  library(lme4)
  library(psychometric)
  RSD <- c()
  Y <- RSD_function_input[[1]]
  condition_variables <- RSD_function_input[[2]]
  I <- condition_variables[[1]]
  K <- condition_variables[[2]]
  D <- condition_variables[[3]]
  J <- condition_variables[[5]]
  # If you have a unique item covariate structure, substitute D (above)
  # with the correct covariate structure, and remove the following code chunk.
  covariates <- matrix(0,nrow=I,ncol=D)
  items_per_group <- I/(D+1)
  covariates <- matrix(nrow=I,ncol=D)
  for (covariate in 1:D){
    covariates[,covariate] <- c(rep(0,items_per_group*covariate),rep(1,items_per_group),
                                rep(0,items_per_group*(D-covariate)))
  }
  # Obtains item discriminations.
  discrimination <- unname(unlist(item.exam(Y, discrim=TRUE)$Discrimination))
  parameter_types <- matrix(0,nrow=I*K,ncol=K-1)
  # Creates a covariate matrix for all item parameter types,
  # to be used with the linear regression later on in this function.
  for (parameter in 2:K){
    parameter_types[((parameter-1)*I+1):(parameter*I),parameter-1] <- 1
  }
  for (threshold in 1:(K-1)){
    Y_threshold <- Y
    Y_threshold[Y_threshold < threshold] <- 0
    Y_threshold[Y_threshold >= threshold] <- 1
    item_means <- colMeans(Y_threshold)
    difficulty <- c()
    parameter_estimates <- discrimination
    for (item in 1:I){
      difficulty[item] <- min(4,max(-4, -qlogis(item_means[item], scale=discrimination[item])))
      # Obtains each item threshold, with a minimum of -4 and maximum of +4.
    }
    parameter_estimates <- c(parameter_estimate, difficulty)
  }
  regression_data <- data.frame(parameter_estimates, parameter_types, covariates)
  # Combines item parameter estimates, parameter covariates, and item covariates for running regression.
  regression_data_names <- c("Estimate")
  for (parameter in 1:(K-1)){
    regression_data_names <- c(regression_data_names, paste('P',parameter,sep=''))
  }
  for (d in 1:D){
    regression_data_names <- c(regression_data_names, paste('X',d,sep=''))
  }
  names(regression_data) <- regression_data_names
  regression_equation <- "Estimate ~ 1"
  for (parameter in 1:(K-1)){
    regression_equation <- paste(regression_equation, paste("+ P",parameter,sep=''),sep=' ')
  }
  for (covariate in 1:D){
    regression_equation <- paste(regression_equation, paste("+ X",covariate,sep=''), sep=' ')
  }
  regression <- lm(formula = as.formula(regression_equation), data = regression_data)
  # Runs linear regression on CTT item parameter estimates.
  regression_summary <- summary(regression)
  RSD <- regression_summary$sigma
  return(RSD)
}

########################################
####### Obtaining RSD estimate #########
########################################

setwd("C:\\Users\\myname\\DataLocation") # Replace this with the directory where the data is located
data <- read.table("data.txt",sep="\t") # Replace "data.txt" with the name of the data file
number_items <- 24
number_categories <- 5
number_covariates <- 5
number_persons <- 100
condition_variables <- list(number_items, number_categories, number_covariates, number_persons)
coefficients <- c(-0.2703644864, 0.8808019178, -0.0001940328, 0.0001785366)
# Based on our simulation study data.
CTT_RSD <- RSD_function(list(data, condition_variables))
RSD_estimate <- coefficients[1] + coefficients[2]*CTT_RSD +
  coefficients[3]*number_items + coefficients[4]*number_persons # Final result.
