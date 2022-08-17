# The following code was designed to obtain MMLE, empirical Bayes, and hierarchical Bayes results.

#################################################
###### Empirical Bayes estimation function: #####
#################################################

EB_estimation_function <- function(EB_results_function_input){
  Y <- EB_results_function_input[[2]]
  condition_variables <- EB_results_function_input[[3]]
  Q <- EB_results_function_input[[4]]
  I <- condition_variables[[1]]
  K <- condition_variables[[2]]
  J <- condition_variables[[3]]
  D <- ncol(Q)
  colnames_Y <- c()
  for (item in 1:I){
    colnames_Y <- c(colnames_Y, paste("Item",item,sep=" "))
  }
  
  colnames(Y) <- colnames_Y
  mirt_Y <- mirt(Y, 1, itemtype='graded', method="EM", SE=TRUE)
  summary_mirt <- coef(mirt_Y, printSE=TRUE, as.data.frame=TRUE)
  if (dim(summary_mirt)[1]==(I*K + 2)){
    mirt_estimates <- matrix(nrow=I,ncol=K)
    mirt_SE <- matrix(nrow=I,ncol=K)
    for (item in 1:I){
      for (parameter in 1:K){
        mirt_estimates[item,parameter] <- summary_mirt[((item-1)*K)+parameter,1]
        mirt_SE[item,parameter] <- summary_mirt[((item-1)*K)+parameter,2]
      }}
    mirt_estimates[,2:K] <- -mirt_estimates[,2:K]
    EB_regression_estimates <- matrix(nrow=(D+1),ncol=K)
    EB_regression_SE <- c()
    for (parameter in 1:K){
      regression_structure <- cbind(mirt_estimates[,parameter],Q)
      colnames(regression_structure) <- c("Parameter",paste("Covariate",seq(1,D),sep='_'))
      regression <- lm(as.formula(paste(colnames(regression_structure)[1],
                                        paste(c(1, colnames(regression_structure)[2:(D+1)]), collapse=" + "), sep=" ~ ")),
                       data = data.frame(regression_structure))
      summary_regression <- summary(regression)
      EB_regression_estimates[,parameter] <- unname(summary_regression$'coefficients'[,1])
      regression_residue <- unname(c(resid(regression)))
      EB_regression_SE[parameter] <- sqrt(sum((regression_residue)^2) / (I - (D+1)))
    }
    EB_estimates <- matrix(nrow=I,ncol=K)
    EB_SE <- matrix(nrow=I,ncol=K)
    for (item in 1:I){
      for (parameter in 1:K){
        EB_estimates[item,parameter] <- (mirt_estimates[item,parameter]*(mirt_SE[item,parameter]^-2) +
                                           sum(Q[item,]*EB_regression_estimates[,parameter])*(EB_regression_SE[parameter]^-2)) / (mirt_SE[item,parameter]^-2 +
                                                                                                                                    EB_regression_SE[parameter]^-2)
        EB_SE[item,parameter] <- sqrt(1/(mirt_SE[item,parameter]^-2 + EB_regression_SE[parameter]^-2))
      }}
    return(list(mirt_estimates, mirt_SE, EB_regression_estimates, EB_regression_SE, EB_estimates, EB_SE))
  } else {
    return(list("Error, MIRT unable to estimate parameters.","Error, MIRT unable to estimate parameters.",
                "Error, MIRT unable to estimate parameters.","Error, MIRT unable to estimate parameters.",
                "Error, MIRT unable to estimate parameters.","Error, MIRT unable to estimate parameters."))
  }}

HB_with_covariates_estimation_function <- function(HB_results_function_input){
  stan_iterations <- HB_results_function_input[[1]]
  Y <- HB_results_function_input[[2]]
  condition_variables <- HB_results_function_input[[3]]
  Q <- HB_results_function_input[[4]] # Now included in input
  Q_type <- HB_results_function_input[[5]]
  I <- condition_variables[[1]]
  K <- condition_variables[[2]]
  J <- condition_variables[[3]]
  D <- ncol(Q)
  Y_long <- c()
  person_long <- c()
  item_long <- c()
  a <- 0
  
  for (person in 1:J){
    for (item in 1:I){
      a <- a + 1
      Y_long[a] <- Y[person,item]
      person_long[a] <- person
      item_long[a] <- item
    }}
  Y_long <- Y_long + 1 # Puts responses on a scale of 1-5 instead of 0-4.
  data_stan <- list(number_categories=K, number_persons=J, number_items=I, Y=Y_long, person=person_long,
                    item=item_long, number_responses=length(Y_long), Q=Q, number_covariates=D)
  if (Q_type=='ME'){
    fit_stan <- sampling(GRM_stan_ME_model, data = data_stan, chains = 1, iter=stan_iterations)
  } else if (Q_type=='NME'){
    fit_stan <- sampling(GRM_stan_NME_model, data = data_stan, chains = 1, iter=stan_iterations)
  }
  summary_stan <- summary(fit_stan)
  summary_stan <- summary_stan$summary
  HB_estimates <- matrix(nrow=I,ncol=K)
  HB_SD <- matrix(nrow=I,ncol=K)
  HB_regression_estimates <- matrix(nrow=(D+1),ncol=K)
  HB_regression_SD <- c()
  start_index <- J + (I*K) + 1
  end_index <- start_index + D
  HB_regression_estimates[,1] <- summary_stan[start_index:end_index, 6]
  HB_regression_SD[1] <- summary_stan[(end_index+1),6]
  start_index <- J + (I*K) + D + 3
  end_index <- start_index + (D+1)*(K-1) - 1
  HB_regression_estimates[,2:K] <- matrix(summary_stan[start_index:end_index, 6],ncol=(K-1),byrow=FALSE)
  HB_regression_SD[2:K] <- summary_stan[((end_index + 1):(end_index + K - 1)),6]
  for (item in 1:I){
    HB_estimates[item,1] <- summary_stan[J+item,6] # Alpha, Median
    HB_SD[item,1] <- summary_stan[J+item,3] # Alpha, SD
    for (threshold in 2:K){
      HB_estimates[item,threshold] <- summary_stan[(J+I+(K-1)*(item-1)+(threshold-1)),6] # Median
      HB_SD[item,threshold] <- summary_stan[(J+I+(K-1)*(item-1)+(threshold-1)),3] # SD
    }}
  return(list(HB_estimates, HB_SD, HB_regression_estimates, HB_regression_SD, summary_stan))
}

############################################################################
###### Hierarchical Bayes without item covariates estimation function: #####
############################################################################

HB_without_covariates_estimation_function <- function(HB_results_function_input){
  stan_iterations <- HB_results_function_input[[1]]
  Y <- HB_results_function_input[[2]]
  condition_variables <- HB_results_function_input[[3]]
  I <- condition_variables[[1]]
  K <- condition_variables[[2]]
  J <- condition_variables[[5]]
  Y_long <- c()
  person_long <- c()
  item_long <- c()
  a <- 0
  for (person in 1:J){
    for (item in 1:I){
      a <- a + 1
      Y_long[a] <- Y[person,item]
      person_long[a] <- person
      item_long[a] <- item
    }}
  Y_long <- Y_long + 1 # Puts responses on a scale of 1-5 instead of 0-4.
  
  data_stan <- list(number_categories=K, number_persons=J, number_items=I, Y=Y_long, person=person_long,
                    item=item_long, number_responses=length(Y_long))
  fit_stan <- sampling(GRM_stan_without_covariates_model, data = data_stan, chains = 1, iter=stan_iterations)
  summary_stan <- summary(fit_stan)
  summary_stan <- summary_stan$summary
  HB_estimates <- matrix(nrow=I,ncol=K)
  HB_SD <- matrix(nrow=I,ncol=K)
  HB_regression_estimates <- c()
  HB_regression_SD <- c()
  start_index <- J + (I*K) + 1
  HB_regression_estimates[1] <- summary_stan[start_index, 6]
  HB_regression_SD[1] <- summary_stan[(start_index + 1),6]
  HB_regression_estimates[2:K] <- summary_stan[(start_index + 2):(start_index + K),6]
  HB_regression_SD[2:K] <- summary_stan[(start_index + K + 1):(start_index + 2*K - 1),6]
  for (item in 1:I){
    HB_estimates[item,1] <- summary_stan[J+item,6] # Alpha, Median
    HB_SD[item,1] <- summary_stan[J+item,3] # Alpha, SD
    for (threshold in 2:K){
      HB_estimates[item,threshold] <- summary_stan[(J+I+(K-1)*(item-1)+(threshold-1)),6] # Median
      HB_SD[item,threshold] <- summary_stan[(J+I+(K-1)*(item-1)+(threshold-1)),3] # SD
    }}
  return(list(HB_estimates, HB_SD, HB_regression_estimates, HB_regression_SD, summary_stan))
}

####################################
###### GRM Stan Model Code #########
####################################
library(mirt)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
GRM_ME_Stan <- "
data{
int<lower=0> number_categories;
int<lower=0> number_persons;
int<lower=0> number_items;
int<lower=0> number_responses;
int<lower=1, upper=number_categories> Y[number_responses];
int<lower=1, upper=number_persons> person[number_responses];
int<lower=1, upper=number_items> item[number_responses];
int<lower=0> number_covariates;
int<lower=0, upper=1> Q[number_items, number_covariates];
}
parameters{
vector[number_persons] theta; //latent variable
real<lower=0> alpha[number_items]; //item discrimination
ordered[number_categories - 1] beta[number_items]; //category difficulty
vector[1 + number_covariates] gamma_alpha;
real<lower=0> phi_alpha;
vector[1 + number_covariates] gamma_beta[number_categories - 1];
real<lower=0> phi_beta[number_categories - 1];
}
model{
for (j in 1:number_persons){
theta[j] ~ normal(0,1);
}
phi_alpha ~ cauchy(0,10);
for (k in 1:number_categories-1){
phi_beta[k] ~ cauchy(0,10);
}
}}
for (a in 1:number_responses){ // For each data point
Y[a] ~ ordered_logistic(theta[person[a]]*alpha[item[a]], beta[item[a]]);
}
}
"

GRM_without_covariates_Stan <- "
data{
int<lower=0> number_categories;
int<lower=0> number_persons;
int<lower=0> number_items;
int<lower=0> number_responses;
int<lower=1, upper=number_categories> Y[number_responses];
int<lower=1, upper=number_persons> person[number_responses];
int<lower=1, upper=number_items> item[number_responses];
}
parameters{
vector[number_persons] theta; //latent variable
real<lower=0> alpha[number_items]; //item discrimination
ordered[number_categories - 1] beta[number_items]; //category difficulty
real<lower=0> gamma_alpha;
real<lower=0> phi_alpha;
real gamma_beta[number_categories - 1];
real<lower=0> phi_beta[number_categories - 1];
}
model{
for (j in 1:number_persons){
theta[j] ~ normal(0,1);
}
phi_alpha ~ cauchy(0,10);
for (k in 1:number_categories-1){
phi_beta[k] ~ cauchy(0,10);
}
for (i in 1:number_items){
alpha[i] ~ normal(gamma_alpha, phi_alpha);
}
for (i in 1:number_items){
for (k in 1:(number_categories-1)){
beta[i,k] ~ normal(gamma_beta[k], phi_beta[k]);
}}
gamma_alpha ~ normal(0,10);
for (k in 1:(number_categories-1)){
gamma_beta[k] ~ normal(0,10);
}
for (a in 1:number_responses){ // For each data point
Y[a] ~ ordered_logistic(theta[person[a]]*alpha[item[a]], beta[item[a]]);
}
}
"
GRM_stan_ME_model <- stan_model(model_code = GRM_ME_Stan)
GRM_stan_NME_model <- stan_model(model_code = GRM_NME_Stan)
GRM_stan_without_covariates_model <- stan_model(model_code = GRM_without_covariates_Stan)

##########################################################################
###### Obtaining empirical Bayes & hierarchical Bayes estimates #########
##########################################################################
stan_iterations <- 2000 # Number of iterations stan will run for.
Y <- Y_matrix # Input the data (Y).
I <- 24 # Number of items
K <- 5 # Number of categeories
J <- 2000 # Number of persons
condition_variables <- list(I, K, J)

Q <- Q_matrix # Input the Q-matrix
Q_type <- "NME" # Either "ME" (for mutually exclusive) or "NME" (for non-mutually exclusive)
estimation_input <- list(stan_iterations, Y, condition_variables, Q, Q_type)
EB_estimation_output <- EB_estimation_function(estimation_input)
# Step 1: EB_estimation_output[[1]]: Maximum likelihood estimates
# EB_estimation_output[[2]]: Maximum likelihood SEs
# Step 2: EB_estimation_output[[3]]: Regression estimates
# EB_estimation_output[[3]]: Regression SEs
# Step 3: EB_estimation_output[[5]]: Empirical Bayes estimates
# EB_estimation_output[[6]]: Empirical Bayes SEs
HB_with_covariates_estimation_output <- HB_with_covariates_estimation_function(estimation_input)
# HB_with_covariates_estimation_output[[1]]: Hierarchical Bayes estimates
# HB_with_covariates_estimation_output[[2]]: Hierarchical Bayes SDs
# HB_with_covariates_estimation_output[[3]]: Hierarchical Bayes regression estimates
# HB_with_covariates_estimation_output[[4]]: Hierarchical Bayes RSD estimates
# HB_with_covariates_estimation_output[[5]]: Full Stan output
HB_without_covariates_estimation_output <- HB_without_covariates_estimation_function(estimation_input)