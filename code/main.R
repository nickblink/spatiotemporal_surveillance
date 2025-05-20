library(MASS)
library(Matrix)
library(dplyr)
library(lubridate)
library(ggplot2)
library(rstan)

source('code/functions.R')
rstan_options(auto_write = TRUE)

#### Make the parameters for DGP, missingness, and model fitting ####
# Weinberger-Fulcher DGP Arguments
WF_DGP_arguments <- list(district_sizes = c(4, 6, 10), # size of the districts
                      R = 2, # number of simulated data sets
                      seed = 1, # random seed
                      end_date = '2020-01-01', # last date of data created.
                      b0_mean = 6, # mean value b0 is sampled from.
                      b1_mean = -0.25, # mean value b1 is sampled from.
                      type = 'WF', # type of DGP: 'WF', 'freqGLM', or 'CAR'.
                      family = 'negbin' # family of DGP
)

# freqGLM DGP
freqGLM_DGP_arguments <- list(district_sizes = c(4, 6, 10), # size of the districts
                          R = 2, # number of simulated data sets
                          seed = 1, # random seed
                          end_date = '2020-01-01', # last date of data created.
                          b0_mean = 5.5, # mean value b0 is sampled from.
                          b1_mean = -0.25, # mean value b1 is sampled from.
                          type = 'freqGLM', # type of DGP: 'WF', 'freqGLM', or 'CAR'.
                          family = 'negbin', # family of DGP
                          rho = 0.2, # spatial correlation of freqGLM process in DGP.
                          alpha = 0.2) # temporal correlation of freqGLM process in DGP.

# CAR DGP
CAR_DGP_arguments <- list(district_sizes = c(4, 6, 10), # size of the districts
                          R = 2, # number of simulated data sets
                          seed = 1, # random seed
                          end_date = '2020-01-01', # last date of data created.
                          b0_mean = 6, # mean value b0 is sampled from.
                          b1_mean = -0.25, # mean value b1 is sampled from.
                          type = 'CAR', # type of DGP: 'WF', 'freqGLM', or 'CAR'.
                          family = 'negbin', # family of DGP
                          rho = 0.3, # spatial correlation of CAR process in DGP.
                          alpha = 0.3, # temporal correlation of CAR process in DGP.
                          tau2 = 0.25) # variance of CAR process in DGP

# parameters for missingness
missingness_mcar <- list(missingness = 'mcar',
                         p = 0.2)
missingness_mar <- list(missingness = 'mar', 
                        p = 0.2,
                        alpha_MAR = 0.7,
                        rho_MAR = 0.7,
                        tau2_MAR = 4)
missingness_mnar <- list(missingness = 'mnar', 
                         p = 0.2,
                         gamma = 1)

# parameters used in model fitting. These are set very low right now to ensure that the code runs properly. To conduct the experiments, these need to be set to R_PI = 200 and burnin = 1000, n.sample = 2000.
model_list <- list(WF = list(model = 'WF',
                             params = list(R_PI = 10)),
                   freqGLM = list(model = 'freqGLM',
                                  params = list(R_PI = 10)),
                   CAR = list(model = 'CAR',
                              params = list(burnin = 10,
                                            n.sample = 20)))

#### Run through 9 experiments for the 3 DGPs and 3 missingness patterns ####
# whether to run the code in parallel.
parallelize <- F

# storing the full results.
full_results <- list()
iter = 1

# cycle through each experiment.
for (DGP in c('WF_DGP_arguments', 'freqGLM_DGP_arguments', 'CAR_DGP_arguments')) {
  for (missingness in c('missingness_mcar', 'missingness_mar', 'missingness_mnar')) {
    print(sprintf('---Running experiment for DGP %s and %s---', DGP, missingness))
    
    # Get the actual objects
    DGP_args <- get(DGP)
    missingness_args <- get(missingness)
    
    # Combine the lists (preserving keys)
    params <- c(DGP_args, missingness_args)

    # simulate the data
    data_list <- do.call(simulate_data, params)
    
    # initialize the error catcher for each run
    errors <- list(freqGLM = data.frame(i = NULL, error = NULL),
                   WF = data.frame(i = NULL, error = NULL),
                   CAR = data.frame(i = NULL, error = NULL))
    
    # run the models for each simulation dataset
    system.time({
      if(parallelize){
        library(doRNG)
        library(doParallel)
        fitted_list <- foreach(i=1:params$R) %dorng% 
          one_run(data_list,i, model_list)
      }else{
        fitted_list <- lapply(1:params$R, 
                              function(i) {
                                one_run(data_list,i, model_list)
                              })
        
      }
      
    })
    
    # store the results.
    full_results[[iter]] <- list(DGP = DGP,
                                 missingness = missingness,
                                 data = data_list,
                                 fitted_list = fitted_list)
    iter <- iter + 1
    
  }
}

#### Compute the results for a single simulation ####
# this shows an example of how to use the calculate_metrics function

sim1_list <- lapply(full_results[[1]]$fitted_list, function(xx) xx$df_miss)

results <- calculate_metrics(sim1_list, methods = c("y_pred_WF_negbin", "y_pred_freqGLMepi_negbin", "y_CARstan"), family = 'negbin') 

# (remember if you didn't increase the number of simulations (R), the R_PI or MCMC parameters the results will be bad.)
results