library(MASS)
#library(CARBayesST)
library(Matrix)
library(dplyr)
library(lubridate)
library(ggplot2)
library(rstan)

source('code/functions.R')
rstan_options(auto_write = TRUE)

# Weinberger-Fulcher Quasi-poisson DGP Arguments
WF_QP_DGP_arguments <- list(district_sizes = c(4, 6, 10), # size of the districts
                      R = 10, # number of simulated data sets
                      seed = 1, # random seed
                      end_date = '2020-12-01', # last date of data created.
                      b0_mean = 6, # mean value b0 is sampled from.
                      b1_mean = -0.25, # mean value b1 is sampled from.
                      type = 'WF', # type of DGP: 'WF', 'freqGLM', or 'CAR'.
                      family = 'quasipoisson', # family of DGP
                      theta = 4) # dispersion parameter for QP DGP

# freqGLM DGP
freqGLM_DGP_arguments <- list(district_sizes = c(4, 6, 10), # size of the districts
                          R = 10, # number of simulated data sets
                          seed = 1, # random seed
                          end_date = '2020-12-01', # last date of data created.
                          b0_mean = 6, # mean value b0 is sampled from.
                          b1_mean = -0.25, # mean value b1 is sampled from.
                          type = 'freqGLM', # type of DGP: 'WF', 'freqGLM', or 'CAR'.
                          family = 'poisson', # family of DGP
                          rho = 0.3, # spatial correlation of freqGLM process in DGP.
                          alpha = 0.3) # temporal correlation of freqGLM process in DGP.

# CAR DGP
CAR_DGP_arguments <- list(district_sizes = c(4, 6, 10), # size of the districts
                          R = 10, # number of simulated data sets
                          seed = 1, # random seed
                          end_date = '2020-12-01', # last date of data created.
                          b0_mean = 6, # mean value b0 is sampled from.
                          b1_mean = -0.25, # mean value b1 is sampled from.
                          type = 'CAR', # type of DGP: 'WF', 'freqGLM', or 'CAR'.
                          family = 'poisson', # family of DGP
                          rho = 0.3, # spatial correlation of CAR process in DGP.
                          alpha = 0.3, # temporal correlation of CAR process in DGP.
                          tau2 = 1) # variance of CAR process in DGP

# Simulate the data
lst <- do.call(simulate_data, freqGLM_DGP_arguments)

# initialize the error catcher for each run
errors <- list(freqGLM = data.frame(i = NULL, error = NULL),
               WF = data.frame(i = NULL, error = NULL),
               CAR = data.frame(i = NULL, error = NULL))

}


# function to run all models for a specific dataset
one_run <- function(lst, i, models = c('freq', 'WF', 'CAR'), WF_params = list(R_PI = 200), freqGLM_params = list(R_PI = 200), MCMC_params = list(burnin.stan = 1000, n.sample.stan = 2000, burnin.CARBayesST = 5000, n.sample.CARBayesST = 10000)){
  
  set.seed(i)
  
  t0 <- Sys.time()
  print(sprintf('i = %i',i))
  df = lst$df_list[[i]]
  
  set.seed(i)
   # add in missingness
  if(params[['missingness']] == 'mcar'){
    df_miss = MCAR_sim(df, p = params[['p']], by_facility = T)
  }else if(params[['missingness']] == 'mar'){
    df_miss = MAR_spatiotemporal_sim(df, p = params[['p']], rho = params[['rho_MAR']], alpha = params[['alpha_MAR']], tau2 = params[['tau2_MAR']])
  }else{
    df_miss <- MNAR_sim(df, p = params[['p']], direction = 'upper', gamma = params[['gamma']], by_facility = T)
  }
  
  # initializing the return list
  return_list <- list(df_miss = df_miss, district_df = lst$district_list[[i]], errors = errors, timing = list())
  rm(df_miss)
  
  set.seed(i)
  # run the freqGLM_epi complete case analysis
  if('freq' %in% models){
    t1 = Sys.time()
    print('running freqGLM_epi')
    return_list <- tryCatch({
      freqGLMepi_list = freqGLMepi_CCA(return_list[['df_miss']], R_PI = freqGLM_params[['R_PI']], verbose = F)
      return_list[['df_miss']] <- freqGLMepi_list$df
      return_list[['district_df']] <- merge(return_list[['district_df']], freqGLMepi_list$district_df, by = c('district','date'))
      return_list[['freqGLM_params']] <- freqGLMepi_list$param_results
      return_list
    }, error = function(e){
      return_list[['errors']][['freqGLM']] <- rbind(return_list[['errors']][['freqGLM']], data.frame(i, error = e[[1]]))
      return_list
    })
    return_list[['timing']][['freq']] <- as.numeric(difftime(Sys.time(), t1, units = 'm'))
  }
  
  set.seed(i)
  # run the WF complete case analysis model
  if('WF' %in% models){
    t1 = Sys.time()
    print('running WF CCA')
    return_list <- tryCatch({
      res <- WF_CCA(return_list[['df_miss']], col = "y", family = 'poisson', R_PI = WF_params[['R_PI']])
      return_list[['df_miss']] <- res$df
      return_list[['district_df']] <- merge(return_list[['district_df']], res$district_df, by = c('district','date'))
        res$district_df
      return_list[['WF_betas']] <- res$betas
      return_list
    }, error = function(e){
      return_list[['errors']][['WF']] <- rbind(return_list[['errors']][['WF']], data.frame(i = i, error = e[[1]]))
      return_list
    })
    return_list[['timing']][['WF']] <- as.numeric(difftime(Sys.time(), t1, units = 'm'))
  }

  set.seed(i)
  if('CAR' %in% models){
    t1 = Sys.time()
    print('running CARBayes with stan')
    return_list <- tryCatch({
      res <- CARBayes_wrapper(return_list[['df_miss']], burnin = MCMC_params[['burnin.stan']], n.sample = MCMC_params[['n.sample.stan']], prediction_sample = T, model = 'facility_fixed', predict_start_date = '2016-01-01', MCMC_sampler = 'stan')
      
      # rename the columns
      colnames(res$df) <- gsub('y_pred_CAR', 'y_CARstan', colnames(res$df))
      colnames(res$district_df) <- gsub('y_pred_CAR', 'y_CARstan', colnames(res$district_df))
      # update the results list
      return_list[['df_miss']] <- res$df
      return_list[['district_df']] <- merge(return_list[['district_df']], res$district_df, by = c('district','date'))
      return_list[['CARstan_summary']] <- res$CARstan_summary
      return_list
    }, error = function(e){
      print(e)
      return_list[['errors']][['CAR']] <- rbind(return_list[['errors']][['CAR']], data.frame(i, error = e[[1]]))
      return_list
    })
    return_list[['timing']][['CAR']] <- as.numeric(difftime(Sys.time(), t1, units = 'm'))
  }

  print(sprintf('i = %i; time = %f minutes', i, difftime(Sys.time(), t0, units = 'm')))
  
  # check I got the correct result names
  # outcome_name_checker(return_list, models = models)
  
  return(return_list)
}

# run the models for each simulation dataset
system.time({
  #imputed_list <- foreach(i=seq) %dorng% one_run(lst, i, models = c('freq', 'WF', 'CARstan'))
  imputed_list <- foreach(i=1:R_new) %dorng% one_run(lst, i, models = c('freq', 'WF', 'CARstan'))
})

# res <- one_run(lst, 1, freqGLM_params = list(R_PI = 5), MCMC_params = list(burnin.stan = 20, n.sample.stan = 50, burnin.CARBayesST = 100, n.sample.CARBayesST = 200))

true_betas <- lst$betas

save(imputed_list, seq, params, arguments, true_betas, file = results_file)

