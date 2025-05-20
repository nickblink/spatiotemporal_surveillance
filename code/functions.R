
#### Model fitting wrapper ####
# function to run all models for a specific dataset
one_run <- function(lst, i, model_list = NULL){
  
  t0 <- Sys.time()
  print(sprintf('i = %i',i))
  df = lst$df_list[[i]]
  
  set.seed(i)
  
  # add in missingness to the data
  if(params[['missingness']] == 'mcar'){
    df_miss = MCAR_sim(df, p = params[['p']], by_facility = T)
  }else if(params[['missingness']] == 'mar'){
    df_miss = MAR_spatiotemporal_sim(df, p = params[['p']], rho = params[['rho_MAR']], alpha = params[['alpha_MAR']], tau2 = params[['tau2_MAR']])
  }else{
    df_miss <- MNAR_sim(df, p = params[['p']], direction = 'upper', gamma = params[['gamma']], by_facility = T)
  }
  
  # initialize errors 
  errors <- list()
  for(j in 1:length(model_list)){
    errors[[model_list[[j]]$model]] <- data.frame(i = NULL, error = NULL)
  }
  
  # initializing the return list
  return_list <- list(df_miss = df_miss, district_df = lst$district_list[[i]], errors = errors, timing = list())
  
  # cycle through models
  for(sub_lst in model_list){
    t1 = Sys.time()
    model = sub_lst$model
    print(sprintf('------(%s)%s-------', i, model))
    params = sub_lst$params
    if(model == 'WF'){
      model_function <- function(df){
        tmp <- WF_CCA(df, col = "y", family = 'negbin', R_PI = model_list[['WF']]$params$R_PI)
        return(tmp)
      }
    }else if(model == 'freqGLM'){
      model_function <- function(df){
        tmp <- freqGLMepi_CCA(df, R_PI = model_list[['freqGLM']]$params$R_PI, verbose = F, family = 'negbin') 
        return(tmp)
      }
    }else if(model == 'CAR'){
      model_function <- function(df){
        tmp <- CAR_wrapper(df, burnin = model_list[['CAR']]$params$burnin, n.sample = model_list[['CAR']]$params$n.sample, prediction_sample = F, predict_start_date = '2016-01-01', MCMC_sampler = 'stan', family = 'negbin', use_fitted_phi = T)
        return(tmp)
      }
    }
    
    # run the model
    res <- tryCatch({
      model_function(return_list[['df_miss']])
    }, error = function(e){
      e[[1]]
    })
    
    # store the results
    if(class(res) == 'list'){
      return_list[['df_miss']] <- res$df
      return_list[['district_df']] <- merge(return_list[['district_df']], res$district_df, by = c('district','date'))
      if(model == 'WF'){
        return_list[['WF_betas']] <- res$betas
        return_list[['WF_NB_overdisp']] <- res$overdisp
      }else if(model == 'freqGLM'){
        return_list[['freqGLM_params']] <- res$param_results
      }else if(model == 'CAR'){
        return_list[['CAR_summary']] <- res$CARstan_summary
      }
    }else{
      return_list[['errors']][[model]] <- rbind(return_list[['errors']][[model]], data.frame(i = i, error = res[[1]]))
    }
    
    # keep track of timing
    return_list[['timing']][[model]] <- as.numeric(difftime(Sys.time(), t1, units = 'm'))
  }
  
  print(sprintf('i = %i; time = %f minutes', i, difftime(Sys.time(), t0, units = 'm')))
  
  return(return_list)
}
#### Helper Functions ####
### function to add periodic terms to data
# df: Input data frame.
# period: period length in months.
add_periodic_cov <- function(df, period = 12){
  df = df %>%
    dplyr::mutate(year = lubridate::year(date) - min(lubridate::year(date)) + 1,
                  month = month(date),
                  cos1 = cos(2*1*pi*month/period),
                  sin1 = sin(2*1*pi*month/period),
                  cos2 = cos(2*2*pi*month/period),
                  sin2 = sin(2*2*pi*month/period),
                  cos3 = cos(2*3*pi*month/period),
                  sin3 = sin(2*3*pi*month/period))
  return(df)
}

### compute the sum of the values at all the neighbors to a given facility
# df: Input data frame.
# target_col: name of outcome column.
# lag: Temporal lag of neighbors in sum.
# scale_by_num_neighbors: If true, then the weights in each column are scaled to sum to one. 
# W: The adjacency matrix. If Null, W is calculate within this function.
add_neighbors <- function(df, target_col = 'y', lag = 1, scale_by_num_neighbors = F, W = NULL){
  if(lag == 0){
    print('WARNING: not doing a lag of 1 on the neighbors, so not using the same model as in the papers')
  }
  
  # remove the neighbor column
  if('y.neighbors' %in% colnames(df)){
    df$y.neighbors = NULL
  }
  
  # get the adjacency matrix
  if(is.null(W)){
    W <- make_district_adjacency(df, scale_by_num_neighbors)
  }
  
  # get the counts for each facility 
  y.counts <- df %>% 
    dplyr::select(date, facility, UQ(target_col)) %>%
    arrange(facility) %>%
    tidyr::spread(facility,get(target_col)) %>% 
    arrange(date)
  
  # check that the column names match
  if(!identical(colnames(W), colnames(y.counts)[-1])){
    # if the colnames do match but are in the wrong order, reorder them
    if(length(setdiff(colnames(W), colnames(y.counts)[-1])) == 0 &
       length(setdiff(colnames(y.counts)[-1], colnames(W))) == 0){
      
      matid = match(colnames(y.counts)[-1], colnames(W))
      y.counts = y.counts[,c(1, matid + 1)]
      
      # for error checking
      if(!identical(colnames(W), colnames(y.counts)[-1])){browser()}
    }else{
      stop('Adjacency and y matrix column names dont match')
    }
    ind = which(colnames(W) != colnames(y.counts)[-1])
  }
  
  # shift y.counts by lag
  if(lag > 0){
    tmp <- y.counts
    tmp[1:lag,-1] <- NA
    tmp[(lag+1):nrow(tmp),-1] <- y.counts[1:(nrow(y.counts) - lag),-1]
    y.counts <- tmp
  }
  
  # merge back into original data frame
  tmp = cbind(y.counts[,'date',drop = F], as.data.frame(as.matrix(y.counts[,-1])%*%W)) %>% 
    tidyr::gather(facility, y.neighbors, -date)
  if(is.factor(df$facility)){
    tmp$facility = factor(tmp$facility, levels = levels(df$facility))  
  }
  df = merge(df, tmp, by = c('date','facility'))
  
  return(df)
}

### add autoregressive terms to the data.
# df: Input data frame.
# target_col: name of outcome column.
# num_terms: 
add_autoregressive <- function(df, target_col = 'y'){
  
  if(!(target_col %in% colnames(df))){
    stop('need a target column to autoregress')
  }
  
  tmp <- lapply(unique(df$facility), function(xx) {
    tt <- df %>% filter(facility == xx) %>% arrange(date)
    tt[,'y.AR1'] = c(NA, tt[1:(nrow(tt) - 1), target_col, drop = T])
    return(tt)
  })
  
  df <- do.call('rbind',tmp)
  return(df)
}


### Randomly sample from a quasipoisson. Do this by sampling from a negative binomial with dispersion parameters chosen to emulate the quasipoisson.
# n: number of samples.
# mu: mean of the points.
# theta: dispersion parameter.
rqpois <- function(n, mu, theta) {
  rnbinom(n = n, mu = mu, size = mu/(theta-1))
}


### make the adjacency matrix according to all facilities in a district being neighbors.
# df: data frame containing a list of facilities and districts.
# scale_by_num_neighbors: If true, then the weights in each column are scaled to sum to one. 
make_district_adjacency <- function(df, scale_by_num_neighbors = F, include_no_neighbors = F){
  
  # get unique districts and facilities
  D2 = df %>% dplyr::select(district, facility) %>% 
    distinct() %>%
    arrange(facility)
  facs <- D2$facility
  
  # make the adjacency matrix
  W = full_join(D2, D2, by = 'district', relationship = 'many-to-many') %>%
    filter(facility.x != facility.y) %>%
    dplyr::select(-district) %>%
    igraph::graph_from_data_frame() %>%
    igraph::as_adjacency_matrix() %>%
    as.matrix()
  
  # add in rows and columns for places with no neighbors
  if(include_no_neighbors){
    if(nrow(W) != length(facs)){
      W_star <- matrix(0, nrow = length(facs), ncol = length(facs))
      rownames(W_star) <- colnames(W_star) <- facs
      matid <- match(rownames(W), rownames(W_star))
      W_star[matid, matid] <- W
      W <- W_star
    }
  }
  
  # scale the neighbor sum by the number of neighbors
  if(scale_by_num_neighbors){
    W = apply(W, 2, function(xx) xx/sum(xx))
  }
  
  return(W)
}



### Make the "W2" matrix, as I am calling it, which is diag(W1) - W. I.e. this is the negative adjacency matrix with the total number of neighbors for each facility on the diagonal.
# df: Input data frame with columns "facility" and "district"
make_district_W2_matrix <- function(df){
  # create the list of matching facilities
  D2 = df %>% dplyr::select(district, facility) %>% distinct()
  pairs = full_join(D2, D2, by = 'district', relationship = 'many-to-many') %>%
    filter(facility.x != facility.y) %>%
    dplyr::select(-district)
  
  # get unique facilities
  facilities = unique(df$facility) %>% sort()
  
  # initialize the W2 matrix (note that what I am calling W2 here is what the original paper calls diag(W1) - W)
  W2 = matrix(0, nrow = length(facilities), ncol = length(facilities))
  colnames(W2) = facilities
  rownames(W2) = facilities
  
  for(i in 1:length(facilities)){
    f1 = facilities[i]
    
    # get the pairs with this facility
    tmp = pairs %>% filter(facility.x == f1)
    
    if(nrow(tmp) == 0){
      print('havent done single-district facilities yet')
      browser()
    }else{
      # put -1 where there are pairs of facilities
      matid = match(tmp$facility.y, facilities)
      W2[i,matid] = -1
    }
    
    # get the number of neighbors this facility has
    W2[i,i] = pairs %>% filter(facility.x == f1) %>% nrow()
    # print(sprintf('%s: NUM neighbors = %s', f1, W2[i,i]))
  }
  
  return(W2)
}


# get the district and facility list from a data frame
get_district_facilities <- function(df){
  tt = unique(df[,c('district','facility')])
  res_lst <- NULL
  for(i in 1:nrow(tt)){
    res_lst[[as.character(tt[i,1])]] <- c(res_lst[[as.character(tt[i,1])]], as.character(tt[i,2]))
  }
  return(res_lst)
}

#### Simulation Functions ####

### Function to initialize the data frame for data simulation.
# district_sizes: number of facilities in the districts.
# start_date: start date for data.
# end_date: end date for data.
initialize_df <- function(district_sizes, start_date = '2016-01-01', end_date = '2019-12-01', ...){
  facilities = unlist(lapply(1:length(district_sizes), function(xx) {paste0(toupper(letters[xx]), sprintf('%03d',1:district_sizes[xx]))}))
  
  dates = seq(as.Date(start_date), as.Date(end_date), by = 'month')
  
  df = expand.grid(facilities, dates, stringsAsFactors = T)
  colnames(df) = c('facility','date')
  
  df$district = substr(df$facility, 1, 1) 
  
  df = df %>%
    dplyr::select(date, facility, district) %>%
    arrange(facility, date)
  return(df)
}

### Function to sample betas for simulated data.
# facilities: vector of facility names.
# b0_mean: mean value b0 is sampled from.
# b0_sd: standard deviation b0 is sampled from.
# b1_mean: mean value b1 is sampled from.
# b1_sd: standard deviation b1 is sampled from.
sample_betas = function(facilities, b0_mean = 4.3, b0_sd = 1, b1_mean = -0.25, b1_sd = 0.25, ...){
  betas = matrix(0, nrow = length(facilities), ncol = 8)
  
  if(length(b0_mean) == 1){
    betas[,1] = rnorm(mean = b0_mean, sd = b0_sd, n = nrow(betas))
  }else{
    u = sample(b0_mean, nrow(betas), replace = T)
    betas[, 1] = rnorm(mean = u, sd = b0_sd, n = nrow(betas))
  }
  
  betas[,2] = rnorm(b1_mean, b1_sd, n = nrow(betas))
  
  for(j in 3:8){
    betas[,j] = rnorm(0, 0.15, n = nrow(betas))
  }
  
  # name the rows and columns.
  rownames(betas) = facilities
  colnames(betas) = c('intercept', 'year', 'cos1', 'sin1', 'cos2', 'sin2', 'cos3', 'sin3')

  return(betas)
}

### Sample the negative binomial thetas for simulated data.
# facilities: vector of facility names.
# DGP_theta_shape and DGP_theta_rate: thetas are sampled from a Gamma(shape, rate) distribution.
sample_thetas <- function(facilities, DGP_theta_shape = 2.5, DGP_theta_rate = 1/3, ...){
  thetas <- rgamma(n = length(facilities), shape = DGP_theta_shape, rate = DGP_theta_rate)
  
  names(thetas) <- facilities
  
  return(thetas)
}


### Function to simulate the data for a variety of situations
# district_sizes: number of facilities in the districts.
# R: number of simulated data sets.
# seed: random seed for reproducibility.
# type: Model to simulate data. Can be 'WF', 'freqGLM', or 'CAR'.
# family: Family of simulated data. Can be 'poisson' or 'quasipoisson' or 'negative binomial'.
### Function to simulate the data for a variety of situations
simulate_data <- function(district_sizes, R = 1, empirical_betas = F, seed = 10, type = 'WF', family = 'negbin', ...){
  
  # set seed so the betas are always the same (the seed input is used later for simulating the data on top of these betas)
  set.seed(10)
  
  # set up data frame
  df = initialize_df(district_sizes, ...)
  
  # make periodic covariates
  df = add_periodic_cov(df)
  
  # get all facility and district names
  facilities = unique(df$facility)
  districts = unique(df$district)
  
  # get all dates
  dates = unique(df$date) %>% sort()
  
  # sample betas
  if(empirical_betas){
    tmp = sample_real_betas(facilities)
    betas = tmp$betas
    dispersion_vec = tmp$dispersion
  }else{
    betas = sample_betas(facilities, ...)  
  }
  
  # set up the data generating function based off the family.
  if(family == 'poisson'){
    DGP_function = function(n, mu){
      y <- rpois(n, mu)
    }
  }else if(family == 'quasipoisson'){
    DGP_function <- function(n, mu){
      y <- rqpois(n, mu, theta = list(...)$theta)
    }
  }else if(family == 'negbin'){
    dispersion_vec <- sample_thetas(facilities, ...)
    matid = match(df$facility, names(dispersion_vec))
    df$theta <- dispersion_vec[matid]
    # make the DGP below since it depends on a different dispersion parameter for each facility.
  }
  
  # initialize list of data frames
  df_lst = list()
  district_lst = list()
  
  # set seed for the data generation
  set.seed(seed)
  
  # simulate the data according to the DGP
  if(type == 'WF'){
    # make R sampled sets of data
    for(i in 1:R){
      
      # simulate values given the betas
      tmp_lst = lapply(facilities, function(xx){
        tmp = df %>% filter(facility == xx)
        
        # keep the 1 for intercepts
        X = tmp %>% 
          mutate(intercept = 1) %>%
          dplyr::select(intercept, year, cos1, sin1, cos2, sin2, cos3, sin3)
        
        # error checking
        if(!identical(colnames(betas), colnames(X))){
          print(colnamaes(betas))
          print(colnames(X))
          stop(sprintf('colnames of betas not equal to X: %s, %s',paste(colnames(betas), collapse = ';'), paste(colnames(X), collapse = ';') ))
        }
        
        # make the 8x1 beta vector for this facility
        beta_f = t(betas[xx,,drop = F])
        
        # get mean prediction from linear model
        mu = as.matrix(X)%*%beta_f
        tmp$y_exp = exp(mu)[,1]
        
        # get Poisson or quasipoisson variance
        if(family == 'poisson'){
          tmp$y_var <- tmp$y_exp
        }else if(family == 'quasipoisson'){
          tmp$y_var <- list(...)$theta*tmp$y_exp
        }else if(family == 'negbin'){
          if(!exists('dispersion_vec')){
            dispersion_parameter = list(...)$dispersion
          }else{
            dispersion_parameter = dispersion_vec[[xx]]
          }
          tmp$y_var <- tmp$y_exp + tmp$y_exp^2/dispersion_parameter
          DGP_function <- function(n, mu){
            y <- rnbinom(n = n, mu = mu, size = dispersion_parameter)
          }
          
        }else{
          stop('input a proper family')
        }
        
        # simluate random values
        tmp$y = DGP_function(length(mu), exp(mu))
        
        return(tmp)
      })
      
      # combine values into one data frame
      df = do.call('rbind', tmp_lst)
      df_lst[[i]] = df
      
      # group the results by district
      district <- df %>%
        group_by(district, date) %>%
        summarize(y_exp = sum(y_exp),
                  y_var = sum(y_var),
                  y = sum(y),
                  y_true = sum(y), .groups = 'drop')
      district_lst[[i]] = district
      
    }
    
  }else if(type == 'CAR'){
    
    rho = list(...)$rho
    alpha = list(...)$alpha
    tau2 = list(...)$tau2
    
    # add in the mean effects because these are the same for all simulations
    df = do.call('rbind', lapply(facilities, function(xx){
      tmp = df %>% filter(facility == xx)
      
      # keep the 1 for intercepts
      X = tmp %>% 
        mutate(intercept = 1) %>%
        dplyr::select(intercept, year, cos1, sin1, cos2, sin2, cos3, sin3)
      
      # error checking
      if(!identical(colnames(betas), colnames(X))){
        browser()
      }
      
      # make the 8x1 beta vector for this facility
      beta_f = t(betas[xx,,drop = F])
      
      # get mean prediction from linear model
      tmp$mu = (as.matrix(X)%*%beta_f)[,1]
      
      return(tmp)
    }))
    
    # make the spatio-temporal precision matrix
    Q = make_precision_mat(df, rho = rho)
    
    tryCatch({
      V = tau2*solve(Q)
    }, error = function(e){
      print(e)
      print('havent dealt with non-invertible precision matrices yet')
      browser()
    })
    
    # checking ordering of facilities matches
    if(!identical(colnames(V), as.character(facilities))){
      browser()
      stop('names of covariances and facilities dont match')
    }
    
    # add in the marginal variance to original df
    dV = diag(V)
    matid = match(df$facility, names(dV))
    df$sigma2_marginal = dV[matid]
    
    # make R sampled sets of data
    df_lst = lapply(1:R, function(r){
      ### get the spatio-temporal random effects
      # initialize phi
      phi = matrix(0, nrow = length(dates), ncol = length(facilities))
      colnames(phi) = facilities
      
      # first time step
      phi[1,] = MASS::mvrnorm(n = 1, mu = rep(0, nrow(V)), Sigma = V)
      
      # cycle through other time steps, using auto-correlated priors
      for(i in 2:length(dates)){
        phi[i,] = MASS::mvrnorm(n = 1, mu = alpha*phi[i-1,], Sigma = V)
      }
      
      # convert to matching format
      phi = as.data.frame(phi)
      phi$date = dates
      phi_df = tidyr::gather(phi, facility, phi, setdiff(colnames(phi), c('facility','date'))) %>% 
        arrange(facility, date)
      
      # add in the previous phi values
      phi_df$phi_previous <- c(0,phi_df$phi[1:(nrow(phi_df) - 1)])
      phi_df$phi_previous[phi_df$date == min(phi_df$date)] <- 0
      
      # merge the phi values into the original data frame
      df = merge(df, phi_df, by = c('date','facility'))
      
      # calculate the expected value of the Poisson log-normal
      df$mu_marginal = df$mu + alpha*df$phi_previous
      df$y_exp = exp(df$mu_marginal + df$sigma2_marginal/2)
      
      # variance from Poisson or quasi-poisson log-normal
      if(family == 'poisson'){
        df$y_var = df$y_exp + (exp(df$sigma2_marginal) - 1)*exp(2*df$mu_marginal + df$sigma2_marginal)
      }else if(family == 'quasipoisson'){
        df$y_var = list(...)$theta*df$y_exp + (exp(df$sigma2_marginal) - 1)*exp(2*df$mu_marginal + df$sigma2_marginal)
      }else{
        df$y_var = df$y_exp + exp(2*df$mu_marginal + 2*df$sigma2_marginal)/df$theta + (exp(df$sigma2_marginal) - 1)*exp(2*df$mu_marginal + df$sigma2_marginal)
      }
      
      # simulate the observed values
      if(family %in% c('poisson','quasipoisson')){
        df$y <- DGP_function(nrow(df), exp(df$mu + df$phi))
      }else if(family == 'negbin'){
        df$y <- rnbinom(n = nrow(df), mu = exp(df$mu + df$phi), size = df$theta)
      }else{
        stop('please input family = poisson, quasipoisson, or negbin')
      }
      
      return(df)
    })
    
    # cycle through each created data frame
    district_lst = lapply(df_lst, function(df){
      # cycle through each district
      district <- do.call('rbind', lapply(districts, function(d){
        # filter data to this district
        df2 <- df %>% filter(district == d)
        facs = unique(df2$facility)
        
        # get covariance by this district
        V_d = V[facs, facs]
        V_exp = exp(V_d) - 1
        ind_cov = upper.tri(V_exp) + lower.tri(V_exp)
        
        # if(family == 'quasipoisson'){
        #   warning('cant do district level quasipoisson variance for CAR yet. Havent coded it.')
        # }
        
        district_df = do.call('rbind', lapply(1:length(dates), function(i_date){
          df3 <- df2 %>% filter(date == dates[i_date])
          covs = (df3$y_exp%*%t(df3$y_exp))*V_exp
          cov = sum(covs*ind_cov)
          
          tmp_df = data.frame(district = d,
                              date = dates[i_date],
                              y_exp = sum(df3$y_exp),
                              y_var_ind = sum(df3$y_var),
                              y_cov = cov,
                              y_var = sum(df3$y_var) + cov,
                              y = sum(df3$y),
                              y_true = sum(df3$y))
          return(tmp_df)
        }
        ))
        return(district_df)
      }))
      return(district)
    })
    
  }else if(type == 'freqGLM'){
    rho = list(...)$rho
    alpha = list(...)$alpha
    
    # get the adjacency matrix
    W <- make_district_adjacency(df, scale_by_num_neighbors = T)
    
    # add in the mean effects because these are the same for all simulations
    df = do.call('rbind', lapply(facilities, function(xx){
      tmp = df %>% filter(facility == xx)
      
      # keep the 1 for intercepts
      X = tmp %>% 
        mutate(intercept = 1) %>%
        dplyr::select(intercept, year, cos1, sin1, cos2, sin2, cos3, sin3)
      
      # error checking
      if(!identical(colnames(betas), colnames(X))){
        browser()
      }
      
      # make the 8x1 beta vector for this facility
      beta_f = t(betas[xx,,drop = F])
      
      # get mean prediction from linear model
      tmp$mu = (as.matrix(X)%*%beta_f)[,1]
      
      return(tmp)
    }))
    
    df$y_exp = df$y_var = df$y = NA
    
    # make R sampled sets of data
    df_lst = lapply(1:R, function(r){
      
      ### Get the first time point values
      ind = which(df$date == min(dates))
      
      # the adjustment accounts for the fact that there aren't additive auto-regressive and spatial terms at the first time point.
      # this adjustment comes from the sum of a geometric series (since this is roughly the effect that the spatial and autoregressive terms approach as we increase the time series)
      adjustment = 1 + rho/(1-rho) + alpha/(1-alpha)
      if(family == 'poisson'){
        df$y_exp[ind] = df$y_var[ind] = adjustment*exp(df$mu[ind])
        
        df$y[ind] = DGP_function(length(ind), df$y_exp[ind])
      }else if(family == 'negbin'){
        df$y_exp[ind] = adjustment*exp(df$mu[ind])
        df$y_var[ind] = df$y_exp[ind] + df$y_exp[ind]^2/df$theta[ind]
        
        DGP_function <- function(n, mu, theta){
          y <- rnbinom(n = n, mu = mu, size = theta)
        }
        
        df$y[ind] = DGP_function(length(ind), df$y_exp[ind], theta = df$theta[ind])
      }else{
        stop('only coded poisson or negative binomial (negbin) families for freqGLM.')
      }
      
      # update the neighbors and auto-regressive terms
      df <- add_autoregressive(df, 'y') %>%
        add_neighbors(., 'y', scale_by_num_neighbors = T, W = W)
      
      ### remaining time points
      for(d in dates[-1]){
        # get the subset of dates
        ind = which(df$date == d)
        
        # get the mean estimates for these dates
        df$y_exp[ind] <- exp(df$mu[ind]) + 
          alpha*df$y.AR1[ind] + rho*df$y.neighbors[ind]
        
        # get Poisson or quasipoisson variance
        if(family == 'poisson'){
          df$y_var[ind] <-df$y_exp[ind]
          
          # predict!
          df$y[ind] = DGP_function(length(ind), df$y_exp[ind])
        }else if(family == 'quasipoisson'){
          df$y_var[ind] <- list(...)$theta*df$y_exp[ind]
          
          # predict!
          df$y[ind] = DGP_function(length(ind), df$y_exp[ind])
        }else if(family == 'negbin'){
          df$y_var[ind] = df$y_exp[ind] + df$y_exp[ind]^2/df$theta[ind]
          
          # predict
          df$y[ind] = DGP_function(length(ind), df$y_exp[ind], theta = df$theta[ind])
        }else{
          stop('input a proper family')
        }
        
        # update the neighbors and auto-regressive terms
        df <- add_autoregressive(df, 'y') %>%
          add_neighbors(., 'y', scale_by_num_neighbors = T, W = W)
      }
      
      df
    })
    
    # group the results by district
    district_lst <- lapply(df_lst, function(df){
      district <- df %>%
        group_by(district, date) %>%
        summarize(y_exp = sum(y_exp),
                  y_var = sum(y_var),
                  y = sum(y),
                  y_true = sum(y), .groups = 'drop')
      district
    })
  }else{
    stop('please input a proper type')
  }
  
  # make list of values to return
  res_lst = list(df_list = df_lst, district_list = district_lst, betas = betas)
  
  if(family == 'negbin'){
    res_lst$thetas = dispersion_vec
  }
  return(res_lst)
}


#### Function to simulate data under the freqGLM_epi framework
##
## district_sizes = number of facilities in the districts
## R = number of simulated data sets
## lambda = autoregressive term
## phi = neighbor term
simulate_data_freqGLM_epi <- function(district_sizes, R = 1, lambda = -2, phi = -2, num_iters = 10, scale_by_num_neighbors = F, seed = 10, start_date = '2016-01-01', end_date = '2019-12-01', b0_mean = 7, b1_mean = -0.2, b1_sd = 0.2){
  
  warning('counting all districts as neighbors')
  print(sprintf('lambda: exp(%0.2f) = %0.2f; phi: exp(%0.2f) = %0.2f', lambda, exp(lambda), phi, exp(phi)))
  
  # set up data frame
  df = initialize_df(district_sizes, start_date = start_date, end_date = end_date)
  
  # make periodic covariates
  df = add_periodic_cov(df)
  
  # get all facility names
  facilities = unique(df$facility) %>% sort()
  
  # get all dates
  dates = unique(df$date) %>% sort()
  #n = length(dates)
  
  # set random seed and sample betas
  set.seed(seed)
  betas = sample_betas(facilities, b0_mean = b0_mean, b1_mean = b1_mean, b1_sd = b1_sd)
  
  ### get the seasonal effects (since these don't change across the simulation samples)
  df = do.call('rbind', lapply(facilities, function(xx){
    tmp = df %>% filter(facility == xx)
    
    # keep the 1 for intercepts
    X = tmp %>% 
      mutate(intercept = 1) %>%
      dplyr::select(intercept, year, cos1, sin1, cos2, sin2, cos3, sin3)
    
    # error checking
    if(!identical(colnames(betas), colnames(X))){
      browser()
    }
    
    # make the 8x1 beta vector for this facility
    beta_f = t(betas[xx,,drop = F])
    
    # get mean prediction from linear model
    tmp$mu_seasonal = (as.matrix(X)%*%beta_f)[,1]
    
    return(tmp)
  }))
  
  # initial sampling
  df$y_seasonal = rpois(n = nrow(df), lambda = exp(df$mu_seasonal))
  
  # initialize list of final data frame
  df_lst = list()
  
  for(r in 1:R){
    
    # to store the predicted values at successive iterations
    y_pred_list = list()
    
    tmp = df %>%
      mutate(y = y_seasonal)
    
    for(i in 1:num_iters){
      # add the neighbors and auto-regressive
      tmp = add_autoregressive(tmp, 'y') %>%
        add_neighbors(., 'y', scale_by_num_neighbors = scale_by_num_neighbors)
      
      # only resampling after the first time point
      ind = which(!is.na(tmp$y.AR1))
      
      # getting the mean according to the model and resampling
      mean_y = exp(lambda)*tmp$y.AR1 + exp(phi)*tmp$y.neighbors + exp(tmp$mu_seasonal)
      tmp$y[ind] = rpois(n = length(ind), lambda = mean_y[ind])
      
      # storing the results
      y_pred_list[[i]] = tmp$y
      
    }
    
    if(mean(y_pred_list[[num_iters]]) > 5*mean(y_pred_list[[1]])){
      print('there might be some divergence of the estimated values')
      browser()
    }
    
    df_lst[[r]] = tmp
  }
  
  params_true = as.data.frame(t(betas))
  rownames(params_true) = paste0('B', rownames(params_true))
  params_true = rbind(t(data.frame(By.AR1 = rep(lambda, ncol(params_true)), By.neighbors = rep(phi, ncol(params_true)), row.names = colnames(params_true))), params_true)
  
  
  # return it!
  res_lst = list(df_list = df_lst, betas = betas, lambda = lambda, phi = phi, params = params_true)
  return(res_lst)
}

MCAR_sim <- function(df, p, by_facility = F, max_missing_date = '2019-12-01'){
  # save the true y value
  df$y_true = df$y
  
  # split the data into a test/hold-out set and the training set
  df_test <- df %>%
    filter(date > max_missing_date)
  df <- df %>%
    filter(date <= max_missing_date)
  
  # delete y values to add in missingness
  if(by_facility){
    num_impute = round(p*nrow(df)/length(unique(df$facility)))
    df = do.call('rbind', lapply(unique(df$facility), function(xx){
      tmp = df %>% filter(facility == xx)
      tmp$y[sample(nrow(tmp), num_impute)] <- NA
      return(tmp)
    }))
  }else{
    num_impute = round(p*nrow(df))
    df$y[sample(nrow(df), num_impute)] <- NA
  }
  
  # combine the hold-out/non-missing data with the training data with missingness
  df <- rbind(df, df_test)
  
  return(df)
}

MNAR_sim <- function(df, p, direction = NULL, gamma = 1.5, by_facility = T, max_missing_date = '2019-12-01'){
  df$y_true = df$y
  
  # split the data into a test/hold-out set and the training set
  df_test <- df %>%
    filter(date > max_missing_date)
  df <- df %>%
    filter(date <= max_missing_date)
  
  if(by_facility){
    # get the number of points to impute
    num_impute = round(p*nrow(df)/length(unique(df$facility)))
    
    # if removing low and high points
    if(is.null(direction)){
      df = do.call('rbind', lapply(unique(df$facility), function(xx){
        tmp = df %>% filter(facility == xx)
        m = median(tmp$y)
        q = (abs(tmp$y - m))^gamma
        tmp$y[sample(nrow(tmp), num_impute, prob = q)] <- NA
        return(tmp)
      }))
      # if removing high points
    }else if(direction == 'upper'){
      df = do.call('rbind', lapply(unique(df$facility), function(xx){
        tmp = df %>% filter(facility == xx)
        m = min(tmp$y)
        q = (abs(tmp$y - m))^gamma
        tmp$y[sample(nrow(tmp), num_impute, prob = q)] <- NA
        return(tmp)
      }))
      # if removing low points
    }else if(direction == 'lower'){
      df = do.call('rbind', lapply(unique(df$facility), function(xx){
        tmp = df %>% filter(facility == xx)
        m = max(tmp$y)
        q = (abs(tmp$y - m))^gamma
        tmp$y[sample(nrow(tmp), num_impute, prob = q)] <- NA
        return(tmp)
      }))
    }
    
  }else{
    print('havent coded this part - is it necessary?')
    browser()
    #num_impute = round(p*nrow(df))
    #df$y[sample(nrow(df), num_impute)] <- NA
  }
  
  df <- rbind(df[,colnames(df_test)], df_test)
  
  return(df)
}

MAR_spatiotemporal_sim <- function(df, p, rho = 0.3, alpha = 0.3, tau2 = 1, by_facility = T, max_missing_date = '2019-12-01'){
  # make phi for df
  # for all, or for each facility,
  # sample according to expit(phi)
  df$y_true = df$y
  
  # split the data into a test/hold-out set and the training set
  df_test <- df %>%
    filter(date > max_missing_date)
  df <- df %>%
    filter(date <= max_missing_date)
  
  # get all facility names
  facilities = as.character(unique(df$facility))
  
  # get all dates
  dates = unique(df$date) %>% sort()
  
  # make the spatio-temporal precision matrix
  Q = make_precision_mat(df, rho = rho)
  
  tryCatch({
    V = tau2*solve(Q)
  }, error = function(e){
    print(e)
    print('havent dealt with non-invertible precision matrices yet')
    browser()
  })
  
  # checking ordering of facilities matches
  if(!identical(colnames(V), facilities)){
    browser()
    stop('names of covariances and facilities dont match')
  }
  
  phi = matrix(0, nrow = length(dates), ncol = length(facilities))
  colnames(phi) = facilities
  
  # first time step
  phi[1,] = MASS::mvrnorm(n = 1, mu = rep(0, nrow(V)), Sigma = V)
  
  # cycle through other time steps, using auto-correlated priors
  for(i in 2:length(dates)){
    phi[i,] = MASS::mvrnorm(n = 1, mu = alpha*phi[i-1,], Sigma = V)
  }
  
  # convert to matching format
  phi = as.data.frame(phi)
  phi$date = dates
  phi_df = tidyr::gather(phi, facility, phiM, setdiff(colnames(phi), c('facility','date')))
  
  # convert to expit for probability
  phi_df$prob.sample = exp(phi_df$phiM)/(1 + exp(phi_df$phiM))
  
  # merge the phi values into the original data frame
  df = merge(df, phi_df, by = c('date','facility'))
  
  if(by_facility){
    num_impute = round(p*nrow(df)/length(unique(df$facility)))
    df = do.call('rbind', lapply(unique(df$facility), function(xx){
      tmp = df %>% filter(facility == xx)
      tmp$y[sample(nrow(tmp), num_impute, prob = tmp$prob.sample)] <- NA
      return(tmp)
    }))
  }else{
    num_impute = round(p*nrow(df))
    df$y[sample(nrow(df), num_impute, prob = df$prob.sample)] <- NA
  }
  
  # combine the hold-out/non-missing data with the training data with missingness
  df <- rbind(df[,colnames(df_test)], df_test)
  
  return(df)
}

#### WF functions ####
### Weingberger-Fulcher fitting method
# df: data frame to fit model on. 
# col: target column (usually "y").
# group: hierarchical level to run analysis on. 
# family: distribution.
# period: period for adding periodic covariates. 
# R_PI: Number of bootstrap iterations for prediction interval.
# bias_correction_chen: Whether to apply a bias correction for estimating Poisson betas.
# quant_probs: The quant probabilities to generate predictive intervals. 
WF_fit <- function(df, col, group = 'facility', family = 'negbin', period = 12, R_PI = 500, bias_correction_chen = F, quant_probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)){
  
  print(sprintf('WF family: %s', family))
  
  # prep the data with the harmonic functions
  df <- add_periodic_cov(df, period = period)
  
  # pulling unique groups
  uni_group = df %>% pull(get(group)) %>% unique()
  
  # get districts and facilities
  dist_fac <- get_district_facilities(df)
  
  # setting up the formula
  formula_col = as.formula(sprintf("%s ~ year + cos1 + sin1 + cos2 + sin2 + cos3 + sin3", col))
  
  if(family == 'quasipoisson'){
    warning('havent coded PIs in quasipoisson yet')
    
    tmp <- lapply(uni_group, function(xx) {
      tt <- df %>% filter(get(group) == xx)
      mod_col <- glm(formula_col, data = tt, family=quasipoisson)
      tt[,paste0(col, '_pred_WF')] = predict(mod_col, tt, type = 'response')
      return(tt)
    })
    
    #df <- data.table::rbindlist(tmp)
    df <- do.call('rbind',tmp)
    
  }else if(family %in% c('NB','negative binomial','negbin')){
    tmp <- lapply(uni_group, function(xx) {
      tt <- df %>% filter(get(group) == xx)
      #mod_col <- MASS::glm.nb(formula_col, data = tt, control = glm.control(maxit=200,trace = 3))
      
      tt[,paste0(col, '_pred_WF_negbin')] <- tryCatch({
        mod_col <- MASS::glm.nb(formula_col, data = tt)
        predict(mod_col, tt, type = 'response')
      }, error = function(e){
        print(sprintf('glm.nb failed for %s', xx))
        rep(NA, nrow(tt))
      })
      
      # generate prediction intervals.
      if(R_PI > 0){
        # extract information from model fit
        beta_hat <- mod_col$coefficients
        beta_vcov <- vcov(mod_col)
        overdisp <- summary(mod_col)$theta
        
        # store the model results to return
        model_res = list()
        model_res[[as.character(xx)]][['beta_hat']] <- beta_hat
        model_res[[as.character(xx)]][['beta_vcov']] <- beta_vcov
        model_res[[as.character(xx)]][['theta']] <- overdisp
        
        # bootstrap 
        sapply(1:R_PI, function(r){
          beta_boot <- MASS::mvrnorm(1,beta_hat,beta_vcov)
          pred_boot <- (tt %>% 
                          mutate(intercept=1) %>%
                          #dplyr::select(intercept,year,ends_with("1"),ends_with("2"),ends_with("3")) %>%
                          dplyr::select(intercept,year,cos1,sin1,cos2,sin2,cos3,sin3) %>%
                          as.matrix())%*%as.matrix(beta_boot)
          pred_boot_exp <- exp(pred_boot) 
          pred_boot_exp[which(is.na(pred_boot_exp))] <- 1
          x = MASS::rnegbin(n = nrow(tt), mu = pred_boot_exp, theta = overdisp)
          
          return(x)
          
        }) -> sim.boot
        
        # get the quantiles and store them
        fitted_quants = t(apply(sim.boot, 1, function(xx){
          quantile(xx, probs = quant_probs)
        }))
        fitted_quants = as.data.frame(fitted_quants)
        colnames(fitted_quants) = paste0(paste0(col, '_pred_WF_negbin_'), quant_probs)
        
        # merge the quantiles back into data frame
        tt = cbind(tt, fitted_quants)
      }
      tmp_lst = list(tt, sim.boot, model_res)
      return(tmp_lst)
    })
    names(tmp) = uni_group
    
  }else if(family == 'poisson'){
    tmp <- lapply(uni_group, function(xx) {
      tt <- df %>% filter(get(group) == xx)
      mod_col <- glm(formula_col, data = tt, family=poisson)
      
      tt[,paste0(col, '_pred_WF')] = predict(mod_col, tt, type = 'response')
      
      if(R_PI > 0){
        # extract information from model fit
        beta_hat <- mod_col$coefficients
        beta_vcov <- vcov(mod_col)
        
        # if correcting the bias, compute and correct it
        if(bias_correction_chen){
          X = tt %>%
            mutate(intercept = 1) %>%
            dplyr::select(intercept, year, cos1, sin1, cos2, sin2, cos3, sin3) %>%
            as.matrix()
          
          # the not good lin alg way
          tmp = lapply(1:nrow(X), function(i){
            x = X[i,]
            res = tt$y_exp[i]*x%*%t(x)
            return(res)
          })
          
          H1 = -Reduce("+", tmp) / length(tmp)
          Q = solve(H1)
          
          H2 = lapply(1:length(beta_hat), function(p){
            tmp = lapply(1:nrow(X), function(i){
              x = X[i,]
              res = x[p]*tt$y_exp[i]*x%*%t(x)
              return(res)
            })
            
            val = -Reduce("+", tmp) / length(tmp)
            return(val)
          })
          
          H2 = do.call('rbind', H2)
          bias_analytical = 1/(2*nrow(tt))*Q%*%t(H2)%*%c(Q)
          
          # update the beta hats to s
          beta_hat <- beta_hat - bias_analytical  
        }
        
        # store the model results to return
        model_res = list()
        model_res[[as.character(xx)]][['beta_hat']] <- beta_hat
        model_res[[as.character(xx)]][['beta_vcov']] <- beta_vcov
        
        # bootstrap 
        sapply(1:R_PI, function(r){
          beta_boot <- MASS::mvrnorm(1,beta_hat,beta_vcov)
          pred_boot <- (tt %>% 
                          mutate(intercept=1) %>%
                          #dplyr::select(intercept,year,ends_with("1"),ends_with("2"),ends_with("3")) %>%
                          dplyr::select(intercept,year,cos1,sin1,cos2,sin2,cos3,sin3) %>%
                          as.matrix())%*%as.matrix(beta_boot)
          pred_boot_exp <- exp(pred_boot) 
          pred_boot_exp[which(is.na(pred_boot_exp))] <- 1
          x = rpois(n = nrow(tt), pred_boot_exp)
          
          return(x)
          
        }) -> sim.boot
        
        # get the quantiles and store them
        fitted_quants = t(apply(sim.boot, 1, function(xx){
          quantile(xx, probs = quant_probs)
        }))
        fitted_quants = as.data.frame(fitted_quants)
        colnames(fitted_quants) = paste0(paste0(col, '_pred_WF_'), quant_probs)
        
        # merge the quantiles back into data frame
        tt = cbind(tt, fitted_quants)
        
      }
      tmp_lst = list(tt, sim.boot, model_res)
      return(tmp_lst)
    })
    names(tmp) = uni_group
    
  }
  
  # combine the individual facility results into a larger data frame
  df <- do.call('rbind',lapply(tmp,'[[',1))
  
  # take the sum of the bootstrap samples of each facility (this returns an n X R matrix itself)
  sim.full = Reduce('+', lapply(tmp, '[[', 2))
  
  # initialize district results
  district_df = NULL
  
  # district-level analysis
  for(d in names(dist_fac)){
    tt = data.frame(district = d,
                    date = tmp[[1]][[1]]$date)
    facs = dist_fac[[d]]
    
    # get the sums by district: returns n x R data frame
    sum_district = Reduce('+', lapply(facs, function(f){
      tmp[[f]][[2]]
    }))
    
    # get the quantiles and store them
    fitted_quants = t(apply(sum_district, 1, function(xx){
      quantile(xx, probs = quant_probs)
    }))
    fitted_quants = as.data.frame(fitted_quants)
    if(family == 'poisson'){
      colnames(fitted_quants) = paste0(paste0(col, '_pred_WF_'), quant_probs)
    }else if(family == 'negbin'){
      colnames(fitted_quants) = paste0(paste0(col, '_pred_WF_negbin_'), quant_probs)
    }
    
    # merge the quantiles back into data frame
    tt = cbind(tt, fitted_quants)
    
    district_df = rbind(district_df, tt)
  }
  
  betas <- do.call('rbind',lapply(1:length(tmp), function(ii){
    tmp_fac <- tmp[[ii]][[3]]
    tt <- matrix(NA, ncol = 8)
    tt[1,] <- tmp_fac[[1]][[1]]
    rownames(tt) <- names(tmp_fac)
    colnames(tt) <- names(tmp_fac[[1]][[1]])
    return(tt)
  }))
  
  beta_vcovs <- lapply(1:length(tmp), function(ii){
    return(tmp[[ii]][[3]])
  })
  
  res_lst = list(df = df, district_df = district_df, betas = betas, beta_vcovs = beta_vcovs)
  if(family == 'negbin'){
    
    overdisp <- sapply(1:length(tmp), function(ii){
      tmp_fac <- tmp[[ii]][[3]]
      tt <- tmp_fac[[1]][[3]]
      names(tt) <- names(tmp_fac)
      tt
    })
    res_lst[['overdisp']] <- overdisp
  }
  return(res_lst)
}

### run a complete case analysis for the WF method
# df: df to fit on.
# district_df: df of district-level analysis if that's what is done.
# train_end_date: when to fit the model until.
# ...: extra parameters that can be input into WF_fit.
WF_CCA <- function(df, district_df = NULL, train_end_date = '2019-12-01', ...){
  # replace values in the test set with missing ones
  tmp <- df
  tmp$y[tmp$date > train_end_date] <- NA
  
  res <- WF_fit(tmp, ...)
  
  # store the results and return the original y values
  # This is because the y values in 2020 are returned as NA from WF_fit
  tmp <- res$df
  tmp = merge(tmp %>% dplyr::select(-y),
              df %>% dplyr::select(date, facility, y),
              by = c('date','facility'))
  
  res$df <- tmp
  
  return(res)
}

#### freqGLM functions ####

### Compute the model means given the data and parameters.
model.mean.exp <- function(D, params){
  mu = params[1]*D$y.AR1 + # auto-regressive
    params[2]*D$y.neighbors + # neighbors
    exp(params[3] + params[4]*D$year + params[5]*D$cos1 + params[6]*D$sin1 + params[7]*D$cos2 + params[8]*D$sin2 + params[9]*D$cos3 + params[10]*D$sin3) # yearly + seasonal component
  
  return(mu)
}

### Compute the log likelihood given the data D, parameters, and target column.
ll.wrapper.negbin.exp = function(params, D, target_col){
  mu = model.mean.exp(D, params)
  theta = params[11]
  logL_vec <- dnbinom(D[,target_col], mu = mu, size = theta, log = T)
  logL = sum(logL_vec, na.rm = T)
  return(-logL)
}

### the NNLS freqGLMepi model (so lambda and phi are constrained rather than exponentiated)
fit_freqGLMepi_nnls_negbin <- function(df, num_inits = 10, verbose = T, target_col = 'y_imp', init = NULL, BFGS = NULL){
  t0 = Sys.time()
  
  # set up initialization
  if(is.null(init)){
    init = c(0.1, 0.1, rep(0,8), 5)
  }
  
  parnames = c('By.AR1', 'By.neighbors', 'Bintercept', 'Byear', 'Bcos1', 'Bsin1', 'Bcos2', 'Bsin2', 'Bcos3', 'Bsin3', 'theta')
  names(init) = parnames
  
  init_OG = init
  
  # params = nlminb(start = init, objective = ll.wrapper, D = df, target_col = target_col, lower = c(0, 0, rep(-10, 10)), upper = c(1, 1, rep(10, 10)))
  
  tryCatch({
    params = nlminb(start = init, objective = ll.wrapper.negbin.exp, D = df, target_col = target_col, lower = c(0, 0, rep(-10, 8), 0), upper = c(1, 1, rep(10, 8), 100))
  }, error = function(e){
    browser()
  })
  
  # try different initialization values to compare convergence
  if(num_inits > 1){
    for(i in 2:num_inits){
      # init = rnorm(10,0,10*i/num_inits)
      
      init = init_OG + c(0.01*i, 0.01*i, rnorm(8,0,5*i/num_inits), i)
      
      # nelder-mead
      tryCatch({
        params2 = nlminb(start = init, objective = ll.wrapper.negbin.exp, D = df, target_col = target_col, lower = c(0, 0, rep(-10, 8), 0), upper = c(1, 1, rep(10, 8), 100))
      }, error = function(e){
        browser()
      })
      
      if(params2$objective < params$objective & params2$convergence == 0){
        if(verbose){print('using another initialization')}
        params = params2
      }
    }
  }
  
  # for error checking
  if(params$convergence != 0){
    if(verbose){print('didnt converge for one iteration')}
    #browser()
  }
  
  if(verbose){
    print(sprintf('freqGLMepi fit in %s seconds', Sys.time() - t0))
  }
  
  names(params$par) = parnames
  return(params)
}

### Given data with missing values, fits the freqGLMepi model on all data points
#   df: data
#   max_iter: number of iterations of imputation
#   tol: tolerance of the converage of imputed values
#   individual_facility_models: fit each facility separately

#   R_PI: number of bootstrap iterations if doing so
#   quant_probs: the quantiles of the bootstrap to store in the data frame
#   verbose: printing updates
freqGLMepi_CCA = function(df, train_end_date = '2019-12-01', max_iter = 1, tol = 1e-4, individual_facility_models = T, family = 'negbin', R_PI = 100, quant_probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975), verbose = F, optim_init = NULL, scale_by_num_neighbors = T, blocksize = 10, nnls = T){
  # check that we have the right columns
  if(!('y' %in% colnames(df) & 'y_true' %in% colnames(df))){
    warning('make sure the data has y (with NAs) and y_true')
  }
  
  # get outcome col
  if(family == 'poisson'){
    pred_col <- 'y_pred_freqGLMepi'
  }else if(family == 'negbin'){
    pred_col <- 'y_pred_freqGLMepi_negbin'
  }else{
    stop('please input either family = poisson or negbin.')
  }
  
  # get districts and facilities
  dist_fac <- get_district_facilities(df)
  
  # convert to proper format
  train_end_date <- as.Date(train_end_date)
  
  # split the hold out set and train set
  df_test <- df %>%
    filter(date > train_end_date)
  df_test$y_imp <- NA
  df_test[,pred_col] <- NA
  df <- df %>%
    filter(date <= train_end_date)
  
  # get the adjacency matrix
  W <- make_district_adjacency(df, scale_by_num_neighbors)
  
  # set fitting function to be regular (exponentiated spatial and temporal parameters) or non-negative least squares 
  if(nnls){
    model_function <- model.mean.exp
    if(family == 'poisson'){
      fit_function <- fit_freqGLMepi_nnls
    }else if(family == 'negbin'){
      fit_function <- fit_freqGLMepi_nnls_negbin
    }
  }else{
    model_function <- model.mean
    if(family == 'poisson'){
      fit_function <- fit_freqGLMepi
    }else{
      stop('havent coded non-poisson for the non-exponentiated formula.')
    }
  }
  
  y_pred_list = list()
  
  ### Do initial filling of y
  # setting up the formula
  formula_col = as.formula("y ~ year + cos1 + sin1 + cos2 + sin2 + cos3 + sin3")
  
  # the unique facility groups
  uni_group = unique(df$facility)
  
  # run the individual model for each group.
  if(family == 'poisson'){
    tmp <- lapply(uni_group, function(xx) {
      tt <- df %>% filter(facility == xx)
      
      # run the model
      mod_col <- glm(formula_col, family = 'poisson', data = tt)
      
      # update predictions
      tt[,pred_col] = predict(mod_col, tt, type = 'response')
      
      return(tt)
    })
  }else if(family == 'negbin'){
    tmp <- lapply(uni_group, function(xx) {
      tt <- df %>% filter(facility == xx)
      
      # run the model
      mod_col <- MASS::glm.nb(formula_col, data = tt)
      
      # update predictions
      tt[,pred_col] = predict(mod_col, tt, type = 'response')
      
      # save the theta values
      tt$WF_theta = mod_col[['theta']]
      
      return(tt)
    })
  }else{
    stop('please put in a proper family')
  }
  # combine into one data frame
  df = do.call('rbind',tmp)
  
  # filling in missing values by randomly sampling mean prediction from Poisson
  df$y_imp = df$y
  na.ind = which(is.na(df$y))
  
  if(family == 'poisson'){
    df$y_imp[na.ind] <- rpois(n = length(na.ind), df[na.ind,pred_col,drop=T])
  }else if(family == 'negbin'){
    df$y_imp[na.ind] = MASS::rnegbin(n = length(na.ind), mu = df[na.ind,pred_col,drop=T], theta = df$WF_theta[na.ind])
  }
  
  # add the neighbors and auto-regressive
  df = add_autoregressive(df, 'y_imp') %>%
    add_neighbors(., 'y_imp', scale_by_num_neighbors = scale_by_num_neighbors, W = W)
  
  ### Run freqGLM_epidemic model iteratively
  iter = 1
  y_pred_list[[1]] = df[,pred_col]
  prop_diffs = c(1)
  while(prop_diffs[length(prop_diffs)] > tol & iter <= max_iter){
    iter = iter + 1
    
    if(individual_facility_models){
      # run the individual model for each group.
      tmp <- lapply(uni_group, function(xx) {
        # subset data
        tt <- df %>% filter(facility == xx)
        
        # fit the model
        params = fit_function(tt, verbose = verbose, init = optim_init[[xx]])
        
        # update y_pred
        tt[,pred_col] = model_function(tt, params$par)
        
        return(list(df = tt, params = params))
      })
      
      # get matrix of the parameter estimates
      parmat = sapply(tmp, function(xx) xx[[2]]$par)
      
      # naming the parameter columns and rows
      rownames(parmat) = names(tmp[[1]]$params$par)
      colnames(parmat) = uni_group
      
      # get the convergence of each facility
      convergence = sapply(tmp, function(xx) (xx[[2]]$convergence == 0))
      
      # combine into one data frame
      df = do.call('rbind', lapply(tmp, '[[', 1)) %>% arrange(facility, date)
    }else{
      # fit the model
      parmat = fit_function(df)$par
      
      # update y_pred 
      df[,pred_col] = model_function(df, parmat)
      
    }
    
    # store the predictions for this iteration
    y_pred_list[[iter]] = df[,pred_col]
    
    if(length(na.ind) == 0){
      if(verbose){print('only running one iteration because there is no missingness')}
      break
    }
    
    # update y_imp (need to update this code for negbin if using)
    if(max_iter > 1){
      stop('havent coded for max_iter > 1')
      na.ind.2 = intersect(na.ind, which(!is.na(df[,pred_col])))
      df$y_imp[na.ind.2] <- rpois(n = length(na.ind.2), df[na.ind.2,pred_col])
      
      # compare y_imps
      prop_diffs = c(prop_diffs, mean(abs(y_pred_list[[iter]][na.ind] - y_pred_list[[iter-1]][na.ind])/y_pred_list[[iter-1]][na.ind], na.rm = T))
      
      # update the neighbors and auto-regressive
      df = add_autoregressive(df, 'y_imp') %>%
        add_neighbors(., 'y_imp', W = W)
    }
    
    # update
    if(iter %% 10 == 0 & verbose){
      print(iter)
      if(individual_facility_models){
        print(parmat)
      }else{
        print(params$par)
      }
      
    }
  }
  
  if(prop_diffs[length(prop_diffs)] < tol){
    print(sprintf('convergence reached in %s iterations', iter))
  }
  
  ### run the stationary bootstrap
  param_boots <- lapply(1:length(uni_group), function(i) {
    # subset data
    tt <- df %>% filter(facility == uni_group[i])
    tt_test <- df_test %>% filter(facility == uni_group[i])
    
    param_est_function <- function(xx){
      # fit model on bootstrapped data set
      # start the initialization at the values from the main model
      params = fit_function(xx, BFGS = F, num_inits = 1, verbose = verbose, init = parmat[,1])
      
      # predict model on original training data set
      # x = suppressWarnings(rpois(n = nrow(tt), model_function(tt, params$par)))
      # x = suppressWarnings(MASS::rnegbin(n = nrow(tt), mu = model_function(tt, params$par), theta = params$par['theta']))
      
      return(params$par)
    }
    
    # run the stationary bootstrap and return the parameters from each fit
    sim_boot <- boot::tsboot(tt, statistic = param_est_function, R = R_PI, sim = 'geom', l = blocksize)$t
    
    # match the parameter names
    colnames(sim_boot) <- names(fit_function(tt, BFGS = F, num_inits = 1, verbose = verbose, init = parmat[,1])$par)
    
    return(sim_boot)
  })
  names(param_boots) <- as.character(uni_group)
  
  # get the parameter results from the bootstrap
  param_results <- do.call('rbind', lapply(as.character(uni_group), function(fac){
    tt <- param_boots[[as.character(fac)]]
    tmp <- data.frame(t(apply(tt, 2, function(col) {
      #mean(col)
      vec = c(mean(col), quantile(col, probs = c(.025,0.25,0.5,0.75,0.975), na.rm = T))
      names(vec) <- c('mean', paste0('Q_', names(vec)[-1]))
      vec
    })))
    tmp$param = rownames(tmp)
    tmp$facility = fac
    tmp
  }))
  
  param_results <- param_results[,c(8,7,1:6)]
  rownames(param_results) <- NULL
  
  # create long form of estimated parameters
  tmp = parmat %>%
    as.data.frame(.) %>%
    mutate(param = rownames(.)) 
  par_long = NULL
  for(i in 1:ncol(tmp)){
    par_long = rbind(par_long,
                     data.frame(facility = colnames(tmp)[i],
                                param = tmp$param, 
                                full_estimate = tmp[,i]))
  }
  
  # merge them!
  param_results = merge(param_results, par_long, by = c('facility','param'))
  
  # rename columns appropriately
  param_results$param = gsub('y.neighbors','rho',
                             gsub('y.AR1', 'alpha',
                                  gsub('B','',param_results$param)))
  
  # combine the data sets and split by facility
  df_combined <- rbind(df[,colnames(df_test)], df_test) %>%
    add_autoregressive(., 'y_imp') %>%
    add_neighbors(., 'y_imp', scale_by_num_neighbors = scale_by_num_neighbors, W = W) %>%
    arrange(date)
  df_combined$y_pred <- NA
  
  # cycle through all the bootstrap iterations
  pred_boots <- lapply(1:nrow(param_boots[[1]]), function(i) {
    # resetting the data frame
    df_tmp <- df_combined 
    
    # do all the baseline predictions (+ 1)
    for(f in uni_group){
      tt <- df_tmp %>% filter(facility == f, date <= (train_end_date + 31))
      
      if(family == 'poisson'){
        y = suppressWarnings(rpois(n = nrow(tt), model_function(tt, param_boots[[f]][1,]))) 
      }else if(family == 'negbin'){
        y = suppressWarnings(MASS::rnegbin(n = nrow(tt), mu = model_function(tt, param_boots[[f]][1,]), theta = param_boots[[f]][1,'theta']))
      }
      
      # put the results back into the data frame
      df_tmp$y_pred[df_tmp$facility == f][1:length(y)] <- y
    }
    
    # update the neighbors and auto-regressive terms
    df_tmp <- add_autoregressive(df_tmp, 'y_pred') %>%
      add_neighbors(., 'y_pred', scale_by_num_neighbors = scale_by_num_neighbors, W = W)
    
    # get remaining time points
    time_points <- unique(df_tmp$date[is.na(df_tmp$y_pred) & df_tmp$date > train_end_date]) 
    
    # cycle through remaining time points
    for(t in time_points){
      for(f in uni_group){
        tt <- df_tmp %>% filter(facility == f, date == t)
        
        if(family == 'poisson'){
          y = suppressWarnings(rpois(n = nrow(tt), model_function(tt, param_boots[[f]][i,]))) 
        }else if(family == 'negbin'){
          y = suppressWarnings(MASS::rnegbin(n = nrow(tt), mu = model_function(tt, param_boots[[f]][i,]), theta = param_boots[[f]][i,'theta']))
        }
        
        # put the results back into the data frame
        df_tmp$y_pred[df_tmp$facility == f & df_tmp$date == t] <- y
      }
      
      # update the neighbors and auto-regressive terms
      df_tmp <- add_autoregressive(df_tmp, 'y_pred') %>%
        add_neighbors(., 'y_pred', scale_by_num_neighbors = scale_by_num_neighbors, W = W)
    }
    
    return(df_tmp$y_pred)
  })
  
  # combine into a data frame
  pred_boots <- do.call('cbind', pred_boots)
  
  if(family == 'poisson'){
    pred_name <- 'y_pred_freqGLMepi_'
  }else if(family == 'negbin'){
    pred_name <- 'y_pred_freqGLMepi_negbin_'
  }
  
  # get the quantiles and store them
  fitted_quants = t(apply(pred_boots, 1, function(xx){
    quantile(xx, probs = quant_probs, na.rm = T)
  }))
  fitted_quants = as.data.frame(fitted_quants)
  colnames(fitted_quants) = paste0(pred_name, quant_probs)
  
  # combine the final results and return
  df_combined = cbind(df_combined, fitted_quants)
  
  # set the prediction to the median
  df_combined$y_pred_CCA_freqGLMepi <- df_combined$y_pred_CCA_freqGLMepi_0.5
  
  # remove the ambiguous "y_pred" column
  df_combined$y_pred <- NULL
  
  ### Make the district results
  district_df = data.frame(cbind(df_combined[,c('date','district')], pred_boots)) %>% 
    group_by(date, district) %>%
    summarize_all(sum)
  
  district_mat = district_df[,3:ncol(district_df)]
  # get the quantiles and store them
  fitted_quants = t(apply(district_df[,3:ncol(district_df)], 1, function(xx){
    quantile(xx, probs = quant_probs, na.rm = T)
  }))
  fitted_quants = as.data.frame(fitted_quants)
  colnames(fitted_quants) = paste0(pred_name, quant_probs)
  
  district_df = cbind(district_df[,c('date','district')], fitted_quants)
  
  # prep data to return
  return_lst = list(df = df_combined, district_df = district_df, param_results = param_results, convergence = convergence, y_pred_list = y_pred_list, prop_diffs = prop_diffs)
  
  return(return_lst)
}

#### CAR functions ####

### Make the Leroux precision matrix from a data frame.
# df: Input data frame with columns "facility" and "district". 
# rho: The spatial parameter.
make_precision_mat <- function(df, rho, W2 = NULL){
  # check the rho value
  if(rho < 0 | rho >= 1){
    stop('please input a rho in [0,1)')
  }
  
  # get unique facilities
  facilities = unique(df$facility) %>% sort()
  
  # create the <W2 = diag(W1) - W> matrix
  if(is.null(W2)){
    W2 <- make_district_W2_matrix(df)
  }
  
  # make the final Q matrix
  Q = rho*W2 + (1-rho)*diag(rep(1, length(facilities)))
  
  return(Q)
}


# prep the data for fitting the stan model
prep_stan_data_rushworth_sparse <- function(df, formula, theta_shape = NULL, theta_rate = NULL){
  N = nrow(df)
  N_T <- length(unique(df$date))
  N_F <- length(unique(df$facility))
  
  # order the data frame
  df <- df %>% arrange(date, facility)
  
  # W_star = make_district_W2_matrix(df)
  W = make_district_adjacency(df)
  W_n = sum(W)/2
  W2 = make_district_W2_matrix(df)
  
  # get eigenvalues for determinant calculation
  lambda = eigen(W2 - diag(1, nrow(W2)))$values
  
  # make the complete model matrix
  df2 <- df; df2$y[is.na(df2$y)] <- 0
  X = model.matrix(formula, data = df2)
  
  # get the outcome
  y = df$y
  
  # # comparing missingness.
  # if(!identical(as.integer(rownames(X_obs)), which(!is.na(df$y)))){
  #   stop('mismatch of model matrix and df missing rows')
  # }
  
  # missingness data
  N_miss = sum(is.na(y))
  N_obs = sum(!is.na(y))
  ind_miss = which(is.na(y))
  ind_obs = which(!is.na(y))
  y_obs = y[ind_obs]
  
  # compute the priors
  lm_fit <- glm(formula, family = 'poisson', data = df)
  coef_mat <- summary(lm_fit)$coefficients
  prior_mean_beta <- coef_mat[,1]
  sigma_beta = 10*vcov(lm_fit)
  
  # make the stan data frame
  stan_data <- list(
    N = N, # number of observations
    p = ncol(X), # number of variables
    N_F = N_F, # number of facilities
    N_T = N_T, # number of time points
    N_miss = N_miss,
    N_obs = N_obs,
    ind_miss = ind_miss,
    ind_obs = ind_obs,
    X = X, # design matrix
    y_obs = y_obs, # outcome variable 
    mu = prior_mean_beta, # prior mean
    Sigma = sigma_beta, # prior variance
    #W_star = W_star, 
    W = W,
    W_n = W_n,
    I = diag(1.0, N_F),
    lambda = lambda)
  
  if(!is.null(theta_shape) & !is.null(theta_rate)){
    stan_data$theta_shape <- theta_shape
    stan_data$theta_rate <- theta_rate
  }
  
  return(stan_data)
}


### Fit the CAR model.
# df: data frame with facilities, dates, and outcomes.
# col: Outcome column.
# AR: Autoregression (for CARBayes function - not typically used).
CAR_fitting <- function(df, col, AR = 1, return_type = 'all', model = 'facility_fixed', burnin = 1000, n.sample = 2000, prediction_sample = T, thin = 10, prior = 'none', prior_var_scale = 1, prior_mean = NULL, prior_var = NULL, MALA = T, MCMC_sampler = 'stan', family = 'negbin', theta_shape = 5, theta_rate = 1){
  
  # check if this method has already been run
  if(any(grepl('y_CAR_ST', colnames(df)))){
    print('previous CAR Bayes predictions found. Removing them')
    df[,grep('y_CAR_ST', colnames(df))] <- NULL
  }
  
  # checking that we don't have any single districts
  tt = df %>% filter(date == unique(date)[1]) %>% pull(district) %>% table
  if(any(tt == 1)){
    warning('Randomly choosing a district for a facility without that district (this should be fixed later)')
    # replace the single districts with the biggest ones
    for(nn in names(which(tt == 1))){
      df$district = gsub(nn, names(which.max(tt)), df$district)
    }
  }
  
  #districts = df %>% group_by(district) %>% summarize(n = length(unique(facility)))
  
  # create the adjacency matrix
  D2 = df %>% dplyr::select(district, facility) %>% distinct()
  W = full_join(D2, D2, by = 'district', relationship = 'many-to-many') %>%
    filter(facility.x != facility.y) %>%
    dplyr::select(-district) %>%
    igraph::graph_from_data_frame() %>%
    igraph::as_adjacency_matrix() %>%
    as.matrix()
  
  # model formula
  if(model == 'facility_intercept'){
    formula_col = as.formula(sprintf("%s ~ year + cos1 + sin1 + cos2 + sin2 + cos3 + sin3 + facility", col))
  }else if(model == 'facility_fixed'){
    formula_col = as.formula(sprintf("%s ~ facility + facility*year + facility*cos1 + facility*sin1 + facility*cos2 + facility*sin2 + facility*cos3 + facility*sin3", col))
  }else if(model == 'fixed'){
    formula_col = as.formula(sprintf("%s ~ year + cos1 + sin1 + cos2 + sin2 + cos3 + sin3", col))
  }else{
    print('Using user-defined model')
    formula_col = as.formula(model)
  }
  
  if(prior == 'WF'){
    stop('havent implemented priors for CARstan')
    # Get the WF estimates
    lm_fit <- glm(formula_col, family = 'poisson', data = df)
    coef_mat <- summary(lm_fit)$coefficients
    prior_mean_beta <- coef_mat[,1]
    prior_var_beta <- coef_mat[,2]^2*prior_var_scale
  }else if(prior == 'constant'){
    stop('havent implemented priors for CARstan')
    prior_mean_beta = prior_mean
    prior_var_beta = prior_var
  }else if(prior == 'none'){
    prior_mean_beta = NULL
    prior_var_beta = NULL
  }else{
    stop('please input a proper prior value (WF, constant, none)')
  }
  
  # run CAR Bayes
  if(MCMC_sampler == 'CARST'){
    chain1 <- ST.CARar(formula = formula_col, 
                       family = "poisson",
                       data = df, W = W, 
                       prior.mean.beta = prior_mean_beta,
                       prior.var.beta = prior_var_beta,
                       burnin = burnin, 
                       n.sample = n.sample,
                       thin = thin, 
                       AR = AR, 
                       verbose = F,
                       MALA = MALA)
    
    # check that the prior names matched the CAR fitted names
    if(prior == 'WF'){
      if(!(identical(rownames(coef_mat), rownames(chain1$summary.results)[1:160]))){
        stop('error in WF prior and CAR formula name mismatch')
      }
    }
    
    beta_df = as.data.frame(chain1$samples$beta)
    colnames(beta_df) <- setdiff(gsub('\\(|\\)', '', row.names(chain1$summary.results)), c('tau2', 'rho.S','rho.T'))
    
    model_chain <- list(
      fitted_mean = chain1$fitted.values,
      fitted = chain1$samples$fitted,
      beta = beta_df,
      phi = chain1$samples$phi,
      rho = chain1$samples$rho[,'rho.S'],
      alpha = chain1$samples$rho[,'rho.T'],
      tau2 = chain1$samples$tau2,
      CARST_summary = chain1$summary.results
    )
    
  }else if(MCMC_sampler == 'stan'){
    if(family == 'poisson'){
      stan_pars <- c('tau2','rho','alpha','beta')
      stan_data <- prep_stan_data_rushworth_sparse(df, formula_col)
      stan_fit <- stan(file = "code/regression_rushworth_sparse.stan",
                       data = stan_data, 
                       iter = n.sample, 
                       warmup = burnin,
                       chains = 1, 
                       init = '0',
                       cores = 1)
    }else if(family == 'negbin'){
      stan_pars <- c('tau2','rho','alpha','beta','theta')
      stan_data <- prep_stan_data_rushworth_sparse(df, formula_col, theta_shape = theta_shape, theta_rate = theta_rate)
      stan_fit <- stan(file = "code/regression_rushworth_sparse_negbin.stan",
                       data = stan_data, 
                       iter = n.sample, 
                       warmup = burnin,
                       chains = 1, 
                       init = '0',
                       cores = 1)
      # make sure to check the results - is theta ordered properly?
      # browser()
      # test <- grep('theta',names(stan_fit))
      # stan_fit[1,1,test[21:80]]
    }
    
    # browser()
    # out_file <- 'C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results/CAR_stan_fit_real_data_first_eval_date_05122025.RData'
    # save(stan_fit, file = out_file)
    # extract out the important features from the model
    stan_out <- rstan::extract(stan_fit)
    
    # pull out the beta params
    beta_df <- as.data.frame(stan_out$beta)
    model_col_names <- gsub('\\(|\\)', '', colnames(stan_data$X))
    colnames(beta_df) <- model_col_names
    
    # pull out the summary values
    stan_summary = summary(stan_fit, pars = stan_pars)$summary
    rownames(stan_summary)[grep('beta', rownames(stan_summary))] <- model_col_names
    
    model_chain = list(
      fitted_mean = apply(stan_out$y_exp, 2, mean),
      fitted = stan_out$y_exp,
      beta = beta_df,
      phi = stan_out$phi,
      rho = stan_out$rho,
      alpha = stan_out$alpha,
      tau2 = stan_out$tau2,
      CARstan_summary = stan_summary
    )
    
    if(family == 'negbin'){
      model_chain[['theta']] = stan_out$theta
    }
  }else{
    stop('input a proper MCMC sampler')
  }
  
  df[,paste0(col, '_CAR_ST')] = model_chain$fitted_mean
  
  # Sample the fitted values for the posterior predictive distribution
  if(prediction_sample){
    # browser()
    # pull the fitted values and randomly select prediction values based on the Poisson distribution
    tt = model_chain$fitted
    tt = apply(tt, 2, function(xx) rpois(n = length(xx), lambda = xx))
    
    # get the quantiles of fitted values 
    quant_probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)
    fitted_quants = t(apply(tt, 2, function(xx){
      quantile(xx, probs = quant_probs)
    }))
    fitted_quants = as.data.frame(fitted_quants)
    colnames(fitted_quants) = paste0(paste0(col, '_CAR_ST_'), quant_probs)
    
    # ignore the sampling of the Poisson and take mean fitted value quantiles
  }else{
    # get the quantiles of fitted values 
    quant_probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)
    fitted_quants = t(apply(model_chain$fitted, 2, function(xx){
      quantile(xx, probs = quant_probs)
    }))
    fitted_quants = as.data.frame(fitted_quants)
    colnames(fitted_quants) = paste0(paste0(col, '_CAR_ST_'), quant_probs)
  }
  
  df = as.data.frame(cbind(df, fitted_quants))
  
  # mapping matrix from each data point to a date
  dates = sort(unique(df$date))
  mat = as.matrix(sparseMatrix(i = 1:nrow(df), j = match(df$date, dates), x = 1))
  
  # condense the fitted value samples to each date
  date_fits = model_chain$fitted %*% mat
  
  # get the quantiles at the county level
  fitted_quants_C = t(apply(date_fits, 2, function(xx){
    quantile(xx, probs = quant_probs)
  }))
  fitted_quants_C = as.data.frame(fitted_quants_C)
  colnames(fitted_quants_C) = paste0(paste0(col, '_CAR_ST_'), quant_probs)
  
  if(return_type == 'data.frame'){
    return(df)
  }else if(return_type == 'model'){
    return(chain1)
  }else if(return_type == 'all'){
    #lst = list(facility_df = df, county_df = df_county, model_chain = model_chain)
    lst = list(facility_df = df, model_chain = model_chain)
    return(lst)
  }
}


### Wrapper to fit the model on all data using the CAR_fitting function and return results
# df: data frame with facilities, dates, and outcomes.
# R_posterior: A specified number of posterior simulations to run (if you want it to be smaller than the number of returned posterior samples from CAR).
# train_end_date: cutoff date for the training data to fit the model.
# predict_start_date: the starting time point for where predictions. should be run. If null, defaults to all dates after train_end_date.
# col: outcome column.
# quant_probs: quantiles to be returned from prediction samples.
# return_chain: whether to return model chains.
# return_raw_fit: whether to return the full MCMC model fit.
# use_fitted_phi: Whether to use fitted phi values for posterior predictions. This only works if doing one-step ahead predictions. If false, resample phi across time points for posterior predictions.
# family: CAR family to fit on.
# one_prediction_date: If there is only one date to predict on.
CAR_wrapper <- function(df, R_posterior = NULL, train_end_date = '2019-12-01', predict_start_date = NULL, col = 'y', quant_probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975), return_chain = F, return_raw_fit = F, use_fitted_phi = F, model_rename = NULL, family = 'negbin', one_prediction_date = T, ...){
  
  # get districts and facilities
  dist_fac <- get_district_facilities(df)
  facilities <- unique(df$facility)
  
  # get the max date
  max_date <- max(df$date)
  
  if(is.null(predict_start_date)){
    dates = df %>% 
      dplyr::filter(date > train_end_date) %>%
      dplyr::select(date) %>%
      unique() %>%
      .$date
    predict_start_date = min(dates)
    if(one_prediction_date){
      df = df %>%
        filter(date <= predict_start_date)
    }
  }else{
    dates = df %>% 
      dplyr::filter(date >= predict_start_date) %>%
      dplyr::select(date) %>%
      unique() %>%
      .$date
  }
  
  if(predict_start_date > min(dates) & !use_fitted_phi){
    warning(sprintf('Setting phi to 0 at %s in the posterior prediction simulation', predict_start_date))
  }
  
  # make a train set to fit on
  train <- df %>%
    filter(date <= train_end_date) %>%
    arrange(facility, date)
  
  # checking the sorting
  order_check <- order(train$facility, train$date)
  if(!identical(order(order_check), order_check)){
    stop('incorrect ordering of train data frame')
  }
  
  # fit the model!
  res <- CAR_fitting(train, col, family = family, ...)
  
  # return the fit if not doing post-processing
  if(return_raw_fit){
    return(res)
  }
  
  #### the rest is for future model prediction
  betas <- res$model_chain$beta
  phi_fit <- res$model_chain$phi
  rho <- res$model_chain$rho
  alpha <- res$model_chain$alpha
  tau2 <- res$model_chain$tau2
  if(family == 'negbin'){
    thetas <- res$model_chain$theta
    names(thetas) <- facilities
  }
  
  # set the R_posterior
  if(is.null(R_posterior)){
    R_posterior = nrow(betas)
  }else if(R_posterior > nrow(betas)){
    stop('cant sample more posterior predictions than the original model fit returns')
  }else{
    # havent implemented yet
    browser()
  }
  
  # group the betas into their separate facility values
  fac_beta_list <- list()
  beta_ref <- betas[,c("Intercept", "year", "cos1", "sin1", "cos2", "sin2", "cos3", "sin3")]
  for(f in facilities){
    # if this is the reference facility (likely A1 in my simulations)
    if(sum(grepl(f, names(betas))) == 0){
      beta_f <- beta_ref
    }else{
      cols <- paste0('facility', f, c("",":year", ":cos1", ":sin1", ":cos2", ":sin2", ":cos3", ":sin3"))
      beta_f <- betas[,cols] + beta_ref 
      colnames(beta_f) <- c("Intercept", "year", "cos1", "sin1", "cos2", "sin2", "cos3", "sin3")
    }
    fac_beta_list[[f]] <- beta_f
  }
  
  # pull the fixed effects for the simulation
  fixed_effects <- lapply(facilities, function(f){
    tmp = df %>% 
      filter(facility == f,
             date >= predict_start_date) %>%
      arrange(date)
    
    # keep the 1 for intercepts
    X = tmp %>% 
      mutate(Intercept = 1) %>%
      dplyr::select(Intercept, year, cos1, sin1, cos2, sin2, cos3, sin3)
    
    rownames(X) = tmp$date
    
    # check that the ordering is correct
    if(!identical(names(X), names(fac_beta_list[[f]]))){
      browser()
    }
    
    betas <- as.matrix(fac_beta_list[[f]])
    
    # (# posterior fits x # params) x (# params x # of data points)
    mean_sims <- betas%*%t(X)
    return(mean_sims)
  })
  names(fixed_effects) <- facilities
  
  # make the "W2" matrix only once for repeated use
  W2 <- make_district_W2_matrix(df)
  
  # make the spatio-temporal precision matrices
  covar_mats <- lapply(1:R_posterior, function(ii){
    Q = make_precision_mat(df, rho = rho[ii], W2)
    
    tryCatch({
      V = tau2[ii]*solve(Q)
    }, error = function(e){
      print(e)
      print('havent dealt with non-invertible precision matrices yet')
      browser()
    })
    return(V)
  })
  
  # make R sampled sets of data
  phi_lst = lapply(1:R_posterior, function(i){
    ### get the spatio-temporal random effects
    # pull in the phi from the stan fit
    phi_fit_r <- matrix(phi_fit[i,], 
                        byrow = F, 
                        nrow = length(unique(train$date)))
    
    # initialize phi to sample
    phi = matrix(0, nrow = length(dates), ncol = length(facilities))
    colnames(phi) = facilities
    
    if(length(dates) > nrow(phi)+1 | length(dates) > nrow(phi_fit_r) + 1){
      stop('dimension mismatch. Probably due to more than one evaluation point.')
    }
    
    # first time step (at 2016-01-01, not 2020-01-01)
    phi[1,] = MASS::mvrnorm(n = 1, mu = rep(0, ncol(phi)), Sigma = covar_mats[[i]])
    
    # cycle through other time steps, using auto-correlated priors
    if(use_fitted_phi){
      for(t in 2:length(dates)){
        phi[t,] = MASS::mvrnorm(n = 1, mu = alpha[i]*phi_fit_r[t-1,], Sigma = covar_mats[[i]])
      }
    }else{
      for(t in 2:length(dates)){
        phi[t,] = MASS::mvrnorm(n = 1, mu = alpha[i]*phi[t-1,], Sigma = covar_mats[[i]])
      }
    }
    
    # convert to matching format
    phi = as.data.frame(phi)
    phi$date = dates
    phi_df = tidyr::gather(phi, facility, phi, setdiff(colnames(phi), c('facility','date')))
    
    return(phi_df)
  })
  
  # rearrange phi_lst to match the format of fixed effects
  phi_by_fac <- lapply(facilities, function(f){
    sapply(phi_lst, function(xx){
      xx %>% 
        filter(facility == f) %>%
        arrange(date) %>%
        pull(phi)
    })
  })
  names(phi_by_fac) <- facilities
  
  # get mean predictions at each point
  mean_pred <- lapply(facilities, function(f){
    # get mean fits
    res <- t(fixed_effects[[f]]) + phi_by_fac[[f]]
    
    # make sure dimensions line up
    if(nrow(res) != length(dates)){
      browser()
    }
    return(res)
  })
  names(mean_pred) <- facilities
  
  # predict new points using the fixed effects, phi values and poisson distribution
  predicted_vals <- lapply(facilities, function(f){
    # get mean fits
    res <- mean_pred[[f]]
    
    # get poisson sampled values from mean
    if(family == 'poisson'){
      tt <- apply(res, 2, function(xx){
        rpois(length(dates), exp(xx))
      })
    }else if(family == 'negbin'){
      theta <- thetas[[f]]
      tt <- apply(res, 2, function(xx){
        rnbinom(n = length(xx), mu = exp(xx), size = theta)
      })
    }
    
  })
  names(predicted_vals) <- facilities
  
  # get the quantiles of the predictions across the simulations
  fitted_quants = do.call('rbind', lapply(facilities, function(f){
    res <- as.data.frame(t(apply(predicted_vals[[f]], 1, function(xx){
      quantile(xx, probs = quant_probs)
    })))
    colnames(res) <- paste0(paste0(col, '_pred_CAR_'), quant_probs)
    res$facility = f
    res$date = dates
    return(res)
  }))
  
  fitted_quants$y_pred_CAR <- fitted_quants$y_pred_CAR_0.5
  
  # merge the results together
  df <- merge(df, fitted_quants, by = c('date', 'facility'), all = T)
  
  # initialize district results
  district_df = NULL
  
  # district-level analysis
  for(d in names(dist_fac)){
    tt = data.frame(district = d,
                    date = dates)
    facs = dist_fac[[d]]
    
    # get the sums by district: returns n x R data frame
    sum_district = Reduce('+', lapply(facs, function(f){
      predicted_vals[[f]]
    }))
    
    # get the quantiles and store them
    fitted_quants = t(apply(sum_district, 1, function(xx){
      quantile(xx, probs = quant_probs)
    }))
    fitted_quants = as.data.frame(fitted_quants)
    colnames(fitted_quants) = paste0(paste0(col, '_pred_CAR_'), quant_probs)
    
    # merge the quantiles back into data frame
    tt = cbind(tt, fitted_quants)
    
    district_df = rbind(district_df, tt)
  }
  
  res_lst <- list(df = df, district_df = district_df)
  
  if(return_chain){
    res_lst[['model_chain']] <- res$model_chain
  }
  
  if(list(...)$MCMC_sampler == 'stan'){
    res_lst[['CARstan_summary']] <- res$model_chain$CARstan_summary
    colnames(res_lst$df) <- gsub('y_pred_CAR', ifelse(is.null(model_rename), 'y_CARstan', model_rename), colnames(res_lst$df))
    colnames(res_lst$district_df) <- gsub('y_pred_CAR', ifelse(is.null(model_rename), 'y_CARstan', model_rename), colnames(res_lst$district_df))
  }
  
  if(list(...)$MCMC_sampler == 'CARST'){
    res_lst[['CARST_summary']] <- res$model_chain$CARST_summary
    colnames(res_lst$df) <- gsub('y_pred_CAR', ifelse(is.null(model_rename), 'y_CARST', model_rename), colnames(res_lst$df))
    colnames(res_lst$district_df) <- gsub('y_pred_CAR', ifelse(is.null(model_rename), 'y_CARST', model_rename), colnames(res_lst$district_df))
  }
  
  return(res_lst)
}

#### Post-processing functions ####
# process the imputed list for metric calculations
clean_data_list <- function(predicted_list, dates = '2020-01-01',  min_date = NULL, rm_ARna = F, imputed_only = F){
  # filter to only be greater than the specified date
  if(!is.null(min_date)){
    # print(sprintf('only getting metrics with dates on or after %s', min_date))
    predicted_list <- lapply(predicted_list, function(xx){
      xx <- xx %>% dplyr::filter(date >= min_date)
    })
  }
  
  # filter to only be the specified date
  if(!is.null(dates)){
    # print(sprintf('only getting metrics with dates on  %s', dates))
    predicted_list <- lapply(predicted_list, function(xx){
      xx <- xx %>% dplyr::filter(date == dates)
    })
  }
  
  # removing the starting points with NA AR1 values, since these
  if(rm_ARna){
    print('removing the starting points because of NA autoregressive term')
    predicted_list = lapply(predicted_list, function(xx) xx[!is.na(xx$y.AR1),])
  }
  
  # remove all non-missing points
  if(imputed_only){
    predicted_list = lapply(predicted_list, function(xx) xx[is.na(xx$y),])
  }
  
  return(predicted_list)
}

# calculate the metrics across simulations
calculate_metrics <- function(predicted_list, methods = c("y_pred_WF", "y_freqGLMepi", "y_CARBayes"), results_by_point = F, date = '2020-01-01', min_date = NULL, rm_ARna = F, imputed_only = F,  use_point_est = F, k = NULL, district_results = F, family = 'negbin', ...){
  
  if(results_by_point){
    # getting the results for each point across all simulations
    avg_fxn = rowMeans
  }else{
    # getting the results for each simulation across all points
    avg_fxn = colMeans
  }
  
  # process the predicted_list
  predicted_list = clean_data_list(predicted_list, date, min_date, rm_ARna, imputed_only)
  
  ## Get the outcome values across all simulations
  {
    if(district_results){
      # updating the y_true and y missing since we don't need those
      y_true <- do.call('cbind', lapply(predicted_list, function(xx) xx[,'y']))
      rownames(y_true) <- predicted_list[[1]]$district
      # districts = do.call('cbind', lapply(predicted_list, function(xx) xx[,'district']))
      y_missing = matrix(NA, nrow = nrow(y_true), ncol = ncol(y_true))
    }else{
      # get the true values everywhere and at the deleted time points
      y_true = do.call('cbind', lapply(predicted_list, function(xx) xx[,'y_true']))
      y_missing = do.call('cbind', lapply(predicted_list, function(xx) {
        y_true = xx[,'y_true'];
        y_true[!is.na(xx[,'y'])] = NA
        y_true
      }))
    }
    
    # get the expected y values and the variance associated with them
    y_exp = do.call('cbind', lapply(predicted_list, function(xx) xx[,'y_exp']))
    y_var = tryCatch({
      y_var = do.call('cbind', lapply(predicted_list, function(xx) xx[,'y_var']))
    }, error = function(e){
      warning('making var(Y) = E(Y). This DOES NOT hold under the CAR DGP')
      y_var = y_exp
    })
    
    # calculate the outbreak values
    y_outbreak3 <- y_exp + 3*sqrt(y_var)
    y_outbreak5 <- y_exp + 5*sqrt(y_var)
    y_outbreak10 <- y_exp + 10*sqrt(y_var)
    
    # numeric missing matrix
    missing_mat <- apply(y_missing, 2, function(xx) 1 - as.numeric(is.na(xx)))
    missing_mat_NA <- missing_mat; missing_mat_NA[missing_mat_NA == 0] <- NA
    
    # get the number of times each data point was missing across simulations
    num_missing = apply(y_missing, 1, function(xx) sum(!is.na(xx)))
  }
  
  df = NULL
  for(method in methods){
    print(method)
    lower_025 = do.call('cbind', lapply(predicted_list, function(xx) xx[,paste0(method, '_0.025')]))
    upper_975 = do.call('cbind', lapply(predicted_list, function(xx) xx[,paste0(method, '_0.975')]))
    lower_25 = do.call('cbind', lapply(predicted_list, function(xx) xx[,paste0(method, '_0.25')]))
    upper_75 = do.call('cbind', lapply(predicted_list, function(xx) xx[,paste0(method, '_0.75')]))
    
    if(use_point_est){
      point_est = do.call('cbind', lapply(predicted_list, function(xx) xx[,method]))
      #tmp$point_est <- sapply(1:nrow(point_est), function(ii) mean(point_est[ii,]))
      outcome = point_est
    }else{
      median = do.call('cbind', lapply(predicted_list, function(xx) xx[,paste0(method, '_0.5')]))
      #tmp$median <- sapply(1:nrow(median), function(ii) mean(median[ii,]))
      outcome = median
    }
    
    # prep the temporary results data frame
    if(results_by_point){
      if(district_results){
        tmp = predicted_list[[1]][,c('date','district')]
      }else{
        tmp = predicted_list[[1]][,c('date','facility','district')]
      }
      tmp$num_missing = num_missing
      
    }else{
      tmp = data.frame(r = 1:length(predicted_list))
    }
    
    tmp$method = method
    
    # point estimation metrics
    tmp$bias = avg_fxn(sapply(1:ncol(outcome), function(ii) {outcome[,ii] - y_true[,ii]}))
    tmp$relative_bias = avg_fxn(sapply(1:ncol(outcome), function(ii) {(outcome[,ii] - y_true[,ii])/y_exp[,ii]}))
    tmp$absolute_bias = avg_fxn(sapply(1:ncol(outcome), function(ii) {abs(outcome[,ii] - y_true[,ii])}))
    tmp$MAPE = avg_fxn(sapply(1:ncol(outcome), function(ii) {abs(outcome[,ii] - y_true[,ii])/y_true[,ii]}))
    tmp$RMSE = sqrt(avg_fxn(sapply(1:ncol(outcome), function(ii) {(outcome[,ii] - y_true[,ii])^2})))
    
    # coverage metrics and specificity
    tmp$coverage50 = avg_fxn(sapply(1:ncol(lower_25), function(ii) (y_true[,ii] >= lower_25[,ii] & y_true[,ii] <= upper_75[,ii])))
    tmp$coverage95 = avg_fxn(sapply(1:ncol(lower_25), function(ii) (y_true[,ii] >= lower_025[,ii] & y_true[,ii] <= upper_975[,ii])))
    tmp$specificity = avg_fxn(sapply(1:ncol(lower_25), function(ii) (y_true[,ii] <= upper_975[,ii])))
    
    # outbreak detection metrics
    tmp$outbreak_detection3 <- avg_fxn(sapply(1:ncol(y_exp), function(ii) y_outbreak3[,ii] >= upper_975[,ii]))
    tmp$outbreak_detection5 <- avg_fxn(sapply(1:ncol(y_exp), function(ii) y_outbreak5[,ii] >= upper_975[,ii]))
    tmp$outbreak_detection10 <- avg_fxn(sapply(1:ncol(y_exp), function(ii) y_outbreak10[,ii] >= upper_975[,ii]))
    
    # measure of how wide the 95% prediction intervals are
    tmp$interval_width = avg_fxn(upper_975 - lower_025) 
    tmp$prop_interval_width = avg_fxn((upper_975 - lower_025)/y_true)
    
    # update the results
    df = rbind(df, tmp)
  }
  
  #res_lst = list(df = df, num_missing = num_missing)
  return(df)
}