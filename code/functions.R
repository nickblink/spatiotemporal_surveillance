
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
make_district_adjacency <- function(df, scale_by_num_neighbors = F){
  
  # get the adjacency matrix
  D2 = df %>% dplyr::select(district, facility) %>% 
    distinct() %>%
    arrange(facility)
  W = full_join(D2, D2, by = 'district') %>%
    filter(facility.x != facility.y) %>%
    dplyr::select(-district) %>%
    igraph::graph_from_data_frame() %>%
    igraph::as_adjacency_matrix() %>%
    as.matrix()
  
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
  pairs = full_join(D2, D2, by = 'district') %>%
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

### Function to simulate the data for a variety of situations
# district_sizes: number of facilities in the districts.
# R: number of simulated data sets.
# seed: random seed for reproducibility.
# type: Model to simulate data. Can be 'WF', 'freqGLM', or 'CAR'.
# family: Family of simulated data. Can be 'poisson' or 'quasipoisson'.
simulate_data <- function(district_sizes, R = 1, seed = 10, type = 'WF', family = 'poisson', ...){
  
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
  betas = sample_betas(facilities, ...)  
  
  if(family == 'poisson'){
    DGP_function = rpois
  }else if(family == 'quasipoisson'){
    DGP_function <- function(n, mu){
      y <- rqpois(n, mu, theta = list(...)$theta)
    }
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
        tmp$y_exp = exp(mu)
        
        # get Poisson or quasipoisson variance
        if(family == 'poisson'){
          tmp$y_var <- tmp$y_exp
        }else if(family == 'quasipoisson'){
          tmp$y_var <- list(...)$theta*tmp$y_exp
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
    # extracting values.
    rho = list(...)$rho
    alpha = list(...)$alpha
    tau2 = list(...)$tau2
    
    # checking all values are in parameters.
    if(is.null(rho)){stop('please put in a rho_DGP value with CAR DGP')}
    if(is.null(alpha)){stop('please put in a alpha_DGP value with CAR DGP')}
    if(is.null(tau2)){stop('please put in a tau2_DGP value with CAR DGP')}
    
    # add in the mean effects because these are the same for all simulations.
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
    
    # make the spatio-temporal precision matrix.
    Q = make_precision_mat(df, rho = rho)
    
    tryCatch({
      V = tau2*solve(Q)
    }, error = function(e){
      print(e)
      print('havent dealt with non-invertible precision matrices yet')
      browser()
    })
    
    # checking ordering of facilities matches.
    if(!identical(colnames(V), as.character(facilities))){
      browser()
      stop('names of covariances and facilities dont match')
    }
    
    # add in the marginal variance to original df.
    dV = diag(V)
    matid = match(df$facility, names(dV))
    df$sigma2_marginal = dV[matid]
    
    # make R sampled sets of data.
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
        stop('input a proper family')
      }
      
      # simulate the observed values
      df$y = DGP_function(nrow(df), exp(df$mu + df$phi))
      
      return(df)
    })
    
    # cycle through each created data frame.
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
    # extracting values
    rho = list(...)$rho
    alpha = list(...)$alpha
    
    # checking values are in list.
    if(is.null(rho)){stop('please put in a rho_DGP value with freqGLM DGP')}
    if(is.null(alpha)){stop('please put in a alpha_DGP value with freqGLM DGP')}

    # get the adjacency matrix.
    W <- make_district_adjacency(df, scale_by_num_neighbors = T)
    
    # add in the mean effects because these are the same for all simulations.
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
    
    # make R sampled sets of data.
    df_lst = lapply(1:R, function(r){
      
      ### Get the first time point values
      ind = which(df$date == min(dates))
      
      # the adjustment accounts for the fact that there aren't additive auto-regressive and spatial terms at the first time point.
      # this adjustment comes from the sum of a geometric series (since this is roughly the effect that the spatial and autoregressive terms approach as we increase the time series)
      adjustment = 1 + rho/(1-rho) + alpha/(1-alpha)
      df$y_exp[ind] = df$y_var[ind] = adjustment*exp(df$mu[ind])
      df$y[ind] = DGP_function(length(ind), df$y_exp[ind])
      
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
        }else if(family == 'quasipoisson'){
          df$y_var[ind] <- list(...)$theta*df$y_exp[ind]
        }else{
          stop('input a proper family')
        }
        
        # predict!
        df$y[ind] = DGP_function(length(ind), df$y_exp[ind])
        
        # update the neighbors and auto-regressive terms
        df <- add_autoregressive(df, 'y') %>%
          add_neighbors(., 'y', scale_by_num_neighbors = T, W = W)
      }
      # ggplot(df, aes(x = date, y = y)) + 
      #   geom_line() + 
      #   facet_wrap(~facility)
      df
    })
    
    # group the results by district.
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
  
  # make list of values to return.
  res_lst = list(df_list = df_lst, district_list = district_lst, betas = betas)
  return(res_lst)
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
