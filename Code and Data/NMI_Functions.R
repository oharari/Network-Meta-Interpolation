#' Turning IPD into five subgroup analyses
#'
#' Matching the data format of NMI for the IPD dataset with a dichotomous 
#' outcome and binary effect modifiers.
#'
#' @param data An IPD dataset.
#' @param EM_cols A vector of effect modifier column names.
#' @param Trt_col Treatment assignment column name.
#' @param outcome_col Binary outcome column name
#' 
#' @return A data frame consisting of 2*p + 1 rows, the first of which contains 
#' the TE and SE for the entire trial, and the remaining 2*p contain subgroup 
#' analyses with missing EM values for all but one effect modifier.
#'
#' @export
GLMLOGIT = function(data, EM_cols, Trt_col, outcome_col){
  k = length(EM_cols)
  
  helper = function(i){
    dat = data[data[,EM_cols[i]] == 1,]
    dat$Tr = dat[,Trt_col]
    dat$Y = dat[,outcome_col]
    X = apply(dat[,EM_cols], 2, mean) 
    
    U = summary(glm(Y ~ Tr, data = dat, family = "binomial"))
    betase1 = c(X, U[["coefficients"]][2, 1:2])
    
    dat = data[data[,EM_cols[i]] == 0,]
    dat$Tr = dat[,Trt_col]
    dat$Y = dat[,outcome_col]
    X = apply(dat[,EM_cols], 2, mean) 
    
    U = summary(glm(Y ~ Tr, data = dat, family = "binomial"))
    betase2 = c(X, U[["coefficients"]][2, 1:2])
    
    
    rbind(betase1, betase2)
  }
  
  out = do.call(rbind, lapply(as.list(1:k), helper))
  all = summary(glm(Y ~ Tr, data = data, family = "binomial"))
  X = apply(data[,EM_cols], 2, mean)
  out = rbind(c(X, all[["coefficients"]][2, 1:2]), out)
  row.names(out) = c()
  
  Trt = as.data.frame(t(sapply(1:nrow(out), 
                               function(x){
                                 sort(unique(data[,Trt_col]))
                               }
  )
  )
  )
  names(Trt) = paste0('Trt', 1:2)
  
  return(cbind(Trt, out))
}




#' Enriching the NMI AgD table
#'
#' Imputing missing covariate values, using the Best Linear Unbiased Predictor.
#'
#' @param IPD An IPD dataset.
#' @param AgD An NMI-format AgD dataset.
#' @param AgD_EM_cols A vector of effect modifier column names in the AgD.
#' @param IPD_EM_cols A vector of effect modifier column names in the IPD.
#' @param Study_col Study identifier column name for the AgD.
#' @param samp_sizes A vector of sample sizes for the AgD studies
#' @param AgD_Trt_cols a vector of treatment column names for the AgD.
#' @param TE_col Relative treatment effect column name for the AgD.
#' @param IPD_Trt_col Treatment column names for the IPD.
#' @param SE_col Standard error column name for the AgD.
#' 
#' @return The enriched (imputed) data frame, with the summarised IPD data 
#' concatenated at the bottom. 
#'
#' @export
BLUP_impute = function(IPD, AgD, AgD_EM_cols, IPD_EM_cols, Study_col, 
                       samp_sizes, AgD_Trt_cols, TE_col, IPD_Trt_col, SE_col){
  rho = cor(IPD[, IPD_EM_cols])
  n_studs = length(unique(AgD[,Study_col]))
  studies = 1:n_studs
  p = length(AgD_EM_cols)
  
  single_BLUP_impute = function(study){
    n = samp_sizes[study]
    Xbar = AgD[AgD[, Study_col] == study, AgD_EM_cols][1,]
    Sbar = sqrt(Xbar*(1 - Xbar)/n)
    
    missing_mat = function(i){
      a = rep(NA, 2*p + 1)
      a[1] = Xbar[i]
      a[2*i] = 1
      a[2*i + 1] = 0
      
      return(a)
    }
    
    X = sapply(1:p, missing_mat)
    Y = X
    for(i in 1:(2*p)){
      ind = ceiling(i/2)
      for(j in 1:p){
        temp = Sbar[j]/Sbar[ind]*rho[ind,j]*(X[-1,][i,ind] - Xbar[ind]) + Xbar[j]
        temp = unlist(ifelse(Xbar[j] == 1, 1, ifelse(Xbar[j] == 0, 0, temp)))
        Y[-1,][i,j] = max(min(temp, 1), 0)
      }
    }
    
    return(Y)
  }
  
  
  imputed = do.call(rbind, lapply(as.list(studies), single_BLUP_impute))
  out = cbind(AgD[, c(Study_col, AgD_Trt_cols)], 
              imputed, AgD[, c(TE_col, SE_col)])
  names(out)[4:(p+3)] = AgD_EM_cols 
  
  IPD_summ = cbind(Study = max(out$Study) + 1,
                   GLMLOGIT(IPD, IPD_EM_cols, IPD_Trt_col, outcome_col)
  )
  
  names(IPD_summ) = names(out)
  
  
  return(rbind(out, IPD_summ))
}





#' Converting the sigma_hat vector to a Variance-Covariance matrix 
#'
#' @param sigma_hat The estimated variance-covariance vector.
#' 
#' @return A symmetric variance-covariance matrix.
#'
#' @export
sigma_hat_vec_to_mat = function(sigma_hat){
  M = length(sigma_hat)
  K = (sqrt(1 + 8*M) - 1)/2
  
  C = matrix(0, nrow = K, ncol = K)
  t = K + 1
  for(i in 1:(K-1)){
    C[i,(i+1):K] = sigma_hat[t:(t + K - i - 1)]
    t = t + K - i
  }
  
  C = C + t(C)
  diag(C) = sigma_hat[1:K]
  
  return(C)
}






#' NMI interpolation pre-NMA
#'
#' Interpolating treatment effect estimates and standard errors at new 
#' effect modifier values.
#'
#' @param IPD An IPD dataset.
#' @param AgD An NMI-format AgD dataset.
#' @param x_vect a vector consisting of the effect modifier values for ITC.
#' @param AgD_EM_cols A vector of effect modifier column names in the AgD.
#' @param IPD_EM_cols A vector of effect modifier column names in the IPD.
#' @param Study_col Study identifier column name for the AgD.
#' @param samp_sizes A vector of sample sizes for the AgD studies
#' @param AgD_Trt_cols a vector of treatment column names for the AgD.
#' @param TE_col Relative treatment effect column name for the AgD.
#' @param IPD_Trt_col Treatment column names for the IPD.
#' @param SE_col Standard error column name for the AgD.
#' levels at which the ITC will be conducted.
#' 
#' @return a list consisting of three data frames -
#' \item{Imputed}{The imputed AgD - the output of \code{\link{BLUP_impute}}.}
#' \item{Final}{The final NMA-format data set, after interpolation at the new
#' effect modifier levels.}
#' \item{Diagnostics}{Observed and predicted TEs for goodness of interpolation
#' diagnostics.}
#'
#' @export
NMI_interpolation = function(IPD, AgD, x_vect, AgD_EM_cols, IPD_EM_cols, 
                             Study_col, samp_sizes, AgD_Trt_cols, TE_col, 
                             SE_col, IPD_Trt_col){
  imputed = BLUP_impute(IPD, AgD, AgD_EM_cols, IPD_EM_cols, Study_col, 
                        samp_sizes, AgD_Trt_cols, TE_col, IPD_Trt_col, SE_col)
  studies = unique(imputed$Study)
  
  single_study_interpolation = function(study){
    dat = imputed %>% filter(Study == study)
    X = apply(as.matrix(dat[, AgD_EM_cols]), 2, as.numeric)
    x_orig = X[1,]
    m = ncol(X)
    M2 = as.matrix(cbind(1, X^2, 2*X, 
                         apply(combn(1:m, 2), 2, 
                               function(u){2*X[,u[1]]*X[,u[2]]})))
    
    beta_hat = lm(as.numeric(dat$TE) ~ X)$coef
    sigma_hat = c(t(M2)%*%solve(M2%*%t(M2), (dat[,SE_col])^2))
    
    C = sigma_hat_vec_to_mat(sigma_hat)
    eigen_vals = eigen(C)$values
    lambda_min = min(eigen_vals)
    
    if(lambda_min <= 0){
      diag(C) = diag(C) - lambda_min + 1e-6
    }
    
    x_vect_star = c(1, x_vect)
    
    TE = beta_hat%*%x_vect_star
    se = sqrt(t(x_vect_star) %*% C %*% x_vect_star)
    
    TE_orig = dat[, TE_col]
    TE_pred = cbind(1, X)%*%beta_hat
    
    NMI_out = data.frame(Study = study, 
                         Trt1 = unique(dat[,AgD_Trt_cols[1]]),
                         Trt2 = unique(dat[,AgD_Trt_cols[2]]),
                         x = rbind(x_vect),
                         TE = TE, se = se) 
    
    colnames(NMI_out) = gsub('\\.', '', colnames(NMI_out))
    
    Diag_out = data.frame(Study = study,
                          X,
                          TE_orig = TE_orig,
                          TE_pred = TE_pred)
    
    list(NMI_out = NMI_out, Diag_out = Diag_out)
  }
  
  out = lapply(studies, single_study_interpolation)
  
  Final = do.call(rbind, lapply(out, `[[`, 1))
  Diagnostics = do.call(rbind, lapply(out, `[[`, 2))
  
  return(list(Imputed = imputed, Final = Final, Diagnostics = Diagnostics))
}





#' Fixed effect NMA BUGS script for numeric outcomes in a treatment-difference 
#' format
#'
#' Written to a temporary location as a '.bug' file and run from 
#' \code{\link{NMA_run}}.
#'
#' @return None
#'
#' @export
FE_NMA_Normal_Trt_Diffs <- function(){
  # Fixed effects model for multi-arm trials
  for(i in 1:ns){ # LOOP THROUGH STUDIES
    var[i] <- pow(se[i], 2) # calculate variances
    prec[i] <- 1/var[i] # set precisions
    y[i] ~ dnorm(theta[i], prec[i]) # normal likelihood
    theta[i] <- d[t[i,2]] - d[t[i,1]] # model for linear predictor
  }
  # Normal likelihood, identity link
  d[1] <- 0 # treatment effect is zero for reference treatment
  for (k in 2:nt){d[k] ~ dnorm(0, .001)} # vague priors for treatment effects
  #Output
  # pairwise treatment effect for all possible pair-wise comparisons, if nt>2
  for (c in 1:(nt - 1)){
    for (k in (c+1):nt){
      D[c,k] <- (d[k] - d[c])
    }
  }
}





#' Fixed effect NMR BUGS script for continuous outcomes with two effect 
#' modifiers in a treatment-difference format.
#'
#' Written to a temporary location as a '.bug' file and run from 
#' \code{\link{NMA_Meta_Reg_run_2D}}.
#'
#' @return None
#'
#' @export
FE_NMA_Normal_Binary_Covariate_Trt_Diffs_2D <- function(){
  # Normal likelihood, identity link
  # Fixed effects model for 2-arm trials
  for(i in 1:ns){ # LOOP THROUGH STUDIES
    var[i] <- pow(se[i],2) # calculate variances
    prec[i] <- 1/var[i] # set precisions
    y[i]~dnorm(theta[i],prec[i]) # normal likelihood
    theta[i] <- d[t[i,2]] - d[t[i,1]] + (beta1[t[i,2]] - beta1[t[i,1]])*x1[i] + (beta2[t[i,2]] - beta2[t[i,1]])*x2[i]# model for linear predictor
  }
  d[1] <- 0 # treatment effect is zero for reference treatment
  beta1[1] <- 0 # covariate effect is zero for reference treatment
  beta2[1] <- 0 # covariate effect is zero for reference treatment
  
  for (k in 2:nt){
    d[k] ~ dnorm(0,.001) # vague priors for treatment effects
    beta1[k] <- B1
    beta2[k] <- B2
  }
  B1 ~ dnorm(0,.001)
  B2 ~ dnorm(0,.001)
  #Output
  # pairwise treatment effect for all possible pair-wise comparisons, if nt>2
  for(c in 1:(nt-1)){
    for(k in (c+1):nt){
      D[c,k] <- (d[k] - d[c])
    }
  }
}






#' Fixed effect NMR BUGS script for continuous outcomes with three effect 
#' modifiers in a treatment-difference format.
#'
#' Written to a temporary location as a '.bug' file and run from 
#' \code{\link{NMA_Meta_Reg_run_3D}}.
#'
#' @return None
#'
#' @export
FE_NMA_Normal_Binary_Covariate_Trt_Diffs_3D <- function(){
  # Normal likelihood, identity link
  # Fixed effects model for 2-arm trials
  for(i in 1:ns){ # LOOP THROUGH STUDIES
    var[i] <- pow(se[i],2) # calculate variances
    prec[i] <- 1/var[i] # set precisions
    y[i]~dnorm(theta[i],prec[i]) # normal likelihood
    theta[i] <- d[t[i,2]] - d[t[i,1]] + (beta1[t[i,2]] - beta1[t[i,1]])*x1[i] + (beta2[t[i,2]] - beta2[t[i,1]])*x2[i]+ (beta3[t[i,2]] - beta3[t[i,1]])*x3[i]# model for linear predictor
  }
  d[1] <- 0 # treatment effect is zero for reference treatment
  beta1[1] <- 0 # covariate effect is zero for reference treatment
  beta2[1] <- 0 # covariate effect is zero for reference treatment
  beta3[1] <- 0 # covariate effect is zero for reference treatment
  
  for (k in 2:nt){
    d[k] ~ dnorm(0,.001) # vague priors for treatment effects
    beta1[k] <- B1
    beta2[k] <- B2
    beta3[k] <- B3
  }
  B1 ~ dnorm(0,.001)
  B2 ~ dnorm(0,.001)
  B3 ~ dnorm(0,.001)
  #Output
  # pairwise treatment effect for all possible pair-wise comparisons, if nt>2
  for(c in 1:(nt-1)){
    for(k in (c+1):nt){
      D[c,k] <- (d[k] - d[c])
    }
  }
}





#' Running NMA in a treatment-difference format
#'
#' @param dat aggregate level data.
#' @param N_chains number of Markov chains.
#' @param N_iter number of total iterations per chain (including burn in).
#' @param burnin length of burn in, i.e. number of iterations to discard at the 
#' beginning.
#' 
#' @return A \code{\link[R2OpenBUGS]{bugs}} object.
#'
#' @export
NMA_run = function(dat, N_chains, N_iter, burnin){
  dat = dat %>% 
    as_tibble %>% 
    mutate(across(c(Trt1, Trt2), 
                  ~ dplyr::recode(., 'A' = '1', 'B' = '2', 
                                  'C' = '3', 'D' = '4'))
    )
  
  ns = nrow(dat)
  
  se = as.numeric(dat$se)
  y = as.numeric(dat$TE)
  nt = length(unique(dat$Trt2))+1
  t = apply(dat[,2:3], 2, as.numeric) # shift one to right.
  
  data = list(ns=ns, nt=nt, t=t,  y=y, se=se)
  
  d = c(NA, rep(0, nt-1))
  params = c("d", "D")
  inits = function(){list(d=as.vector(d))}
  
  mod = FE_NMA_Normal_Trt_Diffs
  filename = file.path(tempdir(), "FE_NMA_Normal_Trt_Diffs.txt")
  write.model(mod, filename)
  
  sim = R2OpenBUGS::bugs(data, inits,
                         model.file = filename,
                         parameters = params,
                         n.chains = N_chains, n.iter = N_iter,
                         n.burnin = burnin, n.thin = 1)
  
  return(sim)
}




#' Running NMR with two effect modifiers in a treatment-difference format
#'
#' @param dat aggregate level data.
#' @param N_chains number of Markov chains.
#' @param N_iter number of total iterations per chain (including burn in).
#' @param burnin length of burn in, i.e. number of iterations to discard at the 
#' beginning.
#' 
#' @return A \code{\link[R2OpenBUGS]{bugs}} object.
#'
#' @export
NMA_Meta_Reg_run_2D = function(dat, N_chains, N_iter, burnin){
  dat = dat %>% 
    as_tibble %>% 
    mutate(across(c(Trt1, Trt2), 
                  ~ dplyr::recode(., 'A' = '1', 'B' = '2', 
                                  'C' = '3', 'D' = '4'))
    )
  
  ns = nrow(dat)
  
  se = as.numeric(dat$se)
  y = as.numeric(dat$TE)
  nt = length(unique(dat$Trt2))+1
  t = apply(dat[,2:3], 2, as.numeric)# shift one to right.
  x1 = as.numeric(dat$x1)
  x2 = as.numeric(dat$x2)
  
  data = list(ns = ns, nt = nt, t = t, y = y, se = se, 
              x1 = x1, x2 = x2)
  
  d=c(NA, rep(0, nt-1))
  params = c("d", "D", "B1", "B2")
  inits = function(){
    list(d=as.vector(d), B1=0, B2=0)
  }
  
  mod = FE_NMA_Normal_Binary_Covariate_Trt_Diffs_2D
  filename = file.path(tempdir(), 
                       "FE_NMA_Normal_Binary_Covariate_Trt_Diffs_2D.txt")
  write.model(mod, filename)
  
  sim = R2OpenBUGS::bugs(data, inits,
                         model.file = filename,
                         parameters = params,
                         n.chains = N_chains, n.iter = N_iter,
                         n.burnin = burnin, n.thin = 1)
  
  return(sim)
}



#' Running NMR with three effect modifiers in a treatment-difference format
#'
#' @param dat aggregate level data.
#' @param N_chains number of Markov chains.
#' @param N_iter number of total iterations per chain (including burn in).
#' @param burnin length of burn in, i.e. number of iterations to discard at the 
#' beginning.
#' 
#' @return A \code{\link[R2OpenBUGS]{bugs}} object.
#'
#' @export
NMA_Meta_Reg_run_3D = function(dat, N_chains, N_iter, burnin){
  dat = dat %>% 
    as_tibble %>% 
    mutate(across(c(Trt1, Trt2), 
                  ~ dplyr::recode(., 'A' = '1', 'B' = '2', 
                                  'C' = '3', 'D' = '4'))
    )
  
  ns = nrow(dat)
  
  se = as.numeric(dat$se)
  y = as.numeric(dat$TE)
  nt = length(unique(dat$Trt2))+1
  t = apply(dat[,2:3], 2, as.numeric)# shift one to right.
  x1 = as.numeric(dat$x1)
  x2 = as.numeric(dat$x2)
  x3 = as.numeric(dat$x3)
  
  data = list(ns = ns, nt = nt, t = t, y = y, se = se, 
              x1 = x1, x2 = x2, x3 = x3)
  
  d=c(NA, rep(0, nt-1))
  params = c("d", "D", "B1", "B2", "B3")
  inits = function(){
    list(d=as.vector(d), B1=0, B2=0, B3=0)
  }
  
  mod = FE_NMA_Normal_Binary_Covariate_Trt_Diffs_3D
  filename = file.path(tempdir(), 
                       "FE_NMA_Normal_Binary_Covariate_Trt_Diffs_3D.txt")
  write.model(mod, filename)
  
  sim = R2OpenBUGS::bugs(data, inits,
                         model.file = filename,
                         parameters = params,
                         n.chains = N_chains, n.iter = N_iter,
                         n.burnin = burnin, n.thin = 1)
  
  return(sim)
}





#' Running ML-NMR for a dichotomous outcome with two effect modifiers
#'
#' @param ML_NMR_data a list consisting of the following elements - 
#' \item{AgD}{Aggregate-level data.}
#' \item{IPD}{Individual patient-level data}
#' @param N_iter number of total iterations per chain (including burn in).
#' @param N_chains number of Markov chains.
#' @param burnin length of burn in, i.e. number of iterations to discard at the 
#' beginning.
#' @param n_int number of Monte Carlo integration points
#' 
#' @return A \code{\link[multinma]{stan_nma}} object.
#'
#' @export
ML_NMR_Run_2D = function(ML_NMR_data, N_iter, N_chains, burnin, n_int){
  net = combine_network(
    set_ipd(ML_NMR_data$IPD %>%
              mutate(TrC = Tr),
            study = Study,
            trt = Tr,
            trt_class = TrtClass,
            r = Y),
    set_agd_arm(ML_NMR_data$AgD %>%
                  mutate(TrC = Tr),
                study = Study,
                trt = Tr,
                trt_class = TrtClass,
                r = r,
                n = n),
    trt_ref = "A")
  
  net = add_integration(net,
                        x1 = distr(qbern, prob = x1),
                        x2 = distr(qbern, prob = x2),
                        n_int = n_int
  )
  
  fit.FE = nma(net, 
               trt_effects = "fixed",
               link = "logit", 
               likelihood = "bernoulli2",
               regression = ~(x1 + x2)*.trt,
               center = FALSE,
               # class_interactions = "common",
               prior_intercept = normal(scale = 10),
               prior_trt = normal(scale = 10),
               prior_reg = normal(scale = 10),
               init_r = 0.1,
               QR = TRUE,
               iter = N_iter/N_chains,
               warmup = burnin/N_chains,
               chains = N_chains)
  
  return(fit.FE)
}






#' Running ML-NMR for a dichotomous outcome with three effect modifiers
#'
#' @param ML_NMR_data a list consisting of the following elements - 
#' \item{AgD}{Aggregate-level data.}
#' \item{IPD}{Individual patient-level data}
#' @param N_iter number of total iterations per chain (including burn in).
#' @param N_chains number of Markov chains.
#' @param burnin length of burn in, i.e. number of iterations to discard at the 
#' beginning.
#' @param n_int number of Monte Carlo integration points
#' 
#' @return A \code{\link[multinma]{stan_nma}} object.
#'
#' @export
ML_NMR_Run_3D = function(ML_NMR_data, N_iter, N_chains, burnin, n_int){
  net = combine_network(
    set_ipd(ML_NMR_data$IPD %>%
              mutate(TrC = Tr),
            study = Study,
            trt = Tr,
            trt_class = TrtClass,
            r = Y),
    set_agd_arm(ML_NMR_data$AgD %>%
                  mutate(TrC = Tr),
                study = Study,
                trt = Tr,
                trt_class = TrtClass,
                r = r,
                n = n),
    trt_ref = "A")
  
  net = add_integration(net,
                        x1 = distr(qbern, prob = x1),
                        x2 = distr(qbern, prob = x2),
                        x3 = distr(qbern, prob = x3),
                        n_int = n_int
  )
  
  fit.FE = nma(net, 
               trt_effects = "fixed",
               link = "logit", 
               likelihood = "bernoulli2",
               regression = ~(x1 + x2 + x3)*.trt,
               center = FALSE,
               # class_interactions = "common",
               prior_intercept = normal(scale = 10),
               prior_trt = normal(scale = 10),
               prior_reg = normal(scale = 10),
               init_r = 0.1,
               QR = TRUE,
               iter = N_iter/N_chains,
               warmup = burnin/N_chains,
               chains = N_chains)
  
  return(fit.FE)
}




#' Estimating treatment effects from an NMA run
#'
#' @param sim An \code{\link[R2OpenBUGS]{bugs}} object. The output of 
#' \code{\link{NMA_run}}.
#' 
#' @return A data frame of posterior summaries for the treatment effects.
#'
#' @export
NMA_NMI_summary = function(sim){
  cols = grepl('D\\[', colnames(sim$sims.matrix))
  M = sim$sims.matrix[, cols]
  
  out = as.data.frame(
    t(apply(M, 2, function(x){
      c(mean(x), quantile(x, c(.5, .025, .975)))
    })))
  
  out$Parameter = row.names(out)
  names(out)[1] = 'Mean'  
  row.names(out) = c()
  
  return(out)
}





#' Estimating treatment effects from an NMR run with two effect modifiers
#'
#' @param sim An \code{\link[R2OpenBUGS]{bugs}} object. The output of 
#' \code{\link{NMA_Meta_Reg_run_2D}}.
#' @param x_vect a vector of length two, consisting of the effect modifier 
#' levels at which the ITC will be conducted.
#' 
#' @return A data frame of posterior summaries for the treatment effects at x1 
#' and x2.
#'
#' @export
NMA_Metareg_summary_2D = function(sim, x_vect){
  x1 = x_vect[1]
  x2 = x_vect[2]
  x3 = x_vect[3]
  
  cols = grepl('D\\[|B1|B2', colnames(sim$sims.matrix))
  M = sim$sims.matrix[, cols]
  M[,1:3] = apply( M[,1:3], 2,
                   function(col){
                     col + M[,"B1"]*x1 + M[,"B2"]*x2
                   })
  
  out = as.data.frame(
    t(apply(M[,1:(ncol(M)-2)], 2, function(x){
      c(mean(x), quantile(x, c(.5, .025, .975)))
    })))
  
  out$Parameter = row.names(out)
  names(out)[1] = 'Mean' 
  row.names(out) = c()
  
  return(out)
}






#' Estimating treatment effects from an NMR run with three effect modifiers
#'
#' @param sim An \code{\link[R2OpenBUGS]{bugs}} object. The output of 
#' \code{\link{NMA_Meta_Reg_run}}.
#' @param x_vect a vector of length two, consisting of the effect modifier 
#' levels at which the ITC will be conducted.
#' 
#' @return A data frame of posterior summaries for the treatment effects at x1 
#' and x2.
#'
#' @export
NMA_Metareg_summary_3D = function(sim, x_vect){
  x1 = x_vect[1]
  x2 = x_vect[2]
  x3 = x_vect[3]
  
  cols = grepl('D\\[|B1|B2|B3', colnames(sim$sims.matrix))
  M = sim$sims.matrix[, cols]
  M[,1:3] = apply( M[,1:3], 2,
                   function(col){
                     col + M[,"B1"]*x1 + M[,"B2"]*x2 + M[,"B3"]*x3
                   })
  
  out = as.data.frame(
    t(apply(M[,1:(ncol(M)-3)], 2, function(x){
      c(mean(x), quantile(x, c(.5, .025, .975)))
    })))
  
  out$Parameter = row.names(out)
  names(out)[1] = 'Mean' 
  row.names(out) = c()
  
  return(out)
}



#' Estimating treatment effects from an ML-NMR run with two effect modifiers
#'
#' @param n_trts Number of distinct treatments in the network.
#' @param ML_NMR_Fit An \code{\link[multinma]{stan_nma}} object. The output of 
#' \code{\link{ML_NMR_Run}}.
#' @param x_vect a vector of length two, consisting of the effect modifier 
#' levels at which the ITC will be conducted.
#' 
#' @return A data frame of posterior summaries for the treatment effects at 
#' x_vect.
#'
#' @export
ML_NMR_summary_2D = function(n_trts, ML_NMR_Fit, x_vect){
  inds = names(ML_NMR_Fit$stanfit)[grep("beta|^d\\[", 
                                        names(ML_NMR_Fit$stanfit))]
  inds1 = grep("trt", inds)
  inds2 = grep("d", inds)
  inds = c(inds1, inds2)
  
  stan_out = rstan::extract(ML_NMR_Fit$stanfit)
  df = as.data.frame(cbind(stan_out$beta, stan_out$d))[,inds]
  
  df1 = df[,inds1 - 2]
  df2 = df[,inds2 - 2]
  
  x1 = x_vect[1]
  x2 = x_vect[2]
  
  X = rbind(x1*diag(n_trts - 1), x2*diag(n_trts - 1))
  post_samp = as.matrix(df2) + 
    apply(as.matrix(df1)%*%c(x1, x2), 1, sum)
  colnames(post_samp) = paste0("A", LETTERS[2:n_trts])
  
  for(i in 1:(n_trts-2)){
    for(j in (i+1):(n_trts-1)){
      temp = as.matrix(post_samp[,j] - post_samp[,i], ncol = 1)
      colnames(temp) = paste0(LETTERS[2:n_trts][i],
                              LETTERS[2:n_trts][j])
      post_samp = cbind(post_samp, temp)
    }
  }
  
  estimates = t(apply(post_samp, 2, 
                      function(x){
                        c(mean(x), quantile(x, c(.5, .025, .975)))
                      }
  ))
  
  estimates = as.data.frame(cbind(estimates))
  
  temp = t(combn(1:4, 2))
  estimates$Parameter = paste0("D[", temp[,1], ",", temp[,2], "]")
  
  names(estimates)[1] = "Mean"
  row.names(estimates) = c()
  
  return(estimates)
}



#' Estimating treatment effects from an ML-NMR run with three effect modifiers
#'
#' @param n_trts Number of distinct treatments in the network.
#' @param ML_NMR_Fit An \code{\link[multinma]{stan_nma}} object. The output of 
#' \code{\link{ML_NMR_Run}}.
#' @param x_vect a vector of length two, consisting of the effect modifier 
#' levels at which the ITC will be conducted.
#' 
#' @return A data frame of posterior summaries for the treatment effects at 
#' x_vect.
#'
#' @export
ML_NMR_summary_3D = function(n_trts, ML_NMR_Fit, x_vect){
  inds = names(ML_NMR_Fit$stanfit)[grep("beta|^d\\[", 
                                        names(ML_NMR_Fit$stanfit))]
  inds1 = grep("trt", inds)
  inds2 = grep("d", inds)
  inds = c(inds1, inds2)
  
  stan_out = rstan::extract(ML_NMR_Fit$stanfit)
  df = as.data.frame(cbind(stan_out$beta, stan_out$d))[,inds]
  
  df1 = df[,inds1 - 3]
  df2 = df[,inds2 - 3]
  
  x1 = x_vect[1]
  x2 = x_vect[2]
  x3 = x_vect[3]
  
  X = rbind(x1*diag(n_trts - 1), x2*diag(n_trts - 1),  x3*diag(n_trts - 1))
  post_samp = as.matrix(df2) + 
    apply(as.matrix(df1)%*%c(x1, x2, x3), 1, sum)
  colnames(post_samp) = paste0("A", LETTERS[2:n_trts])
  
  for(i in 1:(n_trts-2)){
    for(j in (i+1):(n_trts-1)){
      temp = as.matrix(post_samp[,j] - post_samp[,i], ncol = 1)
      colnames(temp) = paste0(LETTERS[2:n_trts][i],
                              LETTERS[2:n_trts][j])
      post_samp = cbind(post_samp, temp)
    }
  }
  
  estimates = t(apply(post_samp, 2, 
                      function(x){
                        c(mean(x), quantile(x, c(.5, .025, .975)))
                      }
  ))
  
  estimates = as.data.frame(cbind(estimates))
  
  temp = t(combn(1:4, 2))
  estimates$Parameter = paste0("D[", temp[,1], ",", temp[,2], "]")
  
  names(estimates)[1] = "Mean"
  row.names(estimates) = c()
  
  return(estimates)
}





#' Treatment effects and 95% CrIs
#'
#' @param summary_tab A data frame. The output of \code{\link{NMA_NMI_summary}},
#' \code{\link{NMA_Metareg_summary}}, or \code{\link{ML_NMR_summary}}.
#' 
#' @return A single row data frame containing posterior medians and 95% CrIs for
#' all treatment effects.
#'
#' @export
result_table = function(summary_tab){
  out = summary_tab %>% 
    mutate(Est = paste0(
      sprintf('%.2f', `50%`), ' [',
      sprintf('%.2f', `2.5%`), ' ;',
      sprintf('%.2f', `97.5%`), ']')
    ) %>% 
    mutate(Parameter = dplyr::recode(Parameter,
                                     "D[1,2]" = "$d_{\\mathrm{AB}}$",
                                     "D[1,3]" = "$d_{\\mathrm{AC}}$",
                                     "D[1,4]" = "$d_{\\mathrm{AD}}$",
                                     "D[2,3]" = "$d_{\\mathrm{BC}}$",
                                     "D[2,4]" = "$d_{\\mathrm{BD}}$",
                                     "D[3,4]" = "$d_{\\mathrm{CD}}$"
    )) %>% 
    select(Parameter, Est) %>% 
    dcast(NULL ~ Parameter, value.var = 'Est') %>% 
    select(-1)
}




#' Graphical display or results
#'
#' @param NMA_summary The output of \code{\link{NMA_NMI_summary}} for an NMA run.
#' @param NMR_summary The output of \code{\link{NMA_Metareg_summary}}.
#' @param ML_NMR_summ The output of \code{\link{ML_NMR_summary}}.
#' @param NMI_summary The output of \code{\link{NMA_NMI_summary}} for an NMI run.
#' @param trt_effs A vector of the true treatment effect used to generate the 
#' data. Default is NULL.
#' 
#' @return A ggplot2 object.
#'
#' @export
result_forest_plot = function(NMA_summary, NMR_summary, 
                              ML_NMR_summ, NMI_summary,
                              trt_effs = NULL){
  df = rbind(NMA_summary, NMR_summary, ML_NMR_summ, NMI_summary)
  if(!is.null(trt_effs)) df = cbind(df, trt_eff = as.numeric(rep(trt_effs, 4)))
  
  df = df %>% 
    mutate(Parameter = dplyr::recode(Parameter,
                                     "D[1,2]" = "hat(d)[AB]",
                                     "D[1,3]" = "hat(d)[AC]",
                                     "D[1,4]" = "hat(d)[AD]",
                                     "D[2,3]" = "hat(d)[BC]",
                                     "D[2,4]" = "hat(d)[BD]",
                                     "D[3,4]" = "hat(d)[CD]"
    ))
  
  df$Method = factor(rep(c('NMA', 'NMR', 'ML-NMR', 'NMI'), each = 6), 
                     levels = c('NMA', 'NMR', 'ML-NMR', 'NMI'))
  
  p = ggplot(df, aes(Method, `50%`)) + 
    geom_segment(mapping = aes(x = Method, xend = Method,
                               y = `2.5%`, yend = `97.5%`),
                 col = 'royal blue', size = 1) + 
    geom_point(col = 'red', size = 2) 
  
  if(!is.null(trt_effs)){
    p = p + 
      geom_hline(df, mapping = aes(yintercept = trt_eff),
                 col = 'dark green', linetype = 'dashed', size = 1)
  }
  
  p = p + 
    facet_rep_wrap(.~ Parameter, ncol = 1,
                   strip.position = 'left',
                   labeller = label_parsed,
                   repeat.tick.labels = 'all') + 
    coord_flip() + 
    theme(panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_line(colour = "grey92"),
          axis.title.x = element_text(size=12, face = "bold"),
          axis.title.y = element_text(size=12, face = "bold"),
          axis.text.x = element_text(size=10),
          axis.text.y = element_text(size=10),
          strip.text = element_text(size = 11, face = 'bold'),
          strip.text.y.left = element_text(angle=0)) + 
    ylab("Parameter estimate (95% CrI)") 
  
  return(p)
}




#' Goodness of interpolation interactive plot
#'
#' @param NMI_object The output of \code{\link{NMI_interpolation}}
#' 
#' @return An interactive plot displaying the predicted TE estimates from the 
#' NMI algorithm vs. the observed TE estimates.
#'
#' @export
NMI_diagnostic_plotly = function(NMI_object){
  df = NMI_object$Diagnostics
  df$Text = apply(df, 1,
                  function(u){
                    v = u[2:(length(u)-2)]
                    paste0('Study ', u[1], '\n',
                           paste0(names(v), ' = ', 
                                  sprintf('%.1f', 100*v), '%', '\n', 
                                  collapse="")
                    )
                  })
  
  p = ggplot(df, aes(TE_orig, TE_pred, text = Text,
                     fill = as.factor(Study))) + 
    geom_segment(x = min(min(df$TE_orig), min(df$TE_pred)), 
                 y = min(min(df$TE_orig), min(df$TE_pred)),
                 xend = max(max(df$TE_orig), max(df$TE_pred)),
                 yend = max(max(df$TE_orig), max(df$TE_pred)),
                 linetype = "dashed", color = "red") + 
    geom_point(size = 3, shape = 21, col = 'royal blue') + 
    theme(panel.background = element_blank(),
          legend.position = 'none',
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_line(colour = "grey92"),
          axis.title.x = element_text(size=13, face = "bold"),
          axis.title.y = element_text(size=13, face = "bold"),
          axis.text.x = element_text(size=10),
          axis.text.y = element_text(size=10),
          plot.title = element_text(size=14, face = "bold", hjust = .5)) + 
    xlab('Observed treatment effect estimate') + 
    ylab('Predicted treatment effect estimate') + 
    ggtitle('NMI goodness of fit') + 
    scale_fill_brewer(palette = "Set1")
  
  ggplotly(p, tooltip = "text")
}




#' Goodness of interpolation plot
#'
#' @param NMI_object The output of \code{\link{NMI_interpolation}}
#' 
#' @return An ggplot displaying the predicted TE estimates from the 
#' NMI algorithm vs. the observed TE estimates.
#'
#' @export
NMI_diagnostic_plot = function(NMI_object){
  df = NMI_object$Diagnostics
  df$Study = factor(df$Study, levels = sort(unique(df$Study)))
  
  p = ggplot(df, aes(TE_orig, TE_pred, fill = Study)) + 
    geom_segment(x = min(min(df$TE_orig), min(df$TE_pred)), 
                 y = min(min(df$TE_orig), min(df$TE_pred)),
                 xend = max(max(df$TE_orig), max(df$TE_pred)),
                 yend = max(max(df$TE_orig), max(df$TE_pred)),
                 linetype = "dashed",
                 color = "red") + 
    geom_point(col = 'royal blue', size = 4,
               shape = 21, stroke = 1) + 
    theme(panel.background = element_blank(),
          legend.position = 'top',
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_line(colour = "grey92"),
          axis.title.x = element_text(size=13, face = "bold"),
          axis.title.y = element_text(size=13, face = "bold"),
          axis.text.x = element_text(size=10),
          axis.text.y = element_text(size=10),
          plot.title = element_text(size=14, face = "bold", hjust = .5)) + 
    xlab('Observed treatment effect estimate') + 
    ylab('Predicted treatment effect estimate') + 
    ggtitle('NMI goodness of fit') + 
    scale_fill_brewer(palette = "Set1")
  
  return(p)
}


#' Graphical display or results
#'
#' @param NMA_results The output of \code{\link{result_table}} for an NMA run.
#' @param NMR_results The output of \code{\link{result_table}} for an NMR run.
#' @param ML_NMR_results The output of \code{\link{result_table}} for an ML-NMR 
#' run.
#' @param NMI_results The output of \code{\link{result_table}} for an NMI run.
#' 
#' @return A kableExtra table.
#'
#' @export
display_result_table = function(NMA_results, NMR_results, 
                                ML_NMR_results, NMI_results){
  rbind(NMA_results, NMR_results, ML_NMR_results, NMI_results) %>%
    add_column(Method = c('NMA', 'NMR', 'ML-NMR', 'NMI'),
               .before = '$d_{\\mathrm{AB}}$') %>% 
    kableExtra::kable(align = 'c') %>%
    kableExtra::kable_styling(full_width = F, font_size = 14,
                              bootstrap_options = c("striped")) %>%
    kableExtra::row_spec(0, font_size = 14)
}
