library(R2OpenBUGS)
library(dplyr)
library(doParallel)
library(doSNOW)
library(tictoc)
library(tibble)
library(multinma)
library(tidyr)
library(ggplot2)
library(broom)
library(broom.mixed)
library(reshape2)
library(lemon)
library(kableExtra)
library(purrr)
library(plotly)

#Setting working directory to current one
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Sourcing functions
source('NMI_Functions.R')

#*******************************************************************************
#*******************************************************************************
#* 2 effect modifiers example
#*******************************************************************************
#*******************************************************************************

#*******************************************************************************
#* Inputs for this script
#*******************************************************************************

#value of effect modifier x1 to be used in ITC
x1 = .675
#value of effect modifier x2 to be used in ITC
x2 = .475

x_vect = c(x1, x2)

N_chains = 3 #number of MCMC chains
N_iter = 1500 #number of MCMC iterations (in total)
burnin = 510 #number of MCMC samples to discard
n_int = 500 #number of MC integration points for ML-NMR

#AgD and IPD effect modifier column names
AgD_EM_cols = IPD_EM_cols = c('x1', 'x2')

Study_col = 'Study' #study column name
samp_sizes = rep(600, 6) #sample sizes for AgD studies
Trt_cols = c('Trt1', 'Trt2') #AgD treatment column names
TE_col = 'TE' #AgD treatment effect estimate column name
SE_col = 'se' #AgD standard error column name
IPD_Trt_col = 'Tr' #IPD treatment column name
outcome_col = 'Y' #IPD outcome column name
#*******************************************************************************


#*******************************************************************************
#* Reading all datasets
#*******************************************************************************
#reading the data
IPD = read.csv('Example_IPD.csv') #IPD (for all methods)

NMR_AgD = read.csv('Example_AgD_NMR_NMA.csv') #AgD for NMR/NMA
NMI_AgD = read.csv('Example_AgD_NMI.csv') #AgD for NMI
ML_NMR_AgD = read.csv('Example_AgD_ML_NMR.csv')  #AgD for ML-NMR

#estimating the Pearson correlation between x1 and x2
(rho_hat = cor(IPD[grep('x', colnames(IPD))]))
#*******************************************************************************


#*******************************************************************************
#* The network
#*******************************************************************************
ML_NMR_data = list(IPD = IPD, AgD = ML_NMR_AgD)

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


par(mar = c(0,0,0,0))
plot(net, weight_edges = TRUE, weight_nodes = TRUE)
#*******************************************************************************


#*******************************************************************************
#** NMI run
#*******************************************************************************
#Imputing the AgD
NMI_object = NMI_interpolation(IPD, NMI_AgD, x_vect, AgD_EM_cols, IPD_EM_cols, 
                               Study_col, samp_sizes, Trt_cols, TE_col, SE_col, 
                               IPD_Trt_col) 
#Have a look at the imputed dataset
NMI_object$Imputed 
#The data submitted to NMA for ITC
NMI_object$Final 

#NMI interpolation goodness of fit
NMI_diagnostic_plotly(NMI_object)

#Model fitting
NMI_sim = NMA_run(NMI_object$Final, N_chains, N_iter, burnin) 
#Summarising results
NMI_summary = NMA_NMI_summary(NMI_sim) 
#TEs and CrIs
NMI_results = result_table(NMI_summary) 
#*******************************************************************************


#*******************************************************************************
#** NMA run
#*******************************************************************************
#*#Model fitting
NMA_sim = NMA_run(NMR_AgD, N_chains, N_iter, burnin) 
#Summarising results
NMA_summary = NMA_NMI_summary(NMA_sim) 
#TEs and CrIs
NMA_results = result_table(NMA_summary) 
#*******************************************************************************


#*******************************************************************************
#** NMR run
#*******************************************************************************
#*#Model fitting
NMR_sim = NMA_Meta_Reg_run_2D(NMR_AgD, N_chains, N_iter, burnin) 
#Summarising results
NMR_summary = NMA_Metareg_summary_2D(NMR_sim, x_vect)
#TEs and CrIs
NMR_results = result_table(NMR_summary) 
#*******************************************************************************


#*******************************************************************************
#** ML-NMR run
#*******************************************************************************
#Model fitting
ML_NMR_sim = ML_NMR_Run_2D(ML_NMR_data, N_iter, N_chains, burnin, n_int) 
#summarising results
ML_NMR_summ = ML_NMR_summary_2D(n_trts=4, ML_NMR_sim, x_vect) 
#TEs and CrIs
ML_NMR_results = result_table(ML_NMR_summ) 
#*******************************************************************************


#*******************************************************************************
#Displaying results in a table
display_result_table(NMA_results, NMR_results, ML_NMR_results, NMI_results)


#These are the logistic regression coefficients that were used to generate the 
#data used in this demonstration 
beta_AB = c(-1.39, 0, 0, 0.69, 1.00, 0.00)
beta_AC = c(-1.39, 0, 0, 1.00, 1.61, 1.00)
beta_AD = c(-1.39, 0, 0, 1.50, -1.20, -1.00)

#converting the logistic regression coefficients into treatment effects
d = rbind(beta_AB, beta_AC, beta_AD)[,4:6]
trt_effs = d%*%c(1, x1, x2)
trt_effs = c(trt_effs, trt_effs[c(2,3,3)] - trt_effs[c(1,1,2)])

#Plotting a forest plot (dashed green lines are true TEs)
result_forest_plot(NMA_summary, NMR_summary, ML_NMR_summ, NMI_summary, trt_effs)
#*******************************************************************************





#*******************************************************************************
#*******************************************************************************
#* 3 effect modifiers example
#*******************************************************************************
#*******************************************************************************


#*******************************************************************************
#* Inputs for this script
#*******************************************************************************

#value of effect modifier x1 to be used in ITC
x1 = .675
#value of effect modifier x2 to be used in ITC
x2 = .475
#value of effect modifier x3 to be used in ITC
x3 = .4

x_vect = c(x1, x2, x3)


N_chains = 3 #number of MCMC chains
N_iter = 1500 #number of MCMC iterations (in total)
burnin = 510 #number of MCMC samples to discard
n_int = 500 #number of MC integration points for ML-NMR

#AgD and IPD effect modifier column names
AgD_EM_cols = IPD_EM_cols = c('x1', 'x2', 'x3')

Study_col = 'Study' #study column name
samp_sizes = rep(600, 6) #sample sizes for AgD studies
Trt_cols = c('Trt1', 'Trt2') #AgD treatment column names
TE_col = 'TE' #AgD treatment effect estimate column name
SE_col = 'se' #AgD standard error column name
IPD_Trt_col = 'Tr' #IPD treatment column name
outcome_col = 'Y' #IPD outcome column name
#*******************************************************************************


#*******************************************************************************
#* Reading all datasets
#*******************************************************************************
#reading the data
IPD = read.csv('Example_IPD_3D.csv') #IPD (for all methods)

NMR_AgD = read.csv('Example_AgD_NMR_NMA_3D.csv') #AgD for NMR/NMA
NMI_AgD = read.csv('Example_AgD_NMI_3D.csv') #AgD for NMI
ML_NMR_AgD = read.csv('Example_AgD_ML_NMR_3D.csv')  #AgD for ML-NMR

#estimating the Pearson correlation between x1, x2 and x3
(rho_hat = cor(IPD[grep('x', colnames(IPD))]))
#*******************************************************************************


#*******************************************************************************
#* The network
#*******************************************************************************
ML_NMR_data = list(IPD = IPD, AgD = ML_NMR_AgD)

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


par(mar = c(0,0,0,0))
plot(net, weight_edges = TRUE, weight_nodes = TRUE)
#*******************************************************************************


#*******************************************************************************
#** NMI run
#*******************************************************************************
#imputing AgD
NMI_object = NMI_interpolation(IPD, NMI_AgD, x_vect, AgD_EM_cols, IPD_EM_cols, 
                               Study_col, samp_sizes, Trt_cols, TE_col, SE_col, 
                               IPD_Trt_col) 
#Have a look at the imputed dataset
NMI_object$Imputed 
#The data submitted to NMA for ITC
NMI_object$Final

#Interpolation goodness of fit
NMI_diagnostic_plotly(NMI_object)

NMI_sim = NMA_run(NMI_object$Final, N_chains, N_iter, burnin) #Model fitting
NMI_summary = NMA_NMI_summary(NMI_sim) #summarising results
NMI_results = result_table(NMI_summary) #TEs and CrIs
#*******************************************************************************


#*******************************************************************************
#** NMA run
#*******************************************************************************
#Model fitting
NMA_sim = NMA_run(NMR_AgD, N_chains, N_iter, burnin) 
#summarising results
NMA_summary = NMA_NMI_summary(NMA_sim)
#TEs and CrIs
NMA_results = result_table(NMA_summary) 
#*******************************************************************************


#*******************************************************************************
#** NMR run
#*******************************************************************************
#Model fitting
NMR_sim = NMA_Meta_Reg_run_3D(NMR_AgD, N_chains, N_iter, burnin) 
#summarising results
NMR_summary = NMA_Metareg_summary_3D(NMR_sim, x_vect) 
#TEs and CrIs
NMR_results = result_table(NMR_summary) 
#*******************************************************************************


#*******************************************************************************
#** ML-NMR run
#*******************************************************************************
#Model fitting
ML_NMR_sim = ML_NMR_Run_3D(ML_NMR_data, N_iter, N_chains, burnin, n_int) 
#summarising results
ML_NMR_summ = ML_NMR_summary_3D(n_trts=4, ML_NMR_sim, x_vect) 
#TEs and CrIs
ML_NMR_results = result_table(ML_NMR_summ) 
#*******************************************************************************


#*******************************************************************************
#Displaying results in a table
display_result_table(NMA_results, NMR_results, ML_NMR_results, NMI_results)

#These are the logistic regression coefficients that were used to generate the 
#data used in this demonstration 
beta_AB = c(-1.39, 0, 0, 0 ,0.69, 1.00, 0.00, .25)
beta_AC = c(-1.39, 0, 0, 0, 1.00, 1.61, 1.00, 0)
beta_AD = c(-1.39, 0, 0, 0, 1.50, -1.20, -1.00, .25)

#converting the logistic regression coefficients into treatment effects
d = rbind(beta_AB, beta_AC, beta_AD)[,5:8]

trt_effs = d%*%c(1, x1, x2, x3)
trt_effs = c(trt_effs, trt_effs[c(2,3,3)] - trt_effs[c(1,1,2)])

#Plotting a forest plot (dashed green lines are true TEs)
result_forest_plot(NMA_summary, NMR_summary, ML_NMR_summ, NMI_summary, trt_effs)
#*******************************************************************************
