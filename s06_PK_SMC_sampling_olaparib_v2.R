rm(list=ls())
library(ospsuite)
library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(mvtnorm)
library("PerformanceAnalytics")

source("Functions/ABC_SMC_PKsim_fun_v3.R") # Session -> set working directory -> to source file location

### Objective function value calculations

# Calculates goodness-of-fit between simulated and observed plasma/urine concentration data
OFVcalc_fun3 <- function(sim_df, data, amt_urine_exp=0.1, amt_urine_rse=0.4) {
  if (nrow(sim_df)==0){
    return(10000) # Returns 10000 (high penalty) if simulation data is empty
  }
  data$PRED <- approx(x = sim_df$Time, y = sim_df$DV, xout = data$TIME)$y # match the simulated concentration to the exact time points of the observed data
  data$WRES <- (log(1 + data$PRED) - log(1+ data$DV))^2 # WRES - Weighted Residuals (the squared difference between the log-transformed predicted (PRED) and observed (DV) concentrations)
  return(sum(data$WRES)/(nrow(data)+1))
} # PRED - predicted data, DV - observed data 


# Description of the function:
# Extended version that separately handles urine compartment data ->
# Splits simulation data into non-urine and urine components ->
# Calculates plasma/urine concentration fit (same as above) ->
# Adds an additional penalty (5000) if cumulative urine amount doesn't match expected value within specified relative standard error ->
# Returns combined objective function value
OFVcalc_fun3urine <- function(sim_df, data, amt_urine_exp, amt_urine_rse=0.4){
  sim_df1 <- sim_df[!grepl(pattern="Urine", sim_df$compartment),]
  sim_df2 <- sim_df[grepl(pattern="Urine", sim_df$compartment),]
  if (nrow(sim_df1)==0 | nrow(sim_df2)==0 ){
    return(10000)
  }
  data$PRED <- approx(x = sim_df1$Time, y = sim_df1$DV, xout = data$TIME)$y
  data$WRES <- (log(1 + data$PRED) - log(1+ data$DV))^2
  ofv <- sum(data$WRES)/(nrow(data)+1)
  
  if (amt_urine_exp >0){
    max_t <- max(sim_df2$Time)
    amt_urine <-  sim_df2[sim_df2$Time== max_t,'DV']
    ofv <- ofv + ifelse(abs(amt_urine - amt_urine_exp)< (amt_urine_exp*amt_urine_rse), 0, 5000)
  }  
  return(ofv)
}


# Approximate Bayesian Computation Sequential Monte Carlo (ABC-SMC) analysis for pharmacokinetic (PK) modeling of Olaparib (cancer drug)
# Goal - to estimate key physiological parameters for the cancer drug Olaparib by fitting a simulation model to experimental data.
### All iteration calculation
######---- test 1 - data from one dosing only (3 articles) -----#####

pkdat_mlx <- read.csv("Data/Olaparib_PK_data_mlx2.csv") # table with mixed observed data from 4 articles 


pkdat_mlx %>% filter(ID=="300_Plummer_2015") # filter mixed table to obtain table containing observed data from "300_Plummer_2015" article 
pkdat_mlx %>% filter(ID=="300_Rolfo_2020") 
pkdat_mlx %>% filter(ID=="300_Rolfo_2019")

# Defines which PK output variables to simulate - plasma concentration and urine excretion
outputs_v  <- c("Organism|PeripheralVenousBlood|*|Plasma (Peripheral Venous Blood)",
                "Organism|Kidney|Urine|olaparib_compound")

# Olaparib_300mg_BID_Weibull|Neighborhoods|Kidney_pls_Kidney_ur|olaparib_compound|Renal Clearances-Olaparib renal clearance-olaparib_compound|Specific clearance
# A prior distribution is defined for key PK parameters (make a dataframe)
priors_setup_df <- data.frame(ParamPath = c("olaparib_compound|Lipophilicity",
                                            "olaparib_compound-CYP3A4-Reddy et al paper|Vmax",
                                            "olaparib_compound-CYP3A4-Reddy et al paper|Km",
                                            "olaparib_compound|Intestinal permeability (transcellular)",  
                                            "Neighborhoods|Kidney_pls_Kidney_ur|olaparib_compound|Renal Clearances-Olaparib renal clearance-olaparib_compound|Plasma clearance"),
                                            
                                         
                              ParamName= c("logP", "Vmax","Km","Pint","CLren"), 
                              ParamDistr = c("Normal"), # assume normal distribution of the priors
                              # parameters were taken from the literature 
                              ParamTV = c(1.55,   1,   60,  1e-04, 3e-04 ), # typical value = mean
                              ParamCV = c(0.4,   1.0,  0.4,  1.0,   1.0)) # coefficient of variation

# SMC-ABC Iterations run in two phases:

# Phase 1 (Iterations 1-3)
iter_info <- data.frame( Iteration =c(1,2,3), # three sequential iterations
                         dosingID=c("300_Plummer_2015"),
                         Nposterior=c(30), # The number of "good" parameter sets (particles) to accept in each iteration.
                         Epsilon=c(5, 2, 1), # The tolerance threshold. In Iteration 1, any simulation with an OFV < 5 is accepted. This threshold is gradually tightened (5 -> 2 -> 1) in subsequent iterations, forcing the algorithm to find better-fitting parameters (iteration 2: parameters based on the 30 "good" sets from Iteration 1 and only accepts them if the OFV is less than 2).
                         ScCov= c(0.25))

startTime <- Sys.time() 
# SMC_all_iterations function - runs the ABC-SMC algorithm (For each iteration, it uses the accepted parameters from the previous iteration to intelligently propose new ones, getting closer to the target with each step)
SMC_all_iterations(iter_info=iter_info, Resultsfile="Results/Olaparib_PK_posterior",
                   pkdat_mlx=pkdat_mlx, # Proposes parameters from priors
                   outputs_v = outputs_v, 
                   drug_MW=434, drug_in_urine_part=0.15, amt_urine_rse=0.4,
                   pkml_file="PKML/Olaparib_300mg_BID_Weibull_v2.pkml", # Simulates PK profiles using PKSim
                   OFVcalc_fun = OFVcalc_fun3urine, # Compares simulations to data using OFVcalc_fun3urine (сomputes goodness-of-fit between simulated and observed PK data)
                   priors_setup_df=priors_setup_df, # Accepts/rejects parameters based on Epsilon
                   fixed_setup_df= data.frame())
# The algorithm has found 40 particles that meet the strict Epsilon=1 criterion. This set of 40 will be saved in the Olaparib_PK_posterior3.csv file.

endTime <- Sys.time() 
# prints recorded time 
print(endTime - startTime)

posterior3 <- read.csv("Results/Olaparib_PK_posterior3.csv")
PosteriorPlots_fun(posterior3) # Visualize posterior distributions
cor(posterior3[,1:nrow(priors_setup_df)]) # Calculates the correlation between the estimated parameters (e.g., if increasing one parameter has the same effect as decreasing another)
summary(posterior3) # Provides summary statistics (mean, median, etc.) for the posterior distributions

# Visualize PK profile fits
# It takes 10 of the accepted parameter sets, runs simulations with them, and plots the resulting PK profiles on top of the real data. This provides a visual check that the accepted parameters do, in fact, produce simulations that look like the real data.
simPlot_fun(pkdat_mlx, pkml_file ="PKML/Olaparib_300mg_BID_Weibull_v2.pkml", 
            outputs_v = outputs_v,
            priors_setup_df= priors_setup_df, fixed_setup_df= data.frame() ,
            posterior_set=posterior3[1:10,], 
            selID="300_Plummer_2015")

# Phase 2 (Iterations 4-6)
iter_info <- data.frame( Iteration =c(4, 5, 6), 
                         dosingID=c("300_Plummer_2015"),
                         Nposterior=c(50,  50,  200),
                         Epsilon=c( 0.5,   0.3,  0.15),
                         ScCov= c(0.25))

startTime <- Sys.time() 
SMC_all_iterations(iter_info=iter_info, Resultsfile="Results/Olaparib_PK_posterior",
                   pkdat_mlx=pkdat_mlx,
                   outputs_v = outputs_v,
                   drug_MW=434, drug_in_urine_part=0.15, amt_urine_rse=0.4,
                   pkml_file="PKML/Olaparib_300mg_BID_Weibull_v2.pkml",
                   OFVcalc_fun = OFVcalc_fun3urine,
                   priors_setup_df=priors_setup_df,
                   fixed_setup_df= data.frame())

endTime <- Sys.time() 
# prints recorded time 
print(endTime - startTime)

posterior6 <- read.csv("Results/Olaparib_PK_posterior6.csv")
PosteriorPlots_fun(posterior6)
cor(posterior6[,1:nrow(priors_setup_df)])
summary(posterior6)

simPlot_fun(pkdat_mlx, pkml_file ="PKML/Olaparib_300mg_BID_Weibull_v2.pkml", 
            outputs_v = outputs_v,
            priors_setup_df= priors_setup_df, fixed_setup_df= data.frame() ,
            posterior_set=posterior6[sample(size=10,x=c(1:200)),], 
            selID="300_Plummer_2015")
#####
# the final summary of the entire analysis
my_data <- posterior6 %>% select(  priors_setup_df$ParamName ) # prepares the final dataset (posterior6) for plotting.
pdf("./Results/Olaparib_5param_corrplot_femail.pdf") # saves the resulting plot as a high-quality PDF file

# creates a matrix plot
# On the diagonal: The final posterior distribution (as a histogram) for each estimated parameter. This shows the most likely value and its uncertainty
# In the lower triangle: Scatter plots showing the relationship between every pair of parameters.
# In the upper triangle: The calculated correlation coefficient for every pair of parameters.
chart.Correlation(my_data, histogram=T, pch=19)
dev.off()





################### Find 300 sets of parameters with epsilon = 0.15 using sequential Monte-Carlo
priors_setup_df <- data.frame(ParamPath = c("olaparib_compound|Lipophilicity",
                                            "olaparib_compound-CYP3A4-Reddy et al paper|Vmax",
                                            "olaparib_compound-CYP3A4-Reddy et al paper|Km",
                                            "olaparib_compound|Intestinal permeability (transcellular)",  
                                            "Neighborhoods|Kidney_pls_Kidney_ur|olaparib_compound|Renal Clearances-Olaparib renal clearance-olaparib_compound|Plasma clearance"),
                              
                              
                              ParamName= c("logP", "Vmax","Km","Pint","CLren"), 
                              ParamDistr = c("Normal"), # assume normal distribution of the priors
                              ParamTV = c(1.55,   1,   60,  1e-04, 3e-04 ), # typical value = mean
                              ParamCV = c(0.4,   1.0,  0.4,  1.0,   1.0)) # coefficient of variation

# SMC-ABC Iterations run in two phases:

# Iterations 1-6
iter_info <- data.frame( Iteration =c(1,2,3, 4, 5, 6), 
                         dosingID=c("300_Plummer_2015"),
                         Nposterior=c(30, 30, 30, 50, 50, 300), 
                         Epsilon=c(5, 2, 1, 0.5, 0.3, 0.15), 
                         ScCov= c(0.25))

startTime <- Sys.time() 
# SMC_all_iterations function - runs the ABC-SMC algorithm (For each iteration, it uses the accepted parameters from the previous iteration to intelligently propose new ones, getting closer to the target with each step)
SMC_all_iterations(iter_info=iter_info, Resultsfile="Results/Anya/sequentialMC1/Olaparib_PK_posterior",
                   pkdat_mlx=pkdat_mlx, # Proposes parameters from priors
                   outputs_v = outputs_v, 
                   drug_MW=434, drug_in_urine_part=0.15, amt_urine_rse=0.4,
                   pkml_file="PKML/Olaparib_300mg_BID_Weibull_v2.pkml", # Simulates PK profiles using PKSim
                   OFVcalc_fun = OFVcalc_fun3urine, # Compares simulations to data using OFVcalc_fun3urine (сomputes goodness-of-fit between simulated and observed PK data)
                   priors_setup_df=priors_setup_df, # Accepts/rejects parameters based on Epsilon
                   fixed_setup_df= data.frame())


endTime <- Sys.time() 
# prints recorded time 
print(endTime - startTime) # Time difference of  mins

posterior6 <- read.csv("Results/Anya/sequentialMC1/Olaparib_PK_posterior6.csv")
PosteriorPlots_fun(posterior6) # Visualize posterior distributions
cor(posterior6[,1:nrow(priors_setup_df)]) # Calculates the correlation between the estimated parameters (e.g., if increasing one parameter has the same effect as decreasing another)
summary(posterior6) # Provides summary statistics (mean, median, etc.) for the posterior distributions

# Visualize PK profile fits
# It takes 10 of the accepted parameter sets, runs simulations with them, and plots the resulting PK profiles on top of the real data. This provides a visual check that the accepted parameters do, in fact, produce simulations that look like the real data.
simPlot_fun(pkdat_mlx, pkml_file ="PKML/Olaparib_300mg_BID_Weibull_v2.pkml", 
            outputs_v = outputs_v,
            priors_setup_df= priors_setup_df, fixed_setup_df= data.frame() ,
            posterior_set=posterior6[1:10,], 
            selID="300_Plummer_2015")





################### Find 300 sets of parameters with epsilon = 0.15 using "rejection from prior" (only 1st iteration)
iter_info <- data.frame( Iteration =c(1), 
                         dosingID=c("300_Plummer_2015"),
                         Nposterior=c(300), 
                         Epsilon=c(0.15), 
                         ScCov= c(0.25))

startTime <- Sys.time() 
# SMC_all_iterations function - runs the ABC-SMC algorithm (For each iteration, it uses the accepted parameters from the previous iteration to intelligently propose new ones, getting closer to the target with each step)
SMC_all_iterations(iter_info=iter_info, Resultsfile="Results/Anya/RejectionFromPrior1/Olaparib_PK_posterior",
                   pkdat_mlx=pkdat_mlx, # Proposes parameters from priors
                   outputs_v = outputs_v, 
                   drug_MW=434, drug_in_urine_part=0.15, amt_urine_rse=0.4,
                   pkml_file="PKML/Olaparib_300mg_BID_Weibull_v2.pkml", # Simulates PK profiles using PKSim
                   OFVcalc_fun = OFVcalc_fun3urine, # Compares simulations to data using OFVcalc_fun3urine (сomputes goodness-of-fit between simulated and observed PK data)
                   priors_setup_df=priors_setup_df, # Accepts/rejects parameters based on Epsilon
                   fixed_setup_df= data.frame())


endTime <- Sys.time() 
# prints recorded time 
print(endTime - startTime) # Time difference of 42.60823 mins

posterior1 <- read.csv("Results/Anya/RejectionFromPrior1/Olaparib_PK_posterior1.csv")
PosteriorPlots_fun(posterior1) # Visualize posterior distributions
cor(posterior1[,1:nrow(priors_setup_df)]) # Calculates the correlation between the estimated parameters (e.g., if increasing one parameter has the same effect as decreasing another)
summary(posterior1) # Provides summary statistics (mean, median, etc.) for the posterior distributions

# Visualize PK profile fits
# It takes 10 of the accepted parameter sets, runs simulations with them, and plots the resulting PK profiles on top of the real data. This provides a visual check that the accepted parameters do, in fact, produce simulations that look like the real data.
simPlot_fun(pkdat_mlx, pkml_file ="PKML/Olaparib_300mg_BID_Weibull_v2.pkml", 
            outputs_v = outputs_v,
            priors_setup_df= priors_setup_df, fixed_setup_df= data.frame() ,
            posterior_set=posterior1[1:10,], 
            selID="300_Plummer_2015")






################### Create a specific type of function that visualize the comparison of posterior parameter distributions from ABC-SMC algorithm and Rejection from priors algorithm
install.packages("dplyr")
library(dplyr)
install.packages("tidyr")
library(tidyr)
library(ggplot2)

PosteriorPlots_fun2 <- function(posterior1, posterior2, appr1, appr2){
  posterior1v <- posterior1 %>% select(-OFV, -acceptF,-w) %>% # removes the columns named OFV (Objective Function Value), acceptF (acceptance flag), and w (particle weight)
    gather(key="PARname", value="PARvalue") %>% mutate(Method=appr1) # converts the data from a "wide" format to a "long" format
  posterior2v <- posterior2 %>% select(-OFV, -acceptF,-w) %>%
    gather(key="PARname", value="PARvalue") %>% mutate(Method=appr2)
  posterior12v <- rbind( posterior1v,  posterior2v) # stacks the two long-format data frames
  
  ggplot(posterior12v, aes(x=PARvalue, fill=as.factor(Method)))+ geom_density(alpha=0.2)  + theme_bw()+
    facet_wrap(~PARname, ncol = 2, scales = "free")
}

PosteriorPlots_fun2(posterior1=posterior1, posterior2= posterior6, appr1="Rejections from prior", appr2="SMC")
# В принципе,  posterior distributions, полученные разными методами, не сильно отличаются. Но бимодальность для SMC скорее всего артефакт расчета. 
# Попробуем дополнительные итерации, чтобы от него уйти.

# Carry out 3 more iterations (sequential MC) to try not to see bimodal in distribution
iter_info <- data.frame( Iteration =c(7, 8, 9),
                         dosingID=c("300_Plummer_2015"),
                         Nposterior=c(100, 200,  300),
                         Epsilon=c( 0.15,   0.15,  0.15),
                         ScCov= c(0.25))

startTime <- Sys.time()
SMC_all_iterations(iter_info=iter_info, Resultsfile="Results/Anya/sequentialMC1/Olaparib_PK_posterior",
                   pkdat_mlx=pkdat_mlx,
                   outputs_v = outputs_v,
                   drug_MW=434, drug_in_urine_part=0.15, amt_urine_rse=0.4,
                   pkml_file="PKML/Olaparib_300mg_BID_Weibull_v2.pkml",
                   OFVcalc_fun = OFVcalc_fun3urine,
                   priors_setup_df=priors_setup_df,
                   fixed_setup_df= data.frame())

endTime <- Sys.time()
# prints recorded time
print(endTime - startTime)

posterior9 <- read.csv("Results/Anya/sequentialMC1/Olaparib_PK_posterior9.csv")
PosteriorPlots_fun(posterior9)
cor(posterior9[,1:nrow(priors_setup_df)])
summary(posterior9)

PosteriorPlots_fun2(posterior1=posterior1, posterior2= posterior9, appr1="Rejections from prior", appr2="SMC")
# It works!








################### Create a function that calculates a difference score between two posterior distributions from ABC-SMC and Rejection Sampling
calculate_distribution_metric <- function(ref_df, # posterior data frame from Rejection Sampling
                                          comp_df, # posterior data frame from SMC
                                          median_threshold = 0.10, 
                                          iqr_threshold = 0.10,
                                          corr_threshold = 0.50) {
  
  # Isolate only the parameter columns
  # Remove non-parameter columns that the ABC-SMC function adds
  params_ref <- ref_df %>% select(-any_of(c("OFV", "acceptF", "w")))
  params_comp <- comp_df %>% select(-any_of(c("OFV", "acceptF", "w")))
  
  # Get the names of the parameters being compared
  param_names <- names(params_ref)
  
  # Initialize the metric score
  difference_score <- 0
  
  # Create a data frame to store detailed results
  breakdown <- data.frame(
    Parameter = character(),
    Metric = character(),
    Penalty_Added = numeric(),
    Reason = character()
  )
  
  # Compare medians and IQRs for each parameter
  for (param in param_names) {
    # Calculate stats for reference distribution (from rejection sampling)
    median_ref <- median(params_ref[[param]])
    iqr_ref <- IQR(params_ref[[param]])
    
    # Calculate stats for comparison distribution (from SMC)
    median_comp <- median(params_comp[[param]])
    iqr_comp <- IQR(params_comp[[param]])
    
    # Check median difference
    if (median_ref != 0) {
      median_diff <- abs((median_comp - median_ref) / median_ref)
      if (median_diff > median_threshold) {
        difference_score <- difference_score + 1
        breakdown <- rbind(breakdown, data.frame(Parameter=param, Metric="Median", Penalty_Added=1, Reason=paste0("Differs by ", round(median_diff*100,1), "% (> ", median_threshold*100, "%)")))
      } else {
        breakdown <- rbind(breakdown, data.frame(Parameter=param, Metric="Median", Penalty_Added=0, Reason="Within threshold"))
      }
    }
    
    # Check IQR difference
    if (iqr_ref != 0) {
      iqr_diff <- abs((iqr_comp - iqr_ref) / iqr_ref)
      if (iqr_diff > iqr_threshold) {
        difference_score <- difference_score + 1
        breakdown <- rbind(breakdown, data.frame(Parameter=param, Metric="IQR", Penalty_Added=1, Reason=paste0("Differs by ", round(iqr_diff*100,1), "% (> ", iqr_threshold*100, "%)")))
      } else {
        breakdown <- rbind(breakdown, data.frame(Parameter=param, Metric="IQR", Penalty_Added=0, Reason="Within threshold"))
      }
    }
  }
  
  # Compare correlation matrices
  cor_ref <- cor(params_ref)
  cor_comp <- cor(params_comp)
  
  # Iterate through the upper triangle of the correlation matrix to avoid duplicates
  for (i in 1:(nrow(cor_ref) - 1)) {
    for (j in (i + 1):ncol(cor_ref)) {
      
      param_pair_name <- paste(rownames(cor_ref)[i], "-", colnames(cor_ref)[j])
      
      is_strong_ref <- abs(cor_ref[i, j]) > corr_threshold
      is_strong_comp <- abs(cor_comp[i, j]) > corr_threshold
      
      # If one is strong and the other is not, add a penalty
      if (is_strong_ref != is_strong_comp) {
        difference_score <- difference_score + 1
        reason_text <- paste0("Ref_cor=", round(cor_ref[i, j], 2), ", Comp_cor=", round(cor_comp[i, j], 2), ". Mismatch in 'strong' correlation.")
        breakdown <- rbind(breakdown, data.frame(Parameter=param_pair_name, Metric="Correlation", Penalty_Added=1, Reason=reason_text))
      } else {
        breakdown <- rbind(breakdown, data.frame(Parameter=param_pair_name, Metric="Correlation", Penalty_Added=0, Reason="Correlation strength matches"))
      }
    }
  }
  
  # Return the final results
  return(list(
    total_difference_score = difference_score,
    detailed_breakdown = breakdown
  ))
}

# Run the metric function
metric_results <- calculate_distribution_metric(ref_df = posterior1, 
                                                comp_df = posterior9)
# View the results
print(paste("Total Difference Score:", metric_results$total_difference_score))
print("Detailed Breakdown:")
detailed_breakdown <- metric_results$detailed_breakdown
print(detailed_breakdown)







##############################################################################################
##############################   Parameter optimization   ####################################


priors_setup_df <- data.frame(ParamPath = c("olaparib_compound|Lipophilicity",
                                            "olaparib_compound-CYP3A4-Reddy et al paper|Vmax",
                                            "olaparib_compound-CYP3A4-Reddy et al paper|Km",
                                            "olaparib_compound|Intestinal permeability (transcellular)",  
                                            "Neighborhoods|Kidney_pls_Kidney_ur|olaparib_compound|Renal Clearances-Olaparib renal clearance-olaparib_compound|Plasma clearance"),
                              
                              
                              ParamName= c("logP", "Vmax","Km","Pint","CLren"),
                              ParamDistr = c("Normal"),
                              ParamTV = c(1.55,   1,   60,  1e-04, 3e-04 ),
                              ParamCV = c(0.2,   2.0,  0.2,  2.0,   2.0))

################### Run rejection from priors 
iter_info <- data.frame( Iteration =c(1), 
                         dosingID=c("300_Plummer_2015"),
                         Nposterior=c(300),
                         Epsilon=c(0.15),
                         ScCov= c(0.25))

startTime <- Sys.time() 
SMC_all_iterations(iter_info=iter_info, Resultsfile="Results/Anya/CV_2/RejectionFromPrior/Olaparib_PK_PriorRej_posterior",
                   pkdat_mlx=pkdat_mlx,
                   outputs_v = outputs_v,
                   drug_MW=434, drug_in_urine_part=0.15, amt_urine_rse=0.4,
                   pkml_file="PKML/Olaparib_300mg_BID_Weibull_v2.pkml",
                   OFVcalc_fun = OFVcalc_fun3urine,
                   priors_setup_df=priors_setup_df,
                   fixed_setup_df= data.frame())

endTime <- Sys.time() 
# prints recorded time 
print(endTime - startTime) # 32 min on Yuri's laptop

posterior1 <- read.csv("Results/Anya/CV_2/RejectionFromPrior/Olaparib_PK_PriorRej_posterior1.csv")





################### Run ABC_SMC
# 0. Base metrics 

output_dir <- "Results/Anya/CV_2/SequentialMC/0_base_metrics/"

iter_info <- data.frame( Iteration =c(1, 2, 3, 4, 5, 6),
                         dosingID=c("300_Plummer_2015"),
                         Nposterior=c(50, 50, 50, 100, 200, 300),
                         Epsilon=c(0.5, 0.3, 0.25, 0.20, 0.18, 0.15),
                         ScCov= c(0.25))

startTime <- Sys.time()
SMC_all_iterations(iter_info=iter_info, Resultsfile="Results/Anya/CV_2/SequentialMC/0_base_metrics/Olaparib_PK_posterior",
                   pkdat_mlx=pkdat_mlx,
                   outputs_v = outputs_v,
                   drug_MW=434, drug_in_urine_part=0.15, amt_urine_rse=0.4,
                   pkml_file="PKML/Olaparib_300mg_BID_Weibull_v2.pkml",
                   OFVcalc_fun = OFVcalc_fun3urine,
                   priors_setup_df=priors_setup_df,
                   fixed_setup_df= data.frame())

endTime <- Sys.time()
# prints recorded time
time0 <- endTime - startTime # 35.03735 mins
print(time0)

cat("Saving execution time...\n")
time_text <- format(time0) # format() makes it a nice character string
writeLines(
  time_text,
  con = file.path(output_dir, "execution_time.txt")
)

posterior6 <- read.csv("Results/Anya/CV_2/SequentialMC/0_base_metrics/Olaparib_PK_posterior6.csv")

cat("Saving posterior distribution plots...\n")
posterior_dist_plot <- PosteriorPlots_fun(posterior6)
# Now, add a title to the plot object
posterior_dist_plot <- posterior_dist_plot + 
  labs(title = "Posterior parameter distributions (base metrics)")
ggsave(
  filename = file.path(output_dir, "posterior_distributions.pdf"),
  plot = posterior_dist_plot,
  width = 8,
  height = 7
)

cat("Saving correlation matrix...\n")
correlation_matrix <- cor(posterior6[,1:nrow(priors_setup_df)])
write.csv(
  correlation_matrix,
  file = file.path(output_dir, "correlation_matrix.csv")
)

cat("Saving summary statistics...\n")
summary_text <- capture.output(summary(posterior6))
writeLines(
  summary_text,
  con = file.path(output_dir, "posterior_summary.txt")
)


# Visualize PK profile fits
# It takes 10 of the accepted parameter sets, runs simulations with them, and plots the resulting PK profiles on top of the real data. This provides a visual check that the accepted parameters do, in fact, produce simulations that look like the real data.
pk_fit_plot <- simPlot_fun(pkdat_mlx, pkml_file ="PKML/Olaparib_300mg_BID_Weibull_v2.pkml", 
            outputs_v = outputs_v,
            priors_setup_df= priors_setup_df, fixed_setup_df= data.frame() ,
            posterior_set=posterior6[1:10,], 
            selID="300_Plummer_2015")
cat("Saving PK profile fit plot...\n")
ggsave(
  filename = file.path(output_dir, "pk_profile_fit.pdf"),
  plot = pk_fit_plot,
  width = 7,
  height = 5
)

# Run the metric function
metric_results <- calculate_distribution_metric(ref_df = posterior1, 
                                                comp_df = posterior6)
# View the results
print(paste("Total Difference Score:", metric_results$total_difference_score))
print("Detailed Breakdown:")
detailed_breakdown_0 <- metric_results$detailed_breakdown
print(detailed_breakdown_0)

cat("Saving metric breakdown...\n")
write.csv(
  detailed_breakdown_0,
  file = file.path(output_dir, "metric_breakdown.csv"),
  row.names = FALSE # Usually good practice for analysis data frames
)


# Plot with comparison of posterior distributions (rejection from prior vs. SMC)
posterior_dist_plot2 <- PosteriorPlots_fun2(posterior1=posterior1, posterior2= posterior6, appr1="Rejections from prior", appr2="SMC")
cat("Saving posterior distribution plots...\n")
# Now, add a title to the plot object
posterior_dist_plot2 <- posterior_dist_plot2 + 
  labs(title = "Comparison of posterior parameter distributions (base metrics)")
ggsave(
  filename = file.path(output_dir, "posterior_distributions_comparison.pdf"),
  plot = posterior_dist_plot2,
  width = 8,
  height = 7
)




################### 1. Iteration number 
# 1a. N = 4

output_dir <- "Results/Anya/CV_2/SequentialMC/1a_N_of_iterations/"

iter_info <- data.frame( Iteration =c(1, 2, 3, 4),
                         dosingID=c("300_Plummer_2015"),
                         Nposterior=c(50, 100, 200, 300),
                         Epsilon=c(0.5, 0.35, 0.2, 0.15),
                         ScCov= c(0.25))

startTime <- Sys.time()
SMC_all_iterations(iter_info=iter_info, Resultsfile="Results/Anya/CV_2/SequentialMC/1a_N_of_iterations/Olaparib_PK_posterior",
                   pkdat_mlx=pkdat_mlx,
                   outputs_v = outputs_v,
                   drug_MW=434, drug_in_urine_part=0.15, amt_urine_rse=0.4,
                   pkml_file="PKML/Olaparib_300mg_BID_Weibull_v2.pkml",
                   OFVcalc_fun = OFVcalc_fun3urine,
                   priors_setup_df=priors_setup_df,
                   fixed_setup_df= data.frame())

endTime <- Sys.time()
# prints recorded time
time1a <- endTime - startTime # 26.65571 mins
print(time1a)

cat("Saving execution time...\n")
time_text <- format(time1a) # format() makes it a nice character string
writeLines(
  time_text,
  con = file.path(output_dir, "execution_time.txt")
)

posterior4 <- read.csv("Results/Anya/CV_2/SequentialMC/1a_N_of_iterations/Olaparib_PK_posterior4.csv")

cat("Saving posterior distribution plots...\n")
posterior_dist_plot <- PosteriorPlots_fun(posterior4)
# Now, add a title to the plot object
posterior_dist_plot <- posterior_dist_plot + 
  labs(title = "Posterior parameter distributions (N = 4)")
ggsave(
  filename = file.path(output_dir, "posterior_distributions.pdf"),
  plot = posterior_dist_plot,
  width = 8,
  height = 7
)

cat("Saving correlation matrix...\n")
correlation_matrix <- cor(posterior4[,1:nrow(priors_setup_df)])
write.csv(
  correlation_matrix,
  file = file.path(output_dir, "correlation_matrix.csv")
)

cat("Saving summary statistics...\n")
summary_text <- capture.output(summary(posterior4))
writeLines(
  summary_text,
  con = file.path(output_dir, "posterior_summary.txt")
)

# Visualize PK profile fits
# It takes 10 of the accepted parameter sets, runs simulations with them, and plots the resulting PK profiles on top of the real data. This provides a visual check that the accepted parameters do, in fact, produce simulations that look like the real data.
pk_fit_plot <- simPlot_fun(pkdat_mlx, pkml_file ="PKML/Olaparib_300mg_BID_Weibull_v2.pkml", 
            outputs_v = outputs_v,
            priors_setup_df= priors_setup_df, fixed_setup_df= data.frame() ,
            posterior_set=posterior4[1:10,], 
            selID="300_Plummer_2015")
cat("Saving PK profile fit plot...\n")
ggsave(
  filename = file.path(output_dir, "pk_profile_fit.pdf"),
  plot = pk_fit_plot,
  width = 7,
  height = 5
)

# Run the metric function
metric_results <- calculate_distribution_metric(ref_df = posterior1, 
                                                comp_df = posterior4)
# View the results
print(paste("Total Difference Score:", metric_results$total_difference_score))
print("Detailed Breakdown:")
detailed_breakdown_1a <- metric_results$detailed_breakdown
print(detailed_breakdown_1a)

cat("Saving metric breakdown...\n")
write.csv(
  detailed_breakdown_1a,
  file = file.path(output_dir, "metric_breakdown.csv"),
  row.names = FALSE # Usually good practice for analysis data frames
)


# Plot with comparison of posterior distributions (rejection from prior vs. SMC)
posterior_dist_plot2 <- PosteriorPlots_fun2(posterior1=posterior1, posterior2= posterior4, appr1="Rejections from prior", appr2="SMC")
cat("Saving posterior distribution plots...\n")
# Now, add a title to the plot object
posterior_dist_plot2 <- posterior_dist_plot2 + 
  labs(title = "Comparison of posterior parameter distributions (N = 4)")
ggsave(
  filename = file.path(output_dir, "posterior_distributions_comparison.pdf"),
  plot = posterior_dist_plot2,
  width = 8,
  height = 7
)





# 1b. N = 9

output_dir <- "Results/Anya/CV_2/SequentialMC/1b_N_of_iterations/"

iter_info <- data.frame( Iteration =c(1, 2, 3, 4, 5, 6, 7, 8, 9),
                         dosingID=c("300_Plummer_2015"),
                         Nposterior=c(50, 50, 50, 75, 75, 75, 100, 200, 300),
                         Epsilon=c(0.5, 0.4, 0.35, 0.3, 0.25, 0.20, 0.18, 0.16, 0.15),
                         ScCov= c(0.25))

startTime <- Sys.time()
SMC_all_iterations(iter_info=iter_info, Resultsfile="Results/Anya/CV_2/SequentialMC/1b_N_of_iterations/Olaparib_PK_posterior",
                   pkdat_mlx=pkdat_mlx,
                   outputs_v = outputs_v,
                   drug_MW=434, drug_in_urine_part=0.15, amt_urine_rse=0.4,
                   pkml_file="PKML/Olaparib_300mg_BID_Weibull_v2.pkml",
                   OFVcalc_fun = OFVcalc_fun3urine,
                   priors_setup_df=priors_setup_df,
                   fixed_setup_df= data.frame())

endTime <- Sys.time()
# prints recorded time
time1b <- endTime - startTime # 32.41298 mins
print(time1b)

cat("Saving execution time...\n")
time_text <- format(time1b) # format() makes it a nice character string
writeLines(
  time_text,
  con = file.path(output_dir, "execution_time.txt")
)

posterior9 <- read.csv("Results/Anya/CV_2/SequentialMC/1b_N_of_iterations/Olaparib_PK_posterior9.csv")

cat("Saving posterior distribution plots...\n")
posterior_dist_plot <- PosteriorPlots_fun(posterior9)
# Now, add a title to the plot object
posterior_dist_plot <- posterior_dist_plot + 
  labs(title = "Posterior parameter distributions (N = 9)")
ggsave(
  filename = file.path(output_dir, "posterior_distributions.pdf"),
  plot = posterior_dist_plot,
  width = 8,
  height = 7
)

cat("Saving correlation matrix...\n")
correlation_matrix <- cor(posterior9[,1:nrow(priors_setup_df)])
write.csv(
  correlation_matrix,
  file = file.path(output_dir, "correlation_matrix.csv")
)

cat("Saving summary statistics...\n")
summary_text <- capture.output(summary(posterior9))
writeLines(
  summary_text,
  con = file.path(output_dir, "posterior_summary.txt")
)

# Visualize PK profile fits
# It takes 10 of the accepted parameter sets, runs simulations with them, and plots the resulting PK profiles on top of the real data. This provides a visual check that the accepted parameters do, in fact, produce simulations that look like the real data.
pk_fit_plot <- simPlot_fun(pkdat_mlx, pkml_file ="PKML/Olaparib_300mg_BID_Weibull_v2.pkml", 
            outputs_v = outputs_v,
            priors_setup_df= priors_setup_df, fixed_setup_df= data.frame() ,
            posterior_set=posterior9[1:10,], 
            selID="300_Plummer_2015")
cat("Saving PK profile fit plot...\n")
ggsave(
  filename = file.path(output_dir, "pk_profile_fit.pdf"),
  plot = pk_fit_plot,
  width = 7,
  height = 5
)

# Run the metric function
metric_results <- calculate_distribution_metric(ref_df = posterior1, 
                                                comp_df = posterior9)
# View the results
print(paste("Total Difference Score:", metric_results$total_difference_score))
print("Detailed Breakdown:")
detailed_breakdown_1b <- metric_results$detailed_breakdown
print(detailed_breakdown_1b)

cat("Saving metric breakdown...\n")
write.csv(
  detailed_breakdown_1b,
  file = file.path(output_dir, "metric_breakdown.csv"),
  row.names = FALSE # Usually good practice for analysis data frames
)


# Plot with comparison of posterior distributions (rejection from prior vs. SMC)
posterior_dist_plot2 <- PosteriorPlots_fun2(posterior1=posterior1, posterior2= posterior9, appr1="Rejections from prior", appr2="SMC")
cat("Saving posterior distribution plots...\n")
# Now, add a title to the plot object
posterior_dist_plot2 <- posterior_dist_plot2 + 
  labs(title = "Comparison of posterior parameter distributions (N = 9)")
ggsave(
  filename = file.path(output_dir, "posterior_distributions_comparison.pdf"),
  plot = posterior_dist_plot2,
  width = 8,
  height = 7
)





################### 2. Start epsilon 
# 2a. epsilon = 0.3

output_dir <- "Results/Anya/CV_2/SequentialMC/2a_epsilon/"

iter_info <- data.frame( Iteration =c(1, 2, 3, 4, 5, 6),
                         dosingID=c("300_Plummer_2015"),
                         Nposterior=c(50, 50, 50, 100, 200, 300),
                         Epsilon=c(0.3, 0.25, 0.2, 0.18, 0.16, 0.15),
                         ScCov= c(0.25))

startTime <- Sys.time()
SMC_all_iterations(iter_info=iter_info, Resultsfile="Results/Anya/CV_2/SequentialMC/2a_epsilon/Olaparib_PK_posterior",
                   pkdat_mlx=pkdat_mlx,
                   outputs_v = outputs_v,
                   drug_MW=434, drug_in_urine_part=0.15, amt_urine_rse=0.4,
                   pkml_file="PKML/Olaparib_300mg_BID_Weibull_v2.pkml",
                   OFVcalc_fun = OFVcalc_fun3urine,
                   priors_setup_df=priors_setup_df,
                   fixed_setup_df= data.frame())

endTime <- Sys.time()
# prints recorded time
time2a <- endTime - startTime # 24.8635 mins
print(time2a)

cat("Saving execution time...\n")
time_text <- format(time2a) # format() makes it a nice character string
writeLines(
  time_text,
  con = file.path(output_dir, "execution_time.txt")
)

posterior6 <- read.csv("Results/Anya/CV_2/SequentialMC/2a_epsilon/Olaparib_PK_posterior6.csv")

cat("Saving posterior distribution plots...\n")
posterior_dist_plot <- PosteriorPlots_fun(posterior6)
# Now, add a title to the plot object
posterior_dist_plot <- posterior_dist_plot + 
  labs(title = "Posterior parameter distributions (epsilon = 0.3)")
ggsave(
  filename = file.path(output_dir, "posterior_distributions.pdf"),
  plot = posterior_dist_plot,
  width = 8,
  height = 7
)

cat("Saving correlation matrix...\n")
correlation_matrix <- cor(posterior6[,1:nrow(priors_setup_df)])
write.csv(
  correlation_matrix,
  file = file.path(output_dir, "correlation_matrix.csv")
)

cat("Saving summary statistics...\n")
summary_text <- capture.output(summary(posterior6))
writeLines(
  summary_text,
  con = file.path(output_dir, "posterior_summary.txt")
)

# Visualize PK profile fits
# It takes 10 of the accepted parameter sets, runs simulations with them, and plots the resulting PK profiles on top of the real data. This provides a visual check that the accepted parameters do, in fact, produce simulations that look like the real data.
pk_fit_plot <- simPlot_fun(pkdat_mlx, pkml_file ="PKML/Olaparib_300mg_BID_Weibull_v2.pkml", 
                           outputs_v = outputs_v,
                           priors_setup_df= priors_setup_df, fixed_setup_df= data.frame() ,
                           posterior_set=posterior6[1:10,], 
                           selID="300_Plummer_2015")
cat("Saving PK profile fit plot...\n")
ggsave(
  filename = file.path(output_dir, "pk_profile_fit.pdf"),
  plot = pk_fit_plot,
  width = 7,
  height = 5
)

# Run the metric function
metric_results <- calculate_distribution_metric(ref_df = posterior1, 
                                                comp_df = posterior6)
# View the results
print(paste("Total Difference Score:", metric_results$total_difference_score))
print("Detailed Breakdown:")
detailed_breakdown_2a <- metric_results$detailed_breakdown
print(detailed_breakdown_2a)

cat("Saving metric breakdown...\n")
write.csv(
  detailed_breakdown_2a,
  file = file.path(output_dir, "metric_breakdown.csv"),
  row.names = FALSE # Usually good practice for analysis data frames
)

# Plot with comparison of posterior distributions (rejection from prior vs. SMC)
posterior_dist_plot2 <- PosteriorPlots_fun2(posterior1=posterior1, posterior2= posterior6, appr1="Rejections from prior", appr2="SMC")
cat("Saving posterior distribution plots...\n")
# Now, add a title to the plot object
posterior_dist_plot2 <- posterior_dist_plot2 + 
  labs(title = "Comparison of posterior parameter distributions (epsilon = 0.3)")
ggsave(
  filename = file.path(output_dir, "posterior_distributions_comparison.pdf"),
  plot = posterior_dist_plot2,
  width = 8,
  height = 7
)





# 2b. epsilon = 1.0
output_dir <- "Results/Anya/CV_2/SequentialMC/2b_epsilon/"

iter_info <- data.frame( Iteration =c(1, 2, 3, 4, 5, 6),
                         dosingID=c("300_Plummer_2015"),
                         Nposterior=c(50, 50, 50, 100, 200, 300),
                         Epsilon=c(1.0, 0.7, 0.5, 0.3, 0.2, 0.15),
                         ScCov= c(0.25))

startTime <- Sys.time()
SMC_all_iterations(iter_info=iter_info, Resultsfile="Results/Anya/CV_2/SequentialMC/2b_epsilon/Olaparib_PK_posterior",
                   pkdat_mlx=pkdat_mlx,
                   outputs_v = outputs_v,
                   drug_MW=434, drug_in_urine_part=0.15, amt_urine_rse=0.4,
                   pkml_file="PKML/Olaparib_300mg_BID_Weibull_v2.pkml",
                   OFVcalc_fun = OFVcalc_fun3urine,
                   priors_setup_df=priors_setup_df,
                   fixed_setup_df= data.frame())

endTime <- Sys.time()
# prints recorded time
time2b <- endTime - startTime # 24.77847 mins
print(time2b)

cat("Saving execution time...\n")
time_text <- format(time2b) # format() makes it a nice character string
writeLines(
  time_text,
  con = file.path(output_dir, "execution_time.txt")
)

posterior6 <- read.csv("Results/Anya/CV_2/SequentialMC/2b_epsilon/Olaparib_PK_posterior6.csv")

cat("Saving posterior distribution plots...\n")
posterior_dist_plot <- PosteriorPlots_fun(posterior6)
# Now, add a title to the plot object
posterior_dist_plot <- posterior_dist_plot + 
  labs(title = "Posterior parameter distributions (epsilon = 1.0)")
ggsave(
  filename = file.path(output_dir, "posterior_distributions.pdf"),
  plot = posterior_dist_plot,
  width = 8,
  height = 7
)

cat("Saving correlation matrix...\n")
correlation_matrix <- cor(posterior6[,1:nrow(priors_setup_df)])
write.csv(
  correlation_matrix,
  file = file.path(output_dir, "correlation_matrix.csv")
)

cat("Saving summary statistics...\n")
summary_text <- capture.output(summary(posterior6))
writeLines(
  summary_text,
  con = file.path(output_dir, "posterior_summary.txt")
)


# Visualize PK profile fits
# It takes 10 of the accepted parameter sets, runs simulations with them, and plots the resulting PK profiles on top of the real data. This provides a visual check that the accepted parameters do, in fact, produce simulations that look like the real data.
pk_fit_plot <- simPlot_fun(pkdat_mlx, pkml_file ="PKML/Olaparib_300mg_BID_Weibull_v2.pkml", 
                           outputs_v = outputs_v,
                           priors_setup_df= priors_setup_df, fixed_setup_df= data.frame() ,
                           posterior_set=posterior6[1:10,], 
                           selID="300_Plummer_2015")
cat("Saving PK profile fit plot...\n")
ggsave(
  filename = file.path(output_dir, "pk_profile_fit.pdf"),
  plot = pk_fit_plot,
  width = 7,
  height = 5
)

# Run the metric function
metric_results <- calculate_distribution_metric(ref_df = posterior1, 
                                                comp_df = posterior6)
# View the results
print(paste("Total Difference Score:", metric_results$total_difference_score))
print("Detailed Breakdown:")
detailed_breakdown_2b <- metric_results$detailed_breakdown
print(detailed_breakdown_2b)

cat("Saving metric breakdown...\n")
write.csv(
  detailed_breakdown_2b,
  file = file.path(output_dir, "metric_breakdown.csv"),
  row.names = FALSE # Usually good practice for analysis data frames
)


# Plot with comparison of posterior distributions (rejection from prior vs. SMC)
posterior_dist_plot2 <- PosteriorPlots_fun2(posterior1=posterior1, posterior2= posterior6, appr1="Rejections from prior", appr2="SMC")
cat("Saving posterior distribution plots...\n")
# Now, add a title to the plot object
posterior_dist_plot2 <- posterior_dist_plot2 + 
  labs(title = "Comparison of posterior parameter distributions (epsilon = 1.0)")
ggsave(
  filename = file.path(output_dir, "posterior_distributions_comparison.pdf"),
  plot = posterior_dist_plot2,
  width = 8,
  height = 7
)





################### 3. Start N posterior
# 3a. Nposterior = 20
output_dir <- "Results/Anya/CV_2/SequentialMC/3a_N_of_posteriors/"

iter_info <- data.frame( Iteration =c(1, 2, 3, 4, 5, 6),
                         dosingID=c("300_Plummer_2015"),
                         Nposterior=c(20, 50, 50, 100, 200, 300),
                         Epsilon=c(0.5, 0.3, 0.25, 0.20, 0.18, 0.15),
                         ScCov= c(0.25))

startTime <- Sys.time()
SMC_all_iterations(iter_info=iter_info, Resultsfile="Results/Anya/CV_2/SequentialMC/3a_N_of_posteriors/Olaparib_PK_posterior",
                   pkdat_mlx=pkdat_mlx,
                   outputs_v = outputs_v,
                   drug_MW=434, drug_in_urine_part=0.15, amt_urine_rse=0.4,
                   pkml_file="PKML/Olaparib_300mg_BID_Weibull_v2.pkml",
                   OFVcalc_fun = OFVcalc_fun3urine,
                   priors_setup_df=priors_setup_df,
                   fixed_setup_df= data.frame())

endTime <- Sys.time()
# prints recorded time
time3a <- endTime - startTime # 31.74183 mins
print(time3a)

cat("Saving execution time...\n")
time_text <- format(time3a) # format() makes it a nice character string
writeLines(
  time_text,
  con = file.path(output_dir, "execution_time.txt")
)

posterior6 <- read.csv("Results/Anya/CV_2/SequentialMC/3a_N_of_posteriors/Olaparib_PK_posterior6.csv")

cat("Saving posterior distribution plots...\n")
posterior_dist_plot <- PosteriorPlots_fun(posterior6)
# Now, add a title to the plot object
posterior_dist_plot <- posterior_dist_plot + 
  labs(title = "Posterior parameter distributions (Nposterior = 20)")
ggsave(
  filename = file.path(output_dir, "posterior_distributions.pdf"),
  plot = posterior_dist_plot,
  width = 8,
  height = 7
)

cat("Saving correlation matrix...\n")
correlation_matrix <- cor(posterior6[,1:nrow(priors_setup_df)])
write.csv(
  correlation_matrix,
  file = file.path(output_dir, "correlation_matrix.csv")
)

cat("Saving summary statistics...\n")
summary_text <- capture.output(summary(posterior6))
writeLines(
  summary_text,
  con = file.path(output_dir, "posterior_summary.txt")
)


# Visualize PK profile fits
# It takes 10 of the accepted parameter sets, runs simulations with them, and plots the resulting PK profiles on top of the real data. This provides a visual check that the accepted parameters do, in fact, produce simulations that look like the real data.
pk_fit_plot <- simPlot_fun(pkdat_mlx, pkml_file ="PKML/Olaparib_300mg_BID_Weibull_v2.pkml", 
                           outputs_v = outputs_v,
                           priors_setup_df= priors_setup_df, fixed_setup_df= data.frame() ,
                           posterior_set=posterior6[1:10,], 
                           selID="300_Plummer_2015")
cat("Saving PK profile fit plot...\n")
ggsave(
  filename = file.path(output_dir, "pk_profile_fit.pdf"),
  plot = pk_fit_plot,
  width = 7,
  height = 5
)

# Run the metric function
metric_results <- calculate_distribution_metric(ref_df = posterior1, 
                                                comp_df = posterior6)
# View the results
print(paste("Total Difference Score:", metric_results$total_difference_score))
print("Detailed Breakdown:")
detailed_breakdown_3a <- metric_results$detailed_breakdown
print(detailed_breakdown_3a)

cat("Saving metric breakdown...\n")
write.csv(
  detailed_breakdown_3a,
  file = file.path(output_dir, "metric_breakdown.csv"),
  row.names = FALSE # Usually good practice for analysis data frames
)


# Plot with comparison of posterior distributions (rejection from prior vs. SMC)
posterior_dist_plot2 <- PosteriorPlots_fun2(posterior1=posterior1, posterior2= posterior6, appr1="Rejections from prior", appr2="SMC")
cat("Saving posterior distribution plots...\n")
# Now, add a title to the plot object
posterior_dist_plot2 <- posterior_dist_plot2 + 
  labs(title = "Comparison of posterior parameter distributions (Nposterior = 20)")
ggsave(
  filename = file.path(output_dir, "posterior_distributions_comparison.pdf"),
  plot = posterior_dist_plot2,
  width = 8,
  height = 7
)





# 3b. Nposterior = 100
output_dir <- "Results/Anya/CV_2/SequentialMC/3b_N_of_posteriors/"

iter_info <- data.frame( Iteration =c(1, 2, 3, 4, 5, 6),
                         dosingID=c("300_Plummer_2015"),
                         Nposterior=c(100, 100, 100, 150, 200, 300),
                         Epsilon=c(0.5, 0.3, 0.25, 0.20, 0.18, 0.15),
                         ScCov= c(0.25))

startTime <- Sys.time()
SMC_all_iterations(iter_info=iter_info, Resultsfile="Results/Anya/CV_2/SequentialMC/3b_N_of_posteriors/Olaparib_PK_posterior",
                   pkdat_mlx=pkdat_mlx,
                   outputs_v = outputs_v,
                   drug_MW=434, drug_in_urine_part=0.15, amt_urine_rse=0.4,
                   pkml_file="PKML/Olaparib_300mg_BID_Weibull_v2.pkml",
                   OFVcalc_fun = OFVcalc_fun3urine,
                   priors_setup_df=priors_setup_df,
                   fixed_setup_df= data.frame())

endTime <- Sys.time()
# prints recorded time
time3b <- endTime - startTime # 47.2296 mins
print(time3b)

cat("Saving execution time...\n")
time_text <- format(time3b) # format() makes it a nice character string
writeLines(
  time_text,
  con = file.path(output_dir, "execution_time.txt")
)

posterior6 <- read.csv("Results/Anya/CV_2/SequentialMC/3b_N_of_posteriors/Olaparib_PK_posterior6.csv")

cat("Saving posterior distribution plots...\n")
posterior_dist_plot <- PosteriorPlots_fun(posterior6)
# Now, add a title to the plot object
posterior_dist_plot <- posterior_dist_plot + 
  labs(title = "Posterior parameter distributions (Nposterior = 100)")
ggsave(
  filename = file.path(output_dir, "posterior_distributions.pdf"),
  plot = posterior_dist_plot,
  width = 8,
  height = 7
)

cat("Saving correlation matrix...\n")
correlation_matrix <- cor(posterior6[,1:nrow(priors_setup_df)])
write.csv(
  correlation_matrix,
  file = file.path(output_dir, "correlation_matrix.csv")
)

cat("Saving summary statistics...\n")
summary_text <- capture.output(summary(posterior6))
writeLines(
  summary_text,
  con = file.path(output_dir, "posterior_summary.txt")
)


# Visualize PK profile fits
# It takes 10 of the accepted parameter sets, runs simulations with them, and plots the resulting PK profiles on top of the real data. This provides a visual check that the accepted parameters do, in fact, produce simulations that look like the real data.
pk_fit_plot <- simPlot_fun(pkdat_mlx, pkml_file ="PKML/Olaparib_300mg_BID_Weibull_v2.pkml", 
                           outputs_v = outputs_v,
                           priors_setup_df= priors_setup_df, fixed_setup_df= data.frame() ,
                           posterior_set=posterior6[1:10,], 
                           selID="300_Plummer_2015")
cat("Saving PK profile fit plot...\n")
ggsave(
  filename = file.path(output_dir, "pk_profile_fit.pdf"),
  plot = pk_fit_plot,
  width = 7,
  height = 5
)

# Run the metric function
metric_results <- calculate_distribution_metric(ref_df = posterior1, 
                                                comp_df = posterior6)
# View the results
print(paste("Total Difference Score:", metric_results$total_difference_score))
print("Detailed Breakdown:")
detailed_breakdown_3b <- metric_results$detailed_breakdown
print(detailed_breakdown_3b)

cat("Saving metric breakdown...\n")
write.csv(
  detailed_breakdown_3b,
  file = file.path(output_dir, "metric_breakdown.csv"),
  row.names = FALSE # Usually good practice for analysis data frames
)


# Plot with comparison of posterior distributions (rejection from prior vs. SMC)
posterior_dist_plot2 <- PosteriorPlots_fun2(posterior1=posterior1, posterior2= posterior6, appr1="Rejections from prior", appr2="SMC")
cat("Saving posterior distribution plots...\n")
# Now, add a title to the plot object
posterior_dist_plot2 <- posterior_dist_plot2 + 
  labs(title = "Comparison of posterior parameter distributions (Nposterior = 100)")
ggsave(
  filename = file.path(output_dir, "posterior_distributions_comparison.pdf"),
  plot = posterior_dist_plot2,
  width = 8,
  height = 7
)






################### 4. ScCov
# 4a. ScCov = 0.5
output_dir <- "Results/Anya/CV_2/SequentialMC/4a_ScCov/"

iter_info <- data.frame( Iteration =c(1, 2, 3, 4, 5, 6),
                         dosingID=c("300_Plummer_2015"),
                         Nposterior=c(50, 50, 50, 100, 200, 300),
                         Epsilon=c(0.5, 0.3, 0.25, 0.20, 0.18, 0.15),
                         ScCov= c(0.5))

startTime <- Sys.time()
SMC_all_iterations(iter_info=iter_info, Resultsfile="Results/Anya/CV_2/SequentialMC/4a_ScCov/Olaparib_PK_posterior",
                   pkdat_mlx=pkdat_mlx,
                   outputs_v = outputs_v,
                   drug_MW=434, drug_in_urine_part=0.15, amt_urine_rse=0.4,
                   pkml_file="PKML/Olaparib_300mg_BID_Weibull_v2.pkml",
                   OFVcalc_fun = OFVcalc_fun3urine,
                   priors_setup_df=priors_setup_df,
                   fixed_setup_df= data.frame())

endTime <- Sys.time()
# prints recorded time
time4a <- endTime - startTime # 31.45221 mins
print(time4a)

cat("Saving execution time...\n")
time_text <- format(time4a) # format() makes it a nice character string
writeLines(
  time_text,
  con = file.path(output_dir, "execution_time.txt")
)

posterior6 <- read.csv("Results/Anya/CV_2/SequentialMC/4a_ScCov/Olaparib_PK_posterior6.csv")

cat("Saving posterior distribution plots...\n")
posterior_dist_plot <- PosteriorPlots_fun(posterior6)
# Now, add a title to the plot object
posterior_dist_plot <- posterior_dist_plot + 
  labs(title = "Posterior parameter distributions (ScCov = 0.5)")
ggsave(
  filename = file.path(output_dir, "posterior_distributions.pdf"),
  plot = posterior_dist_plot,
  width = 8,
  height = 7
)

cat("Saving correlation matrix...\n")
correlation_matrix <- cor(posterior6[,1:nrow(priors_setup_df)])
write.csv(
  correlation_matrix,
  file = file.path(output_dir, "correlation_matrix.csv")
)

cat("Saving summary statistics...\n")
summary_text <- capture.output(summary(posterior6))
writeLines(
  summary_text,
  con = file.path(output_dir, "posterior_summary.txt")
)


# Visualize PK profile fits
# It takes 10 of the accepted parameter sets, runs simulations with them, and plots the resulting PK profiles on top of the real data. This provides a visual check that the accepted parameters do, in fact, produce simulations that look like the real data.
pk_fit_plot <- simPlot_fun(pkdat_mlx, pkml_file ="PKML/Olaparib_300mg_BID_Weibull_v2.pkml", 
                           outputs_v = outputs_v,
                           priors_setup_df= priors_setup_df, fixed_setup_df= data.frame() ,
                           posterior_set=posterior6[1:10,], 
                           selID="300_Plummer_2015")
cat("Saving PK profile fit plot...\n")
ggsave(
  filename = file.path(output_dir, "pk_profile_fit.pdf"),
  plot = pk_fit_plot,
  width = 7,
  height = 5
)

# Run the metric function
metric_results <- calculate_distribution_metric(ref_df = posterior1, 
                                                comp_df = posterior6)
# View the results
print(paste("Total Difference Score:", metric_results$total_difference_score))
print("Detailed Breakdown:")
detailed_breakdown_4a <- metric_results$detailed_breakdown
print(detailed_breakdown_4a)

cat("Saving metric breakdown...\n")
write.csv(
  detailed_breakdown_4a,
  file = file.path(output_dir, "metric_breakdown.csv"),
  row.names = FALSE # Usually good practice for analysis data frames
)


# Plot with comparison of posterior distributions (rejection from prior vs. SMC)
posterior_dist_plot2 <- PosteriorPlots_fun2(posterior1=posterior1, posterior2= posterior6, appr1="Rejections from prior", appr2="SMC")
cat("Saving posterior distribution plots...\n")
# Now, add a title to the plot object
posterior_dist_plot2 <- posterior_dist_plot2 + 
  labs(title = "Comparison of posterior parameter distributions (ScCov = 0.5)")
ggsave(
  filename = file.path(output_dir, "posterior_distributions_comparison.pdf"),
  plot = posterior_dist_plot2,
  width = 8,
  height = 7
)





# 4b. ScCov = 1.0
output_dir <- "Results/Anya/CV_2/SequentialMC/4b_ScCov/"

iter_info <- data.frame( Iteration =c(1, 2, 3, 4, 5, 6),
                         dosingID=c("300_Plummer_2015"),
                         Nposterior=c(50, 50, 50, 100, 200, 300),
                         Epsilon=c(0.5, 0.3, 0.25, 0.20, 0.18, 0.15),
                         ScCov= c(1.0))

startTime <- Sys.time()
SMC_all_iterations(iter_info=iter_info, Resultsfile="Results/Anya/CV_2/SequentialMC/4b_ScCov/Olaparib_PK_posterior",
                   pkdat_mlx=pkdat_mlx,
                   outputs_v = outputs_v,
                   drug_MW=434, drug_in_urine_part=0.15, amt_urine_rse=0.4,
                   pkml_file="PKML/Olaparib_300mg_BID_Weibull_v2.pkml",
                   OFVcalc_fun = OFVcalc_fun3urine,
                   priors_setup_df=priors_setup_df,
                   fixed_setup_df= data.frame())

endTime <- Sys.time()
# prints recorded time
time4b <- endTime - startTime # 34.84068 mins
print(time4b)

cat("Saving execution time...\n")
time_text <- format(time4b) # format() makes it a nice character string
writeLines(
  time_text,
  con = file.path(output_dir, "execution_time.txt")
)

posterior6 <- read.csv("Results/Anya/CV_2/SequentialMC/4b_ScCov/Olaparib_PK_posterior6.csv")

cat("Saving posterior distribution plots...\n")
posterior_dist_plot <- PosteriorPlots_fun(posterior6)
# Now, add a title to the plot object
posterior_dist_plot <- posterior_dist_plot + 
  labs(title = "Posterior parameter distributions (ScCov = 1.0)")
ggsave(
  filename = file.path(output_dir, "posterior_distributions.pdf"),
  plot = posterior_dist_plot,
  width = 8,
  height = 7
)

cat("Saving correlation matrix...\n")
correlation_matrix <- cor(posterior6[,1:nrow(priors_setup_df)])
write.csv(
  correlation_matrix,
  file = file.path(output_dir, "correlation_matrix.csv")
)

cat("Saving summary statistics...\n")
summary_text <- capture.output(summary(posterior6))
writeLines(
  summary_text,
  con = file.path(output_dir, "posterior_summary.txt")
)


# Visualize PK profile fits
# It takes 10 of the accepted parameter sets, runs simulations with them, and plots the resulting PK profiles on top of the real data. This provides a visual check that the accepted parameters do, in fact, produce simulations that look like the real data.
pk_fit_plot <- simPlot_fun(pkdat_mlx, pkml_file ="PKML/Olaparib_300mg_BID_Weibull_v2.pkml", 
                           outputs_v = outputs_v,
                           priors_setup_df= priors_setup_df, fixed_setup_df= data.frame() ,
                           posterior_set=posterior6[1:10,], 
                           selID="300_Plummer_2015")
cat("Saving PK profile fit plot...\n")
ggsave(
  filename = file.path(output_dir, "pk_profile_fit.pdf"),
  plot = pk_fit_plot,
  width = 7,
  height = 5
)

# Run the metric function
metric_results <- calculate_distribution_metric(ref_df = posterior1, 
                                                comp_df = posterior6)
# View the results
print(paste("Total Difference Score:", metric_results$total_difference_score))
print("Detailed Breakdown:")
detailed_breakdown_4b <- metric_results$detailed_breakdown
print(detailed_breakdown_4b)

cat("Saving metric breakdown...\n")
write.csv(
  detailed_breakdown_4b,
  file = file.path(output_dir, "metric_breakdown.csv"),
  row.names = FALSE # Usually good practice for analysis data frames
)


# Plot with comparison of posterior distributions (rejection from prior vs. SMC)
posterior_dist_plot2 <- PosteriorPlots_fun2(posterior1=posterior1, posterior2= posterior6, appr1="Rejections from prior", appr2="SMC")
cat("Saving posterior distribution plots...\n")
# Now, add a title to the plot object
posterior_dist_plot2 <- posterior_dist_plot2 + 
  labs(title = "Comparison of posterior parameter distributions (ScCov = 1.0)")
ggsave(
  filename = file.path(output_dir, "posterior_distributions_comparison.pdf"),
  plot = posterior_dist_plot2,
  width = 8,
  height = 7
)






################### 5. 400_Gao_2023
# 5a. Base metrics 

pkdat_mlx %>% filter(ID=="400_Gao_2023")

output_dir <- "Results/Anya/CV_2/SequentialMC/5a_Gao_base/"

iter_info <- data.frame( Iteration =c(1, 2, 3, 4, 5, 6),
                         dosingID=c("400_Gao_2023"),
                         Nposterior=c(50, 50, 50, 100, 200, 300),
                         Epsilon=c(0.5, 0.3, 0.25, 0.20, 0.18, 0.15),
                         ScCov= c(0.25))

startTime <- Sys.time()
SMC_all_iterations(iter_info=iter_info, Resultsfile="Results/Anya/CV_2/SequentialMC/5a_Gao_base/Olaparib_PK_posterior",
                   pkdat_mlx=pkdat_mlx,
                   outputs_v = outputs_v,
                   drug_MW=434, drug_in_urine_part=0.15, amt_urine_rse=0.4,
                   pkml_file="PKML/Olaparib_300mg_BID_Weibull_v2.pkml",
                   OFVcalc_fun = OFVcalc_fun3urine,
                   priors_setup_df=priors_setup_df,
                   fixed_setup_df= data.frame())

endTime <- Sys.time()
# prints recorded time
time5a <- endTime - startTime # 49.15525 mins
print(time5a)

cat("Saving execution time...\n")
time_text <- format(time5a) # format() makes it a nice character string
writeLines(
  time_text,
  con = file.path(output_dir, "execution_time.txt")
)

posterior6 <- read.csv("Results/Anya/CV_2/SequentialMC/5a_Gao_base/Olaparib_PK_posterior6.csv")

cat("Saving posterior distribution plots...\n")
posterior_dist_plot <- PosteriorPlots_fun(posterior6)
# Now, add a title to the plot object
posterior_dist_plot <- posterior_dist_plot + 
  labs(title = "Posterior parameter distributions (Gao_base metrics)")
ggsave(
  filename = file.path(output_dir, "posterior_distributions.pdf"),
  plot = posterior_dist_plot,
  width = 8,
  height = 7
)

cat("Saving correlation matrix...\n")
correlation_matrix <- cor(posterior6[,1:nrow(priors_setup_df)])
write.csv(
  correlation_matrix,
  file = file.path(output_dir, "correlation_matrix.csv")
)

cat("Saving summary statistics...\n")
summary_text <- capture.output(summary(posterior6))
writeLines(
  summary_text,
  con = file.path(output_dir, "posterior_summary.txt")
)


# Visualize PK profile fits
# It takes 10 of the accepted parameter sets, runs simulations with them, and plots the resulting PK profiles on top of the real data. This provides a visual check that the accepted parameters do, in fact, produce simulations that look like the real data.
pk_fit_plot <- simPlot_fun(pkdat_mlx, pkml_file ="PKML/Olaparib_300mg_BID_Weibull_v2.pkml", 
                           outputs_v = outputs_v,
                           priors_setup_df= priors_setup_df, fixed_setup_df= data.frame() ,
                           posterior_set=posterior6[1:10,], 
                           selID="400_Gao_2023")
cat("Saving PK profile fit plot...\n")
ggsave(
  filename = file.path(output_dir, "pk_profile_fit.pdf"),
  plot = pk_fit_plot,
  width = 7,
  height = 5
)

# Run the metric function
metric_results <- calculate_distribution_metric(ref_df = posterior1, 
                                                comp_df = posterior6)
# View the results
print(paste("Total Difference Score:", metric_results$total_difference_score))
print("Detailed Breakdown:")
detailed_breakdown_5a <- metric_results$detailed_breakdown
print(detailed_breakdown_5a)

cat("Saving metric breakdown...\n")
write.csv(
  detailed_breakdown_5a,
  file = file.path(output_dir, "metric_breakdown.csv"),
  row.names = FALSE # Usually good practice for analysis data frames
)


# Plot with comparison of posterior distributions (rejection from prior vs. SMC)
posterior_dist_plot2 <- PosteriorPlots_fun2(posterior1=posterior1, posterior2= posterior6, appr1="Rejections from prior", appr2="SMC")
cat("Saving posterior distribution plots...\n")
# Now, add a title to the plot object
posterior_dist_plot2 <- posterior_dist_plot2 + 
  labs(title = "Comparison of posterior parameter distributions (Gao_base metrics)")
ggsave(
  filename = file.path(output_dir, "posterior_distributions_comparison.pdf"),
  plot = posterior_dist_plot2,
  width = 8,
  height = 7
)



# 5b. epsilon = 0.3

output_dir <- "Results/Anya/CV_2/SequentialMC/5b_Gao_epsilon/"

iter_info <- data.frame( Iteration =c(1, 2, 3, 4, 5, 6),
                         dosingID=c("400_Gao_2023"),
                         Nposterior=c(50, 50, 50, 100, 200, 300),
                         Epsilon=c(0.3, 0.25, 0.2, 0.18, 0.16, 0.15),
                         ScCov= c(0.25))

startTime <- Sys.time()
SMC_all_iterations(iter_info=iter_info, Resultsfile="Results/Anya/CV_2/SequentialMC/5b_Gao_epsilon/Olaparib_PK_posterior",
                   pkdat_mlx=pkdat_mlx,
                   outputs_v = outputs_v,
                   drug_MW=434, drug_in_urine_part=0.15, amt_urine_rse=0.4,
                   pkml_file="PKML/Olaparib_300mg_BID_Weibull_v2.pkml",
                   OFVcalc_fun = OFVcalc_fun3urine,
                   priors_setup_df=priors_setup_df,
                   fixed_setup_df= data.frame())

endTime <- Sys.time()
# prints recorded time
time5b <- endTime - startTime # 31.77514 mins
print(time5b)

cat("Saving execution time...\n")
time_text <- format(time5b) # format() makes it a nice character string
writeLines(
  time_text,
  con = file.path(output_dir, "execution_time.txt")
)

posterior6 <- read.csv("Results/Anya/CV_2/SequentialMC/5b_Gao_epsilon/Olaparib_PK_posterior6.csv")

cat("Saving posterior distribution plots...\n")
posterior_dist_plot <- PosteriorPlots_fun(posterior6)
# Now, add a title to the plot object
posterior_dist_plot <- posterior_dist_plot + 
  labs(title = "Posterior parameter distributions (Ga0_epsilon = 0.3)")
ggsave(
  filename = file.path(output_dir, "posterior_distributions.pdf"),
  plot = posterior_dist_plot,
  width = 8,
  height = 7
)

cat("Saving correlation matrix...\n")
correlation_matrix <- cor(posterior6[,1:nrow(priors_setup_df)])
write.csv(
  correlation_matrix,
  file = file.path(output_dir, "correlation_matrix.csv")
)

cat("Saving summary statistics...\n")
summary_text <- capture.output(summary(posterior6))
writeLines(
  summary_text,
  con = file.path(output_dir, "posterior_summary.txt")
)

# Visualize PK profile fits
# It takes 10 of the accepted parameter sets, runs simulations with them, and plots the resulting PK profiles on top of the real data. This provides a visual check that the accepted parameters do, in fact, produce simulations that look like the real data.
pk_fit_plot <- simPlot_fun(pkdat_mlx, pkml_file ="PKML/Olaparib_300mg_BID_Weibull_v2.pkml", 
                           outputs_v = outputs_v,
                           priors_setup_df= priors_setup_df, fixed_setup_df= data.frame() ,
                           posterior_set=posterior6[1:10,], 
                           selID="400_Gao_2023")
cat("Saving PK profile fit plot...\n")
ggsave(
  filename = file.path(output_dir, "pk_profile_fit.pdf"),
  plot = pk_fit_plot,
  width = 7,
  height = 5
)

# Run the metric function
metric_results <- calculate_distribution_metric(ref_df = posterior1, 
                                                comp_df = posterior6)
# View the results
print(paste("Total Difference Score:", metric_results$total_difference_score))
print("Detailed Breakdown:")
detailed_breakdown_5b <- metric_results$detailed_breakdown
print(detailed_breakdown_5b)

cat("Saving metric breakdown...\n")
write.csv(
  detailed_breakdown_5b,
  file = file.path(output_dir, "metric_breakdown.csv"),
  row.names = FALSE # Usually good practice for analysis data frames
)


# Plot with comparison of posterior distributions (rejection from prior vs. SMC)
posterior_dist_plot2 <- PosteriorPlots_fun2(posterior1=posterior1, posterior2= posterior6, appr1="Rejections from prior", appr2="SMC")
cat("Saving posterior distribution plots...\n")
# Now, add a title to the plot object
posterior_dist_plot2 <- posterior_dist_plot2 + 
  labs(title = "Comparison of posterior parameter distributions (Gao_epsilon = 0.3)")
ggsave(
  filename = file.path(output_dir, "posterior_distributions_comparison.pdf"),
  plot = posterior_dist_plot2,
  width = 8,
  height = 7
)






################### Make the summary table with parameter optimization results
exp_num <- c("0", "1a", "1b", "2a", "2b", "3a", "3b", "4a", "4b", "5a", "5b")

observed_data <- c(
  "300_Plummer_2015", "300_Plummer_2015", "300_Plummer_2015", "300_Plummer_2015", 
  "300_Plummer_2015", "300_Plummer_2015", "300_Plummer_2015", "300_Plummer_2015", 
  "300_Plummer_2015", "400_Gao_2023", "400_Gao_2023"
)

tested_params <- c(
  "base metrics (Niterations = 6, Start_epsilon = 0.5, Nposterior = 50, ScCov = 0.25)",
  "Niterations = 4", "Niterations = 9", "Start_epsilon = 0.3", "Start_epsilon = 1.0",
  "Nposterior = 20", "Nposterior = 100", "ScCov = 0.5", "ScCov = 1.0",
  "base metrics", "Start_epsilon = 0.3"
)

similarity_metric <- c(4, 3, 3, 1, 4, 7, 5, 4, 8, 5, 5)

execution_time <- c(
  35.03735, 26.65571, 32.41298, 24.8635, 24.77847, 31.74183, 
  47.2296, 31.45221, 34.84068, 49.15525, 31.77514
)

# 2. Create the data frame.
summary_table <- data.frame(
  `id` = exp_num,
  `Observed data` = observed_data,
  `Tested parameters` = tested_params,
  `Distribution similarity metric` = similarity_metric,
  `Execution time (min)` = execution_time,
  stringsAsFactors = FALSE # Good practice to avoid factors unless needed
)

# 3. Print the table to the console to verify it looks correct.
print(summary_table)

# 4. Define the output path and filename
output_path <- "Results/Anya/CV_2/SequentialMC/"
output_filename <- "final_experiment_summary.csv"

# 5. Save the final table to a CSV file
write.csv(
  summary_table,
  file = file.path(output_path, output_filename),
  row.names = FALSE  # Crucial: Prevents R from adding an extra column with its own row numbers
)


script_path <- rstudioapi::getActiveDocumentContext()$path
full_script_path <- path.expand(script_path)
full_script_path









##############################################################################################
############################    Result Reproduction Test    ##################################



# Let's try to reproduce the results of 0 and 2a starts. We did not use random seed, so the results should be different even if we use the same input parameters
################### 0_2. Base metrics 

output_dir <- "Results/Anya/CV_2/SequentialMC/0_base_metrics_2/"

iter_info <- data.frame( Iteration =c(1, 2, 3, 4, 5, 6),
                         dosingID=c("300_Plummer_2015"),
                         Nposterior=c(50, 50, 50, 100, 200, 300),
                         Epsilon=c(0.5, 0.3, 0.25, 0.20, 0.18, 0.15),
                         ScCov= c(0.25))

startTime <- Sys.time()
SMC_all_iterations(iter_info=iter_info, Resultsfile="Results/Anya/CV_2/SequentialMC/0_base_metrics_2/Olaparib_PK_posterior",
                   pkdat_mlx=pkdat_mlx,
                   outputs_v = outputs_v,
                   drug_MW=434, drug_in_urine_part=0.15, amt_urine_rse=0.4,
                   pkml_file="PKML/Olaparib_300mg_BID_Weibull_v2.pkml",
                   OFVcalc_fun = OFVcalc_fun3urine,
                   priors_setup_df=priors_setup_df,
                   fixed_setup_df= data.frame())

endTime <- Sys.time()
# prints recorded time
time0 <- endTime - startTime # 38.01595 mins
print(time0)

cat("Saving execution time...\n")
time_text <- format(time0) # format() makes it a nice character string
writeLines(
  time_text,
  con = file.path(output_dir, "execution_time.txt")
)

posterior6 <- read.csv("Results/Anya/CV_2/SequentialMC/0_base_metrics_2/Olaparib_PK_posterior6.csv")

cat("Saving posterior distribution plots...\n")
posterior_dist_plot <- PosteriorPlots_fun(posterior6)
# Now, add a title to the plot object
posterior_dist_plot <- posterior_dist_plot + 
  labs(title = "Posterior parameter distributions (base metrics)")
ggsave(
  filename = file.path(output_dir, "posterior_distributions.pdf"),
  plot = posterior_dist_plot,
  width = 8,
  height = 7
)

cat("Saving correlation matrix...\n")
correlation_matrix <- cor(posterior6[,1:nrow(priors_setup_df)])
write.csv(
  correlation_matrix,
  file = file.path(output_dir, "correlation_matrix.csv")
)

cat("Saving summary statistics...\n")
summary_text <- capture.output(summary(posterior6))
writeLines(
  summary_text,
  con = file.path(output_dir, "posterior_summary.txt")
)


# Visualize PK profile fits
# It takes 10 of the accepted parameter sets, runs simulations with them, and plots the resulting PK profiles on top of the real data. This provides a visual check that the accepted parameters do, in fact, produce simulations that look like the real data.
pk_fit_plot <- simPlot_fun(pkdat_mlx, pkml_file ="PKML/Olaparib_300mg_BID_Weibull_v2.pkml", 
                           outputs_v = outputs_v,
                           priors_setup_df= priors_setup_df, fixed_setup_df= data.frame() ,
                           posterior_set=posterior6[1:10,], 
                           selID="300_Plummer_2015")
cat("Saving PK profile fit plot...\n")
ggsave(
  filename = file.path(output_dir, "pk_profile_fit.pdf"),
  plot = pk_fit_plot,
  width = 7,
  height = 5
)

# Run the metric function
metric_results <- calculate_distribution_metric(ref_df = posterior1, 
                                                comp_df = posterior6)
# View the results
print(paste("Total Difference Score:", metric_results$total_difference_score))
print("Detailed Breakdown:")
detailed_breakdown_0 <- metric_results$detailed_breakdown
print(detailed_breakdown_0)

cat("Saving metric breakdown...\n")
write.csv(
  detailed_breakdown_0,
  file = file.path(output_dir, "metric_breakdown.csv"),
  row.names = FALSE # Usually good practice for analysis data frames
)


# Plot with comparison of posterior distributions (rejection from prior vs. SMC)
posterior_dist_plot2 <- PosteriorPlots_fun2(posterior1=posterior1, posterior2= posterior6, appr1="Rejections from prior", appr2="SMC")
cat("Saving posterior distribution plots...\n")
# Now, add a title to the plot object
posterior_dist_plot2 <- posterior_dist_plot2 + 
  labs(title = "Comparison of posterior parameter distributions (base metrics)")
ggsave(
  filename = file.path(output_dir, "posterior_distributions_comparison.pdf"),
  plot = posterior_dist_plot2,
  width = 8,
  height = 7
)




################### 0_3. Base metrics 

output_dir <- "Results/Anya/CV_2/SequentialMC/0_base_metrics_3/"

iter_info <- data.frame( Iteration =c(1, 2, 3, 4, 5, 6),
                         dosingID=c("300_Plummer_2015"),
                         Nposterior=c(50, 50, 50, 100, 200, 300),
                         Epsilon=c(0.5, 0.3, 0.25, 0.20, 0.18, 0.15),
                         ScCov= c(0.25))

startTime <- Sys.time()
SMC_all_iterations(iter_info=iter_info, Resultsfile="Results/Anya/CV_2/SequentialMC/0_base_metrics_3/Olaparib_PK_posterior",
                   pkdat_mlx=pkdat_mlx,
                   outputs_v = outputs_v,
                   drug_MW=434, drug_in_urine_part=0.15, amt_urine_rse=0.4,
                   pkml_file="PKML/Olaparib_300mg_BID_Weibull_v2.pkml",
                   OFVcalc_fun = OFVcalc_fun3urine,
                   priors_setup_df=priors_setup_df,
                   fixed_setup_df= data.frame())

endTime <- Sys.time()
# prints recorded time
time0 <- endTime - startTime # 26.72996 mins
print(time0)

cat("Saving execution time...\n")
time_text <- format(time0) # format() makes it a nice character string
writeLines(
  time_text,
  con = file.path(output_dir, "execution_time.txt")
)

posterior6 <- read.csv("Results/Anya/CV_2/SequentialMC/0_base_metrics_3/Olaparib_PK_posterior6.csv")

cat("Saving posterior distribution plots...\n")
posterior_dist_plot <- PosteriorPlots_fun(posterior6)
# Now, add a title to the plot object
posterior_dist_plot <- posterior_dist_plot + 
  labs(title = "Posterior parameter distributions (base metrics)")
ggsave(
  filename = file.path(output_dir, "posterior_distributions.pdf"),
  plot = posterior_dist_plot,
  width = 8,
  height = 7
)

cat("Saving correlation matrix...\n")
correlation_matrix <- cor(posterior6[,1:nrow(priors_setup_df)])
write.csv(
  correlation_matrix,
  file = file.path(output_dir, "correlation_matrix.csv")
)

cat("Saving summary statistics...\n")
summary_text <- capture.output(summary(posterior6))
writeLines(
  summary_text,
  con = file.path(output_dir, "posterior_summary.txt")
)


# Visualize PK profile fits
# It takes 10 of the accepted parameter sets, runs simulations with them, and plots the resulting PK profiles on top of the real data. This provides a visual check that the accepted parameters do, in fact, produce simulations that look like the real data.
pk_fit_plot <- simPlot_fun(pkdat_mlx, pkml_file ="PKML/Olaparib_300mg_BID_Weibull_v2.pkml", 
                           outputs_v = outputs_v,
                           priors_setup_df= priors_setup_df, fixed_setup_df= data.frame() ,
                           posterior_set=posterior6[1:10,], 
                           selID="300_Plummer_2015")
cat("Saving PK profile fit plot...\n")
ggsave(
  filename = file.path(output_dir, "pk_profile_fit.pdf"),
  plot = pk_fit_plot,
  width = 7,
  height = 5
)

# Run the metric function
metric_results <- calculate_distribution_metric(ref_df = posterior1, 
                                                comp_df = posterior6)
# View the results
print(paste("Total Difference Score:", metric_results$total_difference_score))
print("Detailed Breakdown:")
detailed_breakdown_0 <- metric_results$detailed_breakdown
print(detailed_breakdown_0)

cat("Saving metric breakdown...\n")
write.csv(
  detailed_breakdown_0,
  file = file.path(output_dir, "metric_breakdown.csv"),
  row.names = FALSE # Usually good practice for analysis data frames
)


# Plot with comparison of posterior distributions (rejection from prior vs. SMC)
posterior_dist_plot2 <- PosteriorPlots_fun2(posterior1=posterior1, posterior2= posterior6, appr1="Rejections from prior", appr2="SMC")
cat("Saving posterior distribution plots...\n")
# Now, add a title to the plot object
posterior_dist_plot2 <- posterior_dist_plot2 + 
  labs(title = "Comparison of posterior parameter distributions (base metrics)")
ggsave(
  filename = file.path(output_dir, "posterior_distributions_comparison.pdf"),
  plot = posterior_dist_plot2,
  width = 8,
  height = 7
)




################### 0_4. Base metrics 

output_dir <- "Results/Anya/CV_2/SequentialMC/0_base_metrics_4/"

iter_info <- data.frame( Iteration =c(1, 2, 3, 4, 5, 6),
                         dosingID=c("300_Plummer_2015"),
                         Nposterior=c(50, 50, 50, 100, 200, 300),
                         Epsilon=c(0.5, 0.3, 0.25, 0.20, 0.18, 0.15),
                         ScCov= c(0.25))

startTime <- Sys.time()
SMC_all_iterations(iter_info=iter_info, Resultsfile="Results/Anya/CV_2/SequentialMC/0_base_metrics_4/Olaparib_PK_posterior",
                   pkdat_mlx=pkdat_mlx,
                   outputs_v = outputs_v,
                   drug_MW=434, drug_in_urine_part=0.15, amt_urine_rse=0.4,
                   pkml_file="PKML/Olaparib_300mg_BID_Weibull_v2.pkml",
                   OFVcalc_fun = OFVcalc_fun3urine,
                   priors_setup_df=priors_setup_df,
                   fixed_setup_df= data.frame())

endTime <- Sys.time()
# prints recorded time
time0 <- endTime - startTime # 36.48316 mins
print(time0)

cat("Saving execution time...\n")
time_text <- format(time0) # format() makes it a nice character string
writeLines(
  time_text,
  con = file.path(output_dir, "execution_time.txt")
)

posterior6 <- read.csv("Results/Anya/CV_2/SequentialMC/0_base_metrics_4/Olaparib_PK_posterior6.csv")

cat("Saving posterior distribution plots...\n")
posterior_dist_plot <- PosteriorPlots_fun(posterior6)
# Now, add a title to the plot object
posterior_dist_plot <- posterior_dist_plot + 
  labs(title = "Posterior parameter distributions (base metrics)")
ggsave(
  filename = file.path(output_dir, "posterior_distributions.pdf"),
  plot = posterior_dist_plot,
  width = 8,
  height = 7
)

cat("Saving correlation matrix...\n")
correlation_matrix <- cor(posterior6[,1:nrow(priors_setup_df)])
write.csv(
  correlation_matrix,
  file = file.path(output_dir, "correlation_matrix.csv")
)

cat("Saving summary statistics...\n")
summary_text <- capture.output(summary(posterior6))
writeLines(
  summary_text,
  con = file.path(output_dir, "posterior_summary.txt")
)


# Visualize PK profile fits
# It takes 10 of the accepted parameter sets, runs simulations with them, and plots the resulting PK profiles on top of the real data. This provides a visual check that the accepted parameters do, in fact, produce simulations that look like the real data.
pk_fit_plot <- simPlot_fun(pkdat_mlx, pkml_file ="PKML/Olaparib_300mg_BID_Weibull_v2.pkml", 
                           outputs_v = outputs_v,
                           priors_setup_df= priors_setup_df, fixed_setup_df= data.frame() ,
                           posterior_set=posterior6[1:10,], 
                           selID="300_Plummer_2015")
cat("Saving PK profile fit plot...\n")
ggsave(
  filename = file.path(output_dir, "pk_profile_fit.pdf"),
  plot = pk_fit_plot,
  width = 7,
  height = 5
)

# Run the metric function
metric_results <- calculate_distribution_metric(ref_df = posterior1, 
                                                comp_df = posterior6)
# View the results
print(paste("Total Difference Score:", metric_results$total_difference_score))
print("Detailed Breakdown:")
detailed_breakdown_0 <- metric_results$detailed_breakdown
print(detailed_breakdown_0)

cat("Saving metric breakdown...\n")
write.csv(
  detailed_breakdown_0,
  file = file.path(output_dir, "metric_breakdown.csv"),
  row.names = FALSE # Usually good practice for analysis data frames
)


# Plot with comparison of posterior distributions (rejection from prior vs. SMC)
posterior_dist_plot2 <- PosteriorPlots_fun2(posterior1=posterior1, posterior2= posterior6, appr1="Rejections from prior", appr2="SMC")
cat("Saving posterior distribution plots...\n")
# Now, add a title to the plot object
posterior_dist_plot2 <- posterior_dist_plot2 + 
  labs(title = "Comparison of posterior parameter distributions (base metrics)")
ggsave(
  filename = file.path(output_dir, "posterior_distributions_comparison.pdf"),
  plot = posterior_dist_plot2,
  width = 8,
  height = 7
)




################### 2a_2. epsilon = 0.3

output_dir <- "Results/Anya/CV_2/SequentialMC/2a_epsilon_2/"

iter_info <- data.frame( Iteration =c(1, 2, 3, 4, 5, 6),
                         dosingID=c("300_Plummer_2015"),
                         Nposterior=c(50, 50, 50, 100, 200, 300),
                         Epsilon=c(0.3, 0.25, 0.2, 0.18, 0.16, 0.15),
                         ScCov= c(0.25))

startTime <- Sys.time()
SMC_all_iterations(iter_info=iter_info, Resultsfile="Results/Anya/CV_2/SequentialMC/2a_epsilon_2/Olaparib_PK_posterior",
                   pkdat_mlx=pkdat_mlx,
                   outputs_v = outputs_v,
                   drug_MW=434, drug_in_urine_part=0.15, amt_urine_rse=0.4,
                   pkml_file="PKML/Olaparib_300mg_BID_Weibull_v2.pkml",
                   OFVcalc_fun = OFVcalc_fun3urine,
                   priors_setup_df=priors_setup_df,
                   fixed_setup_df= data.frame())

endTime <- Sys.time()
# prints recorded time
time2a <- endTime - startTime # 23.49269 mins
print(time2a)

cat("Saving execution time...\n")
time_text <- format(time2a) # format() makes it a nice character string
writeLines(
  time_text,
  con = file.path(output_dir, "execution_time.txt")
)

posterior6 <- read.csv("Results/Anya/CV_2/SequentialMC/2a_epsilon_2/Olaparib_PK_posterior6.csv")

cat("Saving posterior distribution plots...\n")
posterior_dist_plot <- PosteriorPlots_fun(posterior6)
# Now, add a title to the plot object
posterior_dist_plot <- posterior_dist_plot + 
  labs(title = "Posterior parameter distributions (epsilon = 0.3)")
ggsave(
  filename = file.path(output_dir, "posterior_distributions.pdf"),
  plot = posterior_dist_plot,
  width = 8,
  height = 7
)

cat("Saving correlation matrix...\n")
correlation_matrix <- cor(posterior6[,1:nrow(priors_setup_df)])
write.csv(
  correlation_matrix,
  file = file.path(output_dir, "correlation_matrix.csv")
)

cat("Saving summary statistics...\n")
summary_text <- capture.output(summary(posterior6))
writeLines(
  summary_text,
  con = file.path(output_dir, "posterior_summary.txt")
)

# Visualize PK profile fits
# It takes 10 of the accepted parameter sets, runs simulations with them, and plots the resulting PK profiles on top of the real data. This provides a visual check that the accepted parameters do, in fact, produce simulations that look like the real data.
pk_fit_plot <- simPlot_fun(pkdat_mlx, pkml_file ="PKML/Olaparib_300mg_BID_Weibull_v2.pkml", 
                           outputs_v = outputs_v,
                           priors_setup_df= priors_setup_df, fixed_setup_df= data.frame() ,
                           posterior_set=posterior6[1:10,], 
                           selID="300_Plummer_2015")
cat("Saving PK profile fit plot...\n")
ggsave(
  filename = file.path(output_dir, "pk_profile_fit.pdf"),
  plot = pk_fit_plot,
  width = 7,
  height = 5
)

# Run the metric function
metric_results <- calculate_distribution_metric(ref_df = posterior1, 
                                                comp_df = posterior6)
# View the results
print(paste("Total Difference Score:", metric_results$total_difference_score))
print("Detailed Breakdown:")
detailed_breakdown_2a <- metric_results$detailed_breakdown
print(detailed_breakdown_2a)

cat("Saving metric breakdown...\n")
write.csv(
  detailed_breakdown_2a,
  file = file.path(output_dir, "metric_breakdown.csv"),
  row.names = FALSE # Usually good practice for analysis data frames
)

# Plot with comparison of posterior distributions (rejection from prior vs. SMC)
posterior_dist_plot2 <- PosteriorPlots_fun2(posterior1=posterior1, posterior2= posterior6, appr1="Rejections from prior", appr2="SMC")
cat("Saving posterior distribution plots...\n")
# Now, add a title to the plot object
posterior_dist_plot2 <- posterior_dist_plot2 + 
  labs(title = "Comparison of posterior parameter distributions (epsilon = 0.3)")
ggsave(
  filename = file.path(output_dir, "posterior_distributions_comparison.pdf"),
  plot = posterior_dist_plot2,
  width = 8,
  height = 7
)



################### 2a_3. epsilon = 0.3

output_dir <- "Results/Anya/CV_2/SequentialMC/2a_epsilon_3/"

iter_info <- data.frame( Iteration =c(1, 2, 3, 4, 5, 6),
                         dosingID=c("300_Plummer_2015"),
                         Nposterior=c(50, 50, 50, 100, 200, 300),
                         Epsilon=c(0.3, 0.25, 0.2, 0.18, 0.16, 0.15),
                         ScCov= c(0.25))

startTime <- Sys.time()
SMC_all_iterations(iter_info=iter_info, Resultsfile="Results/Anya/CV_2/SequentialMC/2a_epsilon_3/Olaparib_PK_posterior",
                   pkdat_mlx=pkdat_mlx,
                   outputs_v = outputs_v,
                   drug_MW=434, drug_in_urine_part=0.15, amt_urine_rse=0.4,
                   pkml_file="PKML/Olaparib_300mg_BID_Weibull_v2.pkml",
                   OFVcalc_fun = OFVcalc_fun3urine,
                   priors_setup_df=priors_setup_df,
                   fixed_setup_df= data.frame())

endTime <- Sys.time()
# prints recorded time
time2a <- endTime - startTime # 23.75158 mins
print(time2a)

cat("Saving execution time...\n")
time_text <- format(time2a) # format() makes it a nice character string
writeLines(
  time_text,
  con = file.path(output_dir, "execution_time.txt")
)

posterior6 <- read.csv("Results/Anya/CV_2/SequentialMC/2a_epsilon_3/Olaparib_PK_posterior6.csv")

cat("Saving posterior distribution plots...\n")
posterior_dist_plot <- PosteriorPlots_fun(posterior6)
# Now, add a title to the plot object
posterior_dist_plot <- posterior_dist_plot + 
  labs(title = "Posterior parameter distributions (epsilon = 0.3)")
ggsave(
  filename = file.path(output_dir, "posterior_distributions.pdf"),
  plot = posterior_dist_plot,
  width = 8,
  height = 7
)

cat("Saving correlation matrix...\n")
correlation_matrix <- cor(posterior6[,1:nrow(priors_setup_df)])
write.csv(
  correlation_matrix,
  file = file.path(output_dir, "correlation_matrix.csv")
)

cat("Saving summary statistics...\n")
summary_text <- capture.output(summary(posterior6))
writeLines(
  summary_text,
  con = file.path(output_dir, "posterior_summary.txt")
)

# Visualize PK profile fits
# It takes 10 of the accepted parameter sets, runs simulations with them, and plots the resulting PK profiles on top of the real data. This provides a visual check that the accepted parameters do, in fact, produce simulations that look like the real data.
pk_fit_plot <- simPlot_fun(pkdat_mlx, pkml_file ="PKML/Olaparib_300mg_BID_Weibull_v2.pkml", 
                           outputs_v = outputs_v,
                           priors_setup_df= priors_setup_df, fixed_setup_df= data.frame() ,
                           posterior_set=posterior6[1:10,], 
                           selID="300_Plummer_2015")
cat("Saving PK profile fit plot...\n")
ggsave(
  filename = file.path(output_dir, "pk_profile_fit.pdf"),
  plot = pk_fit_plot,
  width = 7,
  height = 5
)

# Run the metric function
metric_results <- calculate_distribution_metric(ref_df = posterior1, 
                                                comp_df = posterior6)

# View the results
print(paste("Total Difference Score:", metric_results$total_difference_score))
print("Detailed Breakdown:")
detailed_breakdown_2a <- metric_results$detailed_breakdown
print(detailed_breakdown_2a)

cat("Saving metric breakdown...\n")
write.csv(
  detailed_breakdown_2a,
  file = file.path(output_dir, "metric_breakdown.csv"),
  row.names = FALSE # Usually good practice for analysis data frames
)

# Plot with comparison of posterior distributions (rejection from prior vs. SMC)
posterior_dist_plot2 <- PosteriorPlots_fun2(posterior1=posterior1, posterior2= posterior6, appr1="Rejections from prior", appr2="SMC")
cat("Saving posterior distribution plots...\n")
# Now, add a title to the plot object
posterior_dist_plot2 <- posterior_dist_plot2 + 
  labs(title = "Comparison of posterior parameter distributions (epsilon = 0.3)")
ggsave(
  filename = file.path(output_dir, "posterior_distributions_comparison.pdf"),
  plot = posterior_dist_plot2,
  width = 8,
  height = 7
)




################### 2a_4. epsilon = 0.3

output_dir <- "Results/Anya/CV_2/SequentialMC/2a_epsilon_4/"

iter_info <- data.frame( Iteration =c(1, 2, 3, 4, 5, 6),
                         dosingID=c("300_Plummer_2015"),
                         Nposterior=c(50, 50, 50, 100, 200, 300),
                         Epsilon=c(0.3, 0.25, 0.2, 0.18, 0.16, 0.15),
                         ScCov= c(0.25))

startTime <- Sys.time()
SMC_all_iterations(iter_info=iter_info, Resultsfile="Results/Anya/CV_2/SequentialMC/2a_epsilon_4/Olaparib_PK_posterior",
                   pkdat_mlx=pkdat_mlx,
                   outputs_v = outputs_v,
                   drug_MW=434, drug_in_urine_part=0.15, amt_urine_rse=0.4,
                   pkml_file="PKML/Olaparib_300mg_BID_Weibull_v2.pkml",
                   OFVcalc_fun = OFVcalc_fun3urine,
                   priors_setup_df=priors_setup_df,
                   fixed_setup_df= data.frame())

endTime <- Sys.time()
# prints recorded time
time2a <- endTime - startTime # 49.55859 mins
print(time2a)

cat("Saving execution time...\n")
time_text <- format(time2a) # format() makes it a nice character string
writeLines(
  time_text,
  con = file.path(output_dir, "execution_time.txt")
)

posterior6 <- read.csv("Results/Anya/CV_2/SequentialMC/2a_epsilon_4/Olaparib_PK_posterior6.csv")

cat("Saving posterior distribution plots...\n")
posterior_dist_plot <- PosteriorPlots_fun(posterior6)
# Now, add a title to the plot object
posterior_dist_plot <- posterior_dist_plot + 
  labs(title = "Posterior parameter distributions (epsilon = 0.3)")
ggsave(
  filename = file.path(output_dir, "posterior_distributions.pdf"),
  plot = posterior_dist_plot,
  width = 8,
  height = 7
)

cat("Saving correlation matrix...\n")
correlation_matrix <- cor(posterior6[,1:nrow(priors_setup_df)])
write.csv(
  correlation_matrix,
  file = file.path(output_dir, "correlation_matrix.csv")
)

cat("Saving summary statistics...\n")
summary_text <- capture.output(summary(posterior6))
writeLines(
  summary_text,
  con = file.path(output_dir, "posterior_summary.txt")
)

# Visualize PK profile fits
# It takes 10 of the accepted parameter sets, runs simulations with them, and plots the resulting PK profiles on top of the real data. This provides a visual check that the accepted parameters do, in fact, produce simulations that look like the real data.
pk_fit_plot <- simPlot_fun(pkdat_mlx, pkml_file ="PKML/Olaparib_300mg_BID_Weibull_v2.pkml", 
                           outputs_v = outputs_v,
                           priors_setup_df= priors_setup_df, fixed_setup_df= data.frame() ,
                           posterior_set=posterior6[1:10,], 
                           selID="300_Plummer_2015")
cat("Saving PK profile fit plot...\n")
ggsave(
  filename = file.path(output_dir, "pk_profile_fit.pdf"),
  plot = pk_fit_plot,
  width = 7,
  height = 5
)

# Run the metric function
metric_results <- calculate_distribution_metric(ref_df = posterior1, 
                                                comp_df = posterior6)

# View the results
print(paste("Total Difference Score:", metric_results$total_difference_score))
print("Detailed Breakdown:")
detailed_breakdown_2a <- metric_results$detailed_breakdown
print(detailed_breakdown_2a)

cat("Saving metric breakdown...\n")
write.csv(
  detailed_breakdown_2a,
  file = file.path(output_dir, "metric_breakdown.csv"),
  row.names = FALSE # Usually good practice for analysis data frames
)

# Plot with comparison of posterior distributions (rejection from prior vs. SMC)
posterior_dist_plot2 <- PosteriorPlots_fun2(posterior1=posterior1, posterior2= posterior6, appr1="Rejections from prior", appr2="SMC")
cat("Saving posterior distribution plots...\n")
# Now, add a title to the plot object
posterior_dist_plot2 <- posterior_dist_plot2 + 
  labs(title = "Comparison of posterior parameter distributions (epsilon = 0.3)")
ggsave(
  filename = file.path(output_dir, "posterior_distributions_comparison.pdf"),
  plot = posterior_dist_plot2,
  width = 8,
  height = 7
)





################### Make a reproduction summary table
exp_num <- c("0", "0_2", "0_3", "0_4", "2a", "2a_2", "2a_3", "2a_4")

observed_data <- c(
  "300_Plummer_2015", "300_Plummer_2015", "300_Plummer_2015", "300_Plummer_2015", 
  "300_Plummer_2015", "300_Plummer_2015", "300_Plummer_2015", "300_Plummer_2015"
)

tested_params <- c(
  "base metrics", "base metrics", "base metrics", "base metrics", "Start_epsilon = 0.3",
  "Start_epsilon = 0.3", "Start_epsilon = 0.3", "Start_epsilon = 0.3"
)

similarity_metric <- c(4, 4, 5, 4, 1, 8, 5, 5)

execution_time <- c(
  35.03735, 38.01595, 26.72996, 36.48316, 24.8635, 23.49269, 23.75158, 49.55859
)

# 2. Create the data frame.
reproduction_table <- data.frame(
  `id` = exp_num,
  `Observed data` = observed_data,
  `Tested parameters` = tested_params,
  `Distribution similarity metric` = similarity_metric,
  `Execution time (min)` = execution_time,
  stringsAsFactors = FALSE # Good practice to avoid factors unless needed
)

# 3. Print the table to the console to verify it looks correct.
print(reproduction_table)

# 4. Define the output path and filename
output_path <- "Results/Anya/CV_2/SequentialMC/"
output_filename <- "reproduction_table.csv"

# 5. Save the final table to a CSV file
write.csv(
  reproduction_table,
  file = file.path(output_path, output_filename),
  row.names = FALSE  # Crucial: Prevents R from adding an extra column with its own row numbers
)
