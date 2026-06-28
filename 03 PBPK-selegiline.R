# Date    : 01/06/2026 ----
# Author  : Yuchen Guo / iguoyuchen@outlook.com ----
# Purpose : PBPK application of [SELEGILINE]
#         : [SELEGILINE] is metabolized by CYP2B6 and CYP2C19


## load packages #----
library(pmxcopula)
library(dplyr)
library(reshape2)
library(ggplot2)
library(tidyr)
library(ggpubr)

## load function #-----
source("functions/Uniform_transform.R")

#########################
# 0.1 Activity data #----
#########################

# Km and Vm data from literature:
# https://pubmed.ncbi.nlm.nih.gov/11602525/
# collected from the recombinant CYP enzymes
# The derived CLint are:
# unit: mL/min/umol # equal to 
# CYP2B6
CL_invitro_2B6_1 <- 2055
CL_invitro_2B6_2 <- 1994

# CYP2C19
CL_invitro_2C19_1 <- 350
CL_invitro_2C19_2 <- 1500

##########################
# 0.2 Abundance data #----
##########################
set.seed(22222)
NN = 26200 
## VP1 mvn ----
cov_matrix <- cov(CYP_log)
VP_mvn <- MASS::mvrnorm(n = NN,
                        mu = colMeans(CYP_log),
                        Sigma = cov_matrix,
                        tol = 1e-6)  %>% 
  as.data.frame() 
VP_mvn <- 2^(VP_mvn) -1

## VP2 transmvn ----
cov_matrix <- cov(CYP_un)
VP_transmvn_Z <- MASS::mvrnorm(n = NN,
                               mu = colMeans(CYP_un),
                               Sigma = cov_matrix,
                               tol = 1e-6) %>% 
  as.data.frame() 
# z-scale simulation data
# require x-scalue simulation data for comparison
### from normal to uniform ----
VP_transmvn_U <- VP_transmvn_Z
for (i in 1:(ncol(VP_transmvn_U))){
  VP_transmvn_U[,i] <- pnorm(VP_transmvn_Z[,i], mean = 0, sd = 1)
}

### from uniform to original ----
VP_transmvn_X <- backuniform(dat_sim = VP_transmvn_U, dat_org = CYP_log)
# x-scale simulation data
VP_transmvn <- 2^(VP_transmvn_X) -1

## VP3 copula ----
copula <- copula_CYP[[4]]
VP_copula_U <- rvinecop(NN, copula) %>% 
  as.data.frame() 
#### 3.3.2.2 back transform
VP_copula_X <- backuniform(dat_sim = VP_copula_U, dat_org = CYP_log)
VP_copula <- 2^(VP_copula_X) -1

# simulation data ----
# vp_mvn
# VP_transmvn
# VP_copula

###########################################################
# 1.0 Gene expression translate to protein expression #----
###########################################################
# data source: dx.doi.org/10.1124/dmd.124.001608
# Anchor enzyme: CYP2C19
# Median protein concentration: 
# 0.306 pmol/million hepatocytes


# 139 million hepatocytes /gram of human liver tissue
# liver weight: 1400 g (70kg, 2%)

# Get the median gene expression level
median(CYP_expression_ins$CYP2C19) # 15.7896

RF <- 0.306/median(CYP_expression_ins$CYP2C19) # 0.01937984
fu <- 0.1

# protein conc: 
# unit: pmol/million hepatocytes

scale_and_extract <- function(VP, RF) {
  VP_df <- as.data.frame(VP)
  VP_df_scaled <- VP_df * RF
  
  VP_df_scaled |>
    dplyr::select(CYP2B6, CYP2C19) |>
    dplyr::filter(CYP2B6 >= 0, CYP2C19 >= 0)
}


VP_obs_p        <- scale_and_extract(CYP_expression_ins, RF)
VP_mvn_p        <- scale_and_extract(VP_mvn, RF)
VP_transmvn_p   <- scale_and_extract(VP_transmvn, RF)
VP_copula_p     <- scale_and_extract(VP_copula, RF)


############################################
# 2.0 Calculation of hepatic clearance #----
############################################

# in vivo intrinsic clearance ----
# Clint mL/min
Q = 1450 # mL/min, https://pubmed.ncbi.nlm.nih.gov/11602525/
Dose = 10 # mg

# unit: mL/min/umol
CL_invitro_2B6_1 <- 2055
CL_invitro_2B6_2 <- 1994
CL_invitro_2C19_1 <- 350
CL_invitro_2C19_2 <- 1500

# 139 million hepatocytes /gram of human liver tissue
# liver weight: 1500 g (70kg, ~2%)

# B:P value
# CLp [plasma clearance] = CL [blood clearance] * B_P [blood plasma ratio]
# Ref 1: B:P ~ 1.34; doi:10.3390/pharmaceutics12100942
# Ref 2: B:P ~ 1.3; derived by Cmax ratio, not steady state; https://link.springer.com/article/10.2165/00003088-199733020-00002
B_P = 1.34


compute_CLint <- function(VP_df, fu, Q, B_P, Dose) {
  VP_df |>
    dplyr::mutate(
      CLint_2B6_1  = CL_invitro_2B6_1  * CYP2B6 * 139 * 1500 * 1e-6,  # mL/min
      CLint_2C19_1 = CL_invitro_2C19_1 * CYP2C19 * 139 * 1500 * 1e-6,
      CLint_2B6_2  = CL_invitro_2B6_2  * CYP2B6 * 139 * 1500 * 1e-6,
      CLint_2C19_2 = CL_invitro_2C19_2 * CYP2C19 * 139 * 1500 * 1e-6
    ) |>
    dplyr::mutate(
      CLhep_2B6_1  = (CLint_2B6_1  * fu * Q) / (fu * CLint_2B6_1  + Q),
      CLhep_2C19_1 = (CLint_2C19_1 * fu * Q) / (fu * CLint_2C19_1 + Q),
      CLhep_2B6_2  = (CLint_2B6_2  * fu * Q) / (fu * CLint_2B6_2  + Q),
      CLhep_2C19_2 = (CLint_2C19_2 * fu * Q) / (fu * CLint_2C19_2 + Q)
    ) |>
    dplyr::mutate(
      CL  = CLhep_2B6_1 + CLhep_2C19_1 + CLhep_2B6_2 + CLhep_2C19_2,
      CLp = CL * B_P,
      AUC = Dose / CLp # mg*min/mL
    )
}


CLint_data_obs      <- compute_CLint(VP_obs_p,      fu, Q, B_P, Dose)
CLint_data_mvn      <- compute_CLint(VP_mvn_p,      fu, Q, B_P, Dose)
CLint_data_transmvn <- compute_CLint(VP_transmvn_p, fu, Q, B_P, Dose)
CLint_data_copula   <- compute_CLint(VP_copula_p,   fu, Q, B_P, Dose)

# Conversion of units

# unit of CL
# mL/min
#  L/h
Factor <- 60*10^-3

# unit of AUC
# mg*min/mL
# convert to ng*h/mL
Factor_AUC <- 10^5 / 6


make_group_df <- function(df, group, Factor, Factor_AUC) {
  df %>%
    mutate(
      CLp  = Factor * CLp,
      AUC  = Factor_AUC * AUC,
      Group = group,
      ratio = CYP2B6/CYP2C19 
    ) %>%
    select(CLp, CYP2B6, CYP2C19, AUC, Group, ratio)
}

sum_table <- bind_rows(
  make_group_df(CLint_data_obs,      "observed",    Factor, Factor_AUC),
  make_group_df(CLint_data_mvn,      "VP_MVN",      Factor, Factor_AUC),
  make_group_df(CLint_data_transmvn, "VP_transMVN", Factor, Factor_AUC),
  make_group_df(CLint_data_copula,   "VP_copula",   Factor, Factor_AUC)
)

facet_data <- sum_table %>% filter(Group != "observed")
facet_data$Group <- factor(facet_data$Group,
                           levels = unique(facet_data$Group)) # simulated data
nof_data <- sum_table %>% filter(Group == "observed") %>% select(-Group) # observed data


#################################
# 3.0 Plotting of donutVPC  #----
#################################

##### data rearrangement ----
sim_data <- facet_data %>% 
  mutate(expression = log(ratio), LogAUC = log(AUC)) %>% 
  select(expression, CLp, LogAUC, Group, CYP2B6, CYP2C19)

# mvn
sim_data_mvn <- sim_data %>% filter(Group == "VP_MVN")
simulation_nr = c(rep(1:99, each =262), rep(100, times = (nrow(sim_data_mvn)-25938)))
sim_data_mvn <- cbind(sim_data_mvn, simulation_nr) # 262*99 = 25938
# transMVN
sim_data_transmvn <- sim_data %>% filter(Group == "VP_transMVN")
simulation_nr = c(rep(1:99, each =262), rep(100, times = (nrow(sim_data_transmvn)-25938)))
sim_data_transmvn <- cbind(sim_data_transmvn, simulation_nr) # 262*99 = 25938
# copula
sim_data_copula <- sim_data %>% filter(Group == "VP_copula")
simulation_nr = c(rep(1:99, each =262), rep(100, times = (nrow(sim_data_copula)-25938)))
sim_data_copula <- cbind(sim_data_copula, simulation_nr) # 262*99 = 25938

obs_data <- nof_data %>% 
  mutate(expression = log(ratio), LogAUC = log(AUC)) %>% 
  select(expression, CLp, LogAUC,CYP2B6, CYP2C19)

## donutVPC #----
VPC_mvn <- donutVPC(sim_data = sim_data_mvn, 
                  obs_data = obs_data, 
                  percentiles = c(5, 50, 95), 
                  sim_nr = 100, 
                  pairs_matrix = matrix(c("expression", "expression", "CLp", "LogAUC"), 2, 2), 
                  conf_band = 95, colors_bands = c("#1B9E77","#A6D8C5"),
                  cores = 10)
VPC_transmvn <- donutVPC(sim_data = sim_data_transmvn, 
                    obs_data = obs_data, 
                    percentiles = c(5, 50, 95), 
                    sim_nr = 100, 
                    pairs_matrix = matrix(c("expression", "expression", "CLp", "LogAUC"), 2, 2), 
                    conf_band = 95, colors_bands = c("#E69F00","#FFC43F"),
                    cores = 10)
VPC_copula <- donutVPC(sim_data = sim_data_copula, 
                         obs_data = obs_data, 
                         percentiles = c(5, 50, 95), 
                         sim_nr = 100, 
                         pairs_matrix = matrix(c("expression", "expression", "CLp", "LogAUC"), 2, 2), 
                         conf_band = 95, colors_bands = c("#7570B3","#BFBCE5"),
                         cores = 10)


########################################
# 4.0 Calculation of statistics  #----
########################################

summarize_metric <- function(df, dataset_name, var, factor = 1) {
  x <- factor * df[[var]]
  
  tibble(
    Dataset = dataset_name,
    Metric  = var,
    Median  = median(x, na.rm = TRUE),
    Mean    = mean(x,   na.rm = TRUE),
    CV      = sd(x,     na.rm = TRUE) / mean(x, na.rm = TRUE) * 100
  )
}

datasets <- list(
  obs_pop      = CLint_data_obs,
  VP_MVN       = CLint_data_mvn,
  VP_transMVN  = CLint_data_transmvn,
  VP_copula    = CLint_data_copula
)

summary_CLp <- bind_rows(lapply(names(datasets), function(nm) {
  summarize_metric(datasets[[nm]], nm, "CLp", factor = Factor)
}))

summary_AUC <- bind_rows(lapply(names(datasets), function(nm) {
  summarize_metric(datasets[[nm]], nm, "AUC", factor = Factor_AUC)
}))

print(summary_CLp)
print(summary_AUC)

