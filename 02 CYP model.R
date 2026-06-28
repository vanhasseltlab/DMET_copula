# Date    : 01/06/2026 ----
# Author  : Yuchen Guo / iguoyuchen@outlook.com ----
# Purpose : Joint distribution modelling for CYP gene family ----
# Note    : This script provides the joint distribution modelling workflow used in the project.
#           The same workflow is applied to all datasets (i.e., UGT, SLC, ABC); only the input data files
#           and a few dataset-specific settings differ----


## load packages #----
library(readxl)
library(dplyr)
library(ggplot2)
library(reshape2)
library(rvinecopulib)
library(pmxcopula)
library(stringr)
library(tidyr)
library(ggpubr)
library(cowplot)
library(ggpubr)
library(sf)
library(kde1d)

## load data #----
# get the name of gene of interest; use CYP family
# https://doi.org/10.1021/acs.molpharmaceut.8b00941
load("data/expression_profile.Rdata")  # 'expression_profile'
# data from the database
# https://gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression
# "Gene TPM by tissue" --> download the liver one
load("GTEx_Analysis_v10_raw_complete.Rdata") # 'Liver_expression'

## load functions #----
source("functions/Uniform_transform.R")

#######################
# 1 DMET selection ----
#######################

## 1.1 collect the data #----
### 1.1.1 based on the DMET publication #----
# https://doi.org/10.1021/acs.molpharmaceut.8b00941
CYP_names <- colnames(expression_profile[[1]]) 

### 1.1.2 based on the coverag of PKSIM #----
# narrow down the range of CYP enzyme
CYP_list <- cbind(CYP_names, c(1 : length(CYP_names)))
CYP_list
# PKSIM selection: 13 DMET for CYP
CYP_ins <- CYP_names[c(1, 2, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 21)] 
length(CYP_ins)
CYP_ins

# Optional: check for non-ASCII characters in gene names (encoding issues)
# stringi::stri_enc_isascii(CYP_ins)

names.use <- colnames(Liver_expression_complt)[(colnames(Liver_expression_complt) %in% CYP_ins)] 
CYP_expression_ins <- Liver_expression_complt[,names.use] 
CYP_expression_ins <- CYP_expression_ins %>%  as.data.frame() %>% mutate(across(everything(), as.numeric))

## 1.2 check data missingness #----
ag <- sapply(CYP_expression_ins, function(x) list(ActualN=nrow(CYP_expression_ins)-sum(is.na(x)),MissingP=100*(sum(is.na(x)))/nrow(CYP_expression_ins),mean=round(mean(x,na.rm=T),2),
                                        sd=round(sd(x,na.rm=T),2),min=round(min(x,na.rm=T),2),max=round(max(x,na.rm=T),2)))
ag1 <- t(ag) %>%  
  data.frame() %>%
  mutate(meanSD=paste(mean,sd,sep="±"), 
         range=paste(min,max,sep="~")) %>%
  mutate(statistic = paste0(meanSD," ", "[",range,"]"))
ag1 # no missing data

##############################
# 2 data transformation  ----
##############################

############## Z scale ###################
## 2.1 log2(n+1) transformed dataset #----
##########################################
CYP_log <- log(CYP_expression_ins + 1, base = 2)

############## U scale ###################
## 2.2 uniform transformed dataset #----
##########################################
CYP_u <- getUniform(CYP_log)

############## U scale ##################
## 2.3 normal transformed dataset   #----
#########################################
# with uniform and normal transformation
# create a empty dataframe
CYP_un <- CYP_u
for (i in 1:ncol(CYP_un)){
  CYP_un[,i] <- qnorm(CYP_u[,i], mean = 0, sd = 1)
}

##########################################
## 2.4 validate the transformation #-----
##########################################
df_o <- melt(CYP_expression_ins)
ggplot(df_o, aes(x = value)) +
  geom_histogram(bins = 30, fill = "skyblue", color = "white") +
  facet_wrap(~variable, scales = "free") +
  theme_minimal()

df_x <- melt(CYP_log)
ggplot(df_x, aes(x = value)) +
  geom_histogram(bins = 30, fill = "skyblue", color = "white") +
  facet_wrap(~variable, scales = "free") +
  theme_minimal()

df_u <- melt(CYP_u)
ggplot(df_u, aes(x = value)) +
  geom_histogram(bins = 30, fill = "skyblue", color = "white") +
  facet_wrap(~variable, scales = "free") +
  theme_minimal()

df_z <- melt(CYP_un)
ggplot(df_z, aes(x = value)) +
  geom_histogram(bins = 30, fill = "skyblue", color = "white") +
  facet_wrap(~variable, scales = "free") +
  theme_minimal()


########################################
# 3 Joint modeling and simulation  #----
########################################

############################################
## 3.1 multivariate Gaussian distribution #----
############################################
set.seed(12345)
cov_matrix <- cov(CYP_log)
CYP_sim_mvn <- MASS::mvrnorm(n = 100*nrow(CYP_log),
                             mu = colMeans(CYP_log),
                             Sigma = cov_matrix,
                             tol = 1e-6)  %>% 
  as.data.frame() %>% 
  mutate(simulation_nr = rep(1:100, each = nrow(CYP_log)))
# x-scale simulation data

##############################################
## 3.2 multivariate Gaussian distribution #----
#########   + double transform   #############
##############################################
cov_matrix <- cov(CYP_un)
CYP_sim_mvn_un <- MASS::mvrnorm(n = 100*nrow(CYP_un),
                                mu = colMeans(CYP_un),
                                Sigma = cov_matrix,
                                tol = 1e-6) %>% 
  as.data.frame() %>% 
  mutate(simulation_nr = rep(1:100, each = nrow(CYP_un)))
# z-scale simulation data
# require x-scalue simulation data for comparison

### 3.2.1 from normal to uniform #----
CYP_sim_mvn_Uscale <- CYP_sim_mvn_un
for (i in 1:(ncol(CYP_sim_mvn_un) - 1)){
  CYP_sim_mvn_Uscale[,i] <- pnorm(CYP_sim_mvn_Uscale[,i], mean = 0, sd = 1)
}

### 3.2.2 from uniform to original #----
CYP_sim_transmvn <- backuniform(dat_sim = CYP_sim_mvn_Uscale, dat_org = CYP_log)
# this is the simulation data on x scale

### 3.2.3 validate the transformation #----
df_z <- melt(CYP_sim_mvn_un[,-14])
ggplot(df_z, aes(x = value)) +
  geom_histogram(bins = 30, fill = "skyblue", color = "white") +
  facet_wrap(~variable, scales = "free") +
  theme_minimal()

df_u <- melt(CYP_sim_mvn_Uscale[,-14])
ggplot(df_u, aes(x = value)) +
  geom_histogram(bins = 30, fill = "skyblue", color = "white") +
  facet_wrap(~variable, scales = "free") +
  theme_minimal()

df_x <- melt(CYP_sim_transmvn[,-14])
ggplot(df_x, aes(x = value)) +
  geom_histogram(bins = 30, fill = "skyblue", color = "white") +
  facet_wrap(~variable, scales = "free") +
  theme_minimal()

###########################################
## 3.3 copula + uniform transform   ##-----
###########################################

### 3.3.1 calculation of vine copula #----
# as the dataset is already uniformly transformed, we could use vinecop() to calculate the density function

#### 3.3.1.1 Full copula #----
copula_1 <- vinecop(data = CYP_u,
                    var_types = rep("c", ncol(CYP_u)),
                    family_set = "parametric",
                    par_method = "mle",
                    trunc_lvl = Inf, # Inf means no truncation
                    selcrit = "aic",
                    cores = 10)
copula_1 
string1 <- c("parametric", "Inf", "Inf", BIC(copula_1), AIC(copula_1))

#### 3.3.1.2 Truncation copula #----
copula_2 <- vinecop(data = CYP_u,
                    var_types = rep("c", ncol(CYP_u)),
                    family_set = "parametric",
                    par_method = "mle",
                    trunc_lvl = NA, # NA indicates that the truncation level should be selected automatically by mBICV()
                    selcrit = "aic",
                    cores = 10)
copula_2  
string2 <- c("parametric", "NA", "Inf", BIC(copula_2), AIC(copula_2))

#### 3.3.1.3 Gaussian vine #----
copula_3 <- vinecop(data = CYP_u,
                    var_types = rep("c", ncol(CYP_u)),
                    family_set = "gaussian",
                    par_method = "mle",
                    trunc_lvl = Inf,
                    selcrit = "aic",
                    cores = 10)
copula_3  
string3 <- c("gaussian", "Inf", "Inf", BIC(copula_3), AIC(copula_3))

#### 3.3.1.4 Pruning copula #----
set.seed(123)
copula_4 <- vinecop(data = CYP_u,
                    var_types = rep("c", ncol(CYP_u)),
                    family_set = "parametric",
                    par_method = "mle",
                    trunc_lvl = Inf, # Inf means no truncation
                    selcrit = "aic",
                    threshold = NA, #  NA indicates that the threshold should be selected automatically by mBICV()
                    cores = 10)
copula_4     
string4 <- c("parametric", "Inf", "NA", BIC(copula_4), AIC(copula_4))

#### 3.3.1.5 Pruning + truncate copula #----
copula_5 <- vinecop(data = CYP_u,
                    var_types = rep("c", ncol(CYP_u)),
                    family_set = "parametric",
                    par_method = "mle",
                    trunc_lvl = NA, # NA indicates that the truncation level should be selected automatically by mBICV()
                    selcrit = "aic",
                    threshold = NA, #  NA indicates that the threshold should be selected automatically by mBICV()
                    cores = 10)
copula_5   
string5 <- c("parametric", "NA", "NA", BIC(copula_5), AIC(copula_5))

#### 3.3.1.6 Gaussian + truncate copula #----
copula_6 <- vinecop(data = CYP_u,
                    var_types = rep("c", ncol(CYP_u)),
                    family_set = "gaussian",
                    par_method = "mle",
                    trunc_lvl = NA, # NA indicates that the truncation level should be selected automatically by mBICV()
                    selcrit = "aic",
                    cores = 10)
copula_6   
string6 <- c("gaussian", "NA", "Inf", BIC(copula_6), AIC(copula_6))

#### 3.3.1.7 Gaussian + Pruning copula #----
copula_7 <- vinecop(data = CYP_u,
                    var_types = rep("c", ncol(CYP_u)),
                    family_set = "gaussian",
                    par_method = "mle",
                    trunc_lvl = Inf, # Inf means no truncation
                    selcrit = "aic",
                    threshold = NA, #  NA indicates that the threshold should be selected automatically by mBICV()
                    cores = 10)
copula_7   

#### 3.3.1.7 Gaussian + Pruning + truncate copula #----
copula_8 <- vinecop(data = CYP_u,
                    var_types = rep("c", ncol(CYP_u)),
                    family_set = "gaussian",
                    par_method = "mle",
                    trunc_lvl = NA, # NA indicates that the truncation level should be selected automatically by mBICV()
                    selcrit = "aic",
                    threshold = NA, #  NA indicates that the threshold should be selected automatically by mBICV()
                    cores = 10)
copula_8 
string8 <- c("gaussian", "NA", "NA", BIC(copula_8), AIC(copula_8))

#### 3.3.1.8 Summarize the model #----
copula_table <- rbind(string1, string2, string3, string4, string5, string6, string7, string8)
colnames(copula_table) <- c("bivariate copula function", "truncation", "prune", "BIC", "AIC")
copula_table <- cbind(copula_table, feature = c(" full model", "truncation", "Gaussian copula", "pruning",
                                                "truncation + pruning", "Gaussian + truncation",
                                                "Gaussian + pruning", "truncation + pruning + Gaussian"))
copula_table <- as.data.frame(copula_table)
copula_table$feature <- factor(copula_table$feature, levels = copula_table$feature)
copula_table$BIC <- round(as.numeric(copula_table$BIC),1)
copula_table$AIC <- round(as.numeric(copula_table$AIC),1)
results_long <- pivot_longer(copula_table, cols = c("AIC", "BIC"), 
                             names_to = "Metric", values_to = "Value")

CYP_copula_metric <- ggplot(results_long, aes(x = feature, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("#EE9964", "#55669B")) +
  labs(title = "CYP family enzyme", 
       x = "Copula Model", 
       y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.title = element_text(hjust = 0.5))

# conclusion:
# copula_4 is selected, for the lowest BIC
copula_final <- copula_4


### 3.3.2 simulate from the selected copula #----
#### 3.3.2.1 simulate from copula function #----
CYP_sim_copula_u <- rvinecop(100*nrow(CYP_u), copula_final) %>% 
  as.data.frame() %>% 
  mutate(simulation_nr = rep(1:100, each = copula_final$nobs))

#### 3.3.2.2 back transform #----
CYP_sim_copula <- backuniform(dat_sim = CYP_sim_copula_u, dat_org = CYP_log)
# x-scale simulation data

##################################################
# 4 Model performance comparision : marginal #----
##################################################
## 4.1 Marginal comparision  #----
### 4.1 Figure 1 #----
#### 4.1.1 density of observation data #----
margin_obs <- list()
for (k in 1: ncol(CYP_log)) {
  range_input <- CYP_log[,k]
  x_grid <- seq(min(range_input), max(range_input), length.out = 512)
  fit <- kde1d(range_input)
  dens_at_x <- dkde1d(x_grid,fit)
  data_sum <- data.frame(x = x_grid, y =  dens_at_x) %>% mutate(DMET = colnames(CYP_log)[k])
  margin_obs[[k]] <- data_sum
}

#### 4.1.2 density of simulation data #----
sim_all <- list(CYP_sim_mvn, CYP_sim_transmvn, CYP_sim_copula)
names(sim_all) <- c("MVN", "transMVN", "copula")

# initialization
margin_sim_sub <- list()
result_sum <- list()
result <- list()

for (m in 1:3) {
  sim_data <- sim_all[[m]]
  
  for (k in 1: ncol(CYP_log))  {
    range_input <- CYP_log[,k]
    CYP_sim_sub <- sim_data[,c(colnames(CYP_log)[k],"simulation_nr")]
    
    for (j in 1:100){
      density_input <- CYP_sim_sub %>% 
        filter(simulation_nr == j) %>% 
        dplyr::select(-simulation_nr)
      
      fit <- kde1d(density_input[[1]])
      x_grid <- seq(min(range_input), max(range_input), length.out = 512)
      dens_vals <- dkde1d(x_grid,fit)
      
      margin_sim_sub[[j]] <- data.frame(x = x_grid, y = dens_vals)
    }

    x_vals <- margin_sim_sub[[1]][, 1]
    y_matrix <- sapply(1:100, function(i) margin_sim_sub[[i]][, 2])
    result_df <- data.frame(x = x_vals, y_matrix)
    colnames(result_df) <- c("x", paste0("y_sim_", 1:100))
    
    q025 <- apply(result_df[, -1], 1, quantile, probs = 0.025)
    q975 <- apply(result_df[, -1], 1, quantile, probs = 0.975)
    
    quantile_df <- data.frame(
      x = result_df$x,
      q025 = q025,
      q975 = q975,
      DMET = colnames(CYP_log)[k],
      VP = names(sim_all)[m]
    )
    
    result_sum[[k]] <- quantile_df
    
  }
  result[[m]] <- do.call(rbind,result_sum)
}
margin_result <- do.call(rbind,result)
margin_obs <- do.call(rbind,margin_obs)

p_all <- ggplot(margin_result, aes(x = x)) +
  geom_ribbon(aes(ymin = q025, ymax = q975, 
                  fill = factor(VP, levels = c("MVN", "transMVN", "copula"))), alpha = 0.4) +
  geom_line(data = margin_obs, aes(x = x, y = y, color = "Observed")) +
  xlab("Expression (DMET)") +
  ylab("Density") +
  facet_wrap(~DMET, scales = 'free',ncol = 5) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 3), 
                     labels = scales::label_number(accuracy = 0.1)) +
  scale_fill_manual(name = "Population", 
                     values = c("MVN" = "#1B9E77", # "#219EBC" 
                                "transMVN" = "#E69F00", # "#D95F02" 
                                "copula" = "#7570B3"),
                     labels = c("MVN" = "VP_MVN", 
                                "transMVN"  = "VP_transMVN", 
                                "copula" = "VP_copula"))+ 
  scale_color_manual(
    name = "", values = c("Observed" = "black"),labels = c("Observed" = "Observed")
  ) +
  theme_bw() +
  guides(
    fill = guide_legend(order = 1, override.aes = list(alpha = 0.4)),
    color = guide_legend(order = 2)
  ) +
  theme(strip.background = element_rect(fill = NA), panel.grid.minor.x = element_blank(),
        strip.text = element_text(size = 8),
        legend.background = element_rect(fill = NA, color = NA),
        legend.spacing.y = unit(-0.8, "cm"),
        axis.text = element_text(size = 8),
        plot.title = element_text(size=12, hjust = 0.5),
        plot.subtitle = element_text(size=10),
        legend.title = element_text(size = 10, face = "plain"))
p_all


### 4.2 Figure 2 #----
mtr_margin_1 <- calc_margin(
  sim_data = CYP_sim_mvn,
  obs_data = CYP_log,
  sim_nr = 100,
  var = NULL,
) %>%  mutate(Population = "MVN")  

mtr_margin_2 <- calc_margin(
  sim_data = CYP_sim_transmvn,
  obs_data = CYP_log,
  sim_nr = 100,
  var = NULL,
) %>%  mutate(Population = "transMVN")  

mtr_margin_3 <- calc_margin(
  sim_data = CYP_sim_copula,
  obs_data = CYP_log,
  sim_nr = 100,
  var = NULL,
) %>%  mutate(Population = "copula")  

marg_metric_dat <- rbind(mtr_margin_1,
                         mtr_margin_2,
                         mtr_margin_3) %>% filter(statistic %in% c("mean","median","sd", "Q5.5%","Q95.95%","min","max"))

marg_metric_dat$statistic <- factor(marg_metric_dat$statistic,
                                    levels = c("mean","median","sd","Q5.5%","Q95.95%","min","max"),
                                    labels = c("mean","median","SD","5th P","95th P","min","max"))
marg_metric_dat$Population <- factor(marg_metric_dat$Population,
                                     levels = c("MVN", "transMVN", "copula"),
                                     labels = c("MVN", "transMVN", "copula"))

CYP_marg_metric <- marg_metric_dat %>% 
  filter(statistic %in% c("median","5th P","95th P", "mean", "SD")) %>% 
  ggplot(aes(y = rel_error, x = variable, color = Population)) +
  geom_vline(xintercept = seq(0.5, 14, by = 1), color = "grey95") +
  geom_boxplot(outlier.shape = NA) + # fill = "white", color = "#F68E60",
  geom_hline(yintercept = c(-0.2, 0.2), linetype = 2, color = "grey65") +
  geom_hline(yintercept = 0, linetype = 1, color = "#747AA9") +
  labs(x = "Covariates", y = "Relative error", color = "Population") +
  scale_x_discrete(guide = guide_axis(angle = 90), expand = expansion(mult = c(0.1, 0.1))) +
  scale_color_manual(values = c("#1B9E77", "#E69F00","#7570B3"))+
  facet_grid(statistic~., scales = "free")+
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"), panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text=element_text(size=6.5),
        axis.title=element_text(size=10),
        strip.text.x = element_text(size=7))

####################################################
# 5 Model performance comparision : dependency #----
####################################################

## 5.1 Calculation of the dependency metrics #----
metric_depend_1 <- calc_dependency(sim_data = CYP_sim_mvn, 
                                   obs_data = CYP_log, 
                                   pairs_matrix = NULL, 
                                   percentile = 95, 
                                   sim_nr = 100,cores = 10)
metric_depend_1 <- metric_depend_1 %>% mutate(Population = "MVN")  

metric_depend_2 <- calc_dependency(sim_data = CYP_sim_transmvn, 
                                   obs_data = CYP_log, 
                                   pairs_matrix = NULL, 
                                   percentile = 95, 
                                   sim_nr = 100,cores = 10)
metric_depend_2 <- metric_depend_2 %>% mutate(Population = "transMVN") 

metric_depend_3 <- calc_dependency(sim_data = CYP_sim_copula, 
                                   obs_data = CYP_log, 
                                   pairs_matrix = NULL, 
                                   percentile = 95, 
                                   sim_nr = 100,cores = 10)
metric_depend_3 <- metric_depend_3 %>% mutate(Population = "copula") 
metric_depend <- rbind(metric_depend_1,
                       metric_depend_2,
                       metric_depend_3) %>% 
  mutate(error = (value - observed)) 
metric_depend$Population <- factor(metric_depend$Population,
                                   levels = c("MVN", "transMVN", "copula"),
                                   labels = c("MVN", "transMVN", "copula"))


# numeric analysis
corr_info_covariate <- metric_depend %>% filter(statistic == "correlation") %>% 
  group_by(var_pair, Population) %>% 
  summarise(median_err = round(median(error),3)) %>% 
  mutate(abs = abs(median_err))
colnames(corr_info_covariate)
corr_info_covariate_pop <- corr_info_covariate %>%
  group_by(Population) %>%
  summarise(median_median_err = median(median_err, na.rm = TRUE),
            median_abs = median(abs, na.rm = TRUE))
corr_info_covariate_pop

overlap_info_covariate <- metric_depend %>% filter(statistic == "overlap") %>% 
  group_by(var_pair, Population) %>% 
  summarise(median = round(median(value),3))

overlap_info_covariate_2 <- overlap_info_covariate %>% 
  group_by(Population) %>% 
  summarise(median_median = round(median(median),3),
            min_median = min(median),
            max_median = max(median)) %>%
  mutate(range_median = paste0( "[",round(min_median,1),",",round(max_median,1),"]"))

overlap_info_covariate
overlap_info_covariate_pop <- overlap_info_covariate %>%
  group_by(Population) %>%
  summarise(median_median = median(median, na.rm = TRUE))
overlap_info_covariate_pop


## 5.2 donutVPC plotting #----
VPC_1 <- donutVPC(sim_data = CYP_sim_mvn, 
                  obs_data = CYP_log, 
                  percentiles = c(5, 50, 95), 
                  sim_nr = 100, 
                  conf_band = 95, colors_bands = c("#1B9E77","#A6D8C5"), 
                  cores = 10)

VPC_2 <- donutVPC(sim_data = CYP_sim_transmvn, 
                  obs_data = CYP_log, 
                  percentiles = c(5, 50, 95), 
                  sim_nr = 100, 
                  conf_band = 95, colors_bands = c("#E69F00","#FFC43F"),
                  cores = 10)

VPC_3 <- donutVPC(sim_data = CYP_sim_copula, 
                  obs_data = CYP_log, 
                  percentiles = c(5, 50, 95), 
                  sim_nr = 100, 
                  conf_band = 95, colors_bands = c("#7570B3","#BFBCE5"),
                  cores = 10)
