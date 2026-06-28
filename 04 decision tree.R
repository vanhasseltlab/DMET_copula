# Date    : 01/06/2026 ----
# Author  : Yuchen Guo / iguoyuchen@outlook.com ----
# Purpose :  Development and design of decision tree ----

#############################################
# 1.0 Load the data #########################
#############################################

# CYP #####
# x scale, CYP_log
load("data/CYP_log.Rdata") 
# Z scale, CYP_un
load("data/CYP_un.Rdata")

# UGT #####
# x scale, UGT_log
load("data/UGT_log.Rdata") 
# Z scale, UGT_un
load("data/UGT_un.Rdata")

# ABC #####
# x scale, ABC_log
load("data/ABC_log.Rdata") 
# Z scale, ABC_un
load("data/ABC_un.Rdata")

# SLC #####
# x scale, SLC_log
load("data/SLC_log.Rdata") 
# Z scale, SLC_un
load("data/SLC_un.Rdata")

# get the name list
paste(colnames(CYP_log), collapse = ", ")
paste(colnames(UGT_log), collapse = ", ")
paste(colnames(ABC_log), collapse = ", ")
paste(colnames(SLC_log), collapse = ", ")


#############################################
# 2.0 Get the MVN tests #####################
#############################################
## approach 1 #----
# install.packages("MVN")
library(MVN)

result <- mvn(data = CYP_log, mvn_test = "mardia", alpha = 0.05)
result$multivariate_normality
result <- mvn(data = CYP_un, mvn_test = "mardia")
result$multivariate_normality

result <- mvn(data = UGT_log, mvn_test = "mardia")
result$multivariate_normality
result <- mvn(data = UGT_un, mvn_test = "mardia")
result$multivariate_normality

result <- mvn(data = ABC_log, mvn_test = "mardia")
result$multivariate_normality
result <- mvn(data = ABC_un, mvn_test = "mardia")
result$multivariate_normality

result <- mvn(data = SLC_log, mvn_test = "mardia")
result$multivariate_normality
result <- mvn(data = SLC_un, mvn_test = "mardia")
result$multivariate_normality

## approach 2 #----
# install.packages("energy")
library(energy)
?mvnorm.etest
# Test for MVN using the energy package
result <- mvnorm.etest(CYP_log, R = 1000)  # R = bootstrap replicates
print(result)
result <- mvnorm.etest(CYP_un, R = 1000) 
print(result)
result <- mvnorm.etest(UGT_log, R = 1000)  
print(result)
result <- mvnorm.etest(UGT_un, R = 1000)  
print(result)
result <- mvnorm.etest(ABC_log, R = 1000) 
print(result)
result <- mvnorm.etest(ABC_un, R = 1000)  
print(result)
result <- mvnorm.etest(SLC_log, R = 1000)  
print(result)
result <- mvnorm.etest(SLC_un, R = 1000)  
print(result)

## approach 3 #----
result <- mvn(data = CYP_log, mvn_test =  "royston")
result$multivariate_normality
result <- mvn(data = CYP_un, mvn_test =  "royston")
result$multivariate_normality

result <- mvn(data = UGT_log, mvn_test =  "royston")
result$multivariate_normality
result <- mvn(data = UGT_un, mvn_test =  "royston")
result$multivariate_normality

result <- mvn(data = ABC_log, mvn_test = "royston")
result$multivariate_normality
result <- mvn(data = ABC_un, mvn_test = "royston")
result$multivariate_normality

result <- mvn(data = SLC_log, mvn_test =  "royston")
result$multivariate_normality
result <- mvn(data = SLC_un, mvn_test = "royston")
result$multivariate_normality

## approach 4 #----
install.packages("MVN")
library(MVN)
### method 1 ------------
library(nvmix)
library(ggplot2)
install.packages("nvmix")
X <- as.matrix(CYP_log)
loc <- colMeans(X)
scale <- cov(X)
qqplot_maha(X, qmix = qmix, df = df, loc = loc, scale = scale)
### method 2 ------------
par(mfrow = c(4, 2), mar = c(4, 4, 2, 1), oma = c(0, 0, 2, 0))
# geom_point(color = "black", fill = "black", alpha = 0.5, shape = 25) +
data_list <- list(CYP_log, CYP_un, UGT_log, UGT_un, ABC_log, ABC_un, SLC_log, SLC_un)
title_list <- list(
  bquote(italic("X")*"-scale CYP"),
  bquote(italic("Z")*"-scale CYP"),
  bquote(italic("X")*"-scale UGT"),
  bquote(italic("Z")*"-scale UGT"),
  bquote(italic("X")*"-scale ABC"),
  bquote(italic("Z")*"-scale ABC"),
  bquote(italic("X")*"-scale SLC"),
  bquote(italic("Z")*"-scale SLC")
)


for (i in seq_along(data_list)) {
  assign(
    paste0("plot_", i),
    multivariate_diagnostic_plot(data_list[[i]], type = "qq") +
      theme_bw() +
      scale_fill_identity() +
      scale_color_identity() +
      guides(color = "none", fill = "none") +
      labs(title = title_list[[i]], subtitle = NULL) +
      theme(plot.title = element_text(hjust = 0.5))
  )
}


library(ggpubr)
Figure <- ggarrange(plot_1,plot_2,plot_3,plot_4,
                      plot_5,plot_6,plot_7,plot_8,
                      ncol = 4, nrow = 2,heights=c(2,2))
Figure
