# Date    : 01/06/2026 ----
# Author  : Yuchen Guo / iguoyuchen@outlook.com ----
# Purpose   : Process and curate the GTEx database for downstream analyses ----

# load packages and functions
library(readxl)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggpubr)
library(rvinecopulib)
library(reshape2)

source("functions/functions.R")
source("functions/plot_distributions.R")
source("functions/overlap_advanced.R")

################################
# 01 Pre-selection of DMET genes ----
################################

## DMET gene selection is based on literature: 
# source 1 
# https://doi.org/10.1021/acs.molpharmaceut.8b00941
# reported data: "data/ExpressionData_23_individuals.xlsx"
load("data/expression_profile.Rdata") # expression_profiles of the 23 individuals, four dataframes
all_colnames <- unique(unlist(lapply(expression_profile, colnames)))

# source 2
# PKSIM

################################
# 02 Preparation of DMET data ----
################################

## download and process GTeX v10 data ----
#  https://gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression

## load meta data ----
SampleAttributes_meta <- read_xlsx(path = "data/GTEx_Analysis_v10_Annotations_SampleAttributesDD.xlsx")
SubjectPhenotype_meta <- read_xlsx(path = "data/GTEx_Analysis_v10_Annotations_SubjectPhenotypesDD.xlsx",na = "")

## load the de-identifier data ----
SampleAttributes <- read.table(url("https://storage.googleapis.com/adult-gtex/annotations/v10/metadata-files/GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt"), quote = "", fill=TRUE, header = T, sep = '\t')
SubjectPhenotype <- read.table(url("https://storage.googleapis.com/adult-gtex/annotations/v10/metadata-files/GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt"),fill=TRUE, header = T, sep = '\t')
save(SampleAttributes, file = "data/GTEx_Analysis_v10_Annotations_SampleAttributesDS_txt.Rdata")
save(SubjectPhenotype, file = "data/GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS_txt.Rdata")
load("data/GTEx_Analysis_v10_Annotations_SampleAttributesDS_txt.Rdata")
load("data/GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS_txt.Rdata")

## load the bulk liver data ----
# "Gene TPM by tissue" --> download the liver one
Liver_expression_v10 <- read.table('data/gene_tpm_v10_liver.gct.gz',header = T,sep = '\t',skip = 2)
Liver_expression_v10 <- t(Liver_expression_v10)
colnames(Liver_expression_v10) = Liver_expression_v10[2,]

Liver_expression_v10 <- Liver_expression_v10[c(-1,-2),]
Liver_expression_v10 <- as.data.frame(Liver_expression_v10)
names.use <- colnames(Liver_expression_v10)[(colnames(Liver_expression_v10) %in% all_colnames)] # 167 IN TOTAL
Liver_expression_ins <- Liver_expression_v10[,names.use]
rm(Liver_expression_v10)

save(Liver_expression_ins, file = "data/GTEx_Analysis_v10_raw.Rdata")

## extract ID info ----
rows <- rownames(Liver_expression_ins) # # Extract colnames from the 3rd to the last
extracted_ids <- sub("^[^.]+\\.([^.]+)\\..*$", "\\1", rows) # Extract the part between the first and second dots
Liver_expression_ins <- cbind(ID= extracted_ids, Liver_expression_ins)

SubjectPhenotype <- SubjectPhenotype %>% 
  mutate(ID = sub(".*-", "", SUBJID)) %>% 
  select(-SUBJID)
Liver_expression_complt <- right_join(SubjectPhenotype, Liver_expression_ins, by = "ID")
Liver_expression_complt <- Liver_expression_complt %>% mutate(across(4:171, as.numeric))

str(Liver_expression_complt)
# save the data
save(Liver_expression_complt, file = "GTEx_Analysis_v10_raw_complete.Rdata")


##################################
# 03 calculate the statistics ----
##################################

# 1 Male, 2 Female 
hist_sex <- ggplot(Liver_expression_complt, aes(x = factor(SEX, levels = c(1, 2), labels = c("Male", "Female")))) +
  geom_bar(fill = "gray", alpha = 0.7) +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  labs(title = "", x = "SEX", y = "Frequency") +
  coord_cartesian(ylim = c(0, 200)) +
  theme_bw()
hist_sex

hist_age <- ggplot(Liver_expression_complt, aes(x = factor(AGE))) +
  geom_bar(fill = "gray", alpha = 0.7) +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  labs(title = "", x = "AGE", y = "Frequency") +
  coord_cartesian(ylim = c(0, 120)) +
  theme_bw()
hist_age

hist_death <- ggplot(Liver_expression_complt, aes(x = factor(DTHHRDY))) +
  geom_bar(fill = "gray", alpha = 0.7) +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  labs(title = "", x = "Death category", y = "Frequency") +
  coord_cartesian(ylim = c(0, 120)) +
  theme_bw()
hist_death

# Specify the columns to summarize
vars_to_summarize <- c("SEX", "AGE", "DTHHRDY")

summary_table <- Liver_expression_complt %>%
  select(all_of(vars_to_summarize)) %>%
  mutate(across(everything(), as.character)) %>% 
  pivot_longer(cols = everything(), names_to = "variable", values_to = "category") %>%
  group_by(variable, category) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(variable) %>%
  mutate(
    N = sum(n),
    percent = round(n / N * 100, 1),
    label = paste0(n, "/", N, " (", percent, "%)")
  )

summary_table

#############################
# 04 Summarize the genes ----
#############################

# Table S1: selected DMET genes
# get the name of gene of interest; use CYP family for example
load("data/expression_profile.Rdata")
paste(colnames(expression_profile[[1]]), collapse = ", ")
