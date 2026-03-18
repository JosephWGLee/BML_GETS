#Reproduction of wet-lab Masigpro Data Analysis 
#Data files are a representation as current research has not currently been published

install.packages("BiocManager")
BiocManager::install("maSigPro")
BiocManager::available()

library(maSigPro)

rm(list = ls())

all_sample_table = read.csv("data/all_sampletable.csv",header=TRUE,row.names = 1)
all_data_table = read.csv("data/all _normalised_counts.csv", header = TRUE, row.names = 1 )
male_sample_table = read.csv("data/male_sampletable.csv",header=TRUE, row.names = 1)
male_data_table = read.csv("data/male_normalised_counts.csv", header = TRUE, row.names = 1, sep = ";")
female_sample_table = read.csv("data/female_sampletable.csv",header=TRUE,row.names = 1)
female_data_table = read.csv("data/female_normalised_counts.csv", header = TRUE, row.names = 1, sep = ";")

#Converting to numeric matrices 

all_data_table <- as.matrix(all_data_table)
male_data_table <- as.matrix(male_data_table)
female_data_table <- as.matrix(female_data_table)

class(all_data_table) <- "numeric"

# Removing all lines in your data with values 0 for all samples

all_data_table <- all_data_table[rowSums(all_data_table) >0, ]
male_data_table <- male_data_table[rowSums(male_data_table) >0, ]
female_data_table <- female_data_table[rowSums(female_data_table) >0, ]

#Masigpro reproduction of wet-lab analysis - 4th degree polynomial regression
#Polynomial regression 
#Making design matrix from tabular file with sample information

design_all <- make.design.matrix(all_sample_table, degree = 4)
design_male <- make.design.matrix(male_sample_table, degree = 4)
design_female <- make.design.matrix(female_sample_table, degree = 4)

fit_all <- p.vector(all_data_table, design_all, Q=0.05, MT.adjust = "BH", min.obs = 15)
fit_male <- p.vector(male_data_table, design_male, Q=0.05, MT.adjust = "BH", min.obs = 15)
fit_female <- p.vector(female_data_table, design_female, Q=0.05, MT.adjust = "BH", min.obs = 15)
write.table(fit_all$SELEC, file = "all_masigpro_significant_genes.csv")
write.table(fit_male$SELEC, file = "male_masigpro_significant_genes.csv")
write.table(fit_female$SELEC, file = "female_masigpro_significant_genes.csv")

all(colnames(all_data_table) %in% rownames(design_all))

head(all_sample_table)
head(all_data_table)
