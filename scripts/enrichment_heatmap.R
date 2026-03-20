if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("enrichplot")
BiocManager::install("clusterProfiler")
BiocManager::install("AnnotationDbi")
BiocManager::install("org.Mm.eg.db")

library(clusterProfiler)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(DOSE)
library(dplyr)
library(tidyr)
library(tibble)
library(RColorBrewer)

install.packages("readr")

setwd("")

#CSV files with the consistently identified genes by GPcounts and masigpro, only GPcounts and only masigpro
#separated by sex

consistent <- read.csv("GO_COMBINED.csv")
female <- read.csv("GO_FEMALE.csv")
male <- read.csv("GO_MALE.csv")

consistent_combined <- consistent["Combined.Gene.Names"]
consistent_masigpro <- consistent["maSigPro.Gene.Names"]
consistent_gpcounts <- consistent["GPcounts.Gene.Names"]
female_combined <- female["Combined.Gene.Names"]
female_masigpro <- female["maSigPro.Gene.Names"]
female_gpcounts <- female["GPcounts.Gene.Names"]
male_combined <- male["Combined.Gene.Names"]
male_masigpro <- male["maSigPro.Gene.Names"]
male_gpcounts <- male["GPcounts.Gene.Names"]

gene_list <- list(GP_Male = as.character(male_gpcounts$X.3),
                  GP_Female = as.character(female_gpcounts$GPcounts.Gene.Names),
                  GP_Combined = as.character(consistent_gpcounts$GPcounts.Gene.Names),
                  MSP_Male = as.character(male_masigpro$X.2),
                  MSP_Female = as.character(female_masigpro$maSigPro.Gene.Names),
                  MSP_Combined = as.character(consistent_masigpro$maSigPro.Gene.Names),
                  C_Combined = as.character(consistent_combined$Combined.Gene.Names),
                  C_GP = as.character(consistent_gpcounts$GPcounts.Gene.Names),
                  C_MSP = as.character(consistent_masigpro$maSigPro.Gene.Names)
                  )


library(clusterProfiler)
library(org.Mm.eg.db)

enrichment <- compareCluster(geneCluster = gene_list,
                             fun= "enrichGO",
                             OrgDb = org.Mm.eg.db,
                             keyType = "SYMBOL",
                             ont = "BP"
                             )
head(enrichment)

library(dplyr)
library(tidyr)
library(tibble)
library(RColorBrewer)

go_df <- as.data.frame(enrichment)

top_terms <- go_df %>%
  group_by(ID) %>%
  summarise(min_p = min(p.adjust)) %>%
  arrange(min_p) %>%
  slice(1:30) %>%
  pull(ID)

filtered <- go_df %>%
  filter(ID %in% top_terms) %>%
  mutate(UniqueID = paste(Description, Cluster, sep = "_"))

desired_order <- c("GP_Male", "GP_Female", "GP_Combined",
                   "MSP_Male", "MSP_Female", "MSP_Combined",
                   "C_GP", "C_MSP", "C_Combined")

heat_data <- filtered %>%
  select(Description, Cluster, p.adjust) %>%
  mutate(logp = -log10(p.adjust)) %>%  # Create the value column
  select(Description, Cluster, logp) %>%  # Only keep the needed columns
  pivot_wider(names_from = Cluster, values_from = logp, values_fill = 0) %>%
  column_to_rownames("Description")

filtered <- filtered %>%
  mutate(UniqueID = paste(Description, Cluster, sep = "_"))

library(pheatmap)

colnames(heat_data) <- c(
  "GPcounts (Male, Unique Genes)",
  "GPcounts (Female, Unique Genes)",
  "GPcounts (Consistent Genes)",
  "maSigPro (Male, Unique Genes)",
  "maSigPro (Female, Unique Genes)",
  "maSigPro (Consistent Genes)",
  "GPcounts & maSigPro (Consistent Genes)",
  "GPcounts (All Genders, Unique Genes)",
  "maSigPro (All Genders, Unique Genes)"
)

pheatmap(
  as.matrix(heat_data),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  main = "Method dependent GO enrichment",
  color = colorRampPalette(brewer.pal(9, "YlOrRd"))(100) # this removes the x-axis (column) labels
)

str(heat_data)
  