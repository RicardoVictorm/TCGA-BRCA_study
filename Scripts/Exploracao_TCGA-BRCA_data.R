# Análise exploratória dos dados de Câncer-RNASeq - Case study: TCGA-BRCA

## Part. 1 - Exploração das informações dos dados do TCGA

# Carregar pacotes e libs para exploração
library(finalfit)
library(dplyr)
library(skimr)

###-----------------------------------------------------------------------------------------------------###

# Converta o objeto S4 em um dataframe padrão do R
expression_data_df <- as.data.frame(colData(expression_data))

# Crie um objeto "skim" a partir do dataframe
data_summary <- skim(expression_data_df)

# Visualize o resumo estatístico
data_summary


###-----------------------------------------------------------------------------------------------------###

## Part. 2 - Create elegant final results tables with finalfit (error)

# Criar tabelas cruzando informações relevantes

explanatory = c("sample_type.factor", "race.factor", "gender.factor")
dependent = 'vital_status.factor'
colon_s %>%
  summary_factorlist(dependent, explanatory,) -> t1
knitr::kable(t1, align=c("l", "l", "r", "r", "r"))


###-----------------------------------------------------------------------------------------------------###

## Part. 3 - Preprocessing to DE and multidimensional scaling (MDS) plot do Glimma

library(TCGAbiolinks)
library(Glimma)
library(edgeR)
library(DESeq2)


# Perform DE analysis

ddsSE <- DESeqDataSet(expression_data, design = ~ paper_IDH.status)
