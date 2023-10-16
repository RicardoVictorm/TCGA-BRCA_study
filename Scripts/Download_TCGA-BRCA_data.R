# TCGAbiolinks package download_data with a Script - Case study: TCGA-BRCA

## Part. 1 - obtenção dos dados do TCGA

###-----------------------------------------------------------------------------------------------------###

# Carregar TCGABiolinks e outros pacotes
library(TCGAbiolinks)
library(SummarizedExperiment)

# Escolha do tipo de dados, como RNA-Seq ou microarray
data_type <- "RNASeq"

# Especificar o tipo de dados de expressão que se deseja, como dados de contagem ou TPM
expression_type <- "count"

# Especificar o tipo de câncer de mama que se deseja, como "Breast Invasive Carcinoma (BRCA)"
cancer_type <- "TCGA-BRCA"

# Obtendo uma lista de projetos disponíveis para o tipo de câncer e tipo de dados especificados
projects <- GDCquery(project = cancer_type, data.category = data_type, data.type = expression_type)

# Download dos dados para um diretório local

query <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "STAR - Counts",
                  sample.type = c("Primary Tumor", "Solid Tissue Normal"))

GDCdownload(query, files.per.chunk = 100)

# Carregue os dados em um dataframe

expression_data  <- GDCprepare(query = query,
                              save = TRUE, 
                              save.filename = "BRCA_Exp.rda")

# obter info. dos subtype
infomation.subtype <- TCGAquery_subtype(tumor = "BRCA")

# obter clinical data
information.clinical <- GDCquery_clinic(project = "TCGA-BRCA",type = "clinical") 

# amostras de "Primary Tumor"
samples.primary.tumour <- expression_data$barcode[expression_data$shortLetterCode == "TP"]

# amostras de "solid tissue normal"
samples.solid.tissue.normal <- expression_data$barcode[expression_data$shortLetterCode == "NT"]

###-----------------------------------------------------------------------------------------------------###

## Part. 2 - Exploração inicial

# Visualize os primeiros registros do dataframe
head(expression_data)

# Informações sobre os genes 
genes.info <- rowRanges(expression_data)
genes.info

# Criar um "datatable" com info. das amostras
datatable(
  as.data.frame(colData(expression_data)), 
  options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
  rownames = FALSE
)

###-----------------------------------------------------------------------------------------------------###

## Part. 3 - uso do TCGAanalyze_Preprocessing (optional)

# Prepare expression matrix with geneID in the rows and samples (barcode) in the columns
BRCA_Matrix <- assay(expression_data,"unstranded")

# For gene expression if you need to see a boxplot correlation and AAIC plot to define outliers you can run
BRCA_RNAseq_CorOutliers <- TCGAanalyze_Preprocessing(expression_data)


###-----------------------------------------------------------------------------------------------------###rm
