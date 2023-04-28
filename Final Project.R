setwd("/Users/matthewgenchev/Desktop/USC/qbio_490_matthewgenchev/qbio_490_final")

dir.create("/Users/matthewgenchev/Desktop/USC/qbio_490_matthewgenchev/qbio_490_final/outputs")
setwd("outputs")

if (!require(BiocManager)){
  library('BiocManager')
}
if (!require(TCGAbiolinks)){
  library('TCGAbiolinks')
}
if (!require(maftools)){
  library('maftools')
}
BiocManager::install("DESeq2")
library(DESeq2)

BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)

if(!require(SummarizedExperiment)) {
  installed.packages("SummarizedExperiment")
}

clin_query <- GDCquery(project = "TCGA-OV", 
                       data.category = "Clinical", 
                       file.type = "xml")
GDCdownload(clin_query)
clinic <- GDCprepare_clinic(clin_query, 
                            clinical.info = "patient")
colnames(clinic)[ colnames(clinic) == "bcr_patient_barcode" ] <- "Tumor_Sample_Barcode"
#write.csv(clinic, "TCGAclinical.csv", row.names=FALSE, quote=FALSE) 

maf_query <- GDCquery(
  project = "TCGA-OV", 
  data.category = "Simple Nucleotide Variation", 
  access = "open",
  data.type = "Masked Somatic Mutation", 
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)
GDCdownload(maf_query)
maf <- GDCprepare(maf_query)

maf_object <- read.maf(maf = maf, 
                       clinicalData = clinic,
                       isTCGA = TRUE)

rna_query <- GDCquery(project ="TCGA-BRCA",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "STAR - Counts")
GDCdownload(rna_query)
rna_se <- GDCprepare(rna_query)

#Creating an oncoplot to find most mutated genes
oncoplot(maf = maf_object,
         top = 20)
ggsave("OVOncoplot.jpg", device = "jpeg")

#Creating a somatic interactions plot to look at 
par(mar = c(5,5,5,5))
CoMutation_plot <- somaticInteractions( maf = maf_object,
                                        top = 20,
                                        pvalue = c(0.05, 0.1),
                                        fontSize = 0.5,
                                        countsFontSize = 0.5)
ggsave("CoMutationPlot.jpg", device = "jpeg")

CoMutation_plot
write.csv(CoMutation_plot, "CoMutationPlot.csv", row.names = FALSE, quote = FALSE)

#CReating relevant dataframes for my Difefrential expression analysis
rna_clinical <- rna_se@colData
rna_clinical <- as.data.frame(rna_clinical)

rna_genes <- rna_se@rowRanges@elementMetadata
rna_genes <- as.data.frame(rna_genes)

rna_counts <- rna_se@assays@data$unstranded
rna_counts <- as.data.frame(rna_counts)

#adding infromative column and row names to data
row.names(rna_clinical) <- rna_clinical$barcode
row.names(rna_genes) <- rna_genes$gene_id
colnames(rna_counts) <- rna_clinical$barcode
rownames(rna_counts) <- rna_genes$gene_id

#Factoring columns that we are going to hold consistent--we know that Pathological stage and age influence transcriptomic profiles
#So we need to hold these consistent during our differential expression analysis
rna_clinical$ajcc_pathologic_stage <-factor(rna_clinical$ajcc_pathologic_stage)

#I want to subset the age into old and young by dividing around the median
age_na_mask <- !is.na(rna_clinical$age_at_index)
ages <- rna_clinical$age_at_index[age_na_mask]
median(ages) #The Median is 58

age_mask <- ifelse(rna_clinical$age_at_index >= 58, "old", "young")
rna_clinical$age_status <- age_mask
rna_clinical$age_status <- factor(rna_clinical$age_status)

#I want to check if vital status expression is different
rna_clinical$vital_status <- factor(rna_clinical$vital_status)

#I need to get rid of any NA values, so I am checking for NAs
sum(is.na(rna_clinical$ajcc_pathologic_stage)) #13 NAs
sum(is.na(rna_clinical$age_status)) #0 NAs
sum(is.na(rna_clinical$vital_status)) #0 NAs

#I am going to get rid of these patients in both rna_clinical and rna_counts
pathologic_mask <- !is.na(rna_clinical$ajcc_pathologic_stage)
rna_clinical <- rna_clinical[pathologic_mask, ]

rna_counts <- rna_counts[ , pathologic_mask]

#Going to pre-process the dataframes to get rid of insignificant (<10) transcript counts
row_sums <- rowSums(rna_counts)
low_counts_mask <- ifelse(row_sums < 10, F, T)
rna_counts <- rna_counts[low_counts_mask, ]
rna_genes <- rna_genes [low_counts_mask, ]

#Cretaing a DESeq DatSet object (dds)
dds <- DESeqDataSetFromMatrix(countData = rna_counts,
                              colData = rna_clinical,
                              design = ~ajcc_pathologic_stage + age_status + vital_status)
dds_obj <- DESeq(dds)

results <- results(dds_obj, format = "DataFrame", contrast = c("vital_status", "Dead", "Alive"))

#Making a dataframe that we can now use to create the volcano plot
modified_results <- data.frame(rna_genes$gene_name, 
                               results@rownames, 
                               results@listData$log2FoldChange, 
                               results@listData$pvalue,
                               results@listData$padj,
                               -log10(results@listData$padj))
names_mask <- ifelse(rna_genes$gene_id %in% results@rownames, T, F)
modified_results <- modified_results[names_mask, ]
colnames(modified_results) <- c("gene_name", "ensembl", "log2_fold_change", "p_value", "p_adjusted", "-log10_p_adjusted")
write.csv(modified_results, "Volcano_Plot_Data.csv", row.names = TRUE)

#Creating the Volcano plot
par(mar=c(4,4,4,4))
EnhancedVolcano(modified_results,
                title = "Gene Expression in Dead vs Alive Patients",
                lab = modified_results$gene_name,
                x = 'log2_fold_change',
                y = 'p_adjusted',
                xlim = c(-5, 5),
                ylim = c(0, 45),
                labSize = 4,
                boxedLabels = TRUE,
                drawConnectors = TRUE,
                max.overlaps = 30)
ggsave("VolcanoPlot.jpg", device = "jpeg")
