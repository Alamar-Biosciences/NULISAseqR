####################################################
##### NULISAseq Data Analysis Demo Script #####
##### 26 November 2024
###################################################
## install NULISAseqR package
# https://github.com/Alamar-Biosciences/NULISAseqR
install.packages('devtools')
devtools::install_github('Alamar-Biosciences/NULISAseqR',
                         ref = 'main'
)

## install other packages
install.packages("tidyverse")
install.packages("readxl")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")

## load packages
library(NULISAseqR)
library(tidyverse)
library(ComplexHeatmap)

###################################################
## read in data
###################################################
# XML files -- Data matrix has raw counts
file <- c('detectability_P1_Tr03.xml', 'detectability_P2_Tr03.xml')
runs <- lapply(file, function(i) loadNULISAseq(file = i))

# Read in xml then perfom normalization
# plate1 <- readNULISAseq(file='detectability_P1_Tr03.xml')
# plate2 <- readNULISAseq(file='detectability_P2_Tr03.xml')
# 
# plate1_IC_normalized <- intraPlateNorm(data_matrix = plate1$Data,
#                                        IC = plate1$IC)
# plate2_IC_normalized <- intraPlateNorm(data_matrix = plate2$Data,
#                                        IC = plate2$IC)
# 
# IPC_normalized <- interPlateNorm(data_list = list(plate1_IC_normalized$normData,
#                                                   plate2_IC_normalized$normData),
#                                  IPC_wells = list(plate1$IPC, plate2$IPC))   
# for reverse curve targets (CRP and APOE in CNS Disease Panle), use transformReverse parameter to specify the targets
# IPC_normalized <- interPlateNorm(data_list = list(plate1_IC_normalized$normData,
#                                                   plate2_IC_normalized$normData),
#                                  IPC_wells = list(plate1$IPC, plate2$IPC),
#                                  transformReverse = c("CRP", "APOE"))  
#
# plate1_NPQ <- IPC_normalized$log2_interNormData[[1]]
# plate2_NPQ <- IPC_normalized$log2_interNormData[[2]]

# CSV files
## first option (with readNULISAseq())
## generate a list
data_all <- readNULISAseq("Alamar_NULISAseq_Detectability_NPQ.csv",
                          file_type = "csv_long")

## second option (with read_csv())
## metadata included
data <- read_csv("Alamar_NULISAseq_Detectability_NPQ.csv") 

###################################################
## transform data into wide format
###################################################
data_mCherry_log2 <-  data %>%  
  mutate(NPQ = as.numeric(NPQ)) %>% 
  filter(!grepl("pooled", SampleType)) %>%  
  select(Target, 
         sampleName_long,
         NPQ) %>%   
  distinct(.keep_all = T) %>%   
  pivot_wider(names_from = sampleName_long, values_from = NPQ) %>%   
  column_to_rownames("Target")

## create metadata dataframe
metadata <- data[, 1:8] %>% 
  filter(SampleType == "plasma") %>% 
  distinct() %>% 
  mutate(disease_type = factor(disease_type, levels = c("normal", "inflam", "cancer", 
                                                        "kidney", "metab", "neuro")))

## read in metadata/sample information if available
# reading metadata
# metadata <- readxl::read_xlsx("Sample_Info.xlsx") %>% 
#   filter(`Sample Type` != 0)
# 
# data %>% left_join(metadata) -> data

## read in detectability data if available
# detect <- read_csv("target_detectability_table.csv")
# 
# ## using detectability >50% cut off
# detect %>%   
#   filter(`Overall: serum (n = 176)` > 50) %>%   
#   pull(...1) -> target_passed_detct_qc_serum

###################################################
# Heatmaps & PCA
###################################################
# function to generate colors for the heatmap annotations
generate_cov_colors <- function(cov, data, set) {
  
  unique_groups <- sort(unique(data[,cov]))
  unique_groups_length <- length(unique_groups)
  
  values_cols <- RColorBrewer::brewer.pal(unique_groups_length, name = set)[1:unique_groups_length]
  names(values_cols) <- unique_groups
  values_cols
}
# function to generate a heatmap using the complex heatmap package
generate_heatmap <- function(data, samples, targets, 
                             annotate_sample_by = NULL,
                             column_split_by = NULL, fontsize = 4,
                             column_split_num = NULL,
                             name = NULL,
                             ...
) {
  
  scaled_data <- t(scale(t(data[targets, samples]),
                         center=TRUE, scale=TRUE))
  
  h1 <- ComplexHeatmap::Heatmap(
    name = if(is.null(name)) "scale" else name,
    scaled_data,
    top_annotation =  if(!is.null(annotate_sample_by)) HeatmapAnnotation(df = sample_annotation[samples,annotate_sample_by, drop = F],
                                                                         col = heatmap_annotation_colors),
    row_names_gp = gpar(fontsize = fontsize),
    column_names_gp = gpar(fontsize = fontsize),
    column_split = if(all(!is.null(column_split_by), is.null(column_split_num))) sample_annotation[samples,column_split_by, drop = F] else column_split_num,
    ...
  )
  
  draw(h1, background = "transparent")
  
}
# function to run and generate pca biplot
run_pca <- function(data, samples, targets, 
                    annotate_sample_by = NULL, ... ) {
  
  scaled_data <- t(scale(t(data[targets, samples]),
                         center=TRUE, scale=TRUE))
  
  
  pca_res <- PCAtools::pca(scaled_data, 
                           metadata = sample_annotation[colnames(scaled_data),],
                           center = TRUE,
                           scale = TRUE
  )
  
  PCAtools::biplot(pca_res, colby = annotate_sample_by, lab=NULL,
                   hline = 0, 
                   vline = 0,
                   pointSize = 2,
                   legendPosition = 'right', ...) + 
    theme(
      text = element_text(size = 25),
      panel.background = element_rect(fill='white'),
      plot.background = element_rect(fill='transparent', color = NA),
      legend.background = element_rect(fill='transparent'),
      legend.key = element_rect(fill = "transparent", color = NA)
    )
  
}
# function to generate heatmaps and PCA plots
generate_heatmap_and_pca <- function(data_log2, sample_list, targets, 
                                     annotate_by, col_annotate_by,
                                     pca_annotate_by, shape_by = NULL, size = 3, 
                                     col, output_prefix) {
  
  scaled_data <- scale(t(data_log2[targets, sample_list]),
                       center=TRUE, scale=TRUE) %>% as.data.frame() 
  
  columns_with_all_nan <- apply(scaled_data, 2, function(col) all(is.nan(col)))
  targets_with_all_nan <- names(columns_with_all_nan)[columns_with_all_nan]
  targets_used <- targets[!(targets %in% targets_with_all_nan)] 
  
  pdf(paste0(output_prefix, "_Heatmap.pdf"), height = 9, width = 10)
  generate_heatmap(
    data_log2,
    samples = sample_list,
    targets = targets_used,
    annotate_sample_by = annotate_by,
    name = "Z-Score",
    fontsize = size,
    row_title_rot = 0,
    clustering_method_columns = "ward.D2",
    row_split = 2,
    cluster_column_slices = FALSE
  )
  dev.off()
  
  pdf(paste0(output_prefix, "_Heatmap_split.pdf"), height = 11, width = 13)
  generate_heatmap(
    data_log2,
    samples = sample_list,
    targets = targets_used,
    annotate_sample_by = annotate_by,
    column_split_by = col_annotate_by,
    name = "Z-Score",
    fontsize = size,
    row_title_rot = 0,
    clustering_method_columns = "ward.D2",
    row_split = 2,
    cluster_column_slices = FALSE
  )
  dev.off()
  
  plot1 <- run_pca(
    data_log2,
    samples = sample_list,
    targets = targets_used,
    annotate_sample_by = pca_annotate_by,
    shape = shape_by,
    encircle = TRUE,
    encircleFill = TRUE
  ) +
    scale_fill_brewer(palette = col) +
    scale_color_brewer(palette = col)
  
  pdf(paste0(output_prefix, "_PCA_plot.pdf"), width = 9, height = 7)
  print(plot1)
  dev.off()
  
  tiff(paste0(output_prefix, "_PCA_plot.tiff"), width = 900, height = 700)
  print(plot1)
  dev.off()
  
  return (list(targets_used = targets_used))
}

# define sample list
metadata %>% pull(sampleName_long) -> full_sample_list   # all samples

# define target list
data %>% pull(Target) %>% unique() -> full_target_list   # all targets

# prepare sample annotation dataframe for the heatmap function
metadata |> 
  filter(sampleName_long %in% full_sample_list) |> 
  column_to_rownames("sampleName_long") -> sample_annotation

sample_annotation$disease_type <- factor(sample_annotation$disease_type, levels = c("normal", "inflam", "cancer", 
                                                                                    "kidney", "metab", "neuro"))

# generate sample annotation colors for the heatmap function
different_sets <- rep_len(c("Set1", "Set2", "Set1", "Set3"), length.out = length(colnames(sample_annotation)))

heatmap_annotation_colors <- mapply(generate_cov_colors, cov = colnames(sample_annotation), 
                                    set = different_sets, MoreArgs = list(data = sample_annotation))

# All samples Heatmap & PCA
plasma_targets_used <- generate_heatmap_and_pca(data_log2 = data_mCherry_log2, 
                                                sample_list = full_sample_list, 
                                                targets = rownames(data_mCherry_log2), # using all targets
                                                annotate_by = c("disease_type", "sex"), 
                                                col_annotate_by = "disease_type", 
                                                pca_annotate_by = "disease_type", 
                                                shape_by = "sex",
                                                col = "Set1", 
                                                output_prefix = "./All_samples")

###################################################
## Differential Expression Analysis
###################################################
lmTest <- lmNULISAseq(data = as.matrix(data_mCherry_log2[full_target_list, full_sample_list]),
                      sampleInfo = metadata %>% filter(sampleName_long %in% full_sample_list),
                      sampleName_var = "sampleName_long",
                      modelFormula = "disease_type + sex")

## another function: lmerNULISAseq for linear mixed effect model

## NOTE: lmNULISAseq with a single binary covariate is similar to t.test(var.equal=TRUE)
x <- rnorm(20, mean=2)
y <- rnorm(20, mean=0, sd=2)
group <- c(rep(1, length(x)), rep(0, length(y)))
t.test(x, y, var.equal = TRUE)
summary(lm(c(x,y) ~ group))

# generate volcano plot with FDR
## comparing Inflam vs Normal
pdf('./disease_inflam_vs_normal_volcano_plot_fdr.pdf', 
    width=5, height=5)
volcanoPlot(coefs = lmTest$modelStats$disease_typeinflam_coef,
            p_vals = lmTest$modelStats$disease_typeinflam_pval_FDR,
            target_labels = lmTest$modelStats$target)
dev.off()

## comparing Cancer vs Normal
pdf('./disease_cancer_vs_normal_volcano_plot_fdr.pdf', 
    width=5, height=5)
volcanoPlot(coefs = lmTest$modelStats$disease_typecancer_coef,
            p_vals = lmTest$modelStats$disease_typecancer_pval_FDR,
            target_labels = lmTest$modelStats$target)
dev.off()

###################################################
## paired samples or longitudinal data
## use linear mixed effect model
###################################################
# lmTest_pair <- lmerNULISAseq(data = as.matrix(data_mCherry_log2[full_target_list, full_sample_list]),
#                       sampleInfo = metadata %>% filter(sampleName_long %in% full_sample_list),
#                       sampleName_var = "sampleName_long",
#                       modelFormula_fixed = "disease_type + Time", # note Time is not in the current metadata, this is hypothetical example
#                       modelFormula_random = "(1|`SubjectID`)") # note SubjectID is not in the current metadata, this is hypothetical example
# #                   ## change SubjectID to the column name of the column denoting subject number in metadata 

