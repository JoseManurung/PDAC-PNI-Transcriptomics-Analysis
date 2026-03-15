# =================================================================================
# INTEGRATIVE TRANSCRIPTOMIC PROFILING OF PERINEURAL INVASION (PNI) SIGNATURES
# IN PANCREATIC DUCTAL ADENOCARCINOMA (PDAC): A MULTI-CONTRAST BIOINFORMATICS STUDY
# =================================================================================
#
# Author         : Yosia Jose Rasdiva Manurung
# Affiliation    : Diponegoro University (UNDIP), Indonesia
# Project Date   : Started on Sunday, 22 February 2026. Completed on Saturday, 08 March 2026.
# Last Update    : Sunday, 15 March 2026 (Added detailed Venn intersection extraction).
# License        : Copyright (c) 2026 [Yosia Jose Rasdiva Manurung]. All rights reserved.
#
# ---------------------------------------------------------------------------------
# UPDATE LOG:
# 12 March, 2026: Refined UMAP stability with set.seed() for reproducibility.
# 15 March, 2026: Implemented Part T.3 for automated extraction of annotated 
#                 Venn intersection data into professional CSV reports.
# ---------------------------------------------------------------------------------
#
# ------------------------------------------------------------------------------
# DATASET INFORMATION
# ------------------------------------------------------------------------------
# GEO accession : GSE102238
# Title         : Gene expression signatures associated with PNI in PDAC
# Database      : Gene Expression Omnibus (GEO, NCBI)
# Platform      : Agilent-052909 CBC_lncRNAmRNA_V3 (Probe Name version) (GPL19072)
# Organism      : Homo sapiens (Human)
# Disease       : Pancreatic Ductal Adenocarcinoma (PDAC)
# Data Size     : ~15.2 MB
# Samples       : 50 pairs (Tumor & Adjacent Normal) | n=100 total
# Status        : Published (Aug 04, 2017) | Last Updated (Jul 25, 2021)
# Citation      : Yang MW, Tao LY, Jiang YS, Yang JY et al. Perineural Invasion 
#                 Reprograms the Immune Microenvironment through Cholinergic 
#                 Signaling in Pancreatic Ductal Adenocarcinoma. Cancer Res 2020
#                 May 15;80(10):1991-2003. PMID: 32098780. 
#
# ------------------------------------------------------------------------------
# RESEARCH OBJECTIVES & COMPARISON GROUPS
# ------------------------------------------------------------------------------
# Primary Goal : Identify DEGs across 6 PNI-related clinical contrasts
#
# Contrasts Explored:
#   1. PNI Effect (Tumor)   : PNI-positive Tumor  vs. PNI-negative Tumor
#   2. PNI Effect (Normal)  : PNI-positive Normal vs. PNI-negative Normal
#   3. Tumor vs. Normal (+) : PNI-positive Tumor  vs. PNI-positive Normal
#   4. Tumor vs. Normal (-) : PNI-negative Tumor  vs. PNI-negative Normal
#   5. Extreme Contrast     : PNI-positive Tumor  vs. PNI-negative Normal
#   6. Reverse Contrast     : PNI-negative Tumor  vs. PNI-positive Normal
#
# ==============================================================================

# ==================================================
# PART A. ENVIRONMENT SETUP & PROJECT INITIALIZATION
# ==================================================
# -----------------------------------------------
# A.1. SYSTEM DEPENDENCIES & PACKAGE INSTALLATION
# -----------------------------------------------
# A.1.1. Initialize BiocManager
# Checks for the presence of BiocManager and installs it if missing
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# A.1.2. Install Core Bioconductor Suites
# Deploys essential libraries for GEO data retrieval, linear modeling, and gene annotation
BiocManager::install(c("GEOquery", "limma", "biomaRt","clusterProfiler", "org.Hs.eg.db", "enrichplot"), ask = FALSE, update = FALSE)

# A.1.3. Install CRAN Visualization and Utility Packages
# Downloads standard R packages for heatmaps, tidy data processing, and dimensional reduction
install.packages(c("pheatmap", "ggplot2", "dplyr", "umap", "tidyr", "ggVennDiagram"))

# --------------------------------
# A.2. GLOBAL SYSTEM CONFIGURATION
# --------------------------------
# A.2.1. Clean Global Environment
# Removes all existing objects from the workspace to prevent variable contamination
rm(list = ls())

# A.2.2. Disable Scientific Notation
# Forces R to display full decimal values (e.g., 0.00001) instead of exponential (1e-05)
options(scipen = 999)

# A.2.3. Set String Conversion Behavior
# Prevents R from automatically converting character strings into factor levels
options(stringsAsFactors = FALSE)

# A.2.4. Configure Graphical Device Background
# Ensures that saved plots maintain a consistent white background regardless of OS
options(bitmapType = 'cairo')

# ---------------------------------
# A.3. PROJECT DIRECTORY MANAGEMENT
# ---------------------------------
# A.3.1. Define Project Directories
# Creates a structured folder hierarchy to separate raw data, plots, and results
folders <- c("data", "plots", "results")

# A.3.2. Automated Directory Initialization
# Checks for existing directories and creates them if they do not exist
for (folder in folders) {
  if (!dir.exists(folder)) {
    dir.create(folder)
  }
}

# A.3.3. Verify Working Environment
# Prints the current working directory to ensure the project is correctly localized
getwd()

# ---------------------------------
# A.4. LOAD BIOINFORMATIC LIBRARIES
# ---------------------------------
# A.4.1. Load Required Analysis Libraries
# Imports the complete suite of Bioconductor and CRAN packages for the pipeline
library(GEOquery)        # For data acquisition from Gene Expression Omnibus (GEO)
library(limma)           # For linear modeling and differential expression analysis
library(biomaRt)         # For genome-wide annotation and biological ID mapping
library(tidyr)           # For data restructuring and tidying operations
library(dplyr)           # For efficient data transformation and manipulation
library(ggplot2)         # For high-performance statistical visualizations
library(umap)            # For non-linear dimensionality reduction and clustering
library(pheatmap)        # For generating publication-quality clustered heatmaps
library(grid)            # For layout control and clearing graphics device pages
library(ggVennDiagram)   # For visualizing gene intersections across multiple contrasts
library(clusterProfiler) # For automated functional enrichment (GO/KEGG) analysis
library(org.Hs.eg.db)    # For genome-wide annotation of Human (Homo sapiens)
library(enrichplot)      # For advanced visualization of enrichment results (Dotplots & Barplots)

# ===========================================================
# PART B. DATA ACQUISITION FROM GENE EXPRESSION OMNIBUS (GEO)
# ===========================================================
# B.1. Fetch dataset from GEO remote repository
# Downloads the ExpressionSet object including both phenotypic metadata and GPL annotations
gset <- getGEO("GSE102238", GSEMatrix = TRUE, AnnotGPL = TRUE, destdir = "data/")[[1]]

# B.2. Extract primary expression matrix
# Isolates the raw signal intensity data into a matrix for downstream statistical processing
ex <- exprs(gset)

# ======================================
# PART C. EXPRESSION DATA PRE-PROCESSING
# ======================================
# C.1. Evaluate logarithmic transformation requirements
# Performs a quantile distribution analysis to determine if the data is in linear or log scale
qx <- as.numeric(quantile(ex, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE))
LogTransform <- (qx[5] > 100)

# C.2. Execute log2 transformation conditionally
# Normalizes the intensity range and handles non-positive values to ensure data stability
if (LogTransform) {
  ex[ex <= 0] <- NA
  ex <- log2(ex)
}

# C.3. Export Normalized Full Expression Matrix
# Persists the complete log2-transformed dataset for all probes
write.csv(ex, "results/Normalized_Expression_Matrix_Full.csv", row.names = TRUE)

# ====================================================
# PART D. SAMPLE GROUPING & CLINICAL METADATA CURATION
# ====================================================
# D.1. Extract sample phenotype metadata
# Retrieves the phenoData (pData) slot containing clinical information for all samples
metadata <- pData(gset)

# D.2. Define PNI-based experimental groups from 'source_name_ch1'
# Implements regex-based string matching to classify samples based on PNI status and tissue type
group_labels <- pData(gset)$source_name_ch1
clean_groups <- ifelse(grepl("with PNI.*tumor", group_labels, ignore.case = TRUE), "PNI_Pos_Tumor",
                       ifelse(grepl("without PNI.*tumor", group_labels, ignore.case = TRUE), "PNI_Neg_Tumor",
                              ifelse(grepl("with PNI.*normal", group_labels, ignore.case = TRUE), "PNI_Pos_Normal",
                                     ifelse(grepl("without PNI.*normal", group_labels, ignore.case = TRUE), "PNI_Neg_Normal", NA))))

# D.3. Assign categorical factors to ExpressionSet
# Converts classifications into a factor with defined levels and attaches them to the gset object
gset$group <- factor(clean_groups,
                     levels = c("PNI_Pos_Tumor",
                                "PNI_Neg_Tumor",
                                "PNI_Pos_Normal",
                                "PNI_Neg_Normal"))

# D.4. Review sample distribution
# Generates a frequency table to confirm the balanced representation of clinical cohorts
table(gset$group)
levels(gset$group)

# D.5. Export Cleaned Clinical Metadata
# Saves the sample mapping and grouping for clinical reporting
write.csv(metadata, "results/Sample_Metadata_Cleaned.csv", row.names = TRUE)

# ==================================
# PART E. DESIGN MATRIX CONSTRUCTION
# ==================================
# E.1. Define experimental group levels
# Standardizes clinical group names by removing spaces and illegal characters
groups <- factor(gset$group)
levels(groups) <- make.names(levels(groups))

# E.2. Construct the Design Matrix
# Creates a model matrix without an intercept to represent each group independently
design <- model.matrix(~0 + groups)
colnames(design) <- levels(groups)

# E.3. Validate matrix columns
# Confirms the successful generation of the 4 distinct clinical condition columns
colnames(design)

# E.4. Define Contrast Matrix for linear modeling
# Specifies the 6 primary comparisons used to evaluate PNI and tissue effects
contrast_matrix <- makeContrasts(
  
  # E.4.1. Tumor vs Tumor (Effect of PNI on cancer cells)
  # Evaluates the transcriptional impact of PNI specifically within malignant tissue
  PNI_Effect_Tumor = PNI_Pos_Tumor - PNI_Neg_Tumor,
  
  # E.4.2. Normal vs Normal (Effect of PNI on healthy tissue)
  # Evaluates the transcriptional impact of PNI within non-malignant healthy tissue
  PNI_Effect_Normal = PNI_Pos_Normal - PNI_Neg_Normal,
  
  # E.4.3. Tumor vs Normal (In positive PNI conditions)
  # Compares malignant vs healthy tissue under positive PNI conditions
  Tumor_vs_Normal_Pos = PNI_Pos_Tumor - PNI_Pos_Normal,
  
  # E.4.4. Tumor vs Normal (In Negative PNI conditions)
  # Compares malignant vs healthy tissue under negative PNI conditions
  Tumor_vs_Normal_Neg = PNI_Neg_Tumor - PNI_Neg_Normal,
  
  # E.4.5. Extreme_Contrast (Positive Tumor vs Negative Normal)
  # Direct comparison between the most aggressive (PNI+ Tumor) and baseline (PNI- Normal) states
  Extreme_Contrast = PNI_Pos_Tumor - PNI_Neg_Normal,
  
  # E.4.6. Reverse_Contrast (Negative Tumor vs Positive Normal)
  # Direct comparison between the baseline tumor (PNI- Tumor) and (PNI+ normal) states
  Reverse_Contrast = PNI_Neg_Tumor - PNI_Pos_Normal,
  levels = design
)

# ===========================================================
# PART F. EXPERIMENTAL DESIGN VALIDATION & METADATA REPORTING
# ===========================================================
# F.1. Define metadata for automated reporting
# Creates a structured list to store comparison labels and group keys
comparisons <- list(
  list(pos = "PNI_Pos_Tumor",  neg = "PNI_Neg_Tumor",  label = "PNI Effect in Tumor Tissue"),
  list(pos = "PNI_Pos_Normal", neg = "PNI_Neg_Normal", label = "PNI Effect in Normal Tissue"),
  list(pos = "PNI_Pos_Tumor",  neg = "PNI_Pos_Normal", label = "Tumor vs Normal (PNI Positive)"),
  list(pos = "PNI_Neg_Tumor",  neg = "PNI_Neg_Normal", label = "Tumor vs Normal (PNI Negative)"),
  list(pos = "PNI_Pos_Tumor",  neg = "PNI_Neg_Normal", label = "Extreme_Contrast (Positive Tumor vs Negative Normal)"),
  list(pos = "PNI_Neg_Tumor",  neg = "PNI_Pos_Normal", label = "Reverse_Contrast (Negative Tumor vs Positive Normal)")
)

# F.2. Initialize summary report
# Prints a header for the automated PNI differential expression analysis log
cat("--- PNI Differential Expression Analysis Summary ---\n")

# F.3. Execute automated metadata extraction loop
# Iterates through the comparison list to log formulas and group identities
for (i in 1:length(comparisons)) {
  
  # F.3.1. Access current list element
  # Retrieves the indexed comparison sub-list for processing
  comp <- comparisons[[i]]
  
  # F.3.2. Extract group identifiers
  # Assigns the positive and negative keys to local variables for report generation
  target_group   <- comp$pos
  baseline_group <- comp$neg
  analysis_name  <- comp$label
  
  # F.3.3. Construct contrast formula string
  # Generates a human-readable representation of the linear contrast math
  contrast_formula <- paste(target_group, "-", baseline_group)
  
  # F.3.4. Print comparison metadata
  # Formats and outputs the analysis details to the R console
  cat(paste0("\n[Comparison ", i, "]\n"))
  cat(paste0("Analysis Type    : ", analysis_name, "\n"))
  cat(paste0("Contrast Formula : ", contrast_formula, "\n"))
  cat(paste0("Target (+)       : ", target_group, "\n"))
  cat(paste0("Baseline (-)     : ", baseline_group, "\n"))
  cat("---------------------------------------------------\n")
}

# =========================================================
# PART G. DIFFERENTIAL EXPRESSION ANALYSIS (LIMMA PIPELINE)
# =========================================================
# G.1. Build linear model
# Estimates the fold change and standard error by fitting a linear model to each gene
fit <- lmFit(ex, design)

# G.2. Apply contrasts to the model
# Computes estimated coefficients and standard errors for the specified comparisons
fit2 <- contrasts.fit(fit, contrast_matrix)

# G.3. Perform Empirical Bayes moderation
# Stabilizes variance estimates by borrowing information across all genes (shrinkage)
fit2 <- eBayes(fit2, proportion = 0.01)

# G.4. Verify contrast generation
# Confirms the dimensionality of the model to ensure all 6 contrasts were processed
ncol(fit2$coefficients)
ncol(contrast_matrix)

# ==================================================================
# PART H. INTEGRATION OF DIFFERENTIAL STATISTICS FOR SIX COMPARISONS
# ==================================================================
# H.1. Initialize integrative storage list
# Creates a temporary container to aggregate outputs while preventing memory fragmentation
contrast_results <- list()

# H.2. Extract and stratify statistical metrics for all 6 contrasts
# Iteratively captures differential profiles to ensure each comparison is uniquely identifiable
for (i in 1:6) {
  
  # H.2.1. Retrieve rank-based topTable results
  # Extracts full feature sets with FDR-adjusted p-values to control multiple testing errors
  tmp_tt <- topTable(fit2, coef = i, number = Inf, adjust.method = "fdr")
  
  # H.2.2. Isolate core metrics and apply contrast-specific labeling
  # Renames columns dynamically to prevent naming collisions during the global merge
  comp_name <- colnames(fit2)[i]
  
  tmp_tt <- tmp_tt[, c("logFC", "B", "adj.P.Val")]
  colnames(tmp_tt) <- c(paste0("logFC_", comp_name), 
                        paste0("B_", comp_name), 
                        paste0("adjP_", comp_name))
  
  # H.2.3. Persist ProbeID as an explicit relational key
  # Promotes row identifiers to a dedicated column to facilitate high-performance joining
  tmp_tt$ProbeID <- rownames(tmp_tt)
  contrast_results[[i]] <- tmp_tt
}

# H.3. Execute multi-way synthesis via functional Reduce
# Consolidates all 6 disparate datasets into a single unified master statistics matrix
combined_stats <- Reduce(function(x, y) merge(x, y, by = "ProbeID", all = TRUE), contrast_results)

# H.4. Export Integrated Statistical Matrix
# Saves the combined results of all 6 clinical contrasts into a single master file
write.csv(combined_stats, "results/Combined_Statistics_6_Contrasts_Raw.csv", row.names = FALSE)

# ===============================================================================
# PART I. GENE ANNOTATION – AGILENT GPL19072 (Sequence Recovery & Symbol Mapping)
# ===============================================================================
# I.1. Initialize biomaRt Connection
# Establishes a remote connection to the Ensembl Human Genome database (GRCh38)
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# I.2. Extract Probe IDs
# Retrieves unique probe identifiers from the differential expression statistics
probe_ids <- combined_stats$ProbeID

# I.3. Fetch Gene Annotations via biomaRt
# Maps Agilent v2 probe IDs to official HGNC external gene names and descriptions
gene_annotations <- getBM(
  attributes = c("agilent_wholegenome_4x44k_v2", "external_gene_name", "description"),
  filters    = "agilent_wholegenome_4x44k_v2",
  values     = probe_ids,
  mart       = mart
)

# I.4. Merge Statistics with Platform Metadata
# Integrates raw statistical results with feature data stored in the ExpressionSet
stats_with_metadata <- merge(
  combined_stats,     
  fData(gset),        
  by.x = "ProbeID",
  by.y = "ID",
  all.x = TRUE
)

# I.5. Final Integration of Annotations
# Combines the statistical metadata with the mapped Gene Symbols from biomaRt
topTable_Master <- merge(
  stats_with_metadata, 
  gene_annotations,    
  by.x = "ProbeID", 
  by.y = "agilent_wholegenome_4x44k_v2", 
  all.x = TRUE
)

# I.6. Column Sanitization & Formatting
# Standardizes column nomenclature and reorders key variables for better readability
topTable_Master <- topTable_Master %>%
  dplyr::rename(
    Symbol = external_gene_name,
    Gene_Description = description
  ) %>%
  dplyr::relocate(
    ProbeID,
    Symbol,
    Gene_Description,
    SEQUENCE,
    contains("logFC"),
    contains("adjP")
  )

# I.7. Generate Clean Dataset
# Filters out uninformative probes that lack an official Gene Symbol assignment
topTable_Master_Clean <- topTable_Master %>%
  filter(!is.na(Symbol) & Symbol != "")

# I.8. Export Annotated Master Datasets
# Persists the comprehensive annotated statistical table to the results directory
write.csv(topTable_Master, "results/Master_Annotated_Results_Full.csv", row.names = FALSE)

# I.9. Export Filtered Biological Dataset
# Saves the finalized list containing only probes with validated Gene Symbols
write.csv(topTable_Master_Clean, "results/Master_Results_Final_Cleaned.csv", row.names = FALSE)

# =====================================================================================
# PART J. GLOBAL VISUALIZATION PARAMETERS & CLINICAL COLOR MAPPING FOR ANALYTICAL PLOTS
# =====================================================================================
# J.1. Define custom color palette for clinical stratification
# Maps phenotype-specific colors (Dark Red: PNI+ Tumor, Light Red: PNI- Tumor, 
# Dark Green: PNI+ Normal, Light Green: PNI- Normal) to PNI status
my_pni_colors <- c("#8B0000", "#FF8282", "#006400", "#90EE90")

# J.2. Map group-specific colors to each individual sample
# Ensures visual identification of tumor (red) vs. normal (green) samples based on PNI status
sample_colors <- my_pni_colors[as.numeric(gset$group)]

# ======================================================================
# PART K. TIDY DATA PREPARATION (Specific to Boxplot & Density Analysis)
# ======================================================================
# K.1. Wide-to-Long Matrix Transformation
# Reshaping the expression matrix 'ex' into a tidy data frame for ggplot2
expr_long <- data.frame(
  Expression = as.vector(ex),
  SampleID   = rep(colnames(ex), each = nrow(ex)),
  Group      = rep(gset$group, each = nrow(ex))
)

# K.2. Numeric Integrity & Type Coercion
# Enforcing numeric data types to prevent factor artifacts during plotting
expr_long$Expression <- as.numeric(as.character(expr_long$Expression))

# K.3. Data Sanitization: Non-Finite Value Exclusion
# This ensures the density estimation algorithm executes without errors
expr_long_clean <- expr_long[is.finite(expr_long$Expression), ]

# K.4. Preprocessing Audit Log
# Reporting the total number of records excluded during the cleaning process
cat("Data Sanitization Complete. Rows removed:", nrow(expr_long) - nrow(expr_long_clean), "\n")

# K.5. Export Tidy Expression Data for Distribution Analysis
# Persists the long-format expression data for reproducibility in group-wise distribution audits
write.csv(expr_long_clean, "results/Tidy_Expression_Data_Long_Format_for_Boxplot_and_Density.csv", row.names = FALSE)

# K.6. Summary Statistics Export
# Generates and saves a descriptive summary of expression levels per clinical group
expression_summary <- expr_long_clean %>%
  dplyr::group_by(Group) %>%
  dplyr::summarise(
    Mean = mean(Expression),
    Median = median(Expression),
    SD = sd(Expression),
    Min = min(Expression),
    Max = max(Expression),
    Total_Observations = n()
  )

# K.7. Export Group-Wise Statistical Metrics
# Persists finalized descriptive statistics to CSV for methodological validation and reporting
write.csv(expression_summary, "results/Group_Expression_Statistics_Summary.csv", row.names = FALSE)

# ================================================
# PART L. BOXPLOT OF EXPRESSION VALUE DISTRIBUTION 
# ================================================
# L.1. Advanced Boxplot Construction with Customized Aesthetics
# Utilizing 'stat_boxplot' to integrate error bars (staples) and refining line geometry
boxplot <- ggplot(expr_long_clean, aes(x = SampleID, y = Expression, fill = Group)) +
  stat_boxplot(geom = "errorbar", width = 0.5, linewidth = 0.5) +  # Add whisker terminals
  geom_boxplot(outlier.shape = NA, linewidth = 0.5) +              # Hide outliers for clarity         
  
  # L.2. Color Mapping and Axis Calibration
  # Applying manual color scales and enforcing log2 expression range constraints  
  scale_fill_manual(values = my_pni_colors) +
  scale_y_continuous(breaks = seq(3, 15, 3)) +    # Defined intervals for precise reading
  coord_cartesian(ylim = c(3, 15)) +              # Hard-limit view without data removal
  
  # L.3. Thematic Refinement and Label Optimization
  # Implementing publication-ready aesthetics with high-density X-axis label rotation
  theme_bw(base_size = 14) +                  # Menambahkan bingkai kotak
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, 
                               size = 7, color = "black"), # Optimized for 100 samples
    axis.text.y = element_text(color = "black", size = 11),
    axis.ticks = element_line(color = "black", linewidth = 0.8), # Pronounced tick marks
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    panel.grid.major.x = element_blank(),    # Remove vertical grid lines for visual focus
    panel.grid.minor = element_blank(),
    legend.position = "right"
  ) +
  
  # L.4. Annotation and Title Specification
  # Defining descriptive labels for metadata groups and measurement units
  labs(
    title = "Distribution of Expression Values of 100 Samples",
    x = "Sample ID (GSM)",
    y = "Expression Value (log2)",
    fill = "Clinicopathological Group"
  )

# L.5. Graphical Rendering and Visual Audit
# Displaying the processed visualization within the active graphics device
print(boxplot)

# L.6. Automated Figure Export with High-Density Calibration
# Saving as 300 DPI PNG with expanded width to ensure label legibility
ggsave(
  filename = "plots/Boxplot_PDAC_100.png", 
  plot = boxplot, 
  width = 15, 
  height = 8, 
  dpi = 300
)

# ================================================
# PART M. DENSITY PLOT OF GENE EXPRESSION PROFILES
# ================================================
# M.1. Generate Density Plot with ggplot2
# Visualizing the global normalization and range of gene expression profiles
density_plot <- ggplot(expr_long_clean, aes(x = Expression, color = Group)) +
  geom_density(linewidth = 1.2) +             # Enhanced line for publication clarity
  scale_color_manual(values = my_pni_colors) + 
  scale_x_continuous(
    limits = c(0, 19),                        # Force X-axis to start from 0 up to 18
    breaks = seq(0, 18, 3),                   # Linear intervals every 3 units
    expand = expansion(mult = c(0, 0.05))     # Remove gap between Y-axis and 0
  ) +
  theme_bw(base_size = 14) +                  # Boxed theme to match Boxplot aesthetics
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    legend.position = "right",
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = "black"),
    axis.ticks = element_line(color = "black", linewidth = 0.8), # Pronounced tick marks
    axis.line = element_blank()                # Clean box border from theme_bw
  ) +
  labs(
    title = "Distribution of Gene Expression Profiles of 100 Samples",
    x = "Expression Value (log2)",
    y = "Density",
    color = "Clinicopathological Group"
  )

# M.2. Graphical Rendering
# Displaying the final distribution plot within the active R graphics device
print(density_plot)

# M.3. Automated Figure Export
# Saving the visualization as a 300 DPI PNG file in the designated plots folder
ggsave(
  filename = "plots/Density_PDAC_100.png", 
  plot = density_plot, 
  width = 15, 
  height = 8, 
  dpi = 300
)

# ==================================================
# PART N. UMAP (NON-LINEAR DIMENSIONALITY REDUCTION)
# ==================================================
# N.1. Prepare the input matrix
# UMAP requires observations (samples) as rows and features (genes) as columns
umap_input <- t(ex)

# N.2. Handle missing values
# UMAP algorithms cannot process NA values; filtering for complete cases only
ex_clean <- exprs(gset)
ex_clean <- ex_clean[complete.cases(ex_clean), ]

# N.3. Export Filtered Matrix for Dimensionality Reduction
# Saves the 'ex_clean' matrix with all NA values removed
write.csv(ex_clean, "results/Normalized_Expression_Matrix_Clean_for_UMAP.csv", row.names = TRUE)

# N.4. Transpose the cleaned expression matrix
# Ensuring the matrix orientation is correct for the UMAP function
umap_input <- t(ex_clean)

# N.5. Initialize Random Number Generator
# Ensures that stochastic processes like UMAP or t-SNE yield identical results every run
set.seed(123)

# N.6. Execute UMAP algorithm
# Projecting high-dimensional gene expression data into 2D space
umap_result <- umap(umap_input)

# N.7. Construct a structured data frame for visualization
# Extracting the 2D layout coordinates and appending clinical group metadata
umap_df <- data.frame(
  UMAP1 = umap_result$layout[, 1],
  UMAP2 = umap_result$layout[, 2],
  Group = pData(gset)$group
)

# N.8. Optional: Export Cleaned Input Matrix for UMAP
# Saves the transposed, NA-filtered expression matrix used as the UMAP algorithm input
write.csv(as.data.frame(umap_input), "results/UMAP_Input_Matrix_Transposed.csv", row.names = TRUE)

# N.9. Export UMAP Dimensionality Reduction Coordinates
# Persists the 2D layout coordinates alongside clinical metadata for reproducibility
write.csv(umap_df, "results/UMAP_Coordinates_and_Metadata.csv", row.names = TRUE)

# N.10. Generate UMAP Visualization
# Visualizing sample clustering patterns based on global transcriptional profiles
umap_plot <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = my_pni_colors) + # Consistent clinical color mapping
  theme_bw(base_size = 14) +                  # Professional boxed layout
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = "black"),
    axis.ticks = element_line(color = "black", linewidth = 0.8), # Pronounced ticks
    legend.position = "right"
  ) +
  labs(
    title = "UMAP Plot of Gene Expression Profiles of 100 Samples",
    x = "UMAP 1",
    y = "UMAP 2",
    color = "Clinicopathological Group"
  )

# N.11. Graphics Device Verification
# Rendering the visualization to audit sample clustering and group separation
print(umap_plot)

# N.12. Automated Publication-Quality Export
# Saving with naming convention optimized for consistency
ggsave(
  filename = "plots/UMAP_PDAC_100.png", 
  plot = umap_plot, 
  width = 15, 
  height = 8, 
  dpi = 300
)

# ==========================================================================
# PART O. STATISTICAL DATA EXTRACTION FOR VOLCANO & HEATMAP PLOT PREPARATION
# ==========================================================================
# O.1. Function for topTable Extraction
# Automates the retrieval of full statistical results and converts rownames to ProbeID
process_tt <- function(fit_obj, coef_idx) {
  tt <- topTable(fit_obj, coef = coef_idx, number = Inf, adjust.method = "fdr")
  tt$ProbeID <- rownames(tt)
  return(tt)
}

# O.2. Extract Comprehensive Statistical Tables
# Retrieves the top-ranked genes for each of the 6 clinical comparisons
tt_1 <- process_tt(fit2, 1) # PNI Effect in Tumor
tt_2 <- process_tt(fit2, 2) # PNI Effect in Normal
tt_3 <- process_tt(fit2, 3) # Tumor vs Normal (PNI_Pos)
tt_4 <- process_tt(fit2, 4) # Tumor vs Normal (PNI_Neg)
tt_5 <- process_tt(fit2, 5) # Extreme Contrast (PNI_Pos_Tumor vs. PNI_Neg_Normal)
tt_6 <- process_tt(fit2, 6) # Reverse Contrast (PNI_Neg_Tumor vs. PNI_Pos_Normal)

# O.3. Define the Robust Annotation and Cleaning Pipeline
# Merges Limma results with platform metadata and resolves probe redundancy
clean_and_annotate <- function(fit_obj, coef_idx, annotation_df) {
  
  # O.3.1. Raw Statistical Extraction
  # Pulls complete differential expression data for the specified contrast
  tt_data <- topTable(fit_obj, coef = coef_idx, number = Inf, adjust.method = "fdr")
  
  # O.3.2. Identifier Synchronization
  # Stores matrix row names in a new column to facilitate relational merging
  tt_data$MatchID <- rownames(tt_data)
  
  # O.3.3. Relational Database Merging
  # Joins statistical data with the Master Annotation table using ProbeID keys
  annotated <- merge(tt_data, annotation_df, by.x = "MatchID", by.y = "ProbeID", all.x = TRUE)
  
  # O.3.4. Biological Feature Filtering
  # Removes unmapped transcripts lacking official Gene Symbols to ensure relevance
  cleaned <- annotated[!is.na(annotated$Symbol) & annotated$Symbol != "", ]
  
  # O.3.5. Duplicate Probe Resolution
  # Retains only the most significant probe per unique Symbol to prevent bias
  cleaned <- cleaned[order(cleaned$adj.P.Val), ] 
  cleaned <- cleaned[!duplicated(cleaned$Symbol), ] 
  
  return(cleaned)
}

# O.4. Prepare Annotated Datasets for Volcano and Heatmap Plots
# Merges statistical results with Gene Symbols for visualization
volcano_1 <- clean_and_annotate (fit2, 1, topTable_Master_Clean)
volcano_2 <- clean_and_annotate (fit2, 2, topTable_Master_Clean)
volcano_3 <- clean_and_annotate (fit2, 3, topTable_Master_Clean)
volcano_4 <- clean_and_annotate (fit2, 4, topTable_Master_Clean)
volcano_5 <- clean_and_annotate (fit2, 5, topTable_Master_Clean)
volcano_6 <- clean_and_annotate (fit2, 6, topTable_Master_Clean)

# O.5. Define Contrast-Specific Column Selection Function
# Isolates essential statistics and removes irrelevant contrast columns via pattern matching
filter_volcano <- function(volcano, pattern){
  
  volcano %>%
    dplyr::select(
      MatchID,
      logFC,
      AveExpr,
      t,
      P.Value,
      adj.P.Val,
      B,
      Symbol,
      Gene_Description,
      SPOT_ID,
      dplyr::contains(pattern)
    )
  
}

# O.6. Execute Final Data Sanitization for Individual Volcano Plots
# Generates discrete, high-fidelity dataframes optimized for targeted visualization
volcano_1_Final <- filter_volcano(volcano_1, "PNI_Effect_Tumor")
volcano_2_Final <- filter_volcano(volcano_2, "PNI_Effect_Normal")
volcano_3_Final <- filter_volcano(volcano_3, "Tumor_vs_Normal_Pos")
volcano_4_Final <- filter_volcano(volcano_4, "Tumor_vs_Normal_Neg")
volcano_5_Final <- filter_volcano(volcano_5, "Extreme_Contrast")
volcano_6_Final <- filter_volcano(volcano_6, "Reverse_Contrast")

# O.7. Systematic Export of Annotated Volcano and Heatmap Plot Data
# Persists cleaned and filtered statistical results for all 6 clinical contrasts to CSV
write.csv(volcano_1_Final, "results/Volcano_and_Heatmap_Data_1_PNI_Effect_Tumor.csv", row.names = FALSE)
write.csv(volcano_2_Final, "results/Volcano_and_Heatmap_Data_2_PNI_Effect_Normal.csv", row.names = FALSE)
write.csv(volcano_3_Final, "results/Volcano_and_Heatmap_Data_3_Tumor_vs_Normal_PNI_Pos.csv", row.names = FALSE)
write.csv(volcano_4_Final, "results/Volcano_and_Heatmap_Data_4_Tumor_vs_Normal_PNI_Neg.csv", row.names = FALSE)
write.csv(volcano_5_Final, "results/Volcano_and_Heatmap_Data_5_Extreme_Contrast.csv", row.names = FALSE)
write.csv(volcano_6_Final, "results/Volcano_and_Heatmap_Data_6_Reverse_Contrast.csv", row.names = FALSE)

# O.8. Automated Extraction of Top 50 Upregulated and Downregulated Genes
# Consolidates all six final volcano dataframes into a unified list for iterative batch processing.
volcano_list <- list(volcano_1_Final, volcano_2_Final, volcano_3_Final, 
                     volcano_4_Final, volcano_5_Final, volcano_6_Final)

# O.9. Iterative Feature Selection and Global Environment Assignment
# Executes a robust loop to identify high-confidence features based on logFC and significance thresholds.
for (i in 1:length(volcano_list)) {
  
  # Extracts the current dataframe and dynamically identifies target columns to prevent mapping errors.
  df <- volcano_list[[i]]
  col_logFC <- names(df)[grep("logFC", names(df))][1]
  col_adjP  <- names(df)[grep("adjP", names(df))][1]
  
  # Filters for Upregulated features (logFC > 1) and ranks the top 50 genes by magnitude of change.
  up_logic <- df[[col_logFC]] > 1 & df[[col_adjP]] < 0.05
  up_df <- df[up_logic, ]
  top_up <- up_df[order(-up_df[[col_logFC]]), ] %>% head(50)
  
  # Filters for Downregulated features (logFC < -1) and isolates the top 50 most suppressed genes.
  down_logic <- df[[col_logFC]] < -1 & df[[col_adjP]] < 0.05
  down_df <- df[down_logic, ]
  top_down <- down_df[order(down_df[[col_logFC]]), ] %>% head(50)
  
  # Dynamically assigns the filtered top 50 results to independent objects within the Global Environment.
  assign(paste0("top50_up_", i), top_up, envir = .GlobalEnv)
  assign(paste0("top50_down_", i), top_down, envir = .GlobalEnv)
  
  # Generates localized CSV backups for each contrast to maintain a modular data storage structure.
  write.csv(top_up, paste0("results/Top50_Up_Contrast_", i, ".csv"), row.names = FALSE)
  write.csv(top_down, paste0("results/Top50_Down_Contrast_", i, ".csv"), row.names = FALSE)
  
  # Provides real-time console feedback to verify the successful completion of each contrast iteration.
  print(paste("Successfully Processed Contrast", i))
}

# O.10. Formal Export of Top-Ranked Upregulated Gene Profiles
# Persists the top 50 significantly upregulated genes for all 6 clinical comparisons into the results directory.
write.csv(top50_up_1, "results/Top50_Up_volcano_and_heatmap_1_PNI_Effect_Tumor.csv", row.names = FALSE)
write.csv(top50_up_2, "results/Top50_Up_volcano_and_heatmap_2_PNI_Effect_Normal.csv", row.names = FALSE)
write.csv(top50_up_3, "results/Top50_Up_volcano_and_heatmap_3_Tumor_vs_Normal_PNI_Pos.csv", row.names = FALSE)
write.csv(top50_up_4, "results/Top50_Up_volcano_and_heatmap_4_Tumor_vs_Normal_PNI_Neg.csv", row.names = FALSE)
write.csv(top50_up_5, "results/Top50_Up_volcano_and_heatmap_5_Extreme_Contrast.csv", row.names = FALSE)
write.csv(top50_up_6, "results/Top50_Up_volcano_and_heatmap_6_Reverse_Contrast.csv", row.names = FALSE)

# O.11. Formal Export of Top-Ranked Downregulated Gene Profiles
# Records the top 50 significantly suppressed genes for all 6 clinical comparisons for downstream functional analysis.
write.csv(top50_down_1, "results/Top50_Down_volcano_and_heatmap_1_PNI_Effect_Tumor.csv", row.names = FALSE)
write.csv(top50_down_2, "results/Top50_Down_volcano_and_heatmap_2_PNI_Effect_Normal.csv", row.names = FALSE)
write.csv(top50_down_3, "results/Top50_Down_volcano_and_heatmap_3_Tumor_vs_Normal_PNI_Pos.csv", row.names = FALSE)
write.csv(top50_down_4, "results/Top50_Down_volcano_and_heatmap_4_Tumor_vs_Normal_PNI_Neg.csv", row.names = FALSE)
write.csv(top50_down_5, "results/Top50_Down_volcano_and_heatmap_5_Extreme_Contrast.csv", row.names = FALSE)
write.csv(top50_down_6, "results/Top50_Down_volcano_and_heatmap_6_Reverse_Contrast.csv", row.names = FALSE)

# ==================================
# PART P. VOLCANO PLOT VISUALIZATION
# ==================================
# # P.1. Define Automated Volcano Plot Function
# Categorizes genes based on significance thresholds and generates a ggplot object
make_volcano <- function(tt_data, title_name) {
  
  # P.1.1. Classification of Differential Expression status
  # Defines UP and DOWN regulation using |log2FC| > 1 and Adj.P.Val < 0.05
  tt_data$diffexpressed <- "NO"
  tt_data$diffexpressed[tt_data$logFC > 1 & tt_data$adj.P.Val < 0.05] <- "Up-regulated"
  tt_data$diffexpressed[tt_data$logFC < -1 & tt_data$adj.P.Val < 0.05] <- "Down-regulated"
  
  # P.1.2. Generate Plot
  # Renders the volcano plot with consistent aesthetics
  ggplot(tt_data, aes(x = logFC, y = -log10(adj.P.Val), color = diffexpressed)) +
    geom_point(alpha = 0.4, size = 1.5) +
    scale_color_manual(values = c("Down-regulated" = "blue", "NO" = "grey", "Up-regulated" = "red")) +
    theme_bw(base_size = 14) + 
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      legend.position = "right",
      panel.grid.minor = element_blank(),
      axis.text = element_text(color = "black"),
      axis.ticks = element_line(color = "black", linewidth = 0.8)
    ) +
    labs(
      title = title_name,
      x = "log2 Fold Change",
      y = "-log10 (Adjusted P-value)",
      color = "Expression"
    ) +
    geom_vline(xintercept = c(-1, 1), col = "black", linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), col = "black", linetype = "dashed")
}

# P.2. Generate Volcano Plots for Each Comparison
# Iterates through all 6 annotated datasets to create high-resolution visuals
v1 <- make_volcano (volcano_1_Final, "PNI Effect in Tumor")
v2 <- make_volcano (volcano_2_Final, "PNI Effect in Normal")
v3 <- make_volcano (volcano_3_Final, "Tumor vs Normal (PNI_Pos)") 
v4 <- make_volcano (volcano_4_Final, "Tumor vs Normal (PNI_Neg)")
v5 <- make_volcano (volcano_5_Final, "Extreme Contrast (PNI_Pos_Tumor vs PNI_Neg_Normal)")
v6 <- make_volcano (volcano_6_Final, "Reverse Contrast (PNI_Neg_Tumor vs PNI_Pos_Normal)")

# P.3. Review Visualization Results
# Prints all generated plots to the RStudio plot pane for verification
print(v1)
print(v2)
print(v3)
print(v4)
print(v5)
print(v6)

# P.4. Automated Publication-Quality Export
# Saves all 6 plots as high-DPI PNG files optimized for repository structure
ggsave("plots/Volcano_1_PDAC_pni_tumor.png", plot = v1, width = 15, height = 8, dpi = 300)
ggsave("plots/Volcano_2_PDAC_pni_normal.png", plot = v2, width = 15, height = 8, dpi = 300)
ggsave("plots/Volcano_3_PDAC_tumor_vs_normal_pos.png", plot = v3, width = 15, height = 8, dpi = 300)
ggsave("plots/Volcano_4_PDAC_tumor_vs_normal_neg.png", plot = v4, width = 15, height = 8, dpi = 300)
ggsave("plots/Volcano_5_PDAC_extreme_contrast.png", plot = v5, width = 15, height = 8, dpi = 300)
ggsave("plots/Volcano_6_PDAC_reverse_contrast.png", plot = v6, width = 15, height = 8, dpi = 300)

# ===============================================
# PART Q. GLOBAL HEATMAP PROFILING OF 100 SAMPLES
# ===============================================
# ----------------------------------
# Q.1. GLOBAL ANOVA DATA PREPARATION
# ----------------------------------
# Q.1.1. Construct ANOVA Design Matrix
# Defines the linear model for all 4 clinical groups to identify global variance
design_anova <- model.matrix(~0 + gset$group)
colnames(design_anova) <- levels(gset$group)

# Q.1.2. Fit Linear Model and Empirical Bayes
# Performs global F-test to detect genes with significant differences across any group
fit_anova <- lmFit(ex, design_anova)
fit_anova <- eBayes(fit_anova)

# Q.1.3. Extract Top 50 ANOVA Genes
# Retrieves genes with the lowest adjusted P-values across the entire 4-group cohort
anova_tt <- topTable(fit_anova, number = Inf, adjust.method = "fdr")
anova_tt$MatchID <- rownames(anova_tt)

# Q.1.4. Merge with Annotation and Clean
# Ensures only probes with valid Gene Symbols are used for the Global Heatmap
anova_clean <- merge(anova_tt, topTable_Master_Clean[, c("ProbeID", "Symbol", "Gene_Description")], 
                     by.x = "MatchID", by.y = "ProbeID", all.x = TRUE)
anova_clean <- anova_clean[!is.na(anova_clean$Symbol) & anova_clean$Symbol != "", ]
anova_clean <- anova_clean[order(anova_clean$adj.P.Val), ]
global_top50_anova <- head(anova_clean[!duplicated(anova_clean$Symbol), ], 50)

# Q.1.5. Export Global Heatmap Source Data
# Saves the top 50 ANOVA-significant genes used for the 100-sample global profiling
write.csv(global_top50_anova, "results/Global_Heatmap_Top50_ANOVA_Genes.csv", row.names = FALSE)

# ------------------------------------------------------------
# PART Q.2. GLOBAL HEATMAP VISUALIZATION: MULTI-GROUP OVERVIEW
# ------------------------------------------------------------
# Q.2.1. Prepare Global Expression Matrix
# Extracts data for all 100 samples using the top 50 ANOVA-significant genes
mat_global <- ex[global_top50_anova$MatchID, ]
rownames(mat_global) <- make.unique(as.character(global_top50_anova$Symbol))

# Q.2.2. Export Global Expression Matrix (Z-scores)
# Persists the actual expression values of the top 50 genes across all 100 samples
write.csv(mat_global, "results/Global_Heatmap_Expression_Matrix.csv", row.names = TRUE)

# Q.2.3. Define Comprehensive Sample Annotation
# Maps all 100 samples to their respective clinical groups for visualization
ann_col_global <- data.frame(`Clinicopathological Group` = gset$group, check.names = FALSE)
rownames(ann_col_global) <- colnames(mat_global)

# Q.2.4. Global Project Color Palette
# Enforces color consistency across the entire 4-group clinical spectrum
global_colors <- list(`Clinicopathological Group` = c(
  PNI_Pos_Tumor  = "#8B0000", PNI_Neg_Tumor  = "#FF8282",
  PNI_Pos_Normal = "#006400", PNI_Neg_Normal = "#90EE90"
))

# Q.2.5. Render and Export Global Overview Heatmap
# Visualizes global clustering patterns to demonstrate inter-group variability
h_global <- pheatmap(
  mat_global, 
  scale = "row", 
  annotation_col = ann_col_global, 
  annotation_colors = global_colors,
  color = colorRampPalette(c("blue","white","red"))(100),
  border_color = "black", 
  show_colnames = FALSE, 
  clustering_method = "ward.D2", 
  fontsize_row = 7.5,
  fontsize = 8,
  main = "Top 50 Differentially Expressed Genes of 100 Samples",
  silent = TRUE
)

# Q.2.6. Review Visualizations in Console
# Run this command to inspect heatmap global plot with a clean canvas
grid.newpage(); grid.draw(h_global$gtable)

# Q.2.7. Automated Image Export
# Export to PNG with Publication-Ready Labels
png("plots/H0_Global_ANOVA_Heatmap_PDAC_100.png", width=16, height=10, units="in", res=300)
grid::grid.newpage()
grid::grid.draw(h_global$gtable)
dev.off()

# =================================================================
# PART R. HEATMAP VISUALIZATION OF TOP 50 DEGs ACROSS 6 COMPARISONS
# =================================================================
# R.1. Heatmap Visualization Function
# Renders a high-resolution heatmap with 50 genes and synchronized clinical colors
plot_my_heatmap <- function(volcano, title, filename, groups_to_keep) {
  
  # R.1.1. Select Top 50 Genes
  # Taking the most significant genes that already passed cleaning in Part P
  top50 <- head(volcano[order(volcano$adj.P.Val), ], 50)
  
  # R.1.2. Identify Relevant Samples
  # Filters the expression set to include only specified clinical groups
  relevant_samples <- which(gset$group %in% groups_to_keep)
  
  # R.1.3. Expression Matrix Extraction
  # Retrieving scaled expression values and mapping unique Symbols to rows
  mat <- ex[top50$MatchID, relevant_samples, drop = FALSE]
  rownames(mat) <- make.unique(as.character(top50$Symbol))
  
  # R.1.4. Define Clinical Annotation Data
  # Subsets metadata and drops unused factor levels for a clean legend
  # 1. Subset the groups
  subset_group <- gset$group[relevant_samples]
  subset_group <- droplevels(as.factor(subset_group))
  ann_col <- data.frame(
    `Clinicopathological Group` = subset_group, 
    check.names = FALSE
  )
  rownames(ann_col) <- colnames(mat)
  
  # R.1.5. Apply Project-Wide Color Mapping
  # Filters the master palette to match active groups in the current subset
  master_colors <- c(
    PNI_Pos_Tumor  = "#8B0000",
    PNI_Neg_Tumor  = "#FF8282",
    PNI_Pos_Normal = "#006400",
    PNI_Neg_Normal = "#90EE90"
  )
  
  # R.1.6. Apply Dynamic Legend Filtering
  # Synchronizes annotation colors by filtering the master palette for active groups
  active_groups <- levels(subset_group)
  ann_colors <- list(
    `Clinicopathological Group` = master_colors[names(master_colors) %in% active_groups]
  )
  
  # R.1.7. Render Publication-Quality Heatmap
  # Uses Ward.D2 clustering and row scaling for optimal pattern recognition
  p <- pheatmap(
    mat, 
    scale = "row",                      # Row-wise Z-score normalization
    annotation_col = ann_col, 
    annotation_colors = ann_colors,
    color = colorRampPalette(c("blue","white","red"))(100),
    border_color = "black", 
    show_colnames = FALSE, 
    clustering_method = "ward.D2", 
    fontsize_row = 7.5,
    main = title, 
    silent = TRUE,
  )
  
  # R.1.8. Automated Image Export
  # Saves 300 DPI PNG files to the 'plots' folder for documentation
  png(paste0("plots/", filename), width=16, height=10, units="in", res=300)
  grid::grid.newpage()
  grid::grid.draw(p$gtable)
  dev.off()

  return(p)
}

# R.2. Heatmap Generation for Six Comparisons
# Generates comparative visualizations for PNI effect and Tumor vs Normal states
# H1: PNI Effect in Tumor (Focuses on positive vs negative PNI in malignant tissue)
h1 <- plot_my_heatmap(volcano_1_Final, 
                      "Top 50 Differentially Expressed Genes: PNI Effect in Tumor", 
                      "H1_PDAC_pni_tumor.png", 
                      groups_to_keep = c("PNI_Pos_Tumor", "PNI_Neg_Tumor"))

# H2: PNI Effect in Normal (Controls for PNI signatures in non-malignant tissue)
h2 <- plot_my_heatmap(volcano_2_Final, 
                      "Top 50 Differentially Expressed Genes: PNI Effect in Normal", 
                      "H2_PDAC_pni_normal.png", 
                      groups_to_keep = c("PNI_Pos_Normal", "PNI_Neg_Normal"))

# H3: Tumor vs Normal (PNI Positive) (Identifies malignancy markers in PNI+ samples)
h3 <- plot_my_heatmap(volcano_3_Final, 
                      "Top 50 Differentially Expressed Genes: Tumor vs Normal (PNI_Pos)", 
                      "H3_PDAC_tumor_vs_normal_pos.png", 
                      groups_to_keep = c("PNI_Pos_Tumor", "PNI_Pos_Normal"))

# H4: Tumor vs Normal (PNI Negative) (Identifies malignancy markers in PNI- samples)
h4 <- plot_my_heatmap(volcano_4_Final, 
                      "Top 50 Differentially Expressed Genes: Tumor vs Normal (PNI_Neg)", 
                      "H4_PDAC_tumor_vs_normal_neg.png", 
                      groups_to_keep = c("PNI_Neg_Tumor", "PNI_Neg_Normal"))

# H5: Extreme Contrast (Compares PNI+ Tumor against PNI- Normal baseline)
h5 <- plot_my_heatmap(volcano_5_Final, 
                      "Top 50 Differentially Expressed Genes: Extreme Contrast (PNI_Pos_Tumor vs PNI_Neg_Normal)", 
                      "H5_PDAC_extreme.png", 
                      groups_to_keep = c("PNI_Pos_Tumor", "PNI_Neg_Normal"))

# H6: Reverse Contrast (Compares PNI- Tumor against PNI+ Normal baseline)
h6 <- plot_my_heatmap(volcano_6_Final, 
                      "Top 50 Differentially Expressed Genes: Reverse Contrast (PNI_Neg_Tumor vs PNI_Pos_Normal)", 
                      "H6_PDAC_reverse.png", 
                      groups_to_keep = c("PNI_Neg_Tumor", "PNI_Pos_Normal"))

# R.3. Review Visualizations in Console
# Run these commands one-by-one to inspect each plot with a clean canvas
grid.newpage(); grid.draw(h1$gtable)
grid.newpage(); grid.draw(h2$gtable)
grid.newpage(); grid.draw(h3$gtable)
grid.newpage(); grid.draw(h4$gtable)
grid.newpage(); grid.draw(h5$gtable)
grid.newpage(); grid.draw(h6$gtable)

# ==================================================================
# PART S. GLOBAL SCATTER PLOT PROFILING OF CLINICOPATHOLOGICAL GROUP
# ==================================================================
# -----------------------------------------
# S.1. GLOBAL SCATTER PLOT DATA PREPARATION
# -----------------------------------------
# S.1.1. Validate sample and group alignment
# Ensures that the expression matrix columns match the metadata row order
colnames(ex)
gset$group

# S.1.2. Convert expression matrix to data frame
# Casts the raw matrix into a flexible data frame structure for processing
expr_df <- as.data.frame(ex)

# S.1.3. Add ProbeID column
# Persists row names as an explicit variable for relational merging
expr_df$ProbeID <- rownames(expr_df)

# S.1.4. Transpose data orientation
# Rotates the matrix so that biological samples are treated as observations
expr_t <- as.data.frame(t(ex))

# S.1.5. Append Group metadata
# Integrates clinical status labels into the transposed expression dataset
expr_t$Group <- gset$group

# S.1.6. Calculate Mean Expression per Group
# Computes the arithmetic average of expression values partitioned by cohort
mean_expression <- expr_t %>%
  group_by(Group) %>%
  summarise(across(where(is.numeric), mean))

mean_expression

# S.1.7. Reshape to Long Format
# Transforms the wide summary into a long format to facilitate data joining
mean_long <- mean_expression %>%
  pivot_longer(
    cols = -Group,
    names_to = "ProbeID",
    values_to = "MeanExpression"
  )

# S.1.8. Preview long-format data
# Checks the integrity of the tidied dataset before final annotation
head(mean_long)

# S.1.9. Calculate Standard Deviation per Group
# Generates variance metrics to evaluate data spread across clinical groups
mean_sd_expression <- expr_t %>%
  group_by(Group) %>%
  summarise(across(where(is.numeric),
                   list(mean = mean, sd = sd)))

# S.1.10. Verify column consistency
# Ensures primary keys (ProbeID) are identical in both tables before merging
colnames(mean_long)
colnames(topTable_Master_Clean)

# S.1.11. Merge expression data with Gene Annotations
# Cross-references the calculated means with Symbol and Description data
mean_annotated <- merge(
  mean_long,
  topTable_Master_Clean[, c("ProbeID", "Symbol", "Gene_Description")],
  by = "ProbeID",
  all.x = TRUE
)

# S.1.12. Filter out unannotated records
# Removes entries that lack official Gene Symbols to ensure biological relevance
mean_annotated <- mean_annotated %>%
  filter(!is.na(Symbol) & Symbol != "")

# S.1.13. Check for potential duplicates
# Audits the dataset for redundant Symbol or Group collisions
mean_annotated %>%
  dplyr::count(Symbol, Gene_Description, Group) %>%
  dplyr::filter(n > 1)

# S.1.14. Handle probe-to-gene redundancy
# Resolves cases where multiple probes map to a single Symbol by averaging
mean_long_clean <- mean_annotated %>%
  dplyr::group_by(Symbol, Gene_Description, Group) %>%
  dplyr::summarise(
    MeanExpression = mean(MeanExpression),
    .groups = "drop"
  )

# S.1.15. Pivot to Wide Format
# Creates a final comparative table with clinical groups as distinct columns
mean_wide <- mean_long_clean %>%
  tidyr::pivot_wider(
    names_from = Group,
    values_from = MeanExpression
  )

# S.1.16. Pre-visualization Data Filtering
# Removing non-finite values (NA, Inf, -Inf) to ensure statistical validity
mean_long_clean <- mean_long[is.finite(mean_long$MeanExpression), ]
mean_long_clean <- mean_long %>%
  filter(is.finite(MeanExpression))

# S.1.17. Data Quality Audit Log
# Reporting the total count of filtered records to track data loss
cat("Total data points:", nrow(mean_long), "\n")
cat("Rows removed (NA/Inf):", nrow(mean_long) - nrow(mean_long_clean), "\n")

# S.1.18. Prepare Data for Scatter Plot
# Maps mean and standard deviation variables using regex-based column partitioning
mean_sd_long <- mean_sd_expression %>%
  pivot_longer(
    cols = -Group,
    names_to = c("ProbeID", ".value"),
    names_pattern = "(.*)_(mean|sd)"
  )

# S.1.19. Integrate Statistical Metrics with Gene Annotations
# Merges calculated Mean and SD values with Symbol and Description data for tabular clarity
mean_sd_summary <- merge(
  mean_sd_long,
  topTable_Master_Clean[, c("ProbeID", "Symbol", "Gene_Description")],
  by = "ProbeID",
  all.x = TRUE
)

# S.1.20. Post-Merge Data Sanitization
# Filters out unmapped probes lacking official Gene Symbols to ensure biological relevance
mean_sd_summary <- mean_sd_summary %>%
  filter(!is.na(Symbol) & Symbol != "")

# S.1.21. Automated Export of Expression Stability Data
# Persists the finalized Mean and SD summary to a CSV file for clinical stability audits
write.csv(mean_sd_summary, "results/Mean_SD_Expression_Summary.csv", row.names = FALSE)

# --------------------------------------
# S.2. GLOBAL SCATTER PLOT VISUALIZATION
# --------------------------------------
# S.2.1. Generate Mean vs SD Scatter Plot
# Visualizes expression stability across all clinical cohorts with GAM smoothing
scatter_mean_vs_sd <- ggplot(mean_sd_long, aes(x = mean, y = sd, color = Group)) +
  geom_point(alpha = 0.4, size = 1.5) + # Slightly increased size and opacity for visibility
  geom_smooth(method = "gam", color = "black", linetype = "dashed", linewidth = 0.8) +
  scale_color_manual(values = my_pni_colors) +
  facet_wrap(~Group) +
  labs(
    title = "Gene Expression Stability Profile of 100 Samples: Mean vs Standard Deviation",
    x = "Log2 Mean Expression",
    y = "Standard Deviation (SD)",
    color = "Clinicopathological Group" # Fixed legend title
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12, margin = margin(b = 15)),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10),
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(face = "bold", size = 11),
    axis.line = element_line(color = "black"),
    panel.grid.minor = element_blank()
  ) +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1)))

# S.2.2. Console Visual Verification
# Renders the finalized ggplot object for immediate structural inspection
print(scatter_mean_vs_sd)

# S.2.3. High-Resolution Scientific Export
# Executes automated file saving via ggsave for optimal vector scaling
ggsave(
  filename = "plots/Mean_vs_SD_Scatter_PDAC_100.png", 
  plot = scatter_mean_vs_sd, 
  width = 16, 
  height = 10, 
  units = "in", 
  dpi = 300
)

# ================================================================================
# PART T. IDENTIFICATION OF CORE BIOMARKERS (Intersection of 6 Clinical Contrasts)
# ================================================================================
# ---------------------------------------
# T.1. VENN DIAGRAM PLOT DATA PREPARATION
# ---------------------------------------
# T.1.1. Define Significance Thresholds
# Standardizing on FDR < 0.05 for high-stringency biomarker selection
p_cutoff <- 0.05

# T.1.2. Filter Common DEGs (The Intersection Logic)
# Aggregates Differentially Expressed Genes (DEGs) across all six clinical contrasts
common_degs_raw <- topTable_Master_Clean %>%
  filter(
      adjP_PNI_Effect_Tumor    < p_cutoff |
      adjP_PNI_Effect_Normal   < p_cutoff |
      adjP_Tumor_vs_Normal_Pos < p_cutoff | 
      adjP_Tumor_vs_Normal_Neg < p_cutoff |
      adjP_Extreme_Contrast    < p_cutoff |
      adjP_Reverse_Contrast    < p_cutoff
  )

# T.1.3. Data Synchronization: Mapping Probes to Unique Genes
# Resolves probe redundancy by collapsing multiple entries into single unique Gene Symbols
common_degs <- common_degs_raw %>%
  arrange(adjP_Tumor_vs_Normal_Pos) %>% 
  distinct(Symbol, .keep_all = TRUE)   

# T.1.4. Data Ranking & Mean LogFC Calculation
# Calculate mean logFC across all contrasts and rank by absolute impact
common_degs$mean_logFC <- rowMeans(common_degs[, grep("logFC", names(common_degs))])
common_degs <- common_degs %>% arrange(desc(abs(mean_logFC)))
print(head(common_degs[, c("ProbeID", "Symbol", "mean_logFC")]))

# T.1.5. Identify Directional Subsets
# Segregate biomarkers based on directional consistency for downstream analysis
common_up <- common_degs %>% filter(if_all(contains("logFC"), ~ . > 0))
common_down <- common_degs %>% filter(if_all(contains("logFC"), ~ . < 0))

# T.1.6. Identify Mixed/Inconsistent DEGs
# Capture genes that show mixed polarity across different clinical contrasts
common_mixed <- common_degs %>% 
  filter(!(Symbol %in% common_up$Symbol) & !(Symbol %in% common_down$Symbol))

# T.1.7. Summary Report to Console
# Generates a quantitative audit of core biomarkers for structural verification
cat("\n====================================\n")
cat("   DIRECTIONAL CONSISTENCY REPORT\n")
cat("====================================\n")
cat(paste("Total Unique DEGs (Ranked) :", nrow(common_degs), "\n"))
cat(paste("Consistently UP            :", nrow(common_up), "\n"))
cat(paste("Consistently DOWN          :", nrow(common_down), "\n"))
cat(paste("Mixed/Inconsistent         :", nrow(common_mixed), "\n"))
cat("------------------------------------\n")

# T.1.8. Automated Export of Finalized Datasets
# Persist all biomarker categories to CSV for validation and reporting
write.csv(common_degs, "results/Core_Biomarkers_All_6_Contrasts.csv", row.names = FALSE)
write.csv(common_up, "results/Core_Biomarkers_Consistently_UP.csv", row.names = FALSE)
write.csv(common_down, "results/Core_Biomarkers_Consistently_DOWN.csv", row.names = FALSE)
write.csv(common_mixed, "results/Core_Biomarkers_Mixed_Direction.csv", row.names = FALSE)

# T.1.9. Visualizing the Intersection: 6-set Venn Diagram
# Prepare lists for Venn visualization (Extracting Symbols for each set)
venn_list <- list(
  PNI_Effect_Tumor    = topTable_Master_Clean$Symbol[topTable_Master_Clean$adjP_PNI_Effect_Tumor < p_cutoff],
  PNI_Effect_Normal   = topTable_Master_Clean$Symbol[topTable_Master_Clean$adjP_PNI_Effect_Normal < p_cutoff],
  Tumor_vs_Normal_Pos = topTable_Master_Clean$Symbol[topTable_Master_Clean$adjP_Tumor_vs_Normal_Pos < p_cutoff],
  Tumor_vs_Normal_Neg = topTable_Master_Clean$Symbol[topTable_Master_Clean$adjP_Tumor_vs_Normal_Neg < p_cutoff],
  Extreme_Contrast    = topTable_Master_Clean$Symbol[topTable_Master_Clean$adjP_Extreme_Contrast < p_cutoff],
  Reverse_Contrast    = topTable_Master_Clean$Symbol[topTable_Master_Clean$adjP_Reverse_Contrast < p_cutoff]
)

# ------------------------------------
# T.2. VENN DIAGRAM PLOT VISUALIZATION
# ------------------------------------
# T.2.1. Generate Plot
# Configures a 6-way elliptical Venn diagram to visualize overlapping gene signatures
p_venn <- ggVennDiagram(
  venn_list, 
  shape_id = "601", 
  set_color = c("#FF0000", "#FF7F00", "#FFD700", "#00FF00", "#0000FF", "#8B00FF"),
  label = "count", 
  label_size = 3.5,
  label_alpha = 0,
  set_size = 2.5  
) +
  scale_fill_gradientn(
    colors = c("white", "#FEE0D2", "#FB6A4A"),
    name = "Number of Genes",
    guide = guide_colorbar(
      title.position = "top", 
      title.hjust = 0.5, 
      barwidth = 18, 
      barheight = 1
    )
  ) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 11),
    plot.margin = margin(0.5, 5, 0.5, 5, "cm"), 
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 12, color = "black")
  ) +
  labs(
    title = "Differential Expression Intersection Profile",
    subtitle = "Total Unique DEGs Across All 6 Comparisons: 9,750 Genes",
    caption = "FDR < 0.05 | LogFC based on Clinical Contrasts"
  )

# T.2.2. Console Visual Verification
# Renders the finalized venn diagram for immediate structural inspection
print(p_venn)

# T.2.3. Save the Diagram
# Exports the final visualization at high resolution with a forced white background
ggsave("plots/Venn_6_Comparisons_PDAC.png", 
       p_venn, 
       width = 18, 
       height = 12, 
       dpi = 300,
       bg = "white",           
       limitsize = FALSE)

# ================================================================
# PART U. SYSTEMIC BIOLOGICAL INTERPRETATION: GO & KEGG ENRICHMENT
# ================================================================
# ------------------------------------------------------------------------
# U.1. DATA CURATION & PRE-ENRICHMENT PROCESSING FOR FUNCTIONAL ANNOTATION
# ------------------------------------------------------------------------
# U.1.1. Identify All Columns starting with 'adjP'
# Automatically detects all six clinical contrast columns for multi-set significance filtering
adjP_cols <- grep("^adjP", colnames(topTable_Master_Clean), value = TRUE)

# U.1.2. Select Genes Significant in AT LEAST ONE Contrast
# Filters for genes meeting the FDR < 0.05 threshold in at least one experimental condition
sig_genes_master <- topTable_Master_Clean[apply(topTable_Master_Clean[, adjP_cols], 1, 
                                                function(x) any(x < 0.05)), ]

# U.1.3. Convert Gene Symbols to Entrez ID
# Utilizes bitr to map Gene Symbols to Entrez IDs for database-standardized functional analysis
gene_conv <- bitr(sig_genes_master$Symbol, 
                  fromType = "SYMBOL", 
                  toType   = "ENTREZID", 
                  OrgDb    = "org.Hs.eg.db")

# U.1.4. Resolve Ambiguous Probe Mapping (Critical Data Integrity Step)
# Eliminates database redundancy by selecting a single unique Entrez ID for each mapped symbol
gene_conv_final <- gene_conv %>%
  distinct(SYMBOL, .keep_all = TRUE)

# U.1.5. Display Final Gene Count
# Outputs the verified count of unique genes ready for downstream enrichment processing
cat("Final Unique Genes for Enrichment:", nrow(gene_conv_final), "\n")

# U.1.6. Gene Ontology (GO) Enrichment Analysis
# Focusing on biological processes to understand what these genes are doing
ego <- enrichGO(gene          = gene_conv_final$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                readable      = TRUE)

# U.1.7. KEGG Pathway Enrichment Analysis
# Maps the finalized gene list to KEGG biochemical pathways for metabolic and signaling insights
ekegg <- enrichKEGG(gene         = gene_conv_final$ENTREZID,
                    organism     = 'hsa', # 'hsa' is for Homo sapiens
                    pvalueCutoff = 0.05)

# U.1.8. Systematic Export of Enrichment Results
# Persists full GO and KEGG tabular data to CSV for comprehensive supplementary documentation
write.csv(as.data.frame(ego), "results/GO_Enrichment_Results.csv", row.names = FALSE)
write.csv(as.data.frame(ekegg), "results/KEGG_Pathway_Results.csv", row.names = FALSE)

# -------------------------------------------------------------------------
# U.2. VISUALIZING FUNCTIONAL LANDSCAPES: DOT PLOTS (GO) & BAR PLOTS (KEGG)
# -------------------------------------------------------------------------
# U.2.1. Generate GO Dotplot
# Visualizes the top 15 Biological Processes ranked by gene ratio and significance
p_go <- dotplot(ego, showCategory = 15) + 
  ggtitle("Top 15 GO Enrichment: Biological Processes") +
  theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5)) + 
  guides(size = guide_legend(override.aes = list(shape = 21))) 

# U.2.2. Generate KEGG Barplot
# Displays the top 15 metabolic pathways based on absolute gene counts
p_kegg <- barplot(ekegg, showCategory = 15, x = "Count") + 
  ggtitle("Top 15 KEGG Pathway Enrichment Analysis") +
  theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5)) 

# U.2.3. Console Verification & Display
# Renders both finalized plot objects in the RStudio Plots pane for immediate review
print(p_go)
print(p_kegg)

# U.2.4. High-Resolution Publication Export
# Utilizes ggsave to export high-fidelity 300 DPI images for formal reporting
ggsave("plots/GO_Enrichment_Dotplot_PDAC.png", plot = p_go, width = 10, height = 8, dpi = 300)
ggsave("plots/KEGG_Enrichment_Barplot_PDAC.png", plot = p_kegg, width = 10, height = 8, dpi = 300)

# =======================================
# PART V. FINALIZATION & AUTHORSHIP AUDIT
# =======================================
# V.1. Success Summary with Ownership Claim
# Explicitly states the author and project completion details
cat("\n==========================================================================================\n")
cat(" PROJECT TITLE: INTEGRATIVE TRANSCRIPTOMIC PROFILING OF PERINEURAL INVASION (PNI) SIGNATURES\n")
cat(" IN PANCREATIC DUCTAL ADENOCARCINOMA (PDAC): A MULTI-CONTRAST BIOINFORMATICS STUDY\n")
cat("============================================================================================\n")
cat(" STATUS       : COMPLETED SUCCESSFULLY\n")
cat(" DEVELOPED BY : Yosia Jose Rasdiva Manurung\n")
cat(" AFFILIATION  : Diponegoro University (UNDIP), Indonesia\n")
cat(" COMPLETED ON :", format(Sys.time(), "%A, %d %B %Y - %H:%M:%S"), "\n")
cat("============================================================================================\n")

# ------------------------------------------------------------------------------
# END OF SCRIPT - Unauthorized distribution or usage is strictly prohibited.
# Created by [Y.J.R. Manurung] for Advanced Transcriptomics Research.
# ------------------------------------------------------------------------------
