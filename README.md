## 📜 Copyright & Licensing
> [!CAUTION]
> **Academic Integrity & Copyright Notice**.
> 
> This project, including all R scripts, datasets, and visualizations, is the academic work of **Yosia Jose Rasdiva Manurung**. 
> Unauthorized use, reproduction, or redistribution of this content for commercial purposes or academic submission by others is **strictly prohibited**. For permissions or collaboration, please contact the author.  

# Integrative Transcriptomic Profiling of Perineural Invasion (PNI) Signatures in Pancreatic Ductal Adenocarcinoma (PDAC): a Multi-Contrast Bioinformatics Study
**Author:** Yosia Jose Rasdiva Manurung  
**Affiliation:** Diponegoro University  

---

## 1. Research Objectives & Comparison Groups
The primary goal of this study is to identify Differentially Expressed Genes (DEGs) across **six PNI-related clinical contrasts**. This systematic approach allows me to isolate the specific transcriptomic signatures driven by neural invasion versus those driven by general tumorigenesis.



### Contrasts Explored:
1.  **PNI Effect (Tumor):** PNI-positive Tumor vs. PNI-negative Tumor
2.  **PNI Effect (Normal):** PNI-positive Normal vs. PNI-negative Normal
3.  **Tumor vs. Normal (+):** PNI-positive Tumor vs. PNI-positive Normal
4.  **Tumor vs. Normal (-):** PNI-negative Tumor vs. PNI-negative Normal
5.  **Extreme Contrast:** PNI-positive Tumor vs. PNI-negative Normal
6.  **Reverse Contrast:** PNI-negative Tumor vs. PNI-positive Normal

---

## 2. Project Overview
Pancreatic Ductal Adenocarcinoma (PDAC) progression is heavily influenced by **Perineural Invasion (PNI)**. This study utilizes the **GSE102238** dataset (100 samples) to map the molecular crosstalk between malignant cells and the peripheral nervous system.

## 3. Methodology & Workflow

The analysis was conducted using **R (v4.5.2)**. The integrated pipeline combines data acquisition, rigorous preprocessing, and functional interpretation as follows:

### 3.1. Analysis Pipeline
1.  **Data Acquisition:** Retrieval of raw data via `GEOquery`.
2.  **Normalization:** Log2 transformation and Quantile Normalization to stabilize expression distributions.
3.  **Differential Expression (DEA):** Modeled using the `limma` package across six distinct clinical contrasts.
4.  **Annotation:** Systematic mapping of probes to HGNC Symbols via `biomaRt` and relational merging with platform metadata.
5.  **Visualization:** Generation of high-fidelity plots to assess data distribution and DEG significance:
    * **Boxplot & Density:** `ggplot()` with `stat_boxplot` and `geom_density` to verify normalization.
    * **UMAP:** `umap()` algorithm followed by `geom_point()` to visualize 2D sample clustering.
    * **Volcano Plot:** Custom function `make_volcano()` mapping log2FC vs. -log10 Adjusted P-value.
    * **Heatmap:** `pheatmap()` using Global ANOVA and Ward.D2 clustering for top 50 DEGs.
    * **Scatter Plot:** `ggplot()` with `geom_smooth(method = "gam")` to profile gene expression stability (Mean vs. SD) across clinical cohorts.
    * **Venn Diagram:** `ggVennDiagram()` with a 6-set elliptical layout (`shape_id = "601"`) to identify core biomarkers across all clinical contrasts.
6.  **Functional Enrichment:** Pathway analysis using Gene Ontology (GO) and KEGG to interpret biological significance.

### 3.2. Pipeline Workflow
Below is the visual representation of the analytical steps performed in this project:

| Analysis Step |
| :---: |
| **GEO Dataset** |
| ▼ |
| **Data Preprocessing** |
| ▼ |
| **Normalization** |
| ▼ |
| **Differential Expression (limma)** |
| ▼ |
| **Annotation & Filtering** |
| ▼ |
| **Visualization** |
| ▼ |
| **Biological Interpretation** |

---


## 4. Key Findings
### 4.1. Transcriptomic Stability vs. Eruption
* **PNI-Specific Effects:** Contrasts within the same tissue type (Tumor-only or Normal-only) showed high transcriptomic stability with no significant DEGs.
* **Malignancy Drivers:** Direct comparisons between Tumor and Normal tissues revealed a massive "eruption" of differential expression.
* **Top Biomarkers:** Key upregulated genes identified include **CEACAM5, S100P, CST2,** and **TMPRSS4**.

### 4.2. Functional Pathways
Enrichment analysis highlighted critical pathways involved in:
* **Axon Guidance:** Providing a direct transcriptomic link to the PNI phenotype.
* **IgSF CAM Signaling:** Facilitating cell-adhesion and tumor dissemination.
* **Immune Modulation:** Significant enrichment in leukocyte dynamics and T-cell activation.

## 5. Conclusion
This research confirms that while the malignant phenotype is the primary driver of variance in PDAC, pathways associated with neural signaling (Axon Guidance) are uniquely positioned as potential master regulators of neural infiltration. These findings offer high-value candidates for biomarker discovery and therapeutic targeting.

## 6. References
* **Chen, Z., et al. (2023).** Cancers, 15(5), 1360.
* **Sarantis, P., et al. (2020).** World Journal of Gastrointestinal Oncology, 12(2), 173.
* **Sun, Y., et al. (2024).** Frontiers in Oncology, 14, 1421067.
* **Yang, M. W., et al. (2020).** Cancer Research, 80(10), 1991.

## 📂 Repository Structure
* **`/Dataset`**: Curated expression matrix and metadata from GSE102238.
* **`/Results`**: 
    * `Data_Tables/`: Statistical output of DEGs (CSV/Excel tables).
    * `Plots/`: Visualizations (Boxplot, density plot, UMAP plot, volcano plots, heatmaps, scatter plot, venn diagram, dot plot, and bar plot).
* **`/Script`**: R scripts used for normalization, DEA, and visualization.
---
© 2026 Yosia Jose Rasdiva Manurung. All Rights Reserved.
