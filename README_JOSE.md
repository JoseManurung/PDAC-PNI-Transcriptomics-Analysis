> [!IMPORTANT]
> This code and analysis are part of a research project. Please contact the author before using this script for any publication or commercial purposes.


# Integrative Transcriptomic Profiling of Perineural Invasion (PNI) Signatures in Pancreatic Ductal Adenocarcinoma (PDAC)

**Author:** Yosia Jose Rasdiva Manurung  
**Affiliation:** Diponegoro University  
**Topic:** Multi-Contrast Bioinformatics Study (GSE102238)

---

## 1. Research Objectives & Comparison Groups
The primary goal of this study is to identify Differentially Expressed Genes (DEGs) across **six PNI-related clinical contrasts**. This systematic approach allows us to isolate the specific transcriptomic signatures driven by neural invasion versus those driven by general tumorigenesis.



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

## 3. Methodology
The analysis was conducted using **R (v4.5.2)** with the following pipeline:
* **Data Acquisition:** Retrieval of raw data via `GEOquery`.
* **Normalization:** Log2 transformation and Quantile Normalization to stabilize expression distributions.
* **Differential Expression (DEA):** Modeled using the `limma` package across six distinct clinical contrasts.
* **Annotation:** Systematic mapping of probes to HGNC Symbols via `biomaRt`.
* **Functional Enrichment:** Pathway analysis using Gene Ontology (GO) and KEGG.

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

---
© 2026 Yosia Jose Rasdiva Manurung. All Rights Reserved.