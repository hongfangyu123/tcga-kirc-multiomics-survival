# Input data (not included)

The pipeline expects three tab-separated files in this folder (or in the folder set by the
environment variable `KIRC_DATA_DIR`). They are processed TCGA-KIRC tables derived from the
NCI Genomic Data Commons (open-access tier) and are not redistributed here.

| File | Content |
|---|---|
| `integrated_data.txt` | ~750 MB. One row per molecular feature, one column per tumour. The first two columns are `Variable` (feature ID) and `Type` (`methylation` = Illumina 450K β-value, `mRNA` = log-scale RNA-seq expression); the remaining columns are named by TCGA patient barcode (e.g. `TCGA-AK-3425`). |
| `variable_info.txt` | Feature annotation, row-aligned with `integrated_data.txt`: `Variable`, `VarName` (CpG probe ID or Ensembl gene ID / gene symbol), `chr`, `strand`, `pos`, `type`, `relatedGene` (semicolon-separated), `relatedLocation` (e.g. TSS200, Body, 3'UTR). |
| `sample_info.txt` | One row per patient from the GDC clinical tables, including `patient`, `vital_status`, `days_to_last_follow_up`, `days_to_death`, `age_at_index`, `gender` and `ajcc_pathologic_stage`. |

Cohort used: 309 primary tumours with positive follow-up (102 deaths); 263 of them have AJCC stage
and are used for the survival models (79 deaths).
