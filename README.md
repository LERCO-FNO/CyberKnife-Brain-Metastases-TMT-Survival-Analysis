# TMT Analysis Pipeline

R scripts for statistical analysis of the **Trail Making Test (TMT)** in a clinical cohort, including descriptive analysis, inter-rater reliability, survival analysis, and concept-drift analysis.

The repository contains the analysis scripts used to process TMT measurements obtained from multiple raters and to evaluate their relationship with overall survival.

---

## Repository structure

```text
TMT-analysis/
│
├── R/
│   ├── 01_demographic_and_clinical_data.R
│   ├── 02_TMT_survival_analysis.R
│   ├── 03_TMT_interrater_reliability.R
│   └── 04_TMT_survival_concept_drift.R
│
├── Input_data/
│   └── README.md
│
├── results/
│   ├── figures/
│   └── tables/
│
└── README.md
