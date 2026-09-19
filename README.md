# Temporal Muscle Thickness and Machine Learning-Based Survival Prediction

R scripts accompanying the study:

**Temporal Muscle Thickness and Machine Learning-Based Survival Prediction in Patients with Brain Metastases Treated with CyberKnife Stereotactic Radiotherapy**

This repository contains the R scripts used for data preprocessing, statistical analysis, inter-rater reliability assessment, survival analysis, and machine learning-based survival prediction using temporal muscle thickness (TMT) in patients with brain metastases treated with CyberKnife stereotactic radiotherapy.

The analysis evaluates TMT measurements obtained from multiple independent raters and investigates their association with overall survival. In addition, machine learning approaches are used to explore survival prediction and temporal changes in model behaviour.

---

## Repository structure

```text
Temporal-Muscle-Thickness-Survival-ML/
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
