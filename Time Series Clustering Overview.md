# Dynamic Time-Series Clustering of Immunisation Coverage Using DTW

## Overview

This repository presents a comprehensive analytical workflow for clustering district-level immunisation time-series data using Dynamic Time Warping (DTW). The objective is to identify groups of districts that exhibit similar temporal patterns in vaccination coverage, even when those patterns are misaligned in time.

The analysis focuses on birth-dose vaccines (BCG, OPV0, HepB0) over a two-year period (April 2020 – March 2022), capturing the dynamic evolution of immunisation coverage across districts. The workflow integrates data preprocessing, time-series clustering, statistical validation, and visualization to generate meaningful public health insights.

---

## What is Time-Series Clustering?

Time-series clustering involves grouping entities based on the similarity of their temporal patterns rather than static values. In this context, each district is represented as a time-series of immunisation coverage across months.

Unlike traditional clustering, time-series clustering:
- Preserves temporal ordering  
- Captures trends, fluctuations, and seasonality  
- Groups districts with similar trajectories over time  

---

## Dynamic Time Warping (DTW)

Dynamic Time Warping (DTW) is a distance measure designed specifically for time-series data.

### Why DTW?
- Handles **temporal misalignment** (e.g., delayed peaks)  
- Allows flexible matching between sequences  
- Captures similarity even when patterns are shifted in time  

### Concept:
DTW computes the optimal alignment between two time-series by stretching or compressing segments to minimize distance.

This makes it particularly suitable for:
- Epidemiological trends  
- Vaccination coverage patterns  
- Non-linear temporal dynamics  

---

## Objectives

- Cluster districts based on **temporal immunisation trajectories**  
- Identify **high-risk and low-risk clusters**  
- Capture **non-linear temporal similarities** using DTW  
- Validate clustering results using statistical metrics  

---

## Data Description

The dataset includes:
- Monthly immunisation coverage (%)  
- District-level observations  
- Vaccines:
  - BCG  
  - OPV 0 (Birth Dose)  
  - Hepatitis B0 (Birth Dose)  

Time period:
- April 2020 to March 2022  

---

## Methodology

### 1. Data Preprocessing
- Parsing and standardizing month information  
- Aggregating duplicate district-month entries  
- Filtering the study period  
- Handling missing values using linear interpolation  

---

### 2. Time-Series Transformation
- Data reshaped into wide format  
- Each district represented as a multivariate time-series  
- Missing values imputed for continuity  

---

### 3. DTW-Based Clustering

- Clustering performed using **partitional clustering (PAM-like)**  
- Distance metric: **DTW (Dynamic Time Warping)**  
- Number of clusters: **k = 4**  

Each cluster represents districts with similar vaccination trajectories.

---

### 4. Cluster Interpretation

- Average time-series computed per cluster  
- Clusters labeled based on coverage levels:
  - High coverage (Low risk)  
  - Moderate coverage (Medium risk)  
  - Low coverage (High risk)  

---

## Statistical Validation of Clusters

To ensure robustness and reliability, multiple validation techniques are applied:

---

### 1. Cophenetic Correlation Coefficient (CCC)

Measures how well the clustering structure preserves pairwise distances.

- Compares original distance matrix with hierarchical clustering distances  
- Values closer to **1** indicate better structure preservation  

---

### 2. Bootstrap Cluster Stability (pvclust)

Evaluates the stability of clusters using resampling:

- Generates bootstrap samples of the data  
- Computes cluster consistency across samples  
- Produces **Approximately Unbiased (AU) p-values**

#### Interpretation:
- AU > 0.95 → Very stable cluster  
- AU between 0.8–0.95 → Moderately stable  
- AU < 0.8 → Weak or unstable  

---

### 3. Silhouette Analysis

Measures how well each district fits within its assigned cluster.

- Range: **-1 to +1**  
- Higher values indicate better clustering  

#### Interpretation:
- Close to +1 → Well-clustered  
- Near 0 → Overlapping clusters  
- Negative → Possible misclassification  

Cluster-wise summaries provide:
- Average silhouette width  
- Cluster compactness and separation  

---

## Visualization

The workflow includes multiple visualization techniques:

### 1. Interactive Time-Series Plots
- Vaccine coverage trends per cluster  
- Implemented using Plotly  

### 2. Dimensionality Reduction
- **t-SNE** for non-linear visualization  
- **UMAP** for structure preservation  

These methods help visualize cluster separation in lower dimensions.

---

## Outputs

- Cluster assignments for each district  
- Risk classification (High, Medium, Low)  
- Interactive HTML visualizations  
- Static plots  
- File with district-level risk labels  

---

## Tools and Libraries

- **R Programming Language**
- Key packages:
  - `dtwclust` – time-series clustering  
  - `dplyr`, `tidyr` – data manipulation  
  - `plotly` – interactive visualization  
  - `cluster` – validation metrics  
  - `pvclust` – bootstrap stability  
  - `Rtsne`, `umap` – dimensionality reduction  

---

## Key Insights

- DTW enables detection of **non-linear temporal similarities**  
- Clustering reveals **hidden structure in vaccination trends**  
- Statistical validation ensures **robust and reliable grouping**  
- Results support **targeted public health interventions**  

---

## Future Work

- Spatio-temporal modeling  
- Integration with socio-economic indicators  
- Bayesian time-series clustering approaches  

---
