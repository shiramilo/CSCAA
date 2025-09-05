# Complete Statistical Correlation & Association Analysis of Large Phenotypic Data and Metadata

**Author:** Shira Milo  
**Year:** 2025

---

## Table of Contents
- [Part 1: Correlation & Association Analysis](#part-1-correlation--association-analysis)
  - [Data Preparation](#data-preparation)
  - [Normality Testing](#normality-testing)
  - [Correlation & Association Testing](#correlation--association-testing)
  - [Multiple Testing Correction](#multiple-testing-correction)
  - [Visualization](#visualization)
  - [Output Files](#output-files)
- [Part 2: Post-hoc Analysis](#part-2-post-hoc-analysis)
  - [Input & Setup](#input--setup)
  - [Identification of Significant Comparisons](#identification-of-significant-comparisons)
  - [Post-hoc Tests](#post-hoc-tests)
  - [Results and Outputs](#results-and-outputs)
- [Part 3: Post-hoc Interpretation](#part-3-post-hoc-interpretation)
  - [Input & Setup](#input--setup-1)
  - [Processing Pairwise Comparisons](#processing-pairwise-comparisons)
  - [Output](#output)

---

## Part 1: Correlation & Association Analysis

### Data Preparation
- **Data loading and identification:** identifies categorical and numerical variables (categorical column names should include a `CAT` suffix).
- **Reformatting:** casts categorical variables to the appropriate format.
- **Data integrity:** handles missing and infinite values.

### Normality Testing
- **Numerical variables:** tests for normality.  
  - **Shapiro–Wilk** for small datasets  
  - **Lilliefors** for larger ones  
- **Outputs:** saves distribution plots for visualization (in `data_dist/`).

### Correlation & Association Testing
**Numerical vs. Numerical**
- **Pearson correlation** (if both variables are normally distributed).
- **Spearman correlation** (if at least one variable is not normally distributed).

**Categorical vs. Numerical**
- **Mann–Whitney U test** (for binary categories).
- **Kruskal–Wallis test** (for multi-category comparisons).

**Categorical vs. Categorical**
- **Chi-square test** (for general associations).
- **Fisher’s Exact Test** (if expected frequencies are low).

### Multiple Testing Correction
- Adjusts *p*-values using **False Discovery Rate (FDR)** to reduce false positives.

### Visualization
- Generates a **clustered correlation heatmap** for numerical comparisons.
- Generates **individual boxplots** for categorical vs. numerical associations.

### Output Files
- **Correlation results:** `corr_res/correlation_results_with_fdr.csv`  
- **Normality test results:** `norm_res/normality_results.csv`  
- **Correlation heatmap:** saved in `corr_plot/`  
- **Boxplots:** saved in `boxplots/`

---

## Part 2: Post-hoc Analysis

### Input & Setup
- Loads **significant correlation results** from the primary analysis (`corr_res/correlation_results_with_fdr.csv`).
- Loads the **original dataset** to retrieve the categorical and numerical data.

### Identification of Significant Comparisons
- **Case selection:** iterates over significant associations (*FDR*-adjusted *p* < 0.05).
- **Data evaluation:** extracts independent (metadata) and dependent (phenotype) variables for each comparison.
- **Cleanup:** removes missing values.

### Post-hoc Tests
**Pairwise Chi-square Test (Categorical vs. Categorical)**
- Performed when a **significant Chi-square** association was observed in the main analysis.
- Builds **contingency tables** for each pair of categories.
- Runs **pairwise Chi-square** tests to determine which groups differ.

**Pairwise Mann–Whitney U Test (Binary Categorical vs. Numerical)**
- Applied when the main test was **Mann–Whitney U**.
- Compares two groups using a **non-parametric rank-based** test.

**Dunn’s Test (Multi-group Categorical vs. Numerical)**
- Performed when a **Kruskal–Wallis** test was significant.
- Executes **pairwise comparisons** between all categories using **Bonferroni** correction.

### Results and Outputs
- Stores all post-hoc results in `posthoc_res/posthoc_results.csv`.

---

## Part 3: Post-hoc Interpretation

### Input & Setup
- Loads post-hoc test results from `posthoc_res/posthoc_results.csv`.
- Filters for **significant** (*p* < 0.05) post-hoc comparisons.

### Processing Pairwise Comparisons
- Extracts the **two compared groups** (e.g., “True vs. False”).
- Identifies which group has **higher or lower** phenotype values.
- Determines whether a group is **overrepresented or underrepresented** in a phenotype.

### Output
- Writes processed results to `posthoc_res/posthoc_interpretation.csv`.
