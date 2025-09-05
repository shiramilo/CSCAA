Complete Statistical Correlation & Association Analysis of large phenotypic data and metadata

Shira Milo
2025
 

PART 1: Correlation & Association Analysis
1.1.	Data Preparation
•	Data loading and identification: identifies categorical and numerical variables (categorical column names should include a “CAT” suffix).
•	Reformats categorical variables to the appropriate format.
•	Ensure data integrity: handles missing and infinite values.
1.2.	Normality Testing
•	Numerical variables: Tests if follow a normal distribution.
•	Shapiro-Wilk for small datasets and Lilliefors for larger ones.
•	Saves distribution plots for visualization (in data_dist ).
1.3.	Correlation & Association Testing
•	Numerical vs. Numerical:
o	Pearson correlation (if both variables are normally distributed).
o	Spearman correlation (if at least one variable is not normally distributed).
•	Categorical vs. Numerical:
o	Mann-Whitney U test (for binary categories).
o	Kruskal-Wallis test (for multi-category comparisons).
•	Categorical vs. Categorical:
o	Chi-square test (for general associations).
o	Fisher’s Exact Test (if expected frequencies are low).
1.4.	Multiple Testing Correction
•	Adjusts p-values using False Discovery Rate (FDR) to reduce false positives.
1.5.	Visualization
•	Generates a clustered correlation heatmap for numerical comparisons.
•	Generates individual boxplots for categorical vs. numerical associations.
Output Files
a.	Correlation results: saved in corr_res/correlation_results_with_fdr.csv
b.	Normality test results: saved in norm_res/normality_results.csv
c.	Correlation heatmap: saved in corr_plot folder
d.	Boxplots: saved in boxplots folder


PART 2: Post-hoc Analysis
2.1.	Input & Setup
•	Loads significant correlation results from the primary analysis (corr_res/correlation_results_with_fdr.csv).
•	Loads the original dataset to retrieve the categorical and numerical data.
2.2.	Identification of Significant Comparisons
•	Case selection: iterates over significant associations (FDR-adjusted p-value < 0.05).
•	Data evaluation: extracts independent (metadata) and dependant (phenotype) variables from each comparison.
•	Data clean up: removes missing values.
2.3.	Post-hoc Tests
•	Pairwise Chi-square Test (Categorical vs. Categorical)
o	 Performed when a significant Chi-square test association was observed in the main analysis.
o	Contingency tables are created for each pair of categories.
o	Pairwise Chi-square tests determines which groups differ significantly.
•	Pairwise Mann-Whitney U Test (Binary Categorical vs. Numerical)
o	Applied when the main test was Mann-Whitney U.
o	Compares two groups using a non-parametric rank-based test.
•	Dunn’s Test (Multi-group Categorical vs. Numerical)
o	Performed when a Kruskal-Wallis test was significant.
o	Performs pairwise comparisons between all categories using Bonferroni correction.
2.4.	Results and Outputs
•	Stores all post-hoc results in posthoc_res/posthoc_results.csv.

 
PART 3: Post-hoc Interpretation
3.1.	Input & Setup
•	Loads post-hoc test results from posthoc_res/posthoc_results.csv.
•	Filters for significant (p < 0.05) post-hoc comparisons.
3.2.	Processing Pairwise Comparisons
•	Extracts the two compared groups (e.g., "True vs. False").
•	Identifies which group has higher or lower phenotype values.
•	Determines whether a group is overrepresented or underrepresented in a phenotype.
3.3.	Output
•	Outputs processed results in posthoc_res/posthoc_interpretation.csv.


