### Numerical data
# Pearson correlation if both variables are normally distributed
# Spearman correlation if at least one variable is not normal.

### Categorical data vs. numerical data
# Mann-Whitney U test for binary categorical data
# Kruskal-Wallis test for multi-category comparisons

### Categorical data vs. categorical data
# Chi-square test for general associations
# Fisher’s Exact Test if expected frequencies are low

import pandas as pd
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.diagnostic import lilliefors
import os

# Ensure output directories exist
os.makedirs("data_dist", exist_ok=True)
os.makedirs("norm_res", exist_ok=True)
os.makedirs("corr_res", exist_ok=True)
os.makedirs("violin_plots", exist_ok=True)
os.makedirs("corr_plot", exist_ok=True)

# Load data
data_file_path = "input/all_cat_for_corr.csv"
data = pd.read_csv(data_file_path, encoding="ISO-8859-1")

# Identify categorical and numerical columns
categorical_cols = [col for col in data.columns if col.endswith("CAT")]
boolean_cols = data.select_dtypes(include=['bool']).columns.tolist()
categorical_cols += [col for col in boolean_cols if col not in categorical_cols]

numerical_cols = [col for col in data.columns if col not in categorical_cols and col != "strain_id"]

# Convert categorical data properly (including boolean columns)
data[categorical_cols] = data[categorical_cols].astype("category")

# Convert numerical data properly
data[numerical_cols] = data[numerical_cols].apply(pd.to_numeric, errors='coerce')

# Remove infinite values from numerical data
data.replace([np.inf, -np.inf], np.nan, inplace=True)

# **Step 1: Normality Testing & Distribution Plots**
normality_results = {}
for column in numerical_cols:
    if data[column].nunique() > 1:
        plt.figure()
        sns.histplot(data[column].dropna(), kde=True)
        plt.title(f'Distribution of {column}')
        plt.savefig(f'data_dist/distribution_{column}.png')
        plt.close()

        try:
            sample_size = len(data[column].dropna())
            if sample_size < 500:
                stat, p = stats.shapiro(data[column].dropna())
                test_type = "Shapiro-Wilk"
            else:
                stat, p = lilliefors(data[column].dropna())
                test_type = "Lilliefors"

            normality_results[column] = {'Test': test_type, 'Statistic': stat, 'p-value': p}

        except Exception:
            normality_results[column] = {'Test': 'Error', 'Statistic': 'Error', 'p-value': 'Error'}

pd.DataFrame(normality_results).T.to_csv('norm_res/normality_results.csv')

# **Step 2: Correlations (Num-Num and Cat-Num)**
results = []
p_values = []  
processed_pairs = set()

for phenotype in numerical_cols:
    for metadata in categorical_cols + numerical_cols:
        if phenotype == metadata or (metadata, phenotype) in processed_pairs:
            continue

        processed_pairs.add((phenotype, metadata))

        valid_data = data[[phenotype, metadata]].dropna()
        if valid_data.shape[0] <= 5:
            continue

        try:
            if metadata in numerical_cols:
                normal_p1 = normality_results.get(phenotype, {}).get('p-value', 1)
                normal_p2 = normality_results.get(metadata, {}).get('p-value', 1)

                if normal_p1 > 0.05 and normal_p2 > 0.05:
                    stat, p = stats.pearsonr(valid_data[metadata], valid_data[phenotype])
                    test_type = "Pearson"
                else:
                    stat, p = stats.spearmanr(valid_data[metadata], valid_data[phenotype])
                    test_type = "Spearman"

            else:
                filtered_groups = [valid_data[valid_data[metadata] == cat][phenotype].dropna() 
                                   for cat in valid_data[metadata].unique() 
                                   if len(valid_data[valid_data[metadata] == cat]) > 5]

                if len(filtered_groups) < 2:
                    continue

                if len(filtered_groups) == 2:
                    stat, p = stats.mannwhitneyu(*filtered_groups, alternative="two-sided", method="asymptotic")
                    test_type = "Mann-Whitney U"
                else:
                    stat, p = stats.kruskal(*filtered_groups)
                    test_type = "Kruskal-Wallis"

            if np.isnan(stat) or np.isnan(p):
                continue

            results.append({'Metadata': metadata, 'Phenotype': phenotype, 'Test': test_type, 'Statistic': stat, 'p-value': p})
            p_values.append(p)

        except Exception:
            continue

# **Categorical vs. Categorical (Chi-square & Fisher's Exact)**
for i, cat1 in enumerate(categorical_cols):
    for cat2 in categorical_cols[i+1:]:
        valid_data = data[[cat1, cat2]].dropna()
        if valid_data.shape[0] <= 5:
            continue

        contingency_table = pd.crosstab(valid_data[cat1], valid_data[cat2])
        if contingency_table.shape[0] < 2 or contingency_table.shape[1] < 2:
            continue

        try:
            if (contingency_table.values < 5).any() and contingency_table.size == 4:
                stat, p = stats.fisher_exact(contingency_table)
                test_type = "Fisher's Exact Test"
            else:
                stat, p, _, _ = stats.chi2_contingency(contingency_table)
                test_type = "Chi-square Test"

            results.append({'Metadata': cat1, 'Phenotype': cat2, 'Test': test_type, 'Statistic': stat, 'p-value': p})
            p_values.append(p)

        except Exception:
            continue

# **Apply FDR Correction**
if p_values:
    _, corrected_p_values, _, _ = multipletests(p_values, method='fdr_bh')
    for i, result in enumerate(results):
        result['corrected p-value'] = corrected_p_values[i]
else:
    for result in results:
        result['corrected p-value'] = 'NA'

pd.DataFrame(results).to_csv('corr_res/correlation_results_with_fdr.csv', index=False)

# **Step 3: Clustered Heatmap**
if len(numerical_cols) > 1:
    correlation_matrix = data[numerical_cols].corr(method='spearman')
    clustergrid = sns.clustermap(correlation_matrix, cmap="BrBG", linewidths=0.1, figsize=(60, 50), annot=False, method='average')
    clustergrid.fig.suptitle("Clustered Numerical Correlation Heatmap", fontsize=24)
    plt.savefig("corr_plot/numerical_correlation_clustered.png")
    plt.close()

# **Step 4: Violin plots with significance and log-transform for MIC variables**
for metadata in categorical_cols:
    for phenotype in numerical_cols:
        valid_data = data[[metadata, phenotype]].dropna().copy()
        
        # Log-transform values if phenotype contains 'MIC'
        if "MIC" in phenotype.upper():
            valid_data[phenotype] = np.log10(valid_data[phenotype].replace(0, np.nan).dropna())

        if valid_data.shape[0] > 5 and valid_data[metadata].nunique() in [2, 3, 4, 5]:
            plt.figure(figsize=(12, 9))
            ax = sns.violinplot(
                x=valid_data[metadata].astype(str), 
                y=valid_data[phenotype], 
                hue=valid_data[metadata].astype(str),
                palette="Accent", 
                inner="box"
            )
            plt.xticks(rotation=45)
            plt.title(f'{metadata} vs. {phenotype}' + (' (log10-transformed)' if 'MIC' in phenotype.upper() else ''))

            groups = [valid_data[valid_data[metadata] == cat][phenotype].dropna() 
                      for cat in valid_data[metadata].unique()]

            if len(groups) == 2:
                _, p_value = stats.mannwhitneyu(*groups, method="asymptotic")
            else:
                _, p_value = stats.kruskal(*groups)

            ax.text(0.5, ax.get_ylim()[1] * 0.95, f"p = {p_value:.3g}",
                    ha='center', fontsize=12, color='black')

            plt.savefig(f'violin_plots/{metadata}_vs_{phenotype}.png')
            plt.close()

print("Analysis completed successfully!")
