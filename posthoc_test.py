# Cutoff: FDR < 0.05

# Pairwise Chi-square for significant Chi-Square Test results
# Pairwise Mann-Whitney U for significant Mann-Whitney U results
# Dunn's Test for significant Kruskal-Wallis results

import pandas as pd
import numpy as np
import scipy.stats as stats
import scikit_posthocs as sp
import os

# Ensure output directory exists
os.makedirs("posthoc_res", exist_ok=True)

# Load correlation results
data = pd.read_csv('corr_res/correlation_results_with_fdr.csv')

# Load original dataset (single file input)
original_data = pd.read_csv('input/all_cat_for_corr.csv', encoding='ISO-8859-1')

# Initialize list for storing post-hoc results
posthoc_results = []

# Iterate over significant results
for _, row in data.iterrows():
    metadata, phenotype, test, p, fdr_p = row['Metadata'], row['Phenotype'], row['Test'], row['p-value'], row['corrected p-value']
    
    # Apply post-hoc only if the FDR-adjusted p-value is significant
    if fdr_p < 0.05:
        valid_data = original_data[[phenotype, metadata]].dropna()
        valid_data[metadata] = valid_data[metadata].astype(str)  # <-- FIX HERE

        try:
            if test == 'Chi-square Test':
                contingency_table = pd.crosstab(valid_data[metadata], valid_data[phenotype])
                comparisons_done = set()
                for row_label in contingency_table.index:
                    for col_label in contingency_table.columns:
                        if row_label != col_label:
                            comparison_key = tuple(sorted([str(row_label), str(col_label)]))
                            if comparison_key not in comparisons_done:
                                comparisons_done.add(comparison_key)
                                observed = np.array([
                                    [contingency_table.loc[row_label, col_label], sum(contingency_table.loc[row_label, :]) - contingency_table.loc[row_label, col_label]],
                                    [sum(contingency_table[col_label]) - contingency_table.loc[row_label, col_label],
                                     contingency_table.values.sum() - sum(contingency_table.loc[row_label, :]) - sum(contingency_table[col_label]) + contingency_table.loc[row_label, col_label]]
                                ])
                                chi2, post_p, _, _ = stats.chi2_contingency(observed)
                                group1_count = observed[0].sum()
                                group2_count = observed[1].sum()
                                posthoc_results.append({
                                    'Metadata': metadata,
                                    'Phenotype': phenotype,
                                    'Comparison': f'{row_label} vs {col_label}',
                                    'Post Hoc Test': 'Pairwise Chi-square',
                                    'p-value': post_p,
                                    'Group1_Count': group1_count,
                                    'Group2_Count': group2_count
                                })

            elif test == 'Mann-Whitney U':
                categories = list(valid_data[metadata].unique())
                comparisons_done = set()
                for i in range(len(categories)):
                    for j in range(i + 1, len(categories)):
                        comparison_key = tuple(sorted([categories[i], categories[j]]))
                        if comparison_key not in comparisons_done:
                            comparisons_done.add(comparison_key)
                            group1 = valid_data[valid_data[metadata] == categories[i]][phenotype]
                            group2 = valid_data[valid_data[metadata] == categories[j]][phenotype]
                            stat, post_p = stats.mannwhitneyu(group1, group2, alternative='two-sided')
                            group1_median = np.median(group1)
                            group2_median = np.median(group2)
                            posthoc_results.append({
                                'Metadata': metadata,
                                'Phenotype': phenotype,
                                'Comparison': f'{categories[i]} vs {categories[j]}',
                                'Post Hoc Test': 'Pairwise Mann-Whitney U',
                                'p-value': post_p,
                                'Group1_Median': group1_median,
                                'Group2_Median': group2_median
                            })

            elif test == 'Kruskal-Wallis':
                posthoc = sp.posthoc_dunn(valid_data, val_col=phenotype, group_col=metadata, p_adjust='bonferroni')
                comparisons_done = set()
                for group1 in posthoc.index:
                    for group2 in posthoc.columns:
                        if group1 != group2:
                            comparison_key = tuple(sorted([group1, group2]))
                            if comparison_key not in comparisons_done:
                                comparisons_done.add(comparison_key)
                                group1_vals = valid_data[valid_data[metadata] == group1][phenotype]
                                group2_vals = valid_data[valid_data[metadata] == group2][phenotype]
                                group1_median = np.median(group1_vals)
                                group2_median = np.median(group2_vals)
                                posthoc_results.append({
                                    'Metadata': metadata,
                                    'Phenotype': phenotype,
                                    'Comparison': f'{group1} vs {group2}',
                                    'Post Hoc Test': "Dunn's Test",
                                    'p-value': posthoc.loc[group1, group2],
                                    'Group1_Median': group1_median,
                                    'Group2_Median': group2_median
                                })
        except Exception as e:
            print(f"Skipping post-hoc for {metadata} vs {phenotype} due to error: {e}")

# Save posthoc results
posthoc_df = pd.DataFrame(posthoc_results)
posthoc_df.to_csv('posthoc_res/posthoc_results.csv', index=False)

print(f'Post-hoc analysis completed successfully! {len(posthoc_results)} comparisons saved.')
