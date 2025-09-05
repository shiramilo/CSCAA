import pandas as pd
import numpy as np

# Load significant post-hoc results
posthoc_df = pd.read_csv('posthoc_res/posthoc_results.csv')

# Filter for significant results only
significant_results = posthoc_df[posthoc_df["p-value"] < 0.05].copy()

# Define interpretation function
def interpret_row(row):
    metadata = row['Metadata'].replace("_CAT", "")
    phenotype = row['Phenotype']
    comparison = row['Comparison']
    p_val = row['p-value']
    group1, group2 = comparison.split(' vs ')
    test_type = row['Post Hoc Test']

    # Decide if we deal with medians or counts
    if 'Group1_Median' in row and not pd.isna(row['Group1_Median']):
        median1 = row['Group1_Median']
        median2 = row['Group2_Median']

        if median1 > median2:
            interpretation = (f"'{phenotype}': '{group1}' shows significantly higher values "
                              f"than '{group2}' in '{metadata}' (median: {median1:.3f} vs {median2:.3f}; p={p_val:.2e}).")
            comparison_direction = f"{group1} > {group2}"
            fold_change = abs(median1) / abs(median2) if median2 != 0 else np.inf

        elif median2 > median1:
            interpretation = (f"'{phenotype}': '{group2}' shows significantly higher values "
                              f"than '{group1}' in '{metadata}' (median: {median2:.3f} vs {median1:.3f}; p={p_val:.2e}).")
            comparison_direction = f"{group2} > {group1}"
            fold_change = abs(median2) / abs(median1) if median1 != 0 else np.inf

        else:
            interpretation = (f"'{phenotype}': '{group1}' and '{group2}' in '{metadata}' have identical medians ({median1:.3f}), "
                              f"but differ significantly (p={p_val:.2e}).")
            comparison_direction = f"{group1} = {group2}"
            fold_change = 1.0

    elif 'Group1_Count' in row and not pd.isna(row['Group1_Count']):
        count1 = row['Group1_Count']
        count2 = row['Group2_Count']

        # Adjusted interpretation for categorical vs categorical tests
        interpretation = (f"There is a significant association between '{metadata}' ('{group1}') "
                          f"and '{phenotype}' ('{group2}') (counts: {count1} vs {count2}; p={p_val:.2e}).")
        comparison_direction = f"{metadata}:{group1} ↔ {phenotype}:{group2}"
        fold_change = count1 / count2 if count2 != 0 else np.inf

    else:
        interpretation = "Insufficient data"
        comparison_direction = "N/A"
        fold_change = np.nan

    return pd.Series({
        "Comparison_Direction": comparison_direction,
        "Fold_Change": fold_change,
        "Interpretation": interpretation
    })

# Apply interpretation
interpretation_results = significant_results.apply(interpret_row, axis=1)

# Insert 'Comparison_Direction' and 'Fold_Change' columns after 'Group2_Count'
insert_pos = significant_results.columns.get_loc("Group2_Count") + 1
significant_results.insert(loc=insert_pos, column="Comparison_Direction",
                           value=interpretation_results["Comparison_Direction"])
significant_results.insert(loc=insert_pos + 1, column="Fold_Change",
                           value=interpretation_results["Fold_Change"])

# Add Interpretation as last column
significant_results["Interpretation"] = interpretation_results["Interpretation"]

# Save results
significant_results.to_csv('posthoc_res/posthoc_interpretation.csv', index=False)

print("Post-hoc interpretation file with fold change created successfully!")
