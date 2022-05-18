from scipy.stats import pearsonr, spearmanr, ttest_rel, ttest_ind
from statsmodels.stats.multitest import multipletests
from tabulate import tabulate
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# Read the data from the CSV files.
healthy = pd.read_csv('data/lusc-rsem-fpkm-tcga_paired.txt', sep='\t')
cancer = pd.read_csv('data/lusc-rsem-fpkm-tcga-t_paired.txt', sep='\t')


pd.options.display.max_columns = None

# Shape is a method to get [rows,columns]
noOfRows = healthy.shape[0]
noOfColumns = healthy.shape[1]

# FILTERING
i = 0

while i < noOfRows:
    hZeros = 0
    cZeros = 0
    j = 2
    # print(i)
    while j < noOfColumns:
        if(healthy.iloc[i, j] == 0.00):
            hZeros += 1
        if(cancer.iloc[i, j] == 0.00):
            cZeros += 1
        j += 1
    if((hZeros > ((noOfColumns-2)/2)) or (cZeros > ((noOfColumns-2)/2))):
        healthy = healthy.drop(healthy.index[i])
        cancer = cancer.drop(cancer.index[i])
        noOfRows -= 1
        i -= 1
    i += 1

# PEARSON CORRELATION
noOfRows = healthy.shape[0]

i = 0
highestCC = -2
lowestCC = 2
highestCC_ind = -1
lowestCC_ind = -1
while i < noOfRows:
    # print(i)
    Gene_h = healthy.iloc[i, 2:]
    Gene_c = cancer.iloc[i, 2:]
    pCoeff, p = pearsonr(Gene_h, Gene_c)
    highestCC = max(pCoeff, highestCC)
    if(pCoeff == highestCC):
        highestCC_ind = i
    lowestCC = min(pCoeff, lowestCC)
    if(pCoeff == lowestCC):
        lowestCC_ind = i
    i += 1

# 1) Values
print('Highest CC is:')
print(highestCC)
print('Its name is:')
print(healthy.iloc[highestCC_ind, 0])
print('Lowest CC is:')
print(lowestCC)
print('Its name is:')
print(healthy.iloc[lowestCC_ind, 0])


# 2) Plotting
highest_h = healthy.iloc[highestCC_ind, 2:]
highest_c = cancer.iloc[highestCC_ind, 2:]

highest_hINT = [int(i) for i in highest_h]
highest_cINT = [int(i) for i in highest_c]

lowest_h = healthy.iloc[lowestCC_ind, 2:]
lowest_c = cancer.iloc[lowestCC_ind, 2:]

lowest_hINT = [int(i) for i in lowest_h]
lowest_cINT = [int(i) for i in lowest_c]

plt.scatter(highest_hINT, highest_cINT)
plt.plot(np.unique(highest_hINT), np.poly1d(np.polyfit(
    highest_hINT, highest_cINT, 1))(np.unique(highest_hINT)), color='black')
plt.show()

plt.scatter(lowest_hINT, lowest_cINT)
plt.plot(np.unique(lowest_hINT), np.poly1d(np.polyfit(
    lowest_hINT, lowest_cINT, 1))(np.unique(lowest_hINT)), color='black')
plt.show()


# Hypothesis Testing
list_rel = []
list_ind = []
list_rel_ind = []
list_ind_ind = []

i = 0
while i < noOfRows:
    # print(i)
    Gene_h = healthy.iloc[i, 2:]
    Gene_c = cancer.iloc[i, 2:]
    p_val_rel = ttest_rel(Gene_h, Gene_c).pvalue
    p_val_ind = ttest_ind(Gene_h, Gene_c).pvalue
    list_rel.append(p_val_rel)
    list_rel_ind.append(healthy.iloc[i, 0])
    list_ind.append(p_val_ind)
    list_ind_ind.append(healthy.iloc[i, 0])
    i += 1

# Applying FDR multiple tests correction method
corrected_p_values_rel = multipletests(list_rel, method='fdr_bh')[1]
corrected_p_values_ind = multipletests(list_ind, method='fdr_bh')[1]

# P-Values before and after FDR correction
significance_genes_rel = pd.DataFrame(
    {'Gene_name': list_rel_ind, 'p_values': list_rel, 'p_values_fdr': corrected_p_values_rel})
significance_genes_ind = pd.DataFrame(
    {'Gene_name': list_ind_ind, 'p_values': list_ind, 'p_values_fdr': corrected_p_values_ind})


significance_genes_rel['significance:p_value'] = significance_genes_rel['p_values'].apply(
    lambda x: x < 0.05)
significance_genes_rel['significance:p_value_fdr'] = significance_genes_rel['p_values_fdr'].apply(
    lambda x: x < 0.05)
significance_genes_ind['significance:p_value'] = significance_genes_ind['p_values'].apply(
    lambda x: x < 0.05)
significance_genes_ind['significance:p_value_fdr'] = significance_genes_ind['p_values_fdr'].apply(
    lambda x: x < 0.05)

# print(
#     tabulate(significance_genes_rel.iloc[:10], headers='keys', tablefmt='psql'))
# print(
#     tabulate(significance_genes_ind.iloc[:10], headers='keys', tablefmt='psql'))


# Comparing common and distinct genes
diffrentially_genes_rel_Unchanged = significance_genes_rel[
    significance_genes_rel['significance:p_value_fdr'] == True]
diffrentially_genes_rel_Changed = significance_genes_rel[(significance_genes_rel['significance:p_value_fdr']
                                                          == False) & (significance_genes_rel['significance:p_value'] == True)]
diffrentially_genes_ind_Unchanged = significance_genes_ind[
    significance_genes_ind['significance:p_value_fdr'] == True]
diffrentially_genes_ind_Changed = significance_genes_ind[(significance_genes_ind['significance:p_value_fdr']
                                                          == False) & (significance_genes_ind['significance:p_value'] == True)]

DEG_rel_names_ch = diffrentially_genes_rel_Changed[[
    'Gene_name', 'p_values', 'p_values_fdr']]
DEG_rel_names_un = diffrentially_genes_rel_Unchanged[[
    'Gene_name', 'p_values', 'p_values_fdr']]
DEG_ind_names_ch = diffrentially_genes_ind_Changed[[
    'Gene_name', 'p_values', 'p_values_fdr']]
DEG_ind_names_un = diffrentially_genes_ind_Unchanged[[
    'Gene_name', 'p_values', 'p_values_fdr']]


# print(len(DEG_rel_names_un))
# print(len(DEG_ind_names_un))
# print(len(DEG_rel_names_ch))
# print(len(DEG_ind_names_ch))


print('Common genes - Paired')
print(tabulate(DEG_rel_names_un.iloc[:10], headers='keys', tablefmt='psql'))
print('Common genes - Independent')
print(tabulate(DEG_ind_names_un.iloc[:10], headers='keys', tablefmt='psql'))
print('Distinct genes - Paired')
print(tabulate(DEG_rel_names_ch.iloc[:10], headers='keys', tablefmt='psql'))
print('Distinct genes - Independet')
print(tabulate(DEG_ind_names_ch.iloc[:10], headers='keys', tablefmt='psql'))
