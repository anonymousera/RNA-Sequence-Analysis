#!/usr/bin/env python
# coding: utf-8
"""RNA-seq differential expression analysis: Parkinson's disease vs. control.


Data source: GEO accession GSE206308 (control vs. Parkinson's disease RNA-seq).
Expected inputs: seven tab-separated, headerless two-column files (gene ID,
count) — three control samples and four disease samples — placed under
data/control/ and data/disease/ relative to the project root. See README.md
for setup and reproduction instructions.

Outputs (written to the current working directory):
    genes_diff.csv     - names of differentially expressed genes
    pvals_diff.csv     - p-values for the differentially expressed genes
    log2fc_diff.csv    - log2 fold-change values for the differentially expressed genes
    Volcano plot and hierarchical-clustering heatmap (displayed via matplotlib/seaborn)

"""

import pandas as pd
import numpy as np
from scipy import stats
import seaborn as sns
import matplotlib.pyplot as plt
import cufflinks as cf
import plotly.offline as pyo
import plotly.graph_objs as go

# import rpy2.robjects as robjects
# from rpy2.robjects.packages import importr
# cummeRbund = importr("cummeRbund")


# --- Load data ---------------------------------------------------------

# Load control and disease data into pandas dataframes.
# Place the corresponding GSE206308 sample files under data/control/ and data/disease/.
control1 = pd.read_csv('data/control/c1.txt', sep='\t', header=None)
control2 = pd.read_csv('data/control/c2.txt', sep='\t',header=None)
control3 = pd.read_csv('data/control/c4.txt', sep='\t',header=None)
disease1 = pd.read_csv('data/disease/pd1.txt', sep='\t',header=None)
disease2 = pd.read_csv('data/disease/pd2.txt', sep='\t',header=None)
disease3 = pd.read_csv('data/disease/pd4.txt', sep='\t',header=None)
disease4 = pd.read_csv('data/disease/pd5.txt', sep='\t',header=None)

genes=control1.iloc[:,0]
genes

print(control1.shape)
print(control2.shape)
print(control3.shape)
print(disease4.shape)

control1.rename(columns={1: 'sample 1'}, inplace=True)
control2.rename(columns={1: 'sample 2'}, inplace=True)
control3.rename(columns={1: 'sample 3'}, inplace=True)
disease1.rename(columns={1: 'sample 4'}, inplace=True)
disease2.rename(columns={1: 'sample 5'}, inplace=True)
disease3.rename(columns={1: 'sample 6'}, inplace=True)
disease4.rename(columns={1: 'sample 7'}, inplace=True)

data_c = pd.concat([control1, control2.iloc[:,1],control3.iloc[:,1]], axis=1)
data_c.rename(columns={0: 'Gene Name'}, inplace=True)
# data_c.rename(columns={1: 'sample 1'}, inplace=True)
# data_c.rename(columns={2: 'sample 2'}, inplace=True)
# data_c.rename(columns={3: 'sample 3'}, inplace=True)

print(data_c)

data_cd=pd.concat([control1, control2.iloc[:,1],control3.iloc[:,1],disease1.iloc[:,1],disease2.iloc[:,1],disease3.iloc[:,1],disease4.iloc[:,1] ], axis=1)

print(data_cd)


# --- Normalize (z-score per sample) ------------------------------------

data_norm = data_cd.iloc[:,1:].apply(lambda x: (x - np.mean(x)) / np.std(x), axis=0)

data_norm

#data_norm=pd.concat([control1.iloc[:,0],data_norm],axis=1)


# --- Differential expression (t-test, control vs. disease) -------------

diff_expr_genes = []

for i in data_norm.index:
    c = data_norm.iloc[i,:3]
    d = data_norm.iloc[i,3:]
    t, pval = stats.ttest_ind(c, d)
    if pval < 0.05 and abs(t) > 1:
        diff_expr_genes.append(i)
#         if np.mean(d) - np.mean(c) > 0:
#
#         if np.mean(d) - np.mean(c) < 0:
#             downreg_genes.append(i)
# print(len(upreg_genes))
# print(len(downreg_genes))

diff_expr_genes

genes_diff=genes.iloc[diff_expr_genes]
# genes_upreg=genes.iloc[upreg_genes]
# genes_downreg=genes.iloc[downreg_genes]

genes_diff

# write the dataframe to a CSV file
genes_diff.to_csv('genes_diff.csv', index=False)

# genes_upreg

# for value in genes_upreg:
#     print(value)

# genes_downreg

# for value in genes_downreg:
#     print(value)


# --- Fold change and p-values for all genes -----------------------------

# Calculate log2 fold change and p-values
log2fc = data_norm.iloc[:,3:].mean(axis=1) - data_norm.iloc[:,:3].mean(axis=1)
_, pvals = stats.ttest_ind(data_norm.iloc[:,:3], data_norm.iloc[:,3:], axis=1)

pvals

# pvals_diff=pvals[diff_expr_genes]
pvals_diff=[]
for i in range(len(diff_expr_genes)):
    pvals_diff.append(pvals[diff_expr_genes[i]])


pvals_diff

log2fc_diff=log2fc.iloc[diff_expr_genes]
log2fc_diff


# --- Volcano plot ---------------------------------------------------------

# # Find indices of max and min log2fc values with p-value < 0.05
# idx_max = np.argmax(log2fc)
# idx_min = np.argmin(log2fc[(pvals < 0.05)])

# Plot volcano plot
plt.figure(figsize=(10, 8))
plt.scatter(log2fc, -np.log10(pvals), c=['red' if gene in diff_expr_genes else 'blue' for gene in data_norm.index])
plt.axhline(-np.log10(0.05), color='gray', linestyle='--')
plt.axvline(-1, color='gray', linestyle='--')
plt.axvline(1, color='gray', linestyle='--')
plt.xlabel('log2 fold change')
plt.ylabel('-log10(p-value)')
plt.title('Volcano Plot')
# plt.text(log2fc[idx_max], -np.log10(pvals[idx_max]), "RAD52", ha='center', va='bottom', fontsize=12)
# plt.text(log2fc[idx_min], -np.log10(pvals[idx_min]), "pam52", ha='center', va='bottom', fontsize=12)
plt.show()

type(log2fc)


# --- Split into up-/down-regulated gene sets ------------------------------

upreg_genes = []
downreg_genes = []
for i in range(len(diff_expr_genes)):
    if log2fc[diff_expr_genes[i]]>0:
        upreg_genes.append(diff_expr_genes[i])
    if log2fc[diff_expr_genes[i]]<0:
        downreg_genes.append(diff_expr_genes[i])

genes_upreg=genes.iloc[upreg_genes]
genes_downreg=genes.iloc[downreg_genes]

genes_upreg

genes_downreg

len(pvals)

pvals_diff=pd.DataFrame(pvals_diff)
pvals_diff

# write the dataframe to a CSV file
pvals_diff.to_csv('pvals_diff.csv', index=False)

log2fc_diff=pd.DataFrame(log2fc_diff)
log2fc_diff

# write the dataframe to a CSV file
log2fc_diff.to_csv('log2fc_diff.csv', index=False)


# --- Heatmap of differentially expressed genes ----------------------------

# Create heatmap of differentially expressed genes
data_diff_expr = data_norm.loc[diff_expr_genes]
sns.clustermap(data_diff_expr.iloc[:,:], cmap='coolwarm', z_score=0)
plt.title('Differentially Expressed Genes')
plt.show()
