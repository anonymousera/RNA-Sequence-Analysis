#!/usr/bin/env python
# coding: utf-8

# In[28]:


get_ipython().system('pip install cufflinks')


# In[31]:


get_ipython().system('pip install rpy2')


# In[3]:


import pandas as pd
import numpy as np
from scipy import stats
import seaborn as sns
import matplotlib.pyplot as plt
import cufflinks as cf
import plotly.offline as pyo
import plotly.graph_objs as go


# In[1]:


# import rpy2.robjects as robjects
# from rpy2.robjects.packages import importr
# cummeRbund = importr("cummeRbund")


# In[ ]:





# In[13]:



# Load control and disease data into pandas dataframes
control1 = pd.read_csv(r'C:\Users\erasa\OneDrive\Desktop\semester 6\sbv887\assign\control\c1.txt', sep='\t', header=None)
control2 = pd.read_csv(r'C:\Users\erasa\OneDrive\Desktop\semester 6\sbv887\assign\control\c2.txt', sep='\t',header=None)
control3 = pd.read_csv(r'C:\Users\erasa\OneDrive\Desktop\semester 6\sbv887\assign\control\c4.txt', sep='\t',header=None)
disease1 = pd.read_csv(r'C:\Users\erasa\OneDrive\Desktop\semester 6\sbv887\assign\disease\pd1.txt', sep='\t',header=None)
disease2 = pd.read_csv(r'C:\Users\erasa\OneDrive\Desktop\semester 6\sbv887\assign\disease\pd2.txt', sep='\t',header=None)
disease3 = pd.read_csv(r'C:\Users\erasa\OneDrive\Desktop\semester 6\sbv887\assign\disease\pd4.txt', sep='\t',header=None)
disease4 = pd.read_csv(r'C:\Users\erasa\OneDrive\Desktop\semester 6\sbv887\assign\disease\pd5.txt', sep='\t',header=None)


# In[178]:


genes=control1.iloc[:,0]
genes


# In[14]:


print(control1.shape)
print(control2.shape)
print(control3.shape)
print(disease4.shape)


# In[154]:


control1.rename(columns={1: 'sample 1'}, inplace=True)
control2.rename(columns={1: 'sample 2'}, inplace=True)
control3.rename(columns={1: 'sample 3'}, inplace=True)
disease1.rename(columns={1: 'sample 4'}, inplace=True)
disease2.rename(columns={1: 'sample 5'}, inplace=True)
disease3.rename(columns={1: 'sample 6'}, inplace=True)
disease4.rename(columns={1: 'sample 7'}, inplace=True)


# In[155]:


data_c = pd.concat([control1, control2.iloc[:,1],control3.iloc[:,1]], axis=1)
data_c.rename(columns={0: 'Gene Name'}, inplace=True)
# data_c.rename(columns={1: 'sample 1'}, inplace=True)
# data_c.rename(columns={2: 'sample 2'}, inplace=True)
# data_c.rename(columns={3: 'sample 3'}, inplace=True)


# In[156]:


print(data_c)


# In[158]:


data_cd=pd.concat([control1, control2.iloc[:,1],control3.iloc[:,1],disease1.iloc[:,1],disease2.iloc[:,1],disease3.iloc[:,1],disease4.iloc[:,1] ], axis=1)


# In[159]:


print(data_cd)


# In[160]:


data_norm = data_cd.iloc[:,1:].apply(lambda x: (x - np.mean(x)) / np.std(x), axis=0)


# In[161]:


data_norm


# In[146]:


#data_norm=pd.concat([control1.iloc[:,0],data_norm],axis=1)


# In[229]:


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


# In[240]:


diff_expr_genes


# In[230]:


genes_diff=genes.iloc[diff_expr_genes]
# genes_upreg=genes.iloc[upreg_genes]
# genes_downreg=genes.iloc[downreg_genes]


# In[239]:


genes_diff


# In[271]:


# write the dataframe to a CSV file
genes_diff.to_csv('genes_diff.csv', index=False)


# In[249]:


# genes_upreg


# In[250]:


# for value in genes_upreg:
#     print(value)


# In[251]:


# genes_downreg


# In[252]:


# for value in genes_downreg:
#     print(value)


# In[281]:


# Calculate log2 fold change and p-values
log2fc = data_norm.iloc[:,3:].mean(axis=1) - data_norm.iloc[:,:3].mean(axis=1)
_, pvals = stats.ttest_ind(data_norm.iloc[:,:3], data_norm.iloc[:,3:], axis=1)


# In[297]:


pvals


# In[290]:


# pvals_diff=pvals[diff_expr_genes]
pvals_diff=[]
for i in range(len(diff_expr_genes)):
    pvals_diff.append(pvals[diff_expr_genes[i]])
    


# In[298]:


pvals_diff


# In[288]:


log2fc_diff=log2fc.iloc[diff_expr_genes]
log2fc_diff


# In[300]:


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


# In[ ]:





# In[237]:


type(log2fc)


# In[241]:


upreg_genes = []
downreg_genes = []
for i in range(len(diff_expr_genes)):
    if log2fc[diff_expr_genes[i]]>0:
        upreg_genes.append(diff_expr_genes[i])
    if log2fc[diff_expr_genes[i]]<0:
        downreg_genes.append(diff_expr_genes[i])


# In[242]:


genes_upreg=genes.iloc[upreg_genes]
genes_downreg=genes.iloc[downreg_genes]


# In[243]:


genes_upreg


# In[244]:


genes_downreg


# In[247]:


len(pvals)


# In[277]:





# In[278]:


pvals_diff=pd.DataFrame(pvals_diff)
pvals_diff


# In[272]:


# write the dataframe to a CSV file
pvals_diff.to_csv('pvals_diff.csv', index=False)


# In[269]:


log2fc_diff=pd.DataFrame(log2fc_diff)
log2fc_diff


# In[273]:


# write the dataframe to a CSV file
log2fc_diff.to_csv('log2fc_diff.csv', index=False)


# In[183]:


# Create heatmap of differentially expressed genes
data_diff_expr = data_norm.loc[diff_expr_genes]
sns.clustermap(data_diff_expr.iloc[:,:], cmap='coolwarm', z_score=0)
plt.title('Differentially Expressed Genes')
plt.show()


# In[ ]:





# In[ ]:





# In[ ]:




