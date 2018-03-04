import pandas as pd
import os

# Methylation files is very large, and repeats many columns. Read and save just the beta value columns, and one copy of gene name

df = pd.read_csv("GBM.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.2.txt",delimiter='\t',index_col=0,header=None,dtype=str)
dff = pd.DataFrame(index=df.index)
dff['GeneSymbol'] = df[df.columns[1]]

for i,elem in enumerate(df.loc['Composite Element REF']):
    if (elem=="Beta_value"):
        patient = df.columns[i]
        dff[patient]=df[patient]

# Remove 0th row, which just has "composite element"
dff.drop(dff.index[0], inplace=True)

# -------- Strip out all patients without SVZ status ----------

# Rename patient ids with only relevant info
cols = ['GeneSymbol']
for i in dff.columns[1:]:
    if int(i.split('-')[3][:2])<10: cols.append(i)
dff = dff[cols]
dff.columns = [i[:12].replace('.','-') for i in dff.columns]

# Fill empty gene names with '-'
dff.loc[dff['GeneSymbol'].isnull(),'GeneSymbol']='-'
dff.dropna(inplace=True)

# Compare patients with SVZ status
clusters = pd.read_csv("../data/patient_status.csv",index_col=0)
clusters.sort_values(by="SVZ",inplace=True)

shared_samples = [i for i in clusters.index if i in dff.columns]   
dff = dff[["GeneSymbol",]+shared_samples]

dff.to_csv("methylation.csv")

