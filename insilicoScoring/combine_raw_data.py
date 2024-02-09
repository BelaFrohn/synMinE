import numpy as np
import pandas as pd

df_combined = pd.read_csv('./final_scores.tsv', delimiter='\t', names = ['Name','combined_ccore'])
df_mem = pd.read_csv('./membrane_binding_scores.tsv', delimiter='\t')
df_minD = pd.read_csv('./interaction_scores.tsv', delimiter='\t')
df_dimer = pd.read_csv('./dimer_PAEs.tsv', delimiter='\t')
df_sol = pd.read_csv('./solubilities.tsv', delimiter='\t')


df = pd.merge(df_combined, df_mem,  on= 'Name')
df = pd.merge(df, df_minD,  on= 'Name')
df = pd.merge(df, df_dimer, on = 'Name')
df = pd.merge(df, df_sol, on = 'Name')


df_blinding = pd.read_csv('./blinding_map.tsv', delimiter='\t', names = ['Name', 'Version'])
df = pd.merge(df, df_blinding, on = 'Name', how='outer')

print(df)

df.to_csv('./raw_data_combined.csv')
