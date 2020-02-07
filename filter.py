''' this file filetrs a file based on values of a column from another file'''

import pandas as pd
import csv
import numpy as np

data=pd.read_csv('/data/salomonis2/NCI-R01/TCGA-BREAST-CANCER/Anukana/SurvivalAnalysis/results_confounding_75p.txt ', sep='\t')
df=pd.read_csv('/data/salomonis2/NCI-R01/TCGA-BREAST-CANCER/Anukana/TF/TFlist.txt', sep='\t')
TF_list = df["Symbol"].tolist() #converts rows of the column to a list
data=data[data['Gene_name'].isin(TF_list)] #filters dataset based on the list
#data=data.replace(r'^\s*$', np.nan)
#data=data.fillna(data.mean())
#data=data.round(decimals=9)

data.to_csv('/data/salomonis2/NCI-R01/TCGA-BREAST-CANCER/Anukana/SurvivalAnalysis/results_counfounding_75p_TF.txt', sep='\t')