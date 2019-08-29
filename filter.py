''' this file filetrs a file based on values of column from another file'''

import pandas as pd
import csv
import numpy as np

data=pd.read_csv('/data/salomonis2/NCI-R01/TCGA-BREAST-CANCER/Anukana/TF/Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt', sep='\t')
df=pd.read_csv('/data/salomonis2/NCI-R01/TCGA-BREAST-CANCER/Anukana/TF/TF2list.txt', sep='\t')
TF_list = df["Symbol"].tolist() #converts the column to a list
data=data[data['Symbol'].isin(TF_list)] #filters dataset based on the list
data=data.replace(r'^\s*$', np.nan)
data=data.fillna(data.mean())
data=data.round(decimals=9)

data.to_csv('/data/salomonis2/NCI-R01/TCGA-BREAST-CANCER/Anukana/TF/PSI_filtered.txt', sep='\t')