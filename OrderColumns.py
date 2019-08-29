''' this file orders the columns the same way as another file--sample order selection'''

import pandas as pd
import csv

data=pd.read_csv('/data/salomonis2/NCI-R01/TCGA-BREAST-CANCER/Anukana/TF/exp.TCGA-BRCA-steady-state.txt', sep='\t')
df=pd.read_csv('/data/salomonis2/NCI-R01/TCGA-BREAST-CANCER/Anukana/TF/PSI_filtered_final.txt', sep='\t')
header_list=my_data.columns.tolist() #get the header columns as list
df = df[header_list] #order the column of df same as data DataFrame

print(df.head())

df.to_csv('/data/salomonis2/NCI-R01/TCGA-BREAST-CANCER/Anukana/TF/PSI_filtered_final.txt', sep='\t')