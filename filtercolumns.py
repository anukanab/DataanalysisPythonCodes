''' this file filetrs columns of a file based on columns from another file'''

import pandas as pd
import csv
import numpy as np
import sys

data=pd.read_csv('/data/salomonis2/NCI-R01/TCGA-BREAST-CANCER/Anukana/UO1analysis/tcga_rsem_isopct.txt', sep='\t', index_col=1)
print data.head()
df=pd.read_csv('/data/salomonis2/NCI-R01/TCGA-BREAST-CANCER/Anukana/patientID_tcga.txt', sep='\t', index_col=0)
#print df.head()
my_list=df.columns.values.tolist() #converts column headers to list
#my_list=list(df)

data=data.filter(my_list)

#data=data[data.iloc[:,2:].isin(my_list)]
#print data.shape
data.to_csv('/data/salomonis2/NCI-R01/TCGA-BREAST-CANCER/Anukana/UO1analysis/tcga_rsem_isopct_filtered1.txt', sep='\t')
