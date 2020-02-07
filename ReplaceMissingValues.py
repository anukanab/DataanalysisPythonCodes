''' this file fills in the missing values with NaN and replaces them with average of the column for big files that excel can not open/ handle'''

import pandas as pd
import csv
import numpy as np
data=pd.read_csv('/data/salomonis2/NCI-R01/TCGA-BREAST-CANCER/Anukana/MetaDataAnalysis/methylation/BRCA_HumanMethylation450.betaValue3.txt', sep='\t') 
df=pd.DataFrame(data)
df=df.replace(r'^\s*$', np.nan)
df=df.fillna(df.mean())

df.to_csv('/data/salomonis2/NCI-R01/TCGA-BREAST-CANCER/Anukana/MetaDataAnalysis/methylation/BRCA_HumanMethylation450.betaValue3.txt', sep='\t', index=False)
print("replaced file created")