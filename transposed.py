''' this file transposes columns and rows for big files that excel can not open/ handle'''

import pandas as pd
import csv
data=pd.read_csv('/data/salomonis2/NCI-R01/TCGA-BREAST-CANCER/Anukana/SurvivalAnalysis/PSI_EventAnnotation-75p_modified.txt', sep='\t') 
df=pd.DataFrame(data)
data_CNA_transposed=pd.DataFrame.transpose(df)
data_CNA_transposed.to_csv('/data/salomonis2/NCI-R01/TCGA-BREAST-CANCER/Anukana/SurvivalAnalysis/PSI_EventAnnotation-75p_modified_transposed.txt', sep='\t', header=None)
print("transposed file created")