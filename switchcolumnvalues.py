import pandas as pd
import csv
data=pd.read_csv('/data/salomonis2/NCI-R01/TCGA-BREAST-CANCER/methylation/jhu-usc.edu_BRCA_HumanMethylation450.betaValue.csv', sep=',')
df=pd.read_csv('/data/salomonis2/NCI-R01/TCGA-BREAST-CANCER/methylation/HumanMethylation450_15017482_v1-2_simple.csv', sep=',')
dict=pd.Series(df.UCSC_RefGene_Name.values,index=df.Name).to_dict()
print dict
data.iloc[:,0].replace(dict, inplace=True)
data.to_csv('/data/salomonis2/NCI-R01/TCGA-BREAST-CANCER/methylation/output.txt', sep='\t', index=False)