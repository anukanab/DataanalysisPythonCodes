import pandas as pd
import sys
import csv

df=pd.read_csv('/data/salomonis2/NCI-R01/TCGA-BREAST-CANCER/Anukana/CNA_data_del_amp.txt', sep='\t')
#print df
df1=df['level_0']
#print df1.head()
df['gene_id']= df['Hugo_Symbol'] + ':' + df['CNV']
print df.head()
df2=df[['level_0', 'gene_id', 'gene_id']]
print df2.head()
df2.to_csv('/data/salomonis2/NCI-R01/TCGA-BREAST-CANCER/Anukana/CNA_data_del_amp_ref.txt', sep='\t')
df3= df2.drop(['TCGA-A8-A09C-01', 'TCGA-AC-A5EI-01', 'TCGA-AR-A0U1-01', 'TCGA-BH-A0BS-01', 'TCGA-C8-A9FZ-01', 'TCGA-E2-A14S-01', 'TCGA-LD-A7W5-01'])
df3.to_csv('/data/salomonis2/NCI-R01/TCGA-BREAST-CANCER/Anukana/CNA_data_del_amp_ref_filtered.txt', sep='\t')