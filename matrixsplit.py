import pandas as pd
import numpy as np
import os

OrigMatrix = pd.read_csv('/data/salomonis2/NCI-R01/TCGA-BREAST-CANCER/Anukana/MetaDataAnalysis/CNA/PSI/data_CNA_matrix.txt', sep='\t', index_col=0)

submatrixColLength = len(OrigMatrix.columns)/8

dfs = np.split(OrigMatrix, [submatrixColLength,submatrixColLength*2,submatrixColLength*3,submatrixColLength*4,submatrixColLength*5,submatrixColLength*6,submatrixColLength*7,submatrixColLength*8], axis=1)

i = 0
t = 1
while i < 8:
   #for df in dfs:
     output_path = '/data/salomonis2/NCI-R01/TCGA-BREAST-CANCER/Anukana/MetaDataAnalysis/CNA/PSI/'
     filePre='data_CNA_matrix'
     file_suffix = '.txt'
     fileMid=t 
     File_Out=os.path.join(output_path, filePre + str(fileMid) + file_suffix)
     dfs[i].to_csv(File_Out, sep='\t', index=True)
     i=i+1
     t=t+1

