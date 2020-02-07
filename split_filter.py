''' This script splits columns by a delimeter and then filters the df based on the values of a column that matched to another list'''

import pandas as pd
import numpy as np
import sys

 #Modifies/ splits columns to extract gene names
df=pd.read_csv('/data/salomonis2/NCI-R01/TCGA-BREAST-CANCER/Anukana/CNA/MarkerFinder/CNA-data-del-amp.txt', sep='\t', header=None)
#df[[]]=df.iloc[:,1].str.split(":", expand=True)
df=df.iloc[:,1].apply(lambda x: pd.Series(str(x).split(":")))
print df.head();sys.exit()
df1.iloc[:,2:3]=df[df.iloc[:,2].str.split(":", expand=True)]
print df1.head(); sys.exit()
#df=df.drop(df.columns[2,4],inplace=True)
#header=["samples", "Genes", "Gene_Expression"]
    
df2 = pd.read_csv('/data/salomonis2/NCI-R01/TCGA-BREAST-CANCER/Anukana/CNA/MarkerFinder/SplicingFactorList.txt', sep='\t')
SF_list = df2["Gene Symbol"].tolist()  # reads the file and creates list from the appropriate column
    
print SF_list


df=df[df.iloc[:,1].isin(SF_list)]  # filters out dataset based on the list

print df.head()
print df.shape
    


df.to_csv('/data/salomonis2/NCI-R01/TCGA-BREAST-CANCER/Anukana/CNA/MarkerFinder/CNA-SF-filtered.txt', sep='\t', header=None, index=False)

#if __name__ == '__main__':
    #df='/data/salomonis2/NCI-R01/TCGA-BREAST-CANCER/Anukana/CNA/MarkerFinder/CNA-data-del-amp.txt'
    #SF_genes='/data/salomonis2/NCI-R01/TCGA-BREAST-CANCER/Anukana/\ proteins_splicing_human\ _summary.csv'
    #outputfile= '/data/salomonis2/NCI-R01/TCGA-BREAST-CANCER/Anukana/CNA/MarkerFinder/CNA-SF-filtered.txt'

