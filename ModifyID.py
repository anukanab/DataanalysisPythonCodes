''' this scripts modifies the names in the rows of a perticular column'''

import pandas as pd
import numpy as np
import csv

results_confounding_filtered=pd.read_csv('/data/salomonis2/NCI-R01/TCGA-BREAST-CANCER/Anukana/SurvivalAnalysis/results_counfounding_75p.txt', sep='\t')
Gene_list=results_confounding_filtered.iloc[:,0]
print results_confounding_filtered.shape
#print results_confounding_filtered.head()
print Gene_list

Genes=[]
for Gene in Gene_list:
    Gene_parts=Gene.split(".")
    #print Gene_parts
    #print ID_parts;sys.exit()
    Gene_reformat=('-'.join(Gene_parts[0:4]))
    #print Gene_reformat.head()
    Genes.append(Gene_reformat)
    #print len(Genes)

###Replace full-length barcode with abbreviated sample barcode
results_confounding_filtered.iloc[:,0]=np.array(Genes)
    
results_confounding_filtered.to_csv('/data/salomonis2/NCI-R01/TCGA-BREAST-CANCER/Anukana/SurvivalAnalysis/results_counfounding_75p_renamed.txt', sep='\t')


