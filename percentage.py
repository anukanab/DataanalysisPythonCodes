import pandas as pd
import numpy as np
df=pd.read_csv('/data/salomonis2/NCI-R01/TCGA-BREAST-CANCER/Anukana/UO1analysis/survival/tcga_rsem_isopct_filtered-transposed.txt', sep='\t')


# create a function called divides100
def divides100(x):
    # that, if x is a string,
    if type(x) is str or x==0:
        # just returns it untouched
        return x
    # but, if not, return it divided by 100
    elif x:
        return x / 100
    # and leave everything else
    else:
        return
        
df=df.applymap(divides100)

print df.head()

df.to_csv('/data/salomonis2/NCI-R01/TCGA-BREAST-CANCER/Anukana/UO1analysis/survival/tcga_rsem_isopct_filtered-transposedp.txt', sep='\t')