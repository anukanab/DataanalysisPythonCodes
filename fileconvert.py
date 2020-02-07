import pandas as pd
import csv
data=pd.read_csv('/data/salomonis2/NCI-R01/TCGA-BREAST-CANCER/Anukana/UO1analysis/tcga_rsem_isopct', sep='\t')
data.to_csv('/data/salomonis2/NCI-R01/TCGA-BREAST-CANCER/Anukana/UO1analysis/tcga_rsem_isopct.txt', sep='\t')