'''this file transposes rows and columns of a csv file'''

import csv
from itertools import izip
a = izip(*csv.reader(open("/Volumes/salomonis2/TCGA_Lung/Survival-Gene-Analysis/exp.TCGA-Lung-steady-state-filtered-names.txt", "rb"),delimiter='\t'))
csv.writer(open("/Volumes/salomonis2/TCGA_Lung/Survival-Gene-Analysis/exp.TCGA-Lung-steady-state-filtered-name-transpose.txt", "wb"),delimiter='\t').writerows(a)