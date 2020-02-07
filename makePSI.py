""" Goal of this script is to modify the methylation data file look like a PSI eventannotation file. We can then pass on the file as a PSI file in splice-ICGS to get clusters based on the methylation data."""

import pandas as pd
from pandas import DataFrame
import os
import os.path
import sys
import numpy
import string
from collections import OrderedDict
import csv

methylationdf=pd.read_csv('/data/salomonis2/NCI-R01/TCGA-BREAST-CANCER/Anukana/MetaDataAnalysis/methylation/BRCA_HumanMethylation450.betaValue3head.txt', sep='\t')

#add columns

methylationdf=methylationdf["Symbol","Description","AltExons","dPSI","ClusterID","UID","Coordinates","Eventannotation"]


                            
                            