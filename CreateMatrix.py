"""creating a binary value matrix from a dataframe of subtypes"""

import pandas as pd
from pandas import DataFrame
import os
import os.path
import sys
import numpy
import string
from collections import OrderedDict
import csv

df=pd.read_csv('/data/salomonis2/NCI-R01/Harvard/BRC_PacBio_Seq/metadataanalysis/PAM50/Filtered-SampleIDS-Basal-LUM.txt', sep='\t', header=None)

df1= pd.crosstab(df.iloc[:,0], df.iloc[:,1])

df1.to_csv('/data/salomonis2/NCI-R01/Harvard/BRC_PacBio_Seq/metadataanalysis/PAM50/groups.PAM50Matrix.txt', sep='\t')