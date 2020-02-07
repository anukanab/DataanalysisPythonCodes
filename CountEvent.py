''' this script counts number of events'''

import pandas as pd
import numpy as np
import sys
import csv

df=pd.read_csv('/data/salomonis2/NCI-R01/TCGA-BREAST-CANCER/Anukana/CNA/MarkerFinder091019/CNA-SF-filtered.txt', sep='\t', header=None, low_memory=False)

d=dict ()  #create an empty dictionary
my_list=df.iloc[:, 5].tolist() #convert the desired colum with events into a list
for event in my_list:   # iterate over each event in list
  if event in d:     #check if the event is already in dictionary
    d[event] = d[event] + 1   #increment count of event by one
  else:
    d[event] = 1              #add the word to the dictionary with count 1
    
#Print the contents of dictionary
for key in list(d.keys()):
    print(key, ":", d[key])
    
df=pd.DataFrame.from_dict([d])  # the sqare bracket converts it to a list of value, surpasses the ValueError: If using all scalar values, you must pass an index
df.to_csv('/data/salomonis2/NCI-R01/TCGA-BREAST-CANCER/Anukana/CNA/MarkerFinder091019/CNA-SF-filtered_EventCount.txt', sep='\t')