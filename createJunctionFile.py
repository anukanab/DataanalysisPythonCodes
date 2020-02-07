import pandas as pd
from pandas import DataFrame
import os
import os.path
import sys
import numpy as np
import string
from collections import OrderedDict
import csv

def getcolumns(file1, file2):
 
 
 
 
 
 
 Dict_FL.HCC{}
 Dict_FL.HME1={}
 Dict_FL.MCF7={}
 Dict_FL.MDA231={}
 Dict_FL.ZR75={}

 #open file 
    Subtypes = open(file1,'r')  
    for line in Subtypes:
        line = line.strip()
        data= string.split(line,'\t')
    file_info=line.split
    isoform=file_info[0]
    FL.HCC_counts=file_info[1]
    FL.HME1_counts=file_info[2]
    FL.MCF7_counts=file_info[3]
    FL.MDA231_counts=file_info[4]
    Fl.ZR75_counts=file_info[5]
    Dict_FL.HCC[isoform]=FL.HCC_counts
    Dict_FL.HME1[isoform]=FL.HME1_counts
    Dict_FL.MCF7[isoform]=FL.MCF7_counts
    Dict_FL.MDA231[isoform]=FL.MDA231_counts
    Dict_FL.ZR75[isoform]=Fl.ZR75_counts
 
 file2=pd.read_csv(file2, sep='\t')  
 d=pd.file2()
 
 
  for isoform in Dict_FL.HCC:
   if isoform in d:
    d[FL.HCC]=Fl.HCC_counts
    
   
 

 
    
    