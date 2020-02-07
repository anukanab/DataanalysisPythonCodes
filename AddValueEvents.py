import pandas as pd
import csv
import numpy as np
import sys

df=pd.read_csv('/data/salomonis2/NCI-R01/Harvard/BRC_PacBio_Seq/modifiedtxt/MDA231.txt', sep='\t', low_memory=False)
print df['genomic_start_coord.isoforms_junctions'].head()

#changing the value of one column based on th evalue of in another column 
m = df['strand.isoforms_junctions']== "-"
n = df['strand.isoforms_junctions']== "+"
df['genomic_start_coord.isoforms_junctions'] = np.where(m, df['genomic_start_coord.isoforms_junctions'].astype(int)-1, df['genomic_start_coord.isoforms_junctions'])
df['genomic_end_coord.isoforms_junctions'] = np.where(m, df['genomic_end_coord.isoforms_junctions'].astype(int)+1, df['genomic_end_coord.isoforms_junctions'])
df['genomic_start_coord.isoforms_junctions'] = np.where(n, df['genomic_start_coord.isoforms_junctions'].astype(int)-1, df['genomic_start_coord.isoforms_junctions'])
df['genomic_end_coord.isoforms_junctions'] = np.where(n, df['genomic_end_coord.isoforms_junctions'].astype(int)+1, df['genomic_end_coord.isoforms_junctions'])

print df['genomic_start_coord.isoforms_junctions'].head()

df['uid']=df['chrom.isoforms_junctions'].map(str)+":"+df['junction_number.isoforms_junctions'].map(str)+":"+df['genomic_start_coord.isoforms_junctions'].map(str)+"-"+df['genomic_end_coord.isoforms_junctions'].map(str)+df['strand.isoforms_junctions'].map(str)

df1=df.groupby('uid')['counts'].sum() #adds up all the values for all unique "uid"

df1.to_csv('/data/salomonis2/NCI-R01/Harvard/BRC_PacBio_Seq/modifiedtxt/MDA231_counts_Il.txt', sep='\t', header=True)

df1=pd.read_csv('/data/salomonis2/NCI-R01/Harvard/BRC_PacBio_Seq/modifiedtxt/MDA231_counts_Il.txt', sep='\t')

print df1.head()

new = df1.iloc[:,0].str.split(":", expand=True) #splits the column by ":" delemiter

print new.head()

df1['chrom']=new[0]

df1['name']=new[2].str[:-1] #takes out the last character

df1['notes']=new[1] + ":" + new[2]
#df1['notes']=df1['note'].str[:-1]
df1['strand']=new[2].str[-1] #only prints the last character

print df1.head()

new1=df1['name'].str.split("-", expand=True)

df1['position1']=new1[0]
df1['position2']=new1[1] 
df1['null1']='255,0,0'
df1['null2']=2
df1['exon length']='0,1'
df1['null3']='0,793'

cols=['chrom', 'position1', 'position2', 'notes', 'counts', 'strand', 'position1', 'position2', 'null1', 'null2','exon length', 'null3' ]
df1=df1[cols]
#df1.insert(8, "null", "255,0,0")
#df1.insert(9, "null2", "2")

print df1.head()

df1.to_csv('/data/salomonis2/NCI-R01/Harvard/BRC_PacBio_Seq/beds/MDA231.junctions.bed', sep='\t', header=False, index=None)

