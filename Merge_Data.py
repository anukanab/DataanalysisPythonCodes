import pandas as pd
from pandas import DataFrame
import os
import sys
import numpy as np
import string
import csv
#import numpy as np

def getTumorSamplesAbbreviatedIDs(file_1, file_2):
    matrix_df1=pd.read_csv(file_1, sep='\t')
    #matrix_df1=pd.read_csv(file_1, sep='\t', low_memory=False, header=None)
    #matrix_transposed=matrix_df.transpose()
    #TCGA_sample_dictionary=matrix_transposed.to_dict()
    #healthy_samples=matrix_transposed[matrix_transposed[0].str.contains("-11")]
    #healthy_samples2=matrix_transposed[matrix_transposed[0].str.contains("-11")]
    
    matrix_df1_filtered=matrix_df1[matrix_df1.iloc[:,0].str.endswith("-01")]
    TCGA_list=matrix_df1_filtered.iloc[:,0]
    print TCGA_list
    TCGA_final=[]
    for ID in TCGA_list:
        ID_parts=ID.split("-")
        #print ID_parts;sys.exit()
        TCGA_reformat=('.'.join(ID_parts[0:3]))
            #print TCGA_reformat
        TCGA_final.append(TCGA_reformat)
    
    ###Replace full-length barcode with abbreviated sample barcode
    matrix_df1_filtered.iloc[:,0]=np.array(TCGA_final)
    print matrix_df1_filtered.shape
    matrix_df2=pd.read_csv(file_2, sep='\t')
    print matrix_df2.shape
    merged_df=pd.concat([matrix_df1_filtered,matrix_df2], axis=1,  join='inner')
    print merged_df.shape
    #= pd.merge(splicing_factors_mutations,splicing_factors_CNV,on=['Tumor_Sample_Barcode','Hugo_Symbol', 'Group'])
    path=os.path.split(file_1)[0]
    filename='MergedFile'
    suffix='.txt'
    output_file=os.path.join(path,filename + suffix)
    
    merged_df.to_csv(output_file, sep='\t', index=False)
    
    print "Done"
    
    return TCGA_final

if __name__ == '__main__':
    File_1= '/data/salomonis2/NCI-R01/TCGA-BREAST-CANCER/Anukana/UO1analysis/survival/tcga_rsem_isopct_filtered-transposedp.txt'
    File_2='/data/salomonis2/NCI-R01/TCGA-BREAST-CANCER/Anukana/UO1analysis/survival/clinicaldata.txt'
    OutputFile='/data/salomonis2/NCI-R01/TCGA-BREAST-CANCER/Anukana/UO1analysis/survival/MergedFile.txt'
    TCGA_samples_short_barcode = getTumorSamplesAbbreviatedIDs(File_1, File_2)