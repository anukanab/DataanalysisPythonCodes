Usage = """
    pairwiseAnalysis.py
    Goal: Output a matrix of the ratio of TCGA IDs that match each cluster between two different SMVOutput files. 
    USAGE: python pairwiseAnalysis.py
    
    INPUTS:
    Subtypes_1 - Subtypes_1 is the MergedResult.txt file, within SMVOuputs folder, generated after running Oncosplice. Subtypes_1 is a sparse matrix where TCGA Sample ID's are listed in col 1 and the associated clusters as row 1. 
    Subtypes_2- Subtypes_2 is the MergedResult.txt file, within SMVOuputs folder, generated after running Oncosplice. Subtypes_2 is a sparse matrix where TCGA Sample ID's are listed in col 1 and the associated clusters as row 1. 
    
    OUTPUTS:
    Heatmap - png file containing a 2 heatmaps of the overlap between TCGA IDs in each cluster between Subtypes_1 and Subtypes_2. Heatmap annotations are the ratio of overlap.
             Ratios from the left heatmap are calculated as follows: (# of overlapping TCGA samples between two clusters)/(# of TCGA samples in cluster from Subtypes1).
             Ratios from the right heatmap are calculated as follows: (# of overlapping TCGA samples between two clusters)/(# of TCGA samples in cluster from Subtypes2)
        ***NOTE- Axis title denoting whether clusters are from Subtypes_1 or Subtypes_2 not included in output files.
        OutputFile_1 - Txt file containg cluster similarities ratios in a matrix form. Ratios calculated as follows: (# of overlapping TCGA samples between two clusters)/(# of TCGA samples in cluster from Subtypes1).
        OutputFile_2 - Txt file containg cluster similarities ratios in a matrix form. Ratios calculated as follows: (# of overlapping TCGA samples between two clusters)/(# of TCGA samples in cluster from Subtypes2).
        OutputFile_3 - Excel file containing table within OutputFile_1 and OutputFile_2.
    
    Output files will be saved in the same file path as Subtypes_1 file
    """
Usage_Warning = """ TCGA IDs that are less than 4 units long (e.g. ID less than: TCGA-XX-XXXX-XX) cannot differentiate healthy vs. normal samples. To use anyway, uncomment line 94, 100, and 143."""

import pandas as pd
from pandas import DataFrame
import os
import os.path
import sys
import numpy
import string
from collections import OrderedDict
import csv
import shutil
import seaborn as sns
import matplotlib.pyplot as plt
import xlsxwriter
      
""" Standardize TCGA ID format and identify the overlapping TCGA Sample IDs per file """
def getTCGASampleNamesLength(file1, file2):
    Subtypes1 = open(file1,'r')
    ID_length1=[]
    TCGA_units1=[]
    next(Subtypes1)              ###skip first line
    for line in Subtypes1:
        file_info=line.split()
        TCGA_samples=file_info[0]
        TCGA_units1.append(TCGA_samples.split("-"))
        for ID in TCGA_units1: 
            ID_length1= len(ID)
    Subtypes2 = open(file2,'r')
    TCGA_units2=[]
    next(Subtypes2)              ###skip first line
    for line in Subtypes2:
        file2_info=line.split()
        TCGA_samples2=file2_info[0]
        TCGA_units2.append(TCGA_samples2.split("-"))
        ID_length2= len(TCGA_units2)
        for ID in TCGA_units2: 
            ID_length2= len(ID)
        """ If TCGA IDs are in different formats (e.g. full-length vs. patient ID only), TCGA ID format will be standardized for comparisons. """
        if ID_length1 != ID_length2:          
            ID_length= min(ID_length1, ID_length2)
            if ID_length<4:
                print Usage_Warning
        else: 
            ID_length="SAME"      
    return ID_length

def getTCGAidentifyers(file1, file2):
    Subtypes1 = open(file1,'r')
    TCGA_ID1=[]
    TCGA_samples1=[]
    TCGA_units1=[]
    next(Subtypes1)              ###skip first line
    for line in Subtypes1:
        file_info=line.split()
        TCGA_samples1.append(file_info[0])
        TCGA_info=file_info[0]
        TCGA_units1.append(TCGA_info.split("-"))
    Subtypes2 = open(file2,'r')
    TCGA_ID2=[]
    TCGA_units2=[]
    TCGA_samples2=[]
    next(Subtypes2)              ###skip first line
    for line in Subtypes2:
        file_info=line.split()
        TCGA_samples2.append(file_info[0])
        TCGA_info=file_info[0]
        TCGA_units2.append(TCGA_info.split("-"))
    
    """ If TCGA IDs are in different formats (e.g. full-length vs. patient ID only), standardize format for comparisons."""
    if TCGASampleLength != "SAME":          
        for ID in TCGA_units1:
            TCGA_reformat=('-'.join(ID[0:4]))
            TCGA_final1=TCGA_reformat.strip("A").strip("B")
            #TCGA_final=('-'.join(ID[0:TCGASampleLength])) ###only use if TCGA sample names are cutoff before specifying tumor vs normal sample (e.g. tumor: -01A, healthy: -11A)
            TCGA_ID1.append(TCGA_final1)
        for ID in TCGA_units2:
            TCGA_reformat=('-'.join(ID[0:4]))
            TCGA_final2=TCGA_reformat.strip("A").strip("B")
            TCGA_ID2.append(TCGA_final2)
            #TCGA_ID2.append('-'.join(ID[0:TCGASampleLength])) ###only use if TCGA sample names are cutoff before specifying tumor vs normal sample (e.g. tumor: -01A, healthy: -11A)
    else: 
            TCGA_ID1=TCGA_samples1
            TCGA_ID2=TCGA_samples2
    TCGA_overlap=list(set(TCGA_ID1) & set(TCGA_ID2))
    return TCGA_overlap
        
def getSampleClusterDictionary(filename):
    Subtypes=open(filename, 'r')
    Sample_cluster=[]
    cluster_sample_dictionary = OrderedDict()
    firstRow=True
    TCGA_IDs=[]
    TCGA_units=[]
    TCGA_sample_barcode=[]
    cluster_name=[]
    list=[]
    for line in Subtypes:
        line = line.strip()
        data= string.split(line,'\t')
        """ Get header of file """
        if firstRow:
            clusters=data[1:]
            for cluster in clusters: ### Make cluster names neater
                if "vs" in cluster:
                    cluster_units=cluster.split("_")
                    cluster_name.append(cluster_units[0])
                else: 
                    cluster_name.append(cluster)
            firstRow = False
        else:
            """ get cluster info for each TCGA sample and put in a list"""
            TCGA_sample_barcode.append(data[0])
            TCGA_samples=data[0]
            TCGA_units.append(TCGA_samples.split("-"))
            
            cluster_info=data[1:]
            values=map(int, cluster_info)
            Sample_cluster.append(numpy.array(values))

    if TCGASampleLength != "SAME":
        for unit in TCGA_units:
            TCGA_reformat=('-'.join(unit[0:4]))
            #TCGA_IDs.append('-'.join(ID[0:TCGASampleLength])) ###only use if TCGA sample names are cutoff before specifying tumor vs normal sample (e.g. tumor: -01A, healthy: -11A)
            TCGA_final=TCGA_reformat.strip("A").strip("B")
            TCGA_IDs.append(TCGA_final)
    else:
        TCGA_IDs=TCGA_sample_barcode

    Sample_cluster = numpy.array(Sample_cluster) ### make a numpy object from the list of lists
    Sample_cluster1 = numpy.transpose(Sample_cluster) ###transpose list
    Sample_cluster2 = zip(*Sample_cluster)

    ci=0
    for cluster in cluster_name:
        i=0
        for sample_call in Sample_cluster1[ci]:
            if sample_call == 1:
                sample_name = TCGA_IDs[i]
                #print TCGA_IDs; sys.exit()
                sample_name_list=[]
                if sample_name in TCGA_Overlap:
                    if cluster not in cluster_sample_dictionary:
                        """this cluster has not yet been added to this dictionary - must initialize first """
                        cluster_sample_dictionary[cluster] = [sample_name]
                        #print cluster_sample_dictionary; sys.exit()
                    else:
                        """ Atleast one sample has already been added to this dictionary key already """
                        cluster_sample_dictionary[cluster].append(sample_name)
            i+=1
        ci+=1
    
    cluster_num= len(cluster_name)
    return cluster_sample_dictionary, cluster_num

### def pairwiseAnalysisMatrix(dictionary1, dictionary2, file1, file2, file3):
def pairwiseAnalysisMatrix(dictionary1, dictionary2, file1,file2): 
    pairwise_analysis_ratios1 = OrderedDict() # denominator for ratios == number of TCGA samples per cluster in dictionary1 
    pairwise_analysis_ratios2 = OrderedDict() # denominator for ratios == number of TCGA samples per cluster in dictionary2
    column_title=[]
    row_title=[]
    for cluster1 in dictionary1:
        for cluster2 in dictionary2:
            if cluster2 not in column_title:
                column_title.append(cluster2)
            if cluster1 not in row_title:
                """ this cluster has not yet been added to this dictionary - must initialize first """
                row_title.append(cluster1)
                ### ratios1 = (# of overlapping TCGA samples between two clusters)/(# of TCGA samples in cluster from Subtypes1)
                ratio1=len(list(set(dictionary1[cluster1]) & set(dictionary2[cluster2])))/(len(list(set(dictionary1[cluster1])))* 1.000)
                pairwise_analysis_ratios1[cluster1] = [round(ratio1,2)]
                ### ratios2 = (# of overlapping TCGA samples between two clusters)/(# of TCGA samples in cluster from Subtypes2)
                ratio2=len(list(set(dictionary1[cluster1]) & set(dictionary2[cluster2])))/(len(list(set(dictionary2[cluster2])))* 1.000)
                pairwise_analysis_ratios2[cluster1] = [round(ratio2,2)]
                """ calculate % overlap between clusters - Optional """
                #pairwise_analysis_ratios[cluster1] = [len(list(set(dictionary1[cluster1]) & set(dictionary2[cluster2])))/(len(list(set(dictionary1[cluster1] + dictionary2[cluster2])))* 1.000)]

            else:
                """ Atleast one sample has already been added to this dictionary key already """
                ratio1=len(list(set(dictionary1[cluster1]) & set(dictionary2[cluster2])))/(len(list(set(dictionary1[cluster1])))* 1.000)
                pairwise_analysis_ratios1[cluster1].append(round(ratio1,2))
                ratio2=len(list(set(dictionary1[cluster1]) & set(dictionary2[cluster2])))/(len(list(set(dictionary2[cluster2])))* 1.000)
                pairwise_analysis_ratios2[cluster1].append(round(ratio2,2))
                """ calculate % overlap between clusters - Optional """
                #pairwise_analysis_ratios[cluster1].append(len(list(set(dictionary1[cluster1]) & set(dictionary2[cluster2])))/(len(list(set(dictionary1[cluster1] + dictionary2[cluster2])))* 1.000))
    
    """ write to Outputfile_1""" # can comment out this part and directly move to write to output file 3
    ### Create file path
    file_path=os.path.split(file1)[0]
    file2_path=file2.split('/')[-3]
    Output_dir= os.path.join(file_path,'PairwiseAnalysis_'+file2_path)
    filename1= "ComparisonMatrix_1"
    suffix = '.txt'
    if os.path.exists(Output_dir) == False:
        Output_path= os.mkdir(Output_dir)
        file_out1=os.path.join(Output_dir, filename1 + suffix)
    else:
        file_out1=os.path.join(Output_dir, filename1 + suffix)

    ### Save dataframe as matrix and write 
    matrix1_df=pd.DataFrame.from_dict(pairwise_analysis_ratios1)
    matrix1_transposed=matrix1_df.transpose()
    matrix1_transposed.columns = column_title
    matrix1_transposed.to_csv(file_out1, sep='\t', index=True) 
    
    """ write to Outputfile_2""" # can comment out this part and directly move to write to output file 3
    filename2= "ComparisonMatrix_2"
    file_out2=os.path.join(Output_dir, filename2 +suffix)
    matrix2_df=pd.DataFrame.from_dict(pairwise_analysis_ratios2)
    matrix2_transposed=matrix2_df.transpose()
    matrix2_transposed.columns = column_title
    matrix2_transposed.to_csv(file_out2, sep='\t', index=True)
    
    """ write to Outputfile_3 """
    suffix='.xlsx'
    filename3= "ComparisonMatrix_1-2"
    file_out3=os.path.join(Output_dir, filename3 +suffix)
    writer = pd.ExcelWriter(file_out3,engine='xlsxwriter')
    df_list=[matrix1_transposed,matrix2_transposed]
    row = 0
    ### OutputMatrix1 appears first, followed by OutputMatrix2 underneath
    for dataframe in df_list:
        dataframe.to_excel(writer,sheet_name="OutputMatrix_1-2",startrow=row , startcol=0)   
        row = row + len(dataframe.index) + 1 + 1
    writer.save()
    
    return (matrix1_transposed, matrix2_transposed, Output_dir)

def makeHeatmap (matrix1, matrix2, file1, file2):
    """ Make Axis titles for heatmap. Generated from part of input file path"""
    path1=file1.split("/")
    axis_y=path1[-3]
    path2=file2.split("/")
    axis_x=path2[-3]
    
    if clusters_num2 <= 12:
        print "TRUE"
        """ Make font smaller so figure looks better """
        sns.set(font_scale=.65)
        
        """ Make two heatmaps side by side seperated with .05 of white space""" 
        fig, (ax,ax2) = plt.subplots(ncols=2)
        fig.subplots_adjust(wspace=0.05)
        
        """ Make first heatmap from pairwiseAnalysisMatrixOutput1"""
        sns.heatmap(pairwiseAnalysisMatrixOutput1, cmap="YlGnBu", ax=ax, cbar=False, annot=True, annot_kws={"size": 5})
        cbar1=fig.colorbar(ax.collections[0], ax=ax,location="top", use_gridspec=False, pad=0.02, shrink=.8)
        plt.ylabel(axis_y)
        cbar1.ax.tick_params(labelsize=4)
        
        """ Make second heatmap from pairwiseAnalysisMatrixOutput2""" 
        sns.heatmap(pairwiseAnalysisMatrixOutput2, cmap="YlGnBu", ax=ax2, cbar=False, annot=True, annot_kws={"size": 5})
        cbar2=fig.colorbar(ax2.collections[0], ax=ax2,location="top", use_gridspec=False, pad=0.02, shrink=.8)
        cbar2.ax.tick_params(labelsize=4)
        ax2.yaxis.tick_right() #move y-axis lables to right hand side
        ax2.set_yticklabels(ax2.get_yticklabels(),rotation=0)
        ax2.set_xticklabels(ax2.get_xticklabels(),rotation=90)
        ax2.tick_params(right=False) #remove tick parts on heatmap
        
        """ Add figure title """ 
        plt.suptitle(''.join(axis_x + "_vs_"+ axis_y), fontsize=10, fontweight=0, color='black', style='italic', y=.89)
        
        """ Add x- and y-axis titles """ 
        plt.text(0.5, 13.5, axis_x, ha='center', va='center')
        ax.set_ylabel(axis_y)
    else:
        """ Make font smaller so figure looks better """
        sns.set(font_scale=.3)
        
        """ Make two heatmaps side by side seperated with .05 of white space""" 
        fig, (ax,ax2) = plt.subplots(ncols=2)
        fig.subplots_adjust(wspace=0.03)
        
        """ Make first heatmap from pairwiseAnalysisMatrixOutput1"""
        sns.heatmap(pairwiseAnalysisMatrixOutput1, cmap="YlGnBu", ax=ax, cbar=False, annot=True, annot_kws={"size": 2})
        cbar1=fig.colorbar(ax.collections[0], ax=ax,location="top", use_gridspec=False, pad=0.02, shrink=.6)
        plt.ylabel(axis_y)
        cbar1.ax.tick_params(labelsize=3)
    
        """ Make second heatmap from pairwiseAnalysisMatrixOutput2""" 
        sns.heatmap(pairwiseAnalysisMatrixOutput2, cmap="YlGnBu", ax=ax2, cbar=False, annot=True, annot_kws={"size": 2})
        cbar2=fig.colorbar(ax2.collections[0], ax=ax2,location="top", use_gridspec=False, pad=0.02, shrink=.6)
        cbar2.ax.tick_params(labelsize=3)
        ax2.yaxis.tick_right() #move y-axis lables to right hand side
        ax2.set_yticklabels(ax2.get_yticklabels(),rotation=0)
        ax2.set_xticklabels(ax2.get_xticklabels(),rotation=90)
        ax2.tick_params(right=False) #remove tick parts on heatmap
        
        """ Add figure title """ 
        plt.suptitle(''.join(axis_x + "_vs_"+ axis_y), fontsize=6, fontweight=0, color='black', style='italic', y=.85)
        
        """ Add x- and y-axis titles """ 
        plt.text(0.5, 25, axis_x, ha='center', va='center')
        ax.set_ylabel(axis_y)
    
    """ Save figure """
    ### Create file path
    file_path=os.path.split(file1)[0]
    file2_path=file2.split('/')[-3]
    Output_dir= os.path.join(file_path,'PairwiseAnalysis_'+file2_path)
    filename="ComparisonMatrix"
    suffix = '.png'
    file_out=os.path.join(Output_dir, filename + suffix)
    plt.savefig(file_out, dpi=400, bbox_inches="tight")
    print "Task Complete"
    print "".join("Output files located within "+Output_dir)


if __name__ == '__main__':
    Subtypes_1='/data/salomonis2/NCI-R01/TCGA-BREAST-CANCER/Anukana/GE_ICGS/MergedResultSEBC.txt'
    Subtypes_2='/data/salomonis2/NCI-R01/TCGA-BREAST-CANCER/Anukana/GE_ICGS/GroupsGEBC_VanillaICGSCosine.txt'
    TCGASampleLength = getTCGASampleNamesLength(Subtypes_1, Subtypes_2)
    TCGA_Overlap = getTCGAidentifyers(Subtypes_1, Subtypes_2)
    clustersSampleDictionary1, clusters_num1 = getSampleClusterDictionary(Subtypes_1)
    clustersSampleDictionary2, clusters_num2 = getSampleClusterDictionary(Subtypes_2)
    pairwiseAnalysisMatrixOutput1, pairwiseAnalysisMatrixOutput2, OutputDir= pairwiseAnalysisMatrix(clustersSampleDictionary1, clustersSampleDictionary2, Subtypes_1,Subtypes_2)
    Heatmap= makeHeatmap(pairwiseAnalysisMatrixOutput1, pairwiseAnalysisMatrixOutput2, Subtypes_1, Subtypes_2)