import pandas as pd
from pandas import DataFrame
import os
import fnmatch
import sys
import numpy
import string
import numpy as np
import matplotlib.pyplot as plt
from collections import OrderedDict


def GetCounts(File):
    matrix_df=pd.read_csv(File, sep='\t')
    lines_list=matrix_df.iloc[:,0].tolist()
    counts=len(lines_list)
    
    return counts

def GeneExpression_MetaData(SVMFile, folder_GE, folder_PSI):
    ''' Step 1. Get Cluster names '''
    matrix_df=pd.read_csv(SVMFile, sep='\t', low_memory=False)
    #print matrix_df.head()
    clusters=matrix_df.columns[1:].tolist()
    print "Clusters:"
    print clusters
    
    '''Step 2. Create Empty dictionary For clusters and # of genes per clusters'''
    MetaData_Counts_GE={}
    
    '''Step 3. Read through gene expression metadata analysis files and get gene counts per cluster'''
    
    #logfold_setting='DEGs-LogFold_0.26_adjp'
    #logfold_setting='DEGs-LogFold_1.0_adjp'
    logfold_setting='DEGs-LogFold_0.58_adjp'
    dPSI_setting='Events-dPSI_0.2_adjp'
    print ''.join('logfold setting: '+logfold_setting)
    print ''.join('dPSI setting: '+dPSI_setting)
    
    print "Processing Gene Expression Metadata Files"
    for subdir, dirs, files in os.walk(folder_GE):
        for folder in dirs:
            if folder in clusters:
                print folder
                cluster=folder
                cluster_folder=os.path.join(folder_GE, folder)
                #print cluster_folder
                for path, dirs, files in os.walk(cluster_folder):
                    for folder in dirs:
                        for folder in fnmatch.filter(dirs, logfold_setting):
                            #print folder
                            expression_folder=os.path.join(path, folder)
                            #print expression_folder
                            for path, dirs, files in os.walk(expression_folder):
                                for name in files:
                                    '''Find file that starts with 'GE.' '''
                                    for name in fnmatch.filter(files, 'GE.*'):
                                        #print name; sys.exit()
                                        expression_file=os.path.join(path, name)
                                        #print expression_file
                                        ''' Step 4. Use CountGenes function to get the # of cluster-defined genes. Save to Dictionary.'''
                                        genes_count=GetCounts(expression_file)
                                        MetaData_Counts_GE[cluster]=genes_count
                                    
    print MetaData_Counts_GE
    
    '''Step 4. Repeat Step 2 and Step 3 to get splicing events per cluster'''
    print "Processing Splicing Metadata Files"
    
    MetaData_Counts_PSI={}
    for path, dirs, files in os.walk(folder_PSI):
        #print folder_PSI
        for folder in dirs:
            for folder in fnmatch.filter(dirs, dPSI_setting):
                dPSI_folder=os.path.join(path, folder)
                #print dPSI_folder
                for path, dirs, files in os.walk(dPSI_folder):
                    for name in files:
                        for name in fnmatch.filter(files, 'PSI.*'):
                            #print name
                            if 'Kmeans'in name:
                                PSI_cluster= name.split('_')[0]
                                Kmeans=name.split('_')[1]
                                Kmeans_cluster="".join(PSI_cluster+"_"+Kmeans)
                                cluster_fullname=Kmeans_cluster.split('.')[1]
                            else:
                                PSI_cluster= name.split('_')[0]
                                #print PSI_cluster
                                cluster= PSI_cluster.split('.')[1]
                                if "R1" in cluster:
                                    cluster_fullname=cluster
                                else:
                                    cluster_fullname="".join("R1-"+cluster)
                                #print cluster_fullname
                            if cluster_fullname in clusters:
                                print cluster_fullname
                                PSI_file=os.path.join(path, name)
                                #print PSI_file
                                events_count=GetCounts(PSI_file)
                                MetaData_Counts_PSI[cluster_fullname]=events_count
                                    
    print MetaData_Counts_PSI
    
    '''Step 5. Sort one of the dictionaries by decending order. This will make the graph look nicer '''
    MetaData_Counts_GE_descending = sorted(MetaData_Counts_GE, key=MetaData_Counts_GE.get, reverse=False)
    
    GE_counts=[]
    PSI_counts=[]
    for i in MetaData_Counts_GE_descending:
       GE_counts.append(MetaData_Counts_GE[i])
       PSI_counts.append(MetaData_Counts_PSI[i])
    #''' Check if there are missing clusters'''
    #events_List1_unique=list(set(MetaData_Counts_GE.keys()) - set(MetaData_Counts_PSI.keys()))
    #print events_List1_unique; sys.exit()
    
    '''Step 6. Make barplot '''
    y_pos = np.arange(len(MetaData_Counts_PSI.keys()))
    width=0.4
    
    fig, ax = plt.subplots()
    ax.barh(y_pos, PSI_counts, width, color='blue', label='Splicing Events')
    ax.barh(y_pos + width, GE_counts, width, color='orange', label='DEGs')
    
    ax.set(yticks=y_pos + width, yticklabels=MetaData_Counts_GE_descending, ylim=[2*width - 1, len(MetaData_Counts_PSI)])
    ax.legend()
    plt.rcParams["figure.figsize"] = (20,3)
    plt.show()
    
    #'''Plot only one variable'''
    #y_pos = np.arange(len(MetaData_Counts_PSI))
    #plt.barh(y_pos, MetaData_Counts_PSI.values())
    #plt.yticks(y_pos, MetaData_Counts_PSI.keys())
    
    
    
if __name__ == '__main__':
    import getopt
    if len(sys.argv[1:])<=2:  ### Indicates that there are insufficient number of command-line arguments
        print "Misssing a needed input file"
    else:
        options, remainder = getopt.getopt(sys.argv[1:],'', ['gene=','PSI=', 'i='])
        #print sys.argv[1:]
        for opt, arg in options:
            if opt == '--gene': folder_GE=arg
            elif opt == '--PSI': folder_PSI=arg
            elif opt == '--i': SVMFile=arg
        GeneExpression_MetaData(SVMFile, folder_GE, folder_PSI)         
