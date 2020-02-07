import pandas as pd
from pandas import DataFrame
import os
import sys
import numpy
import string
from collections import OrderedDict
import csv
import numpy as np


def getSampleClusterDictionary(SVMfile):
    Subtypes = open(SVMfile,'r')
    Sample_cluster=[]
    cluster_sample_dictionary = OrderedDict()
    firstRow=True
    TCGA_IDs=[]
    TCGA_units=[]
    TCGA_sample_barcode=[]
    cluster_name=[]
    sample_barcode_length=[]
    #list=[]
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
            sample_barcode_length=len(TCGA_units)
            
            cluster_info=data[1:]
            values=map(int, cluster_info)
            Sample_cluster.append(numpy.array(values))


    if sample_barcode_length != "4":
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
                #sample_name_list=[]
                if cluster not in cluster_sample_dictionary:
                    """this cluster has not yet been added to this dictionary - must initialize first """
                    cluster_sample_dictionary[cluster] = [sample_name]
                        #print cluster_sample_dictionary; sys.exit()
                else:
                    """ Atleast one sample has already been added to this dictionary key already """
                    cluster_sample_dictionary[cluster].append(sample_name)
            i+=1
        ci+=1
        
    file_path=SVMfile
    print "Task Complete: Cluster Sample Dictionary Generated"
    
    return (cluster_sample_dictionary, file_path)
    
def MakeGroupsFiles(Subtypes, MutFile):
    #### Parameters 
    gene="RBM10"
    cluster="R1-V7"
    mut_all="NO"                #Make groups file based on patients who have a mutation vs. all others
    mut_in_cluster="NO"         #Make groups file based on patients who have specified mutation in cluster vs. all others AND patients who don't have mutation in cluster vs. all others
    mut_not_in_cluster="YES"    #Make groups file based on patients who do NOT have specified mutation in cluster vs. all others
    cluster_all="NO"            #Make groups file based on patients who are in specified cluster
    

    clustersSampleDictionary, file_path=getSampleClusterDictionary(Subtypes)
    cluster_IDs=clustersSampleDictionary[cluster]
   
    all_TCGA_patients=[]
    gene_mut_all=[]
    gene_mut_cluster=[]
    non_mut_cluster=[]
    mut_data = open(MutFile,'r')
    for line in mut_data:
        line = line.strip()
        data= string.split(line,'\t')
        TCGA_ID=data[0]
        mutation=data[1]
        all_TCGA_patients.append(TCGA_ID)
        if TCGA_ID in cluster_IDs:
            if gene in data[1]:
                #gene_mut_all.append(TCGA_ID)
                gene_mut_cluster.append(TCGA_ID)
            else:
                non_mut_cluster.append(TCGA_ID)
        if mutation == gene:
            gene_mut_all.append(TCGA_ID)
    TCGA_IDs=set(all_TCGA_patients)
    folder_path=os.path.split(file_path)[0]
    del folder_path.split('/')[2]
    path=os.path.join(folder_path, "Splicing_Events_Subgroups")
    #path=os.path.join(path, "GeneExpression_Subgroups")
    suffix='.txt'
    
    comps_list=[["1","2"]]
    comps_list_mut=[["1","4"]]
    comps_list_wt=[["2","4"]]
    comps_list_mut_other=[["3","4"]]
    comps_list_cluster=[["1","3"]]

    comps_list=pd.DataFrame(comps_list)
    comps1=pd.DataFrame(comps_list_mut)
    comps2=pd.DataFrame(comps_list_wt)
    comps3=pd.DataFrame(comps_list_mut_other)
    comps=pd.DataFrame(comps_list_cluster)
    
    filename_groups= "groups"
    filename_comps= "comps"
    
   
    group1=[] #Breakdown of R1-V7:RBM10mut, R1-V7:RMB10wt, and RBM10mut not in R1-V7
    group2=[] #R1-V7 and RBM10mut not in R1-V7
    group3=[] #R1-V7: RBM10mut vs. RBM10wt
    
    for ID in TCGA_IDs:
        if ID in cluster_IDs:
            group2.append("".join(ID+":"+"1"+":"+ cluster))
            if ID in gene_mut_cluster:
                group1.append("".join(ID+":"+"1"+":"+ cluster+"_"+gene+"_mut"))
                group3.append("".join(ID+":"+"1"+":"+ cluster+"_"+gene+"_mut"))
            else:
                group1.append("".join(ID+":"+"2"+":"+cluster+"_"+gene+"_wt"))
                group3.append("".join(ID+":"+"2"+":"+ cluster+"_"+gene+"_wt"))
        else:
            if ID in gene_mut_all:
                group1.append("".join(ID+":"+"3"+":"+gene+"_mut_other"))
                group2.append("".join(ID+":"+"2"+":"+gene+"_mut_other"))
            else:
                group1.append("".join(ID+":"+"4"+":"+"NULL"))
                group2.append("".join(ID+":"+"3"+":"+"NULL"))
                
    groups1=pd.DataFrame(group1)
    output_groups1=groups1.iloc[:,0].str.split(":", expand =True)
    sorted_dataframe1=output_groups1.sort_values(output_groups1.columns[1], ascending=True)

    groups2=pd.DataFrame(group2)
    output_groups2=groups2.iloc[:,0].str.split(":", expand =True)
    sorted_dataframe2=output_groups2.sort_values(by=output_groups2.columns[1])
    
    groups3=pd.DataFrame(group3)
    output_groups3=groups3.iloc[:,0].str.split(":", expand =True)
    sorted_dataframe3=output_groups3.sort_values(by=output_groups3.columns[1])
    
    #Output_dir=os.path.join(path, gene)
    Output_dir1=os.path.join(path, cluster+"_"+gene+"_mut")
    if os.path.exists(Output_dir1) == False:
            os.mkdir(Output_dir1)
    file_Out_groups1=os.path.join(Output_dir1, filename_groups + suffix)
    sorted_dataframe1.to_csv(file_Out_groups1, sep='\t', header=None, index=False)
    file_Out_comps1=os.path.join(Output_dir1, filename_comps + suffix)
    comps1.to_csv(file_Out_comps1, sep='\t', header=None, index=False)
    
    Output_dir2=os.path.join(path, cluster+"_"+gene+"_wt")
    if os.path.exists(Output_dir2) == False:
            os.mkdir(Output_dir2)
    file_Out_groups2=os.path.join(Output_dir2, filename_groups + suffix)
    sorted_dataframe1.to_csv(file_Out_groups2, sep='\t', header=None, index=False)
    file_Out_comps2=os.path.join(Output_dir2, filename_comps + suffix)
    comps2.to_csv(file_Out_comps2, sep='\t', header=None, index=False)
    
    Output_dir3=os.path.join(path, gene+"_mut_other")
    if os.path.exists(Output_dir3) == False:
            os.mkdir(Output_dir3)
    file_Out_groups3=os.path.join(Output_dir3, filename_groups + suffix)
    sorted_dataframe1.to_csv(file_Out_groups3, sep='\t', header=None, index=False)
    file_Out_comps3=os.path.join(Output_dir3, filename_comps + suffix)
    comps3.to_csv(file_Out_comps3, sep='\t', header=None, index=False)
    
    
    #Output_dir=os.path.join(path, gene)
    Output_dir4=os.path.join(path, cluster)
    if os.path.exists(Output_dir4) == False:
            os.mkdir(Output_dir4)
    file_Out_groups4=os.path.join(Output_dir4, filename_groups + suffix)
    sorted_dataframe2.to_csv(file_Out_groups4, sep='\t', header=None, index=False)
    file_Out_comps4=os.path.join(Output_dir4, filename_comps + suffix)
    comps.to_csv(file_Out_comps4, sep='\t', header=None, index=False)
    
    Output_dir5=os.path.join(path, cluster+"_"+gene+"mut_vs_"+gene+"wt")
    if os.path.exists(Output_dir5) == False:
            os.mkdir(Output_dir5)
    file_Out_groups=os.path.join(Output_dir5, filename_groups + suffix)
    sorted_dataframe3.to_csv(file_Out_groups, sep='\t', header=None, index=False)
    file_Out_comps5=os.path.join(Output_dir5, filename_comps + suffix)
    comps_list.to_csv(file_Out_comps5, sep='\t', header=None, index=False)
    
    print "Meta-data Analysis Files Created"

if __name__ == '__main__':
    #Enrichment_File= '/Users/croz9k/Desktop/Test_Code/MergedResult_SR.txt'
    #MutationFile= '/Users/croz9k/Desktop/Test_Code/TCGA.LUAD.mutect.0458c57f-316c-4a7c-9294-ccd11c97c2f9.DR-10.0.somatic.txt'
    Subtypes='/Volumes/salomonis2/NCI-R01/TCGA_Audrey/clusters_bams/Splicing_Analysis/SVMOutputs/MergedResult.txt'
    MutationFile= '/Volumes/salomonis2/NCI-R01/TCGA_Audrey/clusters_bams/Splicing_Analysis/SVMOutputs/LUAD-PANCANCER-IDsFiltered-Mutation-cBio-Xena-FireBrowse.txt'
    sig_mutations= MakeGroupsFiles(Subtypes, MutationFile)
