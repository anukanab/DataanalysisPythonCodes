import pandas as pd
from pandas import DataFrame
import os
import sys
import string
import numpy as np
import csv

def GetMatchingSignatures(metaData_File):
    metaData_df = pd.read_csv(metaData_File, sep='\t', header=None)
    headers = metaData_df.iloc[0]
    new_df  = pd.DataFrame(metaData_df.values[1:], columns=headers)

    if "UID" in new_df.columns:
        data_list=new_df["UID"].tolist()                        ###This line will get the specific splicing events. This will be used to find instances where the splicing event overlaps between two groups.
        #genes=new_df["UID"].str.split(":", n = 1, expand = True)[0]     ### This line isolates gust the genes that are alternatively spliced.
        #print genes; sys.exit()
    elif "Symbol" in new_df.columns:
        filtered_df=new_df.dropna()
        #genes=filtered_df["Symbol"].tolist()
        #filtered_matrix.apply(lambda filtered_matrix: [TCGA_ID_dict[y] if y in TCGA_ID_dict else y for y in filtered_matrix])
        #df['X'] = np.where(df['Y']>=50, 'yes', 'no')
        filtered_df['LogFold']=np.where(filtered_df['LogFold']<0, 'downregulated', 'upregulated')
        #filtered_df['LogFold'].apply(lambda filtered_df['LogFold']: "downregulated" if x in  filtered_df['LogFold'] < 0 else "upregulated" for x in filtered_df['LogFold'])
        filtered_df['Gene_Expression']=filtered_df["Symbol"] +':'+filtered_df['LogFold']
        data_list=filtered_df['Gene_Expression'].tolist()
        print data_list
    
    #print data_list    drop(columns=['Entrez_Gene_Id'])
    return data_list

def process_genes_events(List, name):
    List_df=pd.DataFrame(List)
    List_df=List_df.iloc[:,0].str.split("|", n=1, expand=True)
    List_df=List_df.iloc[:,0].str.split(":", n=3, expand=True)
    if len(List_df.columns) ==3:
        List_df=List_df.drop(columns=[1])
        List_df.columns=["Genes", "Splicing_Event"]
    else:
        header=["Genes", "Gene_Expression"]
    
    List_df["Group"]=name
    
    return List_df
        
    #print List_df.head(); sys.exit()

def KeepSignificantGenes(List, sig_genes):
    genes_filtered=[]
    for value in List:
        #print value
        if ":" in value:
            gene=value.split(":")[0]
        else:
            gene=value
            #print gene
        for item in sig_genes:
            if item == gene:
                genes_filtered.append(gene)
            
    genes=list(set(genes_filtered)) ###Remove duplicates. Needed if parsing splicing events
    return genes

def AnnotateGenes(Df, sf_genes, oncogenes, TFs):
    Df['Splicing_Factors']=np.where(Df.iloc[:,0].isin(sf_genes), "1", "0")
    #Df.apply(lambda Df: ["1" if y in sf_genes else "0"])
    Df['Cancer_Drivers']=np.where(Df.iloc[:,0].isin(oncogenes), "1", "0")
    #Df.apply(lambda Df: ["1" if y in oncogenes else "0"])
    Df['Transcription_Factors']=np.where(Df.iloc[:,0].isin(TFs), "1", "0")
    #Df.apply(lambda Df: ["1" if y in TFs else "0"])
    #filtered_matrix.apply(lambda filtered_matrix: [TCGA_ID_dict[y] if y in TCGA_ID_dict else y for y in filtered_matrix])
    
    return Df
    
    
def CompareOverlap(File1, File2,splicing_factors_file,oncogenes_file, transcription_factors):
    #Prepare output file paths
    file_path=os.path.split(File1)[0]
    name1=file_path.split('/')[-2]
    file_path2=os.path.split(File2)[0]
    name2=file_path2.split('/')[-2]
    #print name
    filename_events='Overlap-comparison_events'
    filename_all='Overlap-comparison_all'
    #filename_sf= "splicing_factors"
    #filename_onco= "oncogenes"
    suffix = '.txt'
    path='/'.join(file_path.split('/')[0:-2])
    Output_dir=os.path.join(path,'Overlap_'+name1+"_"+name2)
    if os.path.exists(Output_dir) == False:
        os.mkdir(Output_dir)
    file_Out_events=os.path.join(Output_dir, filename_events + suffix)
    file_Out_all=os.path.join(Output_dir, filename_all + suffix)
    
    
    ### Isolate genes and splicing events/gene expression levels from input files
    eventsList_1=GetMatchingSignatures(File1)
    eventsList_2=GetMatchingSignatures(File2)
    #print genesList_1
    #print genesList_2
    
    ### Compare overlapping and unique splicinge/gene expression values. This will match based on the gene impacted AND the event (i.e. splicing event or gene expression)
    events_overlap=list(set(eventsList_1) & set(eventsList_2))
    events_List1_unique=list(set(eventsList_1) - set(eventsList_2))
    events_List2_unique=list(set(eventsList_2) - set(eventsList_1))
    #print genes_overlap
    #print genes_List1_unique
    
    ### Turn the events_list into a dataframe with the gene and alteration (i.e. splicing event or gene expression) in seperat column.
    events_df_1=process_genes_events(eventsList_1, name1)
    events_df_2=process_genes_events(eventsList_2, name2)
    
    merged_events_df=pd.concat([events_df_1, events_df_2], axis=0)
    #print merged_events_df.head()
    if 'Splicing_Event' in merged_events_df.columns:
        merged_events_df['Group'] = np.where(merged_events_df.duplicated(subset=['Genes', 'Splicing_Event'], keep=False), 'Both', merged_events_df["Group"])
    if 'Gene_Expression' in merged_events_df.columns:
        merged_events_df['Group'] = np.where(merged_events_df.duplicated(subset=['Genes', 'Gene_Expression'], keep=False), 'Both', merged_events_df["Group"])
    merged_events_df=merged_events_df.drop_duplicates()
    merged_events_df.groupby(["Genes"])
    #print merged_events_df.head()
    
    ###Compare overlapping and gene expression. This will seperate the gene and alteration.
    
    '''Get splicing, oncogene, and transcription factor data''' 
    splicing_factors_data = pd.read_csv(splicing_factors_file, sep=',')
    splicing_factors_list= splicing_factors_data.iloc[:, 0].tolist()
    oncogenes_data = pd.read_csv(oncogenes_file, sep='\t', low_memory=False)
    oncogenes_list= oncogenes_data.iloc[:,0].tolist()
    tf_data = pd.read_csv(transcription_factors, sep='\t', low_memory=False)
    tf_list= oncogenes_data.iloc[:,0].tolist()
    #print splicing_factors_list
    #rint oncogenes_list

    
    ''' Alternate usage. Print individual files for splicing factors and cancer drivers ''' 
    #file_Out_sp=os.path.join(Output_dir, filename_sf + suffix)
    #file_Out_onco=os.path.join(Output_dir, filename_onco + suffix)

    #print Output_dir
    #filename_sf= "".join(name1+"_"+name2+"_splicing_factors")
    #filename_onco= "".join(name1+"_"+name2+"_oncogenes")
    #file_Out_sp=os.path.join(Output_dir, filename_sf + suffix)
    #file_Out_onco=os.path.join(Output_dir, filename_onco + suffix)
    #print file_Out_sp
    #print file_Out_onco
    
    events_overlap=list(set(eventsList_1) & set(eventsList_2))
    events_List1_unique=list(set(eventsList_1) - set(eventsList_2))
    events_List2_unique=list(set(eventsList_2) - set(eventsList_1))
    
    ''' Make dataframe for overlapping events (splicing or gene expression)'''
    #nest_overlap_unique=[overlap, List1_unique, List2_unique]
    #overlap_unique_df=pd.DataFrame(nest_overlap_unique, ['Overlap', column2, column3]).T
    column1=''.join(name1+'_Only')
    column2=''.join(name2+'_Only')
    overlap_df=pd.DataFrame(events_overlap, columns=['Overlap'])
    List1_unique_df=pd.DataFrame(events_List1_unique, columns=[column1])
    List2_unique_df=pd.DataFrame(events_List2_unique, columns=[column2])
    
    '''
    Make dataframe for overlapping gene+events (events = splicing or gene expression)
    column1=''.join(name1+'_Only')
    column2=''.join(name2+'_Only')
    overlap_df=pd.DataFrame(events_overlap, columns=['Overlap'])
    List1_unique_df=pd.DataFrame(events_List1_unique, columns=[column1])
    List2_unique_df=pd.DataFrame(events_List2_unique, columns=[column2])
    '''
    
    ''' Annotate list '''
    overlap_annotated=AnnotateGenes(overlap_df, splicing_factors_list, oncogenes_list, tf_list)
    List1_annotated=AnnotateGenes(List1_unique_df, splicing_factors_list, oncogenes_list, tf_list)
    List2_annotated=AnnotateGenes(List2_unique_df, splicing_factors_list, oncogenes_list, tf_list)
    
    merged_events_df_annotated=AnnotateGenes(merged_events_df, splicing_factors_list, oncogenes_list, tf_list)
    
    print merged_events_df_annotated.head()

        
    merged_df=pd.concat([overlap_annotated, List1_annotated, List2_annotated], axis=1)
    print merged_df.head(5)
    print file_Out_events
    print file_Out_all
    merged_df.to_csv(file_Out_events, sep='\t', header=True, index=False)
    merged_events_df_annotated.to_csv(file_Out_all, sep='\t', header=True, index=False)


if __name__ == '__main__':
    #metaData_File1='/Volumes/salomonis2/NCI-R01/TCGA_Audrey/clusters_bams/GeneExpression/R1-V7_RBM10_mut/DEGs-LogFold_0.263034405834_rawp/GE.R1-V7_RBM10_mut_vs_NULL.txt'
    #metaData_File2='/Volumes/salomonis2/NCI-R01/TCGA_Audrey/clusters_bams/GeneExpression/R1-V7_RBM10_wt/DEGs-LogFold_0.263034405834_rawp/GE.R1-V7_RBM10_wt_vs_NULL.txt'
    metaData_File1='/Volumes/salomonis2/NCI-R01/TCGA_Audrey/clusters_bams/Splicing_Analysis/SVMOutputs/Splicing_Events_Subgroups/R1-V7_RBM10_mut/Events-dPSI_0.1_adjp/PSI.R1-V7_RBM10_mut_vs_NULL.txt'
    metaData_File2='/Volumes/salomonis2/NCI-R01/TCGA_Audrey/clusters_bams/Splicing_Analysis/SVMOutputs/Splicing_Events_Subgroups/RBM10_mut_other/Events-dPSI_0.1_adjp/PSI.RBM10_mut_other_vs_NULL.txt'
    SplicingFactors='/Volumes/salomonis2/NCI-R01/TCGA_Audrey/clusters_bams/Mutation_data/proteins_splicing_human_summary.csv'
    Oncogenes='/Volumes/salomonis2/NCI-R01/TCGA_Audrey/clusters_bams/Mutation_data/CosmicMutantExportCensus.tsv'
    TranscriptionFactors='/Volumes/salomonis2/NCI-R01/TCGA_Audrey/clusters_bams/Mutation_data/TF_list_MAtt_Weirauch.txt'
    Overlap=CompareOverlap(metaData_File1,metaData_File2, SplicingFactors,Oncogenes,TranscriptionFactors)
    