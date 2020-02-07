import sys,string,os
import pandas as pd

def CreateHg19LiftoverDict(filename1):    #Step 1: Create UID Conversion Dictionary
    firstLine = False
    Dict1={}          #Create empty dictionary
    for line in open(filename1, 'rU').xreadlines():   #the file it is reading has the coordinates as it is and also flipped to account for a change in strand
         if firstLine:
           firstLine = False
           continue
         else:
           data = line.strip()
           values = string.split(data, '\t')
           #print values
           hg38UID=values[3]
           hg19Coordinates=values[4]
           Dict1[hg19Coordinates]=hg38UID #create a dictionary of UID:Coordinates from Hg19 file
    #print dict(Dict1.items()); sys.exit()
    return Dict1

def CreateHg19SignatureDict(folder1):
    firstLine = True
    Dict2={}
    for files in os.listdir(folder1):
        filename = "".join(folder1 + '/' + files)
        if '.txt' in files and 'PSI.' in files:
           for line in open(filename, 'rU').xreadlines():
                data = line.strip()
                values = string.split(data, '\t')
                if firstLine:
                   firstLine = False
                   continue
                else:
                    Dict={}
                    uid19=values[0]
                    coordinates_full=values[7]
                    coordinates_new=coordinates_full.split("|")
                    coordinates=coordinates_new[0]
                    Dict[coordinates]=uid19    #create new dictionary
                    Dict2.update(Dict)         #update existing dictionary with the new k, v pairs
    #print dict(Dict2.items()); sys.exit()               
    return Dict2

def ReplaceHg19LiftoverUID(Dict1, Dict2):
    # Use a dictionary comprehension to collect Dict2 values of shared key
    Dict3 = {key:[Dict2[key], Dict1[key]] for key in Dict1 if key in Dict2} #dictionary comprehension to iterate through Dict1's keys and, if the key is in both Dict1 and Dict2, store the key in Dict3 with the value from Dict1 and Dict2
    Coordinates_hg19 = Dict3.keys()
    #print Coordinates_hg19[:5]
    UID_hg19 = Dict3.values()
    hg19LookupTable = pd.DataFrame(list(zip(Coordinates_hg19, UID_hg19)), columns =['Coordinates_hg19', 'UID_hg19'])
    hg19LookupTable[['UIDhg19','UID']] = pd.DataFrame(hg19LookupTable.UID_hg19.values.tolist(), index= hg19LookupTable.index)
    hg19LookupTable = hg19LookupTable.drop(columns='UID_hg19')
    return hg19LookupTable  

def ChangePSIFileUID(hg19LookupTable, folder):
    FilePath='/Volumes/salomonis2/NCI-R01/hg19-signatures/BRCA-PSI'     #Specificy output path
    for files in os.listdir(folder):                                    #Open output file
        filename = "".join(folder + '/' + files)
        if '.txt' in files and 'PSI.' in files:
            filename_raw=files[:-4]
            filename_out="".join(filename_raw+'-hg19.txt')
            FilePath_Out="".join(FilePath+"/"+filename_out)
            OutputFile=open(FilePath_Out, 'w')
            df=pd.read_csv(filename, sep='\t')
            df=df.merge(hg19LookupTable, on='UID')
            df['UID'] = df['UIDhg19'].values
            df['Coordinates'] = df['Coordinates_hg19'].values
            df = df.iloc[:, :-2]
            df.to_csv(OutputFile, sep='\t', index=None)
            
if __name__ == '__main__':
    
    folder1 = '/Volumes/salomonis2/NCI-R01/hg19-signatures'
    folder = '/Volumes/salomonis2/NCI-R01/TCGA-BREAST-CANCER/Anukana/MetaDataAnalysis/PSI/Events-dPSI_0.1_adjp'
    filename1 = '/Volumes/salomonis2/NCI-R01/TCGA-BREAST-CANCER/TCGA-files-Ens91/bams/AltResults/AlternativeOutput/hg19-hglft_genome.txt'
    #filename2 = '/Volumes/salomonis2/NCI-R01/TCGA-BREAST-CANCER/TCGA-files-Ens91/bams/AltResults/AlternativeOutput/JunctionsCoordinates-hg38o.txt'
    Dict1 = CreateHg19LiftoverDict(filename1)
    Dict2 = CreateHg19SignatureDict(folder1)
    hg19LookupTable = ReplaceHg19LiftoverUID(Dict1, Dict2)
    BRCAchangedcoordinate = ChangePSIFileUID(hg19LookupTable, folder)
    
                
                
                        
                
                

               