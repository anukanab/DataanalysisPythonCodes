#!/usr/bin/env python

import sys,string,os
sys.path.insert(1, os.path.join(sys.path[0], '..')) ### import parent dir dependencies

import export
import numpy as np
import os.path
from numpy import corrcoef, sum, log, arange
from scipy.stats.stats import pearsonr
import traceback
from stats_scripts import statistics
import unique
import UI
import getopt

def filepath(filename):
    fn = unique.filepath(filename)
    return fn

def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data

        
###################################################
############### External Comparison ###############

def matchNeoantigenSequencesInFreeText(neoantigen_file,genes_file,free_text_file,coordinates=None):
    output_dir = export.findParentDir(neoantigen_file)
    eo = export.ExportFile(output_dir+'/corresponding_overlaps2.txt')
    neoantigen_seq = {}
    for line in open(neoantigen_file,'rU').xreadlines():
        data = cleanUpLine(line)
        neoantigen, uid = string.split(data,'\t')
        #neoantigen_seq[neoantigen[:6]]=uid+'__'+neoantigen
        #neoantigen_seq[neoantigen[-6:]]=uid+'__'+neoantigen
        n = 7
        x = [neoantigen[i:i+n] for i in range(0, len(neoantigen), n)]
        for i in x:
            if len(i)>5:
                neoantigen_seq[i]=uid+'__'+neoantigen
    
    if coordinates != None:
        neoantigen_seq = {}
        for line in open(coordinates,'rU').xreadlines():
            data = cleanUpLine(line)
            t = string.split(data,'\t')
            chr = t[0]
            start = t[1]
            stop = t[2]
            neoantigen = t[3]
            neoantigen_seq[start]=neoantigen
            neoantigen_seq[stop]=neoantigen
            neoantigen_seq[str(int(start)+1)]=neoantigen
            neoantigen_seq[str(int(stop)+1)]=neoantigen
            neoantigen_seq[str(int(start)-1)]=neoantigen
            neoantigen_seq[str(int(stop)-1)]=neoantigen
                
    genes={}
    for line in open(genes_file,'rU').xreadlines():
        data = cleanUpLine(line)
        genes[data]=[]
        
    for line in open(free_text_file,'rU').xreadlines():
        data = cleanUpLine(line)
        data = string.replace(data,' ','')
        for gene in genes:
            if gene in data:
                pass
                #eo.write(gene+'\n')
        for neoantigen in neoantigen_seq:
            if neoantigen in data:
                eo.write(neoantigen_seq[neoantigen]+'\t'+neoantigen+'\t'+data+'\n')
                print neoantigen_seq[neoantigen][:30],data
    eo.close()

###################################################
################### SNAF B ########################

def integrateSNAFBPredictions(file_dir,EventAnnotationFile):
    import collections
    
    entries = collections.OrderedDict()
    junction_db = collections.OrderedDict()
    
    def importSNAF(file,filename):
        if 'sr_str3' in filename:
            annotation_type = 'SR-only'
        if 'sr_str4' in filename:
            annotation_type = 'SR-LR-junction-valdation'
        if 'sr_str5' in filename:
            annotation_type = 'SR-LR-protein-valdation'
        if 'lr' in filename:
            annotation_type = 'LR-only'
        if 'insertion' in filename:
            insertion = 'TRUE'
        else:
            insertion = 'FALSE'
        if 'deletion' in filename:
            deletion = 'TRUE'
        else:
            deletion = 'False'
        if 'True' in filename:
            extracellular = 'TRUE'
        else:
            extracellular = 'FALSE'
        
        header = True
        for line in open(file,'rU').xreadlines():
            data = cleanUpLine(line)
            if header:
                header=False
            else:
                Candidate_id,NeoJunction,mode,evidence,mRNA_sequence,peptide_sequence,gene_symbol,cohort_frequency,tumor_specificity_mean,tumor_specificity_mle,validation = string.split(data,'\t')
                junction_db[NeoJunction]=gene_symbol,cohort_frequency,tumor_specificity_mean,tumor_specificity_mle
                try:
                    db = entries[NeoJunction]
                    try: db[annotation_type].append([len(peptide_sequence),peptide_sequence,mRNA_sequence,mode,evidence,insertion,extracellular])
                    except: db[annotation_type] = [[len(peptide_sequence),peptide_sequence,mRNA_sequence,mode,evidence,insertion,extracellular]]
                except:
                    db={}
                    db[annotation_type] = [[len(peptide_sequence),peptide_sequence,mRNA_sequence,mode,evidence,insertion,extracellular]]
                    entries[NeoJunction] = db

    files = UI.read_directory(file_dir)
    for file in files:
        if 'report' in file:
            report = file[:-4]
            importSNAF(file_dir+'/'+file,report)
    """
    importSNAF(file1,'SR-only')
    importSNAF(file2,'SR-LR-junction-valdation')
    importSNAF(file3,'SR-LR-protein-valdation')
    importSNAF(file4,'LR-only')
    """
            
    junction_annotations = importEventAnnotations(EventAnnotationFile)
    
    gene_to_protein={}
    reference_protein_seq = '/Users/bha1ov/Desktop/Code/AltAnalyze/AltDatabase/EnsMart100/Hs/RNASeq/Hs_all-transcript-matches.txt'
    for line in open(reference_protein_seq,'rU').xreadlines():
        data = cleanUpLine(line)
        uid,transcript1,transcript2 = string.split(data,'\t')
        gene = string.split(uid,':')[0]
        try:
            gene_to_protein[gene].append(transcript1)
            gene_to_protein[gene].append(transcript2)
        except:
            gene_to_protein[gene]=[transcript1]
            gene_to_protein[gene].append(transcript2)
        
    protein_seq_db={}
    alt_protein_seq_db={}
    extracellular_domain_db={}
    
    protein_translation_dir = '/Users/bha1ov/Desktop/Code/AltAnalyze/AltDatabase/EnsMart100/ensembl/Hs/Hs_Ensembl_Protein__100_38.tab'
    protein_translation_db={}
    for line in open(protein_translation_dir,'rU').xreadlines():
        data = cleanUpLine(line)
        gene,transcript,protein = string.split(data,'\t')
        protein_translation_db[protein] = gene
    
    reference_protein_dir = '/Users/bha1ov/Desktop/Code/AltAnalyze/AltDatabase/EnsMart100/Hs/SequenceData/output/sequences' ### can include other sequence databases
    files = UI.read_directory(reference_protein_dir)
    for file in files:
        path = reference_protein_dir+'/'+file
        for line in open(path,'rU').xreadlines():
            data = cleanUpLine(line)
            transcript,protein,seq = string.split(data,'\t')
            protein_seq_db[transcript]=protein,seq,len(seq)
            alt_protein_seq_db[transcript]=seq
    
    uniprot_coordinates = '/Users/bha1ov/Desktop/Code/AltAnalyze/AltDatabase/EnsMart100/uniprot/Hs/Hs_FeatureCoordinate.txt'
    for line in open(uniprot_coordinates,'rU').xreadlines():
        data = cleanUpLine(line) #ensembl_prot	aa_start	aa_stop	genomic_start	genomic_stop	name	interpro_id	description
        try: ensembl_protein,aa_start,aa_stop,start,stop,name, interpro_id,domain = string.split(data,'\t')
        except: print data;sys.exit()

        if ensembl_protein in protein_translation_db:
            gene = protein_translation_db[ensembl_protein]
            if 'TOPO_DOM-Extracell' in domain:
                try: extracellular_domain_db[gene].append([int(start),int(stop)])
                except: extracellular_domain_db[gene] = [[int(start),int(stop)]]
            
    for NeoJunction in entries:
        db = entries[NeoJunction]
        for annotation in db:
            ls = db[annotation]
            ls.sort()
            ls.reverse()
            db[annotation]=ls[0][1:] ### longest protein - example
            
    eo = export.ExportFile(file_dir+'/summary/combined-evidence-SNAF-B.txt')
    header = ['NeoJunction','gene_symbol','event-type','coordinates','AltAnalyze ID','EventAnnotation','clusterID','max_dPSI','protein_predictions','cohort_frequency','tumor_specificity_mean','tumor_specificity_mle']
    header += ['ref_prot','ref_prot_len','ref_seq','ExNeoEpitopes_type','ExNeoEpitopes_evidence','ExNeoEpitopes_prot','ExNeoEpitopes_len','ExNeoEpitopes_mRNA','LR_evidence','SR-LR_evidence','SR-LR_match','SR_evidence',
    'NeoProtein','Overlaps with ECD','UniProt TOPO_DOM-Extracellular Genomic Coordinates','insertion-SNAF','extracellular-SNAF']
    eo.write(string.join(header,'\t')+'\n')
    for NeoJunction in entries:
        uid,protein_predictions,max_dPSI,clusterID,coordinates,EventAnnotation = junction_annotations[NeoJunction]
        gene_symbol,cohort_frequency,tumor_specificity_mean,tumor_specificity_mle = junction_db[NeoJunction]
        db = entries[NeoJunction]
        try: sr_peptide_sequence,sr_mRNA_sequence,sr_mode,sr_evidence,sr_insertion,sr_extracellular = db['SR-only']; SR='TRUE'
        except: sr_peptide_sequence=''; sr_mRNA_sequence=''; sr_mode=''; sr_evidence=''; SR=''
        try: srlr1_peptide_sequence,srlr1_mRNA_sequence,srlr1_mode,srlr1_evidence,slr1_insertion,slr1_extracellular = db['SR-LR-junction-valdation']; SRLR1='TRUE'
        except: srlr1_peptide_sequence=''; srlr1_mRNA_sequence=''; srlr1_mode=''; srlr1_evidence=''; SRLR1=''
        try: srlr2_peptide_sequence,srlr2_mRNA_sequence,srlr2_mode,srlr2_evidence,slr2_insertion,slr2_extracellular = db['SR-LR-protein-valdation']; SRLR2='TRUE'
        except: srlr2_peptide_sequence=''; srlr2_mRNA_sequence=''; srlr2_mode=''; srlr2_evidence=''; SRLR2=''
        try: lr_peptide_sequence,lr_mRNA_sequence,lr_mode,lr_evidence,lr_insertion,lr_extracellular = db['LR-only']; LR='TRUE'
        except: lr_peptide_sequence=''; lr_mRNA_sequence=''; lr_mode=''; lr_evidence=''; LR=''
    
        if LR=='TRUE':
            peptide_sequence,mRNA_sequence,mode,evidence,insertion,extracellular = lr_peptide_sequence,lr_mRNA_sequence,lr_mode,lr_evidence,lr_insertion,lr_extracellular
        elif SRLR2=='TRUE':
            peptide_sequence,mRNA_sequence,mode,evidence,insertion,extracellular = srlr2_peptide_sequence,srlr2_mRNA_sequence,srlr2_mode,srlr2_evidence, slr2_insertion,slr2_extracellular
        elif SRLR1=='TRUE':
            peptide_sequence,mRNA_sequence,mode,evidence,insertion,extracellular = srlr1_peptide_sequence,srlr1_mRNA_sequence,srlr1_mode,srlr1_evidence, slr1_insertion,slr1_extracellular
        else:
            peptide_sequence,mRNA_sequence,mode,evidence,insertion,extracellular = sr_peptide_sequence,sr_mRNA_sequence,sr_mode,sr_evidence, sr_insertion,sr_extracellular
        if sr_evidence in protein_seq_db:
            protein,seq,seq_len=protein_seq_db[sr_evidence]
        else:
            protein=''; seq_len = ''; seq = ''
                        
        ### Determine if it is an inclusion exon or not
        junction1, junction2 = string.split(coordinates,'|')
        chr,junction1 = string.split(junction1,':')
        junction2 = string.split(junction2,':')[1]
        junc1_elements = map(int,string.split(junction1,'-'))
        junc2_elements = map(int,string.split(junction2,'-'))
        diff = abs(junc2_elements[0]-junc2_elements[1])-abs(junc1_elements[0]-junc1_elements[1])
        if diff>0: inclusion = 'inclusion'
        else: inclusion = 'exclusion'
        
        ### Does the protein exist in the AltAnalyze protein database
        neoprotein = 'TRUE'
        ensembl = string.split(NeoJunction,':')[0]
        if ensembl in gene_to_protein:
            #if 'ENSG00000101187:E20.1-E21.2' == NeoJunction:
            ens_proteins = gene_to_protein[ensembl]
            for ens_prot in ens_proteins:
                if ens_prot in alt_protein_seq_db:
                    if peptide_sequence in alt_protein_seq_db[ens_prot]:
                        #print peptide_sequence; break
                        neoprotein = 'FALSE'
                        break
        
        ### Does the junction overlap with a domain
        junc_start, junct_stop = junc1_elements
        extracellular_evidence = 'FALSE'
        evidence = ''
        
        if ensembl in extracellular_domain_db:
            #if 'ENSG00000179915:E21.1-E31.2' == NeoJunction:
            for (ex_start1, ex_stop1) in extracellular_domain_db[ensembl]:
                temp = [ex_start1, ex_stop1]
                temp.sort()
                ex_start, ex_stop = temp
                ex_start1, ex_stop1 = str(ex_start1), str(ex_stop1)
                #print [ex_start, ex_stop],junc1_elements
                if junc_start >= ex_start and junc_start <= ex_stop:
                    extracellular_evidence = 'TRUE'
                    evidence = chr+':'+ex_start1+'|'+chr+':'+ex_stop1
                elif junct_stop >= ex_start and junct_stop <= ex_stop:
                    extracellular_evidence = 'TRUE'
                    evidence = chr+':'+ex_start1+'|'+chr+':'+ex_stop1
                #print extracellular_evidence
                                
        event = [NeoJunction,gene_symbol,inclusion,coordinates,uid,EventAnnotation,clusterID,max_dPSI,protein_predictions,cohort_frequency,tumor_specificity_mean,tumor_specificity_mle]
        event += [sr_evidence,str(seq_len),seq,mode,evidence,peptide_sequence,str(len(peptide_sequence)),mRNA_sequence,LR,SRLR2,SRLR1,SR,neoprotein,extracellular_evidence,evidence,insertion,extracellular]
        eo.write(string.join(event,'\t')+'\n')

    eo.close()

###################################################
################### SNAF T ########################

def collapseRankedNeoantigens(neo_file):

    gene_to_protein={}

    reference_protein_seq = '/Users/bha1ov/Desktop/Code/AltAnalyze/AltDatabase/EnsMart100/Hs/RNASeq/Hs_all-transcript-matches.txt'
    for line in open(reference_protein_seq,'rU').xreadlines():
        data = cleanUpLine(line)
        uid,transcript1,transcript2 = string.split(data,'\t')
        gene = string.split(uid,':')[0]
        try:
            gene_to_protein[gene].append(transcript1)
            gene_to_protein[gene].append(transcript2)
        except:
            gene_to_protein[gene]=[transcript1]
            gene_to_protein[gene].append(transcript2)
        
    protein_seq_db={}
    alt_protein_seq_db={}

    reference_protein_dir = '/Users/bha1ov/Desktop/Code/AltAnalyze/AltDatabase/EnsMart100/Hs/SequenceData/output/sequences'
    files = UI.read_directory(reference_protein_dir)
    for file in files:
        path = reference_protein_dir+'/'+file
        for line in open(path,'rU').xreadlines():
            data = cleanUpLine(line)
            transcript,protein,seq = string.split(data,'\t')
            protein_seq_db[transcript]=protein,seq,len(seq)
            alt_protein_seq_db[transcript]=seq
       
    import collections
    neo_full_db = collections.OrderedDict()
    junction_db = collections.OrderedDict()
    inframe_db = {}
    collaped_neojunction_antigens=collections.OrderedDict()
    eo = export.ExportFile(neo_file[:-4]+'_collapsed.txt')
    firstRow=True
    for line in open(neo_file,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if firstRow:
            eo.write(line)
            inframe_ind = t.index('In-Frame')
            firstRow=False
        else:
            neojunction = t[1]
            neoantigen = t[0]
            #if neojunction == 'ENSG00000164175:E3.2_33963931-E4.2': print neoantigen
            try: junction_db[neojunction].append(neoantigen)
            except: junction_db[neojunction] = [neoantigen]
            neo_full_db[neojunction,neoantigen] = t
            if_call = t[inframe_ind]
            if neojunction in inframe_db:
                iF = inframe_db[neojunction]
                iF[neoantigen]=if_call
            else:
                iF = {neoantigen:if_call}
                inframe_db[neojunction]=iF
        
    temp=0
    ### First objective is for a given neoantigen, to find out how many overlapping sets there are
    neojunction_lists = collections.OrderedDict()
    for neojunction in junction_db:
        #if len(junction_db[neojunction])>4:
        strFrag = junction_db[neojunction]
        overlapping_frags = {}
        for frag1 in strFrag:
            for frag2 in strFrag:
                if frag1!=frag2:
                    if frag1[3:7] in frag2: #### Critical line for potential issues - if > 4AA overlap non-specifically, will break "combined"
                        if frag1[3:7] in overlapping_frags:
                            if frag2 not in overlapping_frags[frag1[3:7]]:
                                overlapping_frags[frag1[3:7]].append(frag2)
                        else:
                            overlapping_frags[frag1[3:7]] = [frag1,frag2]
        if len(overlapping_frags)>0:  
            for f in overlapping_frags:
                # dictionary values of overlapping peptides, appended as sets
                try: neojunction_lists[neojunction].append(overlapping_frags[f])
                except: neojunction_lists[neojunction] = [overlapping_frags[f]]
        else:
            for pep in strFrag:
                inframe_call = inframe_db[neojunction][pep]
                try: collaped_neojunction_antigens[neojunction].append([pep,inframe_call])
                except: collaped_neojunction_antigens[neojunction]=[[pep,inframe_call]]

    ### Here we will collapse the overlapping peptides into sets (e.g., forward and reverse strand)
    for neojunction in neojunction_lists:
        #if neojunction == 'ENSG00000205944:E20.1-E21.1':
        data = neojunction_lists[neojunction]
        """
        for i in data:
            if len(i)>8: print i
            """
        strFrag = junction_db[neojunction]
        #https://stackoverflow.com/questions/30917226/collapse-list-of-lists-to-eliminate-redundancy
        result=[]
        ### Combine lists of overlapping peptides together - step 2
        for d in data: ### d is a list of overlapping peptides based on core sequence overlaps - d is redundant with other d's for the junction
            d = set(d)
            matched = [d]
            unmatched = []
            # first divide into matching and non-matching groups
            for g in result:
                if d & g:
                    matched.append(g)
                else:
                    unmatched.append(g)
            # then combine all matching groups into one group
            # while leaving unmatched groups intact

            result = unmatched + [set().union(*matched)]

        for pep_set in result:
            ### collapse each set into a single combined long string
            pep_set = unique.unique(list(pep_set))
            redundant = []
            for pep1 in pep_set:
                for pep2 in pep_set:
                    if pep1 != pep2:
                        if pep1 in pep2:
                            redundant.append(pep1)
                        elif pep2 in pep1:
                            redundant.append(pep2)
            inframe_call = inframe_db[neojunction][pep1]
            redundant = unique.unique(redundant)
            #print 'redundnat',redundant ### remove these since these will slow down the analysis
            pep_set2=[]
            for pep in pep_set:
                if pep not in redundant:
                    pep_set2.append(pep)
            combined = assemble(pep_set2)
            try: 
                collapsed_neo = str(list(combined)[0])
            except: 
                print pep_set2, combined;sys.exit()
            try: collaped_neojunction_antigens[neojunction].append([collapsed_neo,inframe_call])
            except: collaped_neojunction_antigens[neojunction] = [[collapsed_neo,inframe_call]]

    valid_entries=0
    partial_protein_matches=0
    neojunctions_added = []
    for (neojunction,neoantigen) in neo_full_db:
        if neojunction not in neojunctions_added:
            neojunctions_added.append(neojunction)
            t = neo_full_db[(neojunction,neoantigen)]
            
            match = False
            true_peptide_sequence=''
            for peptide_sequence in junction_db[neojunction]: #collaped_neojunction_antigens[neojunction]
                ensembl = string.split(neojunction,'.')[0]
                if ensembl in gene_to_protein:
                    ens_proteins = gene_to_protein[ensembl]
                    for ens_prot in ens_proteins:
                        if ens_prot in alt_protein_seq_db:
                            if peptide_sequence in alt_protein_seq_db[ens_prot]:
                                match = True ### Hence the peptide exists in our refernece protein database - not interesting
                                #print peptide_sequence; break
                                break
                            else:
                                if peptide_sequence[:3] in alt_protein_seq_db[ens_prot]:
                                    true_peptide_sequence = peptide_sequence
                                    #print peptide_sequence, t[:5]; break
                                elif peptide_sequence[-3:] in alt_protein_seq_db[ens_prot]:
                                    true_peptide_sequence = peptide_sequence
            if match == False:
                peps=[]
                for (pep,inframe_call) in collaped_neojunction_antigens[neojunction]:
                    peps.append(pep + '('+inframe_call+')')
                ### In case the above collapsing step failed to include all predicted neoantigens
                for short_pep in junction_db[neojunction]:
                    #if neojunction == 'ENSG00000164175:E3.2_33963931-E4.2':
                    inframe_call = inframe_db[neojunction][short_pep]
                    found=False
                    for (long_pep, frame) in collaped_neojunction_antigens[neojunction]:
                        if short_pep in long_pep: found = True
                    if found == False:
                        peps.append(short_pep + '('+inframe_call+')')

                peps = string.join(peps,'|')
                t = t+[peps]
                eo.write(string.join(t+[true_peptide_sequence],'\t')+'\n')
                valid_entries+=1
                if len(true_peptide_sequence)>0:
                    partial_protein_matches+=1
    eo.close()
    print valid_entries,partial_protein_matches

def assemble(str_list, min=3, max=15):
    #https://stackoverflow.com/questions/52528744/how-to-merge-strings-with-overlapping-characters-in-python
    if len(str_list) < 2:
        return set(str_list)
    output = set()
    string = str_list.pop()
    for i, candidate in enumerate(str_list):
        matches = set()
        if candidate in string:
            matches.add(string)
        elif string in candidate:
            matches.add(candidate)
        for n in range(min, max + 1):
            if candidate[:n] == string[-n:]:
                matches.add(string + candidate[n:])
                #print string, candidate, matches
            if candidate[-n:] == string[:n]:
                matches.add(candidate[:-n] + string)
        #print matches
        for match in matches:
            output.update(assemble(str_list[:i] + str_list[i + 1:] + [match]))
    return output

class PeptideInformation:
    def __init__(self,symbol,coord,junction_read_count,inFrame,tumor_specificity_mean, tumor_specificity_mle, n_sample, present_in_ensembl):
        self.symbol = symbol
        self.coord = coord
        self.junction_read_count = junction_read_count
        self.inFrame = inFrame
        self.tumor_specificity_mean = tumor_specificity_mean
        self.tumor_specificity_mle = tumor_specificity_mle
        self.n_sample = n_sample
        self.hla_binding_db = {}
        self.immunogenicity = {}
        self.samples = []
        self.present_in_ensembl = present_in_ensembl
    def Symbol(self): return self.symbol
    def Coord(self): return self.coord
    def JunctionCount(self): return self.junction_read_count
    def InFrame(self): return self.inFrame
    def TumorSpecificityMean(self): return self.tumor_specificity_mean
    def TumorSpecificityMLE(self): return self.tumor_specificity_mle
    def NumSamples(self): return self.n_sample
    def PresentInEnsembl(self): return self.present_in_ensembl
    def setAntigenData(self,sample,hla,binding_affinity,immunogenicity):
         self.hla_binding_db[hla,sample] = binding_affinity
         self.immunogenicity[hla,sample] = immunogenicity
         self.samples.append(sample)
    def MeanImmunogenicity(self):
        mean = []
        for (hla,sample) in self.immunogenicity:
            mean.append(float(self.immunogenicity[hla,sample]))
        return str(statistics.avg(mean))
    def MeanBinding(self):
        mean = []
        for (hla,sample) in self.hla_binding_db:
            mean.append(float(self.hla_binding_db[hla,sample]))
        return str(statistics.avg(mean))
    def Samples(self):
        return string.join(unique.unique(self.samples),'|')
  
def importEventAnnotations(EventAnnotationFile):
    junction_annotations={}
    for line in open(EventAnnotationFile,'rU').xreadlines():
        data = cleanUpLine(line)
        #Symbol	Description	Examined-Junction	Background-Major-Junction	AltExons	ProteinPredictions	dPSI	ClusterID	UID	Coordinates	EventAnnotation
        t = string.split(data,'\t')
        symbol = t[0]
        uid = t[8]
        junction = t[2]
        protein_predictions = t[5]
        max_dPSI = t[6]
        clusterID = t[7]
        coordinates=t[9]
        EventAnnotation=t[10]
        #junction_id = string.replace(junction,':','.')
        #junction_id = string.replace(junction_id,'-','.')
        junction_annotations[junction] = uid,protein_predictions,symbol,clusterID,coordinates,EventAnnotation
    print len(junction_annotations), 'unique total junctions'
    return junction_annotations

def genericImport(filepath,key=None):
    db={}
    firstRow = True
    for line in open(filepath,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        uid = t[0]
        values = t[1:]
        if key == 'multiple':
            uid = t[0],t[1]
            values = t[2:]
        if firstRow:
            firstRow = False
            if key == 'multiple':
                uid = 'header'
                values = t[2:]
            else:
                uid = 'header'
        db[uid] = values
    return db

def collapseAnnotateNeoB(filepath,customJunctionData=None):
    eo = export.ExportFile(filepath[:-4]+'-updated.txt')

    ### Import custom combined junction annotaitons
    if customJunctionData != None:
        custom_annotations={}
        firstRow = True
        for line in open(customJunctionData,'rU').xreadlines():
            data = cleanUpLine(line)
            t = string.split(data,'\t')
            if firstRow:
                customHeader=t[1:]
                firstRow=False
            uid = t[0]
            values = t[1:]
            neojunctions = string.split(uid,'|')
            for neojunction in neojunctions:
                custom_annotations[neojunction] = values
    else:
        customHeader=[]
        custom_annotations={}

    firstRow = True
    neopeptide_db={}
    for line in open(filepath,'rU').xreadlines():
            data = cleanUpLine(line)
            t = string.split(data,'\t')
            if firstRow:
                header = t
                level_ind = header.index('Level')
                neoEx_ind = header.index('ExNeoEpitopes_prot')
                firstRow=False
            else:
                neojunction = t[0]
                values = t[1:]
                neoEx = t[neoEx_ind]
                level = int(t[level_ind])
                if neojunction in custom_annotations:
                    values+=custom_annotations[neojunction]
                else:
                    values+=['']*len(customHeader)
                ### Collapse all junctinons into common events by annotation stringency
                try: 
                    neopeptide_db[neoEx].append([level,neojunction,values])    
                except:
                    neopeptide_db[neoEx] = [[level,neojunction,values]]

    eo.write(string.join(header+customHeader,'\t')+'\n')
    for neoEx in neopeptide_db:
        neopeptide_db[neoEx].sort()
        collapsed_junctions = []
        for (level,neojunction,values) in neopeptide_db[neoEx]:
            collapsed_junctions.append(neojunction)
        collapsed_junctions = string.join(collapsed_junctions,'|')
        exemplar_junction_values = neopeptide_db[neoEx][0][-1]
        eo.write(string.join([collapsed_junctions]+exemplar_junction_values,'\t')+'\n')
    eo.close()

def updatedSNAFTMerge(file_dir, EventAnnotationFile, survival_file=None, customPeptideData=None, customJunctionData=None):
    eo = export.ExportFile(file_dir+'/summary/collapsed_SNAF-T.txt')
    junction_annotations = importEventAnnotations(EventAnnotationFile)
    
    survival_associations={}
    if survival_file!=None: 
        for line in open(survival_file,'rU').xreadlines():
            data = cleanUpLine(line)
            #gene, uid, LRT, Wald, SlrT, Zscore, n_sample, symbol, prop = string.split(data,'\t')
            uid, pep, junc, n_sample, Zscore, Wald = string.split(data,'\t')
            #survival_associations[uid] = LRT, Wald, SlrT, Zscore, n_sample, prop
            survival_associations[pep,junc] = n_sample, Zscore, Wald 
    
    if customJunctionData != None:
        customJunction_db = genericImport(customJunctionData)
        print len(customJunction_db),'customJunction_db'
    else:
        customJunction_db = {}

    if customPeptideData != None:
        customPeptide_db = genericImport(customPeptideData,key='multiple')
    else:
        customPeptide_db = {}
   
    peptide_db={}
    files = UI.read_directory(file_dir)
    for file in files:
        path = file_dir+'/'+file
        if '.bed' in file or '.txt' in file:
            firstRow=True
            for line in open(path,'rU').xreadlines():
                data = cleanUpLine(line)
                present_in_ensembl = 'UNK'
                try: sample, peptide, junction, junction_read_count, phase, evidence, hla, binding_affinity, immunogenicity, tumor_specificity_mean, tumor_specificity_mle, n_sample, coord, symbol, present_in_ensembl = string.split(data,'\t')
                except: 
                    sample, peptide, junction, junction_read_count, phase, evidence, hla, binding_affinity, immunogenicity, tumor_specificity_mean, tumor_specificity_mle, n_sample, coord, symbol = string.split(data,'\t')
                sample = string.replace(sample,'.bed','')
                if firstRow:
                    firstRow = False
                else:
                    if len(evidence)>3:
                        inFrame = "True"
                    else:
                        inFrame = "UNK"
                    if (peptide,junction) in peptide_db:
                        pi = peptide_db[(peptide,junction)]
                        pi.setAntigenData(sample,hla,binding_affinity,immunogenicity)
                    else:
                        pi = PeptideInformation(symbol,coord,junction_read_count,inFrame,tumor_specificity_mean, tumor_specificity_mle, n_sample, present_in_ensembl)
                        pi.setAntigenData(sample,hla,binding_affinity,immunogenicity)
                        peptide_db[(peptide,junction)] = pi
                        
    print len(peptide_db), 'unique total junction peptides'
    
    header = ['Peptide','Junction','Symbol','Coord','JunctionCount','In-Frame','Tumor Specificity Mean','Tumor Specificity MLE',
              'Num Samples','Mean Binding','Mean Immunogenicity','Samples','UID','ClusterID','Coordinates','EventAnnotation', 'Present-In-Ensembl',
            'Wald', 'Zscore', 'n_sample', 'peptideID']
    if len(customPeptide_db)>0:
        header+=customPeptide_db['header']
    if len(customJunction_db)>0:
        print customJunction_db['header']
        header+=customJunction_db['header']
    eo.write(string.join(header,'\t')+'\n')
    for (peptide,junction) in peptide_db:
        pi = peptide_db[(peptide,junction)]
        uid,protein_predictions,symbol,clusterID,coordinates,annotation = junction_annotations[junction]
        peptide_id = peptide+'.'+string.replace(junction,':','.'); peptide_id = string.replace(peptide_id,'-','.')
        alt_junction = string.replace(junction,':','.')
        alt_junction = string.replace(alt_junction,'-','.')

        if (peptide,junction) in survival_associations:
            n_sample, Zscore, Wald  = survival_associations[(peptide,junction)]
        elif (peptide,alt_junction) in survival_associations:
            n_sample, Zscore, Wald  = survival_associations[(peptide,alt_junction)]
        else:
            LRT=''; Wald='';SlrT='';Zscore='';n_sample='';prop=''
        outs = [peptide,junction,symbol,pi.Coord(),pi.JunctionCount(),pi.InFrame(),pi.TumorSpecificityMean(),pi.TumorSpecificityMLE(),
                pi.NumSamples(),pi.MeanBinding(),pi.MeanImmunogenicity(),pi.Samples(),uid,clusterID,coordinates,annotation, pi.PresentInEnsembl(),
                Wald, Zscore, n_sample, peptide_id]
        if len(customPeptide_db)>0:
            alt_junction = string.replace(junction,'-','.')
            alt_junction = string.replace(alt_junction,':','.')
            if (peptide,alt_junction) in customPeptide_db:
                values = customPeptide_db[peptide,alt_junction]
                outs+=values
            else:
                outs+=['']*len(customPeptide_db['header'])
        if len(customJunction_db)>0:
            if junction in customJunction_db:
                values = customJunction_db[junction]
                outs+=values
            else:
                outs+=['']*len(customJunction_db['header'])

        eo.write(string.join(outs,'\t')+'\n') 
    eo.close()
        
if __name__ == '__main__':
    survival_file = None
    analysis = 'T'
    customPeptideData = None
    customJunctionData = None
    if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
      print "Insufficient arguments";sys.exit()
    else:
        options, remainder = getopt.getopt(sys.argv[1:],'', ['i=','a=','overlap=','survival=','concordance=',
                                                'EventAnnotation=','custom=','customJunction='])
        for opt, arg in options:
            if opt == '--i': input_dir=arg
            elif opt == '--EventAnnotation': EventAnnotationFile = arg
            elif opt == '--a': analysis=arg
            elif opt == '--survival': survival_file = arg
            elif opt == '--dPSI': dPSI_cutoff = float(arg)
            elif opt == '--custom': customPeptideData = arg
            elif opt == '--customJunction': customJunctionData = arg
            elif opt == '--concordance':
                if string.lower(arg) == 'true' or string.lower(arg) == 'yes':
                    concordance = True
                else:
                    concordance = False
                
    if analysis == 'T':
        updatedSNAFTMerge(input_dir, EventAnnotationFile, survival_file=survival_file, customPeptideData=customPeptideData, customJunctionData=customJunctionData)
        sys.exit()
    if analysis == 'B':
        integrateSNAFBPredictions(input_dir,EventAnnotationFile)
        sys.exit()
    elif analysis == 'collapseRankedNeoantigens':
        #pre-consolidated SNAF-T angtigen file: '/Users/saljh8/Dropbox/Manuscripts/InProgress/SNAF/Figures/Patent/SNAF-T.txt'
        collapseRankedNeoantigens(input_dir);sys.exit()
    elif analysis == 'collapseAnnotateNeoB':
        #pre-consolidated SNAF-T angtigen file: '/Users/saljh8/Dropbox/Manuscripts/InProgress/SNAF/Figures/Patent/SNAF-T.txt'
        collapseAnnotateNeoB(input_dir,customJunctionData);sys.exit()
    elif analysis == 'matchNeoantigenSequencesInFreeText':
        genes_file = '/Users/bha1ov/Downloads/genes.txt'
        neoantigen_file = '/Users/bha1ov/Downloads/NeoB.txt'
        free_text_file = '/Users/bha1ov/Downloads/GBM.txt'
        coordinates = '/Users/bha1ov/Downloads/hg19-TCR.txt'
        matchNeoantigenSequencesInFreeText(neoantigen_file,genes_file,free_text_file,coordinates=coordinates);sys.exit()
