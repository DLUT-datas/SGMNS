# -*- coding: utf-8 -*-
"""
Created on Mon Jan 11 15:48:29 2021

@author: EDG Deft
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 09:37:35 2019

@author: dell
"""
import pandas as pd
from rdkit.Chem import AllChem 
from rdkit import Chem
from rdkit import DataStructs
from itertools import islice 
import os
import shutil
from itertools import islice 
#import xlwt
import re
from itertools import islice 
import numpy as np
import numpy
import gensim
from spec2vec import Spec2Vec
from spec2vec import SpectrumDocument
import spec2vec.utils
from matchms.filtering import add_losses
from matchms.filtering import add_parent_mass
from matchms.filtering import default_filters
from matchms.filtering import normalize_intensities
from matchms.filtering import require_minimum_number_of_peaks
from matchms.filtering import select_by_mz
from matchms.importing import load_from_mgf
import copy

print("\nParameter settings for SGMNS")
exe_adr=input("Please enter the pathway of the SGMNS.exe file，e.g. D:\SGMNS\n")
ion_mode=int(input("Ion mode of dataset：positive ion mode=1；negative ion mode=0, e.g. 1\n"))
MS2_mz=0.01
MS2_rt=30
SD_mz_threshold=float(input("Tolerance of the exact mass match/ppm，e.g. 10\n"))
SD_rt_threshold=float(input("Tolerance for the retention time match/%，e.g. 30\n"))
sim2_threshold=float(input("Threshold of the MS/MS spectral similarity，e.g. 0.5\n"))
result_name=str(input("Please enter the name of the result file as you expecte，e.g. result\n"))
bk_tf=int(input("Please decide whether you want to build a new background network or use the default network for the subsequent metabolite annotation：new=1；default=0, e.g. 0\n"))
def RemoveDir(filepath):
    if not os.path.exists(filepath):
        os.mkdir(filepath)
    else:
        shutil.rmtree(filepath)
        os.mkdir(filepath)
    return filepath
if bk_tf==1:
    beta=float(input("Please enter the structure similarity threshold to build the background network，e.g. 0.4\n"))
    smiles=[]
    molecule_id=[]
    if ion_mode==1:
        bk_adr=exe_adr + "\data/Background database_POS.csv"
    else:
        bk_adr=exe_adr + "\data/Background database_NEG.csv"
    with open(bk_adr)as f:
        for line in islice(f,1,None):
            line=line.strip().split(',')
            molecule_id.append(line[1])
            smiles.append(line[2])
    fingerprint_numb=len(smiles) 
    fp=[]
    for j in range(fingerprint_numb):
        smiles1=smiles[j]
        m1=Chem.MolFromSmiles(smiles1)
        fp.append(AllChem.GetMorganFingerprint(m1,2))
    print("\nConstruction of GCMN")
    print("Saving node information...")
    node_adr=RemoveDir(exe_adr + "\data/Nodes")
    for i in range(fingerprint_numb):
        HDMI_ID=molecule_id[i]
        c=node_adr + "//" + str(HDMI_ID) + '.txt' 
        fixedfp=fp[i]
        with open(c,'a') as net:
            for j in range(fingerprint_numb):
                similarity=round(DataStructs.DiceSimilarity(fixedfp,fp[j]),2)
                if similarity>=beta:    
                    net.write(str(j))
                    net.write(' ')
        net.close()

    print("Saving edge information...") 
    edge_adr=RemoveDir(exe_adr + "\data/Edges")
    for i in range(fingerprint_numb):
        HDMI_ID=molecule_id[i]
        c=edge_adr + "//" + str(HDMI_ID) + '.txt'
        fixedfp=fp[i]
        with open(c,'a') as net:
            for j in range(fingerprint_numb):
                similarity=round(DataStructs.DiceSimilarity(fixedfp,fp[j]),2)
                if similarity>=beta:
                    net.write(str(similarity))
                    net.write(' ')
        net.close()

    print("Saving the layout of background network...")
    net_adr=exe_adr + "\data/net.txt"
    with open(net_adr,'a') as net:
        for i in range(fingerprint_numb-1):
            fixedfp=fp[i]
            for j in range(i+1,fingerprint_numb):
                similarity=round(DataStructs.DiceSimilarity(fixedfp,fp[j]),2)
                if similarity>=beta: 
                    net.write(str(i))
                    net.write(' ')
                    net.write(str(j))
                    net.write(' ')
                    net.write(str(similarity))
                    net.write('\n')           
    print("Background network construction is completed.")

def RemoveDir(filepath):
    if not os.path.exists(filepath):
        os.mkdir(filepath)
    else:
        shutil.rmtree(filepath)
        os.mkdir(filepath)

#filter
def window_filter_peaks(peaks, window_size, top_peaks):
    peak_masses_to_keep = []
    for i in range(len(peaks)):
        peak_window=[]
        mass=peaks[i][0]
        for j in range(len(peaks)):
            candidate_mass=peaks[j][0]
            if candidate_mass-mass<-50:
                continue
            if abs(candidate_mass-mass)<=50:
                peak_window.append(peaks[j])
            if candidate_mass-mass>50:
                break
        peak_window=sorted(peak_window, key=lambda peak: peak[1],reverse=True)
        for peak in peak_window[:top_peaks]:
            if mass==peak[0]:
                peak_masses_to_keep.append(peaks[i])
                break
            
    return peak_masses_to_keep
                
def filter_precursor_peaks(peaks, tolerance_to_precursor, mz):
    new_peaks = []
    for peak in peaks:
        if abs(peak[0] - mz) > tolerance_to_precursor:
            new_peaks.append(peak)
    return new_peaks
#spec2vec  
def apply_my_filters(s):
    """This is how a user would typically design his own pre- and post-
    processing pipeline."""
    s = default_filters(s)
    s = add_parent_mass(s)
    s = normalize_intensities(s)
    #s = reduce_to_number_of_peaks(s, n_required=10, ratio_desired=0.5)
    s = select_by_mz(s, mz_from=0, mz_to=1000)
    #s = add_losses(s, loss_mz_from=0.0, loss_mz_to=1000.0)
    s = add_losses(s, loss_mz_from=10.0, loss_mz_to=200.0)
    s = require_minimum_number_of_peaks(s, n_required=0)
    return s

def check_words(model,spectrum):
    document = SpectrumDocument(spectrum)
    idx_not_in_model = [i for i, x in enumerate(document.words) if x not in model.wv.vocab]
    words_in_model = [x for i, x in enumerate(document.words) if i not in idx_not_in_model]
    if len(words_in_model)>0:
        return True
    else:
        return False
    
def check_missing_percentage(document,intensity_weighting_power):
    """Return True if model covers enough of the document words."""
    idx_not_in_model = [i for i, x in enumerate(document.words) if x not in model.wv.vocab]
    words_in_model = [x for i, x in enumerate(document.words) if i not in idx_not_in_model]
    weights_in_model = numpy.asarray([x for i, x in enumerate(document.weights)
                                      if i not in idx_not_in_model]).reshape(len(words_in_model), 1)
    weights_raised = numpy.power(weights_in_model, intensity_weighting_power)
    missing_percentage = 0
    if len(idx_not_in_model) > 0:
        weights_missing = numpy.array([document.weights[i] for i in idx_not_in_model])
        weights_missing_raised = numpy.power(weights_missing, intensity_weighting_power)
        missing_percentage = 100 * weights_missing_raised.sum() / (weights_raised.sum()
                                                                       + weights_missing_raised.sum())
    return missing_percentage

def compute_sim_score_by_spec2vec_filter(model,spec1,spec2):
    #temp = spectrums[MS2_ID_1]

    documents = [SpectrumDocument(s) for s in [spec1, spec2]]
    
    pre_1 = check_missing_percentage(documents[0],1.0/8)
    pre_2 = check_missing_percentage(documents[0],1.0/8)
    if pre_1<=5 and pre_2<= 5:
        spec2vec_similarity = Spec2Vec(model=model, intensity_weighting_power=1.0/8,
                               allowed_missing_percentage=100.0)
        sim_score = spec2vec_similarity.pair(documents[0], documents[1])
    else:
        sim_score = -1

    
    return sim_score

def compute_sim_score_by_spec2vec(model,spec1,spec2):
    #temp = spectrums[MS2_ID_1]
    #spec1 = spectrums[MS2_ID_1]
    #spec2_2 = spectrums[MS2_ID_2]
    
    documents = [SpectrumDocument(s) for s in [spec1, spec2]]
    # documents = [SpectrumDocument(s) for s in [ spec2]]
    spec2vec_similarity = Spec2Vec(model=model, intensity_weighting_power=1.0/8,
                               allowed_missing_percentage=100.0)
    

    sim_score = spec2vec_similarity.pair(documents[0], documents[1])

    
    # query_vector = calc_vector(model,  documents[0], intensity_weighting_power=1.0/8,
    #                                 allowed_missing_percentage=100.0)

    return sim_score
#MS2id
def seekMS2id(weight,rt,MS2weight,MS2rt):
    
    candidate_MS2id=[]
    MS2List=[]
    for i in range(len(weight)):
        candidate_MS2id.append([])
    
    delata_weight=MS2_mz                    
    last_index=0
    
    for index_weight in range(len(weight)):
        knownmz=weight[index_weight]
        knownrt=rt[index_weight]*60
       
        Time=0
        index=last_index
        if MS2weight[index]-knownmz>delata_weight:
            continue
        for index_MS2weight in range(index,len(MS2weight)+1):
            if index_MS2weight==len(MS2weight):
                break
            if MS2weight[index_MS2weight]-knownmz>delata_weight:
                break
            if abs(knownmz-MS2weight[index_MS2weight])<=delata_weight:
                Time+=1
                if abs(knownrt-MS2rt[index_MS2weight])<MS2_rt:             
                    candidate_MS2id[index_weight].append(index_MS2weight)
                if Time==1:
                    last_index=index_MS2weight
                    
        if last_index==len(MS2weight) and knownmz-MS2weight[last_index]>delata_weight:
            break
    
    for i in range(len(candidate_MS2id)):
        
        error=[]
        if len(candidate_MS2id[i]):
            for j in candidate_MS2id[i]:
                error.append(abs(weight[i]-MS2weight[j])+abs(rt[i]-MS2rt[j]))
            MS2List.append(candidate_MS2id[i][error.index(min(error))])
        else:
            MS2List.append(-1)
    
    return MS2List

def MS2toList(mgf_name):
    MS2List=[]
    with open(mgf_name)as f:
        for line in f:
            line=line.strip()
            line=re.split('=| ',line)
            if line[0]=='BEGIN':
                oneMS2List=[0,0]
                spec=[]
            elif line[0]=='TITLE':
                continue
            elif line[0]=='RTINSECONDS':
                oneMS2List[0]=float(line[1])
            elif line[0]=='PEPMASS':
                oneMS2List[1]=float(line[1])
            elif line[0]=='END':
                oneMS2List.append(spec)
                MS2List.append(oneMS2List)
            elif line[0]=='CHARGE':
                continue
            else:
                mz_intensity=[]
                mz_intensity.append(float(line[0]))
                mz_intensity.append(float(line[1]))
                spec.append(mz_intensity)
    return MS2List

#Background database
smiles=[]
molecule_HMDB=[]
molecule_weight=[]
molecule_formula=[]
molecule_name=[]
molecule_rt=[]
print("\nMetabolite annotation")
if ion_mode==1:
    read_background=exe_adr + "\data/Background database_POS.csv"
else:
    read_background=exe_adr + "\data/Background database_NEG.csv"
with open(read_background)as f:
    for line in islice(f,1,None):
        line=line.strip().split(',')
        molecule_HMDB.append(line[1])
        molecule_name.append(line[6:])
        molecule_formula.append(line[5])
        molecule_weight.append(float(line[3])) 
        molecule_rt.append(float(line[4]))#Exact weight
        smiles.append(line[2])
#the MS/MS spectra
print("Reading the MS/MS spectra...")
MS2List=[]
spectrums=[]
if ion_mode==1:
    msms_adr=exe_adr + "\\data\\MSMS\\POS"
else:
    msms_adr=exe_adr + "\\data\\MSMS\\NEG"
msms_files=os.listdir(msms_adr) 
msms_num=len(msms_files) 
for i in range(msms_num):
    mgf_name=msms_adr + '/MSMS' + str(i+1) + '.mgf'
    MS2List1=MS2toList(mgf_name)
    MS2List=MS2List+MS2List1
    spectrums_1 = load_from_mgf(mgf_name)
    spectrums_1 = [apply_my_filters(s) for s in spectrums_1]
    spectrums = spectrums+spectrums_1

if ion_mode==1:
    model_file = exe_adr + "\data\spec2vec model/spec2vec_POS_uniqueinchikey_ratio05_filtered_iter_8.model"
else:
    model_file = exe_adr + "\data\spec2vec model/spec2vec_NEG_uniqueinchikey_ratio05_filtered_iter_32.model"

model = gensim.models.Word2Vec.load(model_file)
MS2List_backups = MS2List.copy()
MS2List_mz_rank=[]

for rank_MS2List in MS2List_backups:
    MS2List_mz_rank.append(rank_MS2List[1])

rank_id = np.argsort(MS2List_mz_rank)   

MS2List=sorted(MS2List,key=lambda x:x[1])
spec2=[]
for i in range(len(MS2List)):
    pm=MS2List[i][1]
    spec=MS2List[i][2]
    spec=filter_precursor_peaks(spec, 17, pm)
    spec=window_filter_peaks(spec, 50, 6)
    spec2.append(spec)
MS2List_rt=[]
MS2List_mz=[]
for one_MS2List in MS2List:
    MS2List_rt.append(one_MS2List[0])
    MS2List_mz.append(one_MS2List[1])

initial_molecule_HMDB=[]    
initial_molecule_id=[]     
initial_molecule_MS2id=[]
initial_molecule_weight=[]
initial_molecule_rt=[]
initial_molecule_level=[]
#the seed metabolites
print("Reading the seed metabolites...")
if ion_mode==1:
    seed_adr=exe_adr + "\data\Seed metabolites_POS.csv"
else:
    seed_adr=exe_adr + "\data\Seed metabolites_NEG.csv"
with open(seed_adr)as f:
    for line in islice(f,1,None):
        line=line.strip().split(',')
        initial_molecule_HMDB.append(line[3])
        initial_molecule_weight.append(float(line[1]))
        initial_molecule_rt.append(float(line[2]))
        initial_molecule_level.append(0)
for molecule in initial_molecule_HMDB:
    try:
        initial_molecule_id.append(molecule_HMDB.index(molecule))
    except:
        initial_molecule_id.append(-1)
#MS2id
initial_molecule_MS2id_1=seekMS2id(initial_molecule_weight,initial_molecule_rt,MS2List_mz,MS2List_rt)
initial_molecule_MS2id=[]
for i in initial_molecule_MS2id_1:
    initial_molecule_MS2id.append([i])
#the peak table    
unknown_weight=[]
unknown_rt=[]
unknown_MS2id=[]
unknown_w=[]
unknown_calibrate_rt=[]
print("Reading the peak table...")
if ion_mode==1:
    read_unknown=exe_adr + "\data\Peak table_POS.csv"
else:
    read_unknown=exe_adr + "\data\Peak table_NEG.csv"
with open(read_unknown)as f:
    for line in islice(f,1,None):
        line=line.strip().split(',')
        unknown_weight.append(float(line[1]))
        unknown_rt.append(float(line[2]))
        unknown_w.append(float(line[3]))
        unknown_calibrate_rt.append(float(line[4]))
        
#MS2id       
unknown_MS2id=seekMS2id(unknown_weight,unknown_rt,MS2List_mz,MS2List_rt)


candidate_HMDB=[]
candidate_Grade=[]
candidate_Sim1=[]
candidate_Sim2=[]
candidate_Name=[]
candidate_Formula=[]
candidate_Adduct=[]
candidate_Seed=[]
candidate_level=[]
candidate_score=[]
for i in range(len(unknown_weight)):
    candidate_HMDB.append([])
    candidate_Grade.append([])
    candidate_Sim1.append([])    
    candidate_Sim2.append([])
    candidate_Name.append([])
    candidate_Formula.append([])
    candidate_Adduct.append([])    
    candidate_Seed.append([])
    candidate_level.append([])
    candidate_score.append([])
node_adr=exe_adr + "\data/Nodes"
edge_adr=exe_adr + "\data/Edges"
for index_adduct in range(1):
    known_molecule_HMDB=initial_molecule_HMDB.copy()   
    known_molecule_id=initial_molecule_id.copy()
    known_molecule_MS2id=initial_molecule_MS2id.copy()
    known_molecule_weight=initial_molecule_weight.copy()
    known_molecule_rt=initial_molecule_rt.copy()
    known_molecule_level=initial_molecule_level.copy()
  
    known_molecule_index=-1
    for one_molecule_id in known_molecule_id:
        if one_molecule_id==-1:
            known_molecule_index+=1
            continue
        else:
            c=node_adr+'//' + str(molecule_HMDB[one_molecule_id]) + '.txt'
            d=edge_adr+'//' + str(molecule_HMDB[one_molecule_id]) + '.txt'
            
            with open(d) as f:
                line=f.readline()
                fp_score=line.strip().split(' ')
                
                
            with open(c) as f:
                line=f.readline()
                one_fps=line.strip().split(' ')
            
            known_molecule_index+=1
            print(known_molecule_index)
            #seed_weight=known_molecule_weight[known_molecule_index]
                        
             
            for i in range(len(one_fps)):
                one_fps[i]=int(one_fps[i])
            repeat=[]
            for each_molecule_id in known_molecule_id:
                if each_molecule_id in one_fps:
                    repeat.append(one_fps.index(each_molecule_id))
            repeat=sorted(repeat,reverse=True)
            for repeat_index in repeat:
                del(fp_score[repeat_index])
                del(one_fps[repeat_index])
                #print(repeat)
                    
            score_index=-1
            for neighbor in one_fps:
                neighbor=int(neighbor)
                score_index+=1
                sim1=float(fp_score[score_index])
                neighbor_weight=molecule_weight[neighbor]
                neighbor_rt=molecule_rt[neighbor]
                candidate_unknown_id=[]
                candidate_unknown_similarity=[] 
                unknown_molecule_count=0
                score_sum=[]
                
                for unknown_molecule_weight in unknown_w: 
                    SD_mz=abs(float(unknown_molecule_weight) - float(neighbor_weight))/neighbor_weight*1000000
                    SD_rt=abs(float(unknown_calibrate_rt[unknown_molecule_count])-float(neighbor_rt))/neighbor_rt*100
                    if SD_mz<=SD_mz_threshold and SD_rt<=SD_rt_threshold :
                        if unknown_MS2id[unknown_molecule_count]!=-1 and known_molecule_MS2id[known_molecule_index][0]!=-1: 
                            
                            multi_similarity=[]
                            for MS2_ID_1 in known_molecule_MS2id[known_molecule_index]:
                                pm1=MS2List[MS2_ID_1][1]
        
                                MS2_ID_2=unknown_MS2id[unknown_molecule_count]
                                pm2=MS2List[MS2_ID_2][1]
                                
                                if check_words(model,spectrums[rank_id[MS2_ID_1]]) and check_words(model,spectrums[rank_id[MS2_ID_2]]):
                                    # print('TRUE')
                                    #MS2_similarity,reported_alignments=Cosine.score_alignment(spec2[MS2_ID_1],spec2[MS2_ID_2],pm1,pm2,0.01,max_charge_consideration=1) #二级质谱对齐
                                    MS2_similarity = compute_sim_score_by_spec2vec(model,spectrums[rank_id[MS2_ID_1]],spectrums[rank_id[MS2_ID_2]])
                                else:
                                    MS2_similarity = 0
                                    print(MS2_similarity)
                                
                                multi_similarity.append(MS2_similarity)
                            similarity=max(multi_similarity)
                            SD_mz_score=1-SD_mz/SD_mz_threshold
                            SD_rt_score=1-SD_rt/SD_rt_threshold
                            sim1_score=sim1
                            sim2_score=similarity
                            #seed_score=SD_mz_score*0.4+SD_rt_score*0.4+sim2_score*0.2
                            score=SD_mz_score*0.2+SD_rt_score*0.2+sim1_score*0.3+sim2_score*0.3
                            
                            if(similarity>=sim2_threshold):  
                                score_sum.append(score)
                                candidate_unknown_id.append(unknown_molecule_count)
                                
                                candidate_HMDB[unknown_molecule_count].append(molecule_HMDB[neighbor])
                                candidate_Grade[unknown_molecule_count].append(1)
                                candidate_Sim1[unknown_molecule_count].append(sim1)
                                candidate_Sim2[unknown_molecule_count].append(similarity)
                                candidate_Name[unknown_molecule_count].append(molecule_name[neighbor])
                                candidate_Formula[unknown_molecule_count].append(molecule_formula[neighbor])
                                #candidate_Adduct[unknown_molecule_count].append(adduct_type[index_adduct])
                                candidate_score[unknown_molecule_count].append(score)
                                
                                HMDB_name=molecule_HMDB[known_molecule_id[known_molecule_index]]
                                candidate_Seed[unknown_molecule_count].append(HMDB_name)
                                level=known_molecule_level[known_molecule_index]+1
                                candidate_level[unknown_molecule_count].append(level)
                                unknown_molecule_count+=1
                            elif (similarity<sim2_threshold ):  #and seed_score<0.6
                                if molecule_HMDB[neighbor] in candidate_HMDB[unknown_molecule_count]:
                                    unknown_molecule_count+=1
                                    continue
                                else:
                                    candidate_HMDB[unknown_molecule_count].append(molecule_HMDB[neighbor])
                                    candidate_Grade[unknown_molecule_count].append(2)
                                    candidate_Sim1[unknown_molecule_count].append(sim1)
                                    candidate_Sim2[unknown_molecule_count].append(similarity)
                                    candidate_score[unknown_molecule_count].append(score)
                                    candidate_Name[unknown_molecule_count].append(molecule_name[neighbor])
                                    candidate_Formula[unknown_molecule_count].append(molecule_formula[neighbor])
                                    #candidate_Adduct[unknown_molecule_count].append(adduct_type[index_adduct])
                
                                    HMDB_name=molecule_HMDB[known_molecule_id[known_molecule_index]]
                                    candidate_Seed[unknown_molecule_count].append(HMDB_name)
                                    
                                    level=known_molecule_level[known_molecule_index]+1
                                    candidate_level[unknown_molecule_count].append(level)
                                    
                                    unknown_molecule_count+=1
                        else:
                            if molecule_HMDB[neighbor] in candidate_HMDB[unknown_molecule_count]:
                                unknown_molecule_count+=1
                                continue
                            else:
                                candidate_HMDB[unknown_molecule_count].append(molecule_HMDB[neighbor])
                                candidate_Grade[unknown_molecule_count].append(3)
                                candidate_Sim1[unknown_molecule_count].append(sim1)
                                candidate_Name[unknown_molecule_count].append(molecule_name[neighbor])
                                candidate_Formula[unknown_molecule_count].append(molecule_formula[neighbor])
                                #candidate_Adduct[unknown_molecule_count].append(adduct_type[index_adduct])
                                
                                HMDB_name=molecule_HMDB[known_molecule_id[known_molecule_index]]
                                candidate_Seed[unknown_molecule_count].append(HMDB_name)
                                
                                level=known_molecule_level[known_molecule_index]+1
                                candidate_level[unknown_molecule_count].append(level)
                                unknown_molecule_count+=1
                            
                    else:
                        unknown_molecule_count+=1
                #如果candidate_unknown_id 非空
                if len(candidate_unknown_id):
                    if neighbor in known_molecule_id:
                        continue
                    else:
                            known_molecule_id.append(neighbor) #将该代谢物加入已知代谢物中
                            new_MS2id=[]
                            for i in candidate_unknown_id:
                                new_MS2id.append(unknown_MS2id[i])
                            known_molecule_MS2id.append(new_MS2id)
                            level=known_molecule_level[known_molecule_index]+1
                            known_molecule_level.append(level)

   
#Outputting result
list_result_one =[]
list_result = []
result_name1 = ['ID','m/z','RT','ID','Name','Grade','Formula','Sim1','Sim2','Adduct','MS2_ID','Match_Seed','Level','Score']

i=-1
for identity_HMDB in candidate_HMDB:
    list_result_one =['-','-','-','-','-','-','-','-','-','-','-','-','-','-']
    
    i+=1
    
    if len(identity_HMDB):
        k=-1
        l=-1
        for one_indentity_HMDB in identity_HMDB:
            k+=1
            
            list_result_one[0]=i
            list_result_one[1]=unknown_weight[i]
            list_result_one[2]=unknown_rt[i]
            list_result_one[3]=candidate_HMDB[i][k]
            list_result_one[4]=candidate_Name[i][k]
            list_result_one[5]=candidate_Grade[i][k]
            
 
            if candidate_Grade[i][k]<=2:
                l+=1
            list_result_one[6]=candidate_Formula[i][k]
            list_result_one[7]=candidate_Sim1[i][k]
            if l>=0:
                list_result_one[8]=round(candidate_Sim2[i][l],2)
                list_result_one[13]=round(candidate_score[i][l],2)
                
            # list_result_one[9]=candidate_Adduct[i][k]
            list_result_one[9]='-'
            list_result_one[10]=unknown_MS2id[i]
            list_result_one[11]=candidate_Seed[i][k]
            list_result_one[12]=candidate_level[i][k]
            list_b = copy.deepcopy(list_result_one)
            list_result.append(list_b)
            
    else:
            
            list_result_one[0]= i
            list_result_one[1]=unknown_weight[i]
            list_result_one[2]=unknown_rt[i]
            # list_result_one[3]=candidate_HMDB[i]
            # list_result_one[4]=candidate_Name[i]
            # list_result_one[5]=candidate_Grade[i]
            # list_result_one[6]=candidate_Formula[i]
            # list_result_one[7]=candidate_Sim1[i]
            # list_result_one[9]=candidate_Adduct[i]
            list_result_one[10]=unknown_MS2id[i]
            # list_result_one[11]=candidate_Seed[i]
            # list_result_one[12]=candidate_level[i]
            list_b = copy.deepcopy(list_result_one)
            list_result.append(list_b)

resut_csv = pd.DataFrame(columns=result_name1, data=list_result)
print("Outputting result...")
if ion_mode==1:
    write_result = exe_adr + "\\" + "result" + "\\" + result_name + "_POS.csv"
else:
    write_result = exe_adr + "\\" + "result" + "\\" + result_name + "_NEG.csv"
resut_csv.to_csv(write_result, index=None)
print("Finished.")
input('Press Enter to exit...')