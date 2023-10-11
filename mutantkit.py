#!/usr/bin/env python

import pandas as pd
import numpy as np
import pickle, os, glob, warnings, re

import matplotlib.pyplot as plt
import seaborn as sns

import math
from sklearn.metrics import mean_squared_error
from scipy.stats import mannwhitneyu, pearsonr, spearmanr

from Bio.PDB import PDBParser
from Bio.PDB.DSSP import *
from Bio.PDB.Polypeptide import three_to_one, one_to_three
p = PDBParser(QUIET=True)
from Bio import SeqIO, PDB

aa = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']


def splitMutation(dataframe, mutation_column):
    df = dataframe.copy()
    df['wt_aa'] = df[mutation_column].apply(lambda x: x[0])
    df['position'] = df[mutation_column].apply(lambda x: x[1:-1]).astype('int64')
    df['mut_aa'] = df[mutation_column].apply(lambda x: x[-1])
    return df

def getDDG(data, y_ddg, outliers=False, outliers_threshold=None):
    
    def remove_outliers(df, y_ddg):
        Q1 = df[y_ddg].quantile(0.25)
        Q3 = df[y_ddg].quantile(0.75)
        IQR = Q3-Q1
        df_trim = df[~((df[y_ddg]<(Q1-1.5*IQR)) | (df[y_ddg]>(Q3+1.5*IQR)))]
        return df_trim
    
    df = data[data[y_ddg].notna()] # remove NaN ddg (not computed yet)
    
    if outliers == True:
        if outliers_threshold != None:
            # trim outliers if |ddg| > outliers_threshold
            if df[y_ddg].max() > outliers_threshold or df[y_ddg].min() < -outliers_threshold: 
                df = remove_outliers(df, y_ddg)
        else:
            df = remove_outliers(df, y_ddg)
    
    return df

def getStat(x,y, prnt=False):
    pcc = pearsonr(x,y)[0]
    sp = spearmanr(x,y)[0]
    p_pcc = pearsonr(x,y)[1]
    se = math.sqrt((1 - pcc**2)/(len(x)-2))
    sp_se = math.sqrt((1 - sp**2)/(len(x)-2))
    mse = mean_squared_error(x,y, squared=False)
    mwu = mannwhitneyu(x,y)[1]
    if prnt == True:
        print(f'PCC = {pcc:.2f}±{se:.2f}, SCC = {sp:.2f}±{sp_se:.2f}, MSE = {mse:.2f}')
    return round(pcc,2), round(se,2), round(sp,2), round(sp_se,2), round(mse,2)
             
ddg_class_2 = lambda x: 'D' if x >= 0 else 'S'    
ddg_class_3 = lambda x: 'D' if x >= 0.5 else ('S' if x <= -0.5 else 'N') 
    
def removeOutliers(df, y_ddg):
    Q1 = df[y_ddg].quantile(0.25)
    Q3 = df[y_ddg].quantile(0.75)
    IQR = Q3-Q1
    df_trim = df[~((df[y_ddg]<(Q1-1.5*IQR)) | (df[y_ddg]>(Q3+1.5*IQR)))]
    return df_trim

def getAntisymmetry(f, r):
    pcc =  getStat(f, r)[0]
    b = np.sum((f + r)) / (2*len(f))
    return round(pcc,2), round(b,2)

def getMutationProps(mutation):
    # does not distinguish between C(SH) and C(SS)
    # charge
    negative = ['D', 'E']
    positive = ['R', 'K', 'H']
    nocharge = [i for i in aa if i not in positive+negative]
    # hydrophobicity
    hphob = ['F','W','I','L','M','V','C','A']
    hphil = ['K','R','E','D','Q','N']
    neutral = ['Y','P','G','T','S','H']
    # polarity
    nonpolar=hphob+['P','G']
    polar = [i for i in aa if i not in nonpolar]
    # size
    XS = ['S', 'G', 'A']
    S = ['N', 'D', 'T', 'P', 'C']
    M = ['E', 'V', 'Q', 'H']
    L = ['M','I','L','K','R']
    XL = ['F', 'Y', 'W']
    
    wt = mutation[0]
    mut = mutation[-1]
    
    def getFeature(wt, mut, feature):
        for label, group in feature.items():
            if wt in group:
                wt_l=label
            if mut in group:
                mut_l=label
        return wt_l+mut_l
        
    chargeF = getFeature(wt, mut, {'-':negative, '+':positive, '0':nocharge})
    hydroF = getFeature(wt, mut, {'-':hphob, '+':hphil, '0':neutral})
    polarF = getFeature(wt, mut, {'n':nonpolar, 'p':polar})
    
    sizes = {1:XS,2:S,3:M,4:L,5:XL}
    for label, group in sizes.items():
        if wt in group:
            wt_size = label
        if mut in group:
            mut_size = label
    sizeF = mut_size - wt_size
    
    change = lambda x: 0 if len(set(x)) == 1 else 1
    chargeChange = change(chargeF)
    hydroChange = change(hydroF)
    polarChange = change(polarF)
    size_changes = {0:0, -1:'SD', -2:'SD', -3:'LD', -4:'LD', 1:'SI', 2:'SI', 3:'LI', 4:'LI'}
    sizeChange = size_changes[sizeF]
    
    return (chargeF,chargeChange), (hydroF,hydroChange), (polarF,polarChange), (sizeF,sizeChange)

def getSeq(pdb, chain):
    '''
    Returns a tuple of SEQRES sequence and ATOM sequence.
    '''
    file = '/home/m.pak/db/pdb/pdb'+pdb.lower()+'.pdb'
    if os.path.isfile(file):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            for i in SeqIO.parse(file, 'pdb-seqres'):
                if i.annotations['chain'] == chain:
                    seqres = str(i.seq)
            for i in SeqIO.parse(file, 'pdb-atom'):
                if i.annotations['chain'] == chain:
                    atom = str(i.seq)
    return seqres, atom

def getResidueIndex(pdb, chain, mutation):
    '''
    Returns index of the mutated residue and None. 
    If the wt residue doesn't match that in the structure, returns index and 0. 
    '''
    match = None
    try:
        filename = '/home/m.pak/db/pdb/pdb'+pdb.lower()+'.pdb'
        if os.path.isfile(filename): 
            structure = p.get_structure(' ', filename)
            structure = [j for j in structure.get_residues() if j.get_full_id()[3][0] == ' ' and j.get_full_id()[2] == chain]
        for i, res in enumerate(structure):
            num = res.get_full_id()[3][1]
            if num == int(mutation[1:-1]):
                residue = three_to_one(res.get_resname())
                if residue != mutation[0]:
                    match = 0
                    #print(pdb, "Mutated amino acid doesn't match. Residue", int(mutation[1:-1]), 'is', residue)
                return i, match
    except:
        return None

def mutations2020(mutations):
    
    '''
    Returns an array 20x20 of mutations counts (x-axis - mutation from. y-axis - mutation to).
    '''
    
    mutations = [i[0]+i[-1] for i in mutations]
    mutations = {pair:mutations.count(pair) for pair in set(mutations)}
    
    aa_pairs = [a+b for idx, a in enumerate(aa) for b in aa]

    aa2020 = {k:v[k] if k in mutations.keys() else 0 for k in aa_pairs}
    aa2020 = np.array(list(aa2020.values())).reshape((20,20))
    return aa2020

def checkPosition(mutation, chain, pdb='', path=''):
    
    '''
    Returns True if the specified WT residue in the specified position matches that in the specified structure.
    Otherwise returns False.
    '''
    
    wt = mutation[0]
    num = int(mutation[1:-1])
    
    def check(file, wt, num, chain):
        structure = p.get_structure(' ', file)
        structure = [j for j in structure.get_residues() if j.get_full_id()[3][0] == ' ' and j.get_full_id()[2] == chain]
        for i, res in enumerate(structure):
            resnum = res.get_full_id()[3][1]
            #print(res.get_full_id())
            if resnum == num:
                residue = three_to_one(res.get_resname())
                if residue == wt:
                    return True
                else:
                    return False
                
    if pdb != '':
        file = '/home/m.pak/db/pdb/pdb'+pdb.lower()+'.pdb'
        if os.path.isfile(file): 
            return check(file, wt, num, chain)
        else:
            print('No file for', pdb)
    elif path != '':
        if os.path.isfile(path): 
            return check(path, wt, num, chain)
        else:
            print('No such file', path)
            
def getRenumberedPositions(mutations, seq=None):
    
    '''
    Returns a pattern of mutated residues according to their positioins.
    If pattern matches the given sequence, returns the start of the pattrn relative to the sequence 
    and renumbered mutations positions. If not - None, None.
    '''
    
    x = {int(m[1:-1]):m[0] for m in mutations}
    x_src = x
    x = {k: v for k, v in sorted(x.items(), key=lambda item: item[0])}
    
    pattern = []
    for i in range(len(x.keys())):
        aa = list(x.values())[i]
        pattern.append(aa)
        if i < len(x.keys())-1:
            dots = list(x.keys())[i+1] - list(x.keys())[i]
            pattern.append('.'*(dots-1))
    pattern = r''.join(pattern)
    #print(pattern)
    
    start, reindex, reindexed = None, None, None
    if seq != None:
        try:
            start = [m.start(0) for m in re.finditer(pattern, seq)][0] + 1
            reindex = min(x.keys()) - start
            reindexed = [i-reindex for i in x_src.keys()]
        except:
            pass
        
    return pattern, start, reindex, reindexed

#getRenumberedPositions(['Q30A', 'G45A', 'L59A', 'E86A', 'G23A', 'E26A'], \
#           'RVGLTEEQKQEIREAFDLFDTDGSGTIDAKELKVAMRALGFEPKKEEIKKMISEIDKDGSGTIDFEEFLTMMTAKM')
