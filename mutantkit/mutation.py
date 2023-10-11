#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import os, warnings, re

from Bio.PDB import PDBParser
p = PDBParser(QUIET=True)
from Bio.PDB.Polypeptide import three_to_one, one_to_three
from Bio import SeqIO

aa = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

################################# MUTATIONS STATISTICS #################################

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

#getMutationProps('D1K')

def mutations2020(mutations):
    
    '''
    Returns an array 20x20 of mutations counts (x-axis - mutation from. y-axis - mutation to).
    '''
    
    mutations = [i[0]+i[-1] for i in mutations]
    mutations = {pair:mutations.count(pair) for pair in set(mutations)}
    
    aa_pairs = [a+b for idx, a in enumerate(aa) for b in aa]

    aa2020 = {k:mutations[k] if k in mutations.keys() else 0 for k in aa_pairs}
    aa2020 = np.array(list(aa2020.values())).reshape((20,20))
    return aa2020

#mutations2020(['A2N', 'E43D', 'W209Y', 'A24N', 'E12A'])

def mutateSequence(mutation, sequence):
    mut_aa = mutation[-1]
    pos = int(mutation[1:-1]) - 1
    return sequence[:pos] + mut_aa + sequence[pos+1:]

################################# MUTATION POSITION RENUMBERING #################################

def getResidueIndex(pdb, chain, mutation):
    '''
    Returns index of the mutated residue and True. 
    If the wt residue doesn't match that in the structure, returns index and False. 
    '''
    match = True
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
                    match = False
                    #print(pdb, "Mutated amino acid doesn't match. Residue", int(mutation[1:-1]), 'is', residue)
                return i, match
    except:
        return None
    
#getResidueIndex('2AMI', 'A', 'A35N')
#getResidueIndex('2AMI', 'A', 'Z35N')

def getResidueNumber(pdb, chain, mutation):
    '''
    Returns residue number as in PDB file given the index of the mutated residue and True. 
    If the wt residue doesn't match that in the structure, returns index and False. 
    '''
    match = True
    position_index = int(mutation[1:-1]) - 1
    try:
        filename = '/home/m.pak/db/pdb/pdb'+pdb.lower()+'.pdb'
        if os.path.isfile(filename): 
            structure = p.get_structure(' ', filename)
            structure = [j for j in structure.get_residues() if j.get_full_id()[3][0] == ' ' and j.get_full_id()[2] == chain]
        for i, res in enumerate(structure):
            if i == position_index:
                position = res.get_full_id()[3][1]
                residue = three_to_one(res.get_resname())
                if residue != mutation[0]:
                    match = False
                    print(pdb, "Mutated amino acid doesn't match. Residue", int(mutation[1:-1]), 'is', residue)
                return position, match
    except:
        return None

def checkPositionInSequence(mutation, sequence):
    '''
    Returns True if WT residue in mutation matches the residue in the mutations position in the sequence. otherwise False.
    '''
    wt = mutation[0]
    pos = int(mutation[1:-1]) - 1
    if sequence[pos] != wt:
        print(f'Mismatch: residue in position {pos+1} is {sequence[pos]}, not {wt}.')
        return False
    else:
        return True

def checkPositionInPDB(mutation, chain, pdb='', path=''): 
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
            
#checkResidueNumber('D48Q', 'A', path='/home/i.vorobiev/ivan/AF_models/1RN1_A/ranked_0.pdb')            

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
