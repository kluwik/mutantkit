#!/usr/bin/env python

import pandas as pd
import numpy as np
import os, warnings, re, subprocess

from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import three_to_one, one_to_three, is_aa
p = PDBParser(QUIET=True)
from Bio import SeqIO
from Bio import Entrez

import urllib.request
from urllib.error import HTTPError
import ast

from mutantkit.config import pdb_file_path, cif_file_path, downloadPdb

def parseUrl(url):
    try:
        data = urllib.request.urlopen(url).readlines()
        data = [i.decode('utf-8') for i in data]
        return data
    except HTTPError:
        return []

################################# PDB STRUCTURE #################################

def getSeq(pdb, chain, file_format='pdb'):
    '''
    Returns a tuple of SEQRES sequence and ATOM sequence.
    '''
    def seqParser(file, file_format, seq, chain):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            for i, record in enumerate(SeqIO.parse(file, f'{file_format}-{seq}')):
                if type(chain) == int: # :> перепроверить, правильно ли определятся entity таким образом, кажется  нет (см. 2dr9); если сильно надо, то пользуемся mapEntitytoChain
                    if i == chain - 1:
                        return str(record.seq)
                else:
                    if record.annotations['chain'] == chain:
                        return str(record.seq)
    # :> этот аргумент вообще нужен? может убрать и сделать так, чтобы при отсутствии пдб использовался циф?
    if file_format == 'pdb':
        file = pdb_file_path(pdb)
    elif file_format == 'cif':
        file = cif_file_path(pdb)
    
    # downloads file if absent or empty
    downloadPdb(pdb, file)
    seqres = seqParser(file, file_format, 'seqres', chain)
    atom = seqParser(file, file_format, 'atom', chain)
    return seqres, atom

def getSeqres(pdb, chain, pdb_seqres_file = '/home/m.pak/db/pdb_seqres.txt'):
    '''
    Extract seqres from pdb_seqres.txt file (https://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz)
    Медленная. 
    '''
    pdbid = f'{pdb}_{chain}'
    line = subprocess.check_output(f'grep -n {pdbid} {pdb_seqres_file}', shell=True).decode('utf-8')
    line = int(line.split(':')[0]) + 1
    seqres = subprocess.check_output(f"sed -n '{line}p' {pdb_seqres_file}", shell=True).decode('utf-8')
    seqres = seqres.strip()
    return seqres

def parsePDBFasta(pdb, chain=None):
    fasta = parseUrl('https://www.rcsb.org/fasta/entry/'+pdb+'/display')
    if fasta != []:
        chains = [i.split('|')[1].replace('Chain', '') for i in fasta if '>' in i]
        p_o = [i.split('|')[2:] for i in fasta if '>' in i]
        p_o = [j.replace('\n','') for i in p_o for j in i]
        protein = [p_o[i] for i in range(0, len(p_o), 2)]
        organism = [' '.join(p_o[i].split()[:-1]) for i in range(1, len(p_o), 2)] # delete numbers in brackets after organism name
        if chain != None:
            index = [i for i, c in enumerate(chains) if chain in c][0]
            return protein[index], organism[index]
        return [(protein[i], organism[i]) for i in range(len(protein))]
    else:
        return None, None

def mapEntitytoChain(pdb, entity=None):
    '''
    Maps entity to chains: takes entity (integer) and returns the corresponding chain. Parses fasta header from PDB website.
    Полноценно не дебажила.
    '''
    fasta = parseUrl('https://www.rcsb.org/fasta/entry/'+pdb+'/display')
    entity_dict = {}
    for i in range(0, len(fasta), 2):
        header = fasta[i].split('|')
        e = header[0].split('_')[-1]
        chain = header[1].replace('Chain ', '').replace('Chains ', '').split(',')[0] # if several chains for one entity, leave the first one
        if 'auth' in chain:
            chain = chain.split('[auth ')[-1].split(']')[0]
        entity_dict[e] = chain
    if entity != None:
        return entity_dict[str(entity)]
    else:
        return entity_dict
    
def getResolution(pdb, download_cif = False):
    '''
    Parse resolution from local cif file. If not available returns None.
    '''
    def parseResFromCif(cif_file):
        with open(cif_file) as f:
            cif_data = f.read()
        if '_refine.ls_d_res_high' in cif_data:
            try: # может не быть (NMR) или бывает точка вместо числа (scattering)
                res = [i for i in cif_data.split('\n') if '_refine.ls_d_res_high' in i][0]
                res = float(res.split()[-1])
                return res
            except:
                print(pdb)
    
    cif_file = cif_file_path(pdb)
    if os.path.isfile(cif_file):
        return parseResFromCif(cif_file)
    else:
        if download_cif == True:
            os.system(f'wget https://files.rcsb.org/download/{pdb.lower()}.cif -O {cif_file}')
            return parseResFromCif(cif_file)    

        
def getAminoacidLigands(pdb, chain=None):
    '''
    Returns a list of tuples with amino acid ligands (ligand name, residue number, chain)
    !!! Biopython 1.81 is required! is_aa() doesn't work properly for older versions.
    '''
    def extract_section(cif_data, section_name):
        hits = [i for i in cif_data if section_name in i] # ideally returns a list with a single item that is a query section
        # but sometimes section name can be mentioned in other comment-like sections, e.g. "audit"
        # thus, several sections are returned: a query and those, where query is mentioned
        if len(hits) > 1: 
            for i in hits: 
                # each line in the query section starts with its name, that's how it can be distinguished from others
                # but firstly, we need to remove 'loop_\n' in the beginning of the section
                # it is present if section is organised as a table
                item = i.replace('loop_\n', '')
                # now we check that the first line of the section starts with its name
                if item.split('\n')[0].startswith(section_name):
                    return i
        elif hits != []:
            return hits[0]
        else:
            return hits

    def parse_section(section, section_name):
        if 'loop_' in section:
            section = section.replace('loop_\n', '').strip().split('\n')
            data = [i for i in section if i.startswith(section_name) == False]
            return [(i.split()[2], i.split()[4], i.split()[-2]) for i in data] # name residue_number auth_chain
        else:
            return section.split('\n')[2].split()[1], section.split('\n')[4].split()[1], section.split('\n')[-2].split()[1]
    
    cif_file = cif_file_path(pdb)
    downloadPdb(pdb, cif_file)
    with open(cif_file) as f:
        cif_data = f.read()
    cif_data = cif_data.split('# \n')
    
    hetnam = extract_section(cif_data, '_pdbx_nonpoly_scheme') # аналогично можно modres = extract_section(cif_data, '_pdbx_struct_mod_residue')
    if hetnam != []:
        het_residues = parse_section(hetnam, '_pdbx_nonpoly_scheme')
        het_residues = [i for i in het_residues if is_aa(i[0], standard=False)]
        het_residues = list(set(het_residues))
        if chain:
            het_residues = [i for i in het_residues if i[-1] == chain]
        if het_residues != []:
            return het_residues
        
################################# CROSS-REFERENCES #################################

def getUniprotTxt(uni):
    url = 'https://rest.uniprot.org/uniprotkb/'+str(uni)+'.txt'
    data = parseUrl(url)
    return data

def getUniprotFasta(uniprot):
    url = 'https://rest.uniprot.org/uniprotkb/'+uniprot+'.fasta'
    data = parseUrl(url)
    seq = data[1:]
    seq = ''.join([i.strip() for i in seq])
    name = organism = protein = gene = None
    name = data[0].split('|')[2].split(' ', 1)[0]
    if 'OS=' in data[0]:
        protein = data[0].split('|')[2].split(' ', 1)[1].split(' OS=')[0]
    if 'OX=' in data[0]:        
        organism = data[0].split(' OS=')[1].split(' OX=')[0]
    if 'GN=' in data[0]:
        gene = data[0].split('GN=')[1].split(' PE=')[0]
    return {'data' : data, 'seq' : seq, 'name' : name, 'organism' : organism, 'protein' : protein, 'gene' : gene}
#getUniprotFasta('P00044')

def getUniprotData(uniprot):
    url = 'https://rest.uniprot.org/uniprotkb/search/?query=(accession:'+uniprot+')'
    data = parseUrl(url)
    data = data[0].replace(':false', ':False').replace(':true', ':True').replace(':null', ':None')
    data = ast.literal_eval(data)['results'][0]
    seq = data['sequence']['value']
    name = organism = protein = gene = None
    if 'uniProtkbId' in data.keys():
        name = data['uniProtkbId']
    if 'organism' in data.keys():
        organism = data['organism']['scientificName']
    if 'proteinDescription' in data.keys():
        protein = data['proteinDescription']['recommendedName']['fullName']['value']
    if 'genes' in data.keys():
        gene = data['genes'][0]['geneName']['value']
    # :> переписать, чтобы возвращался объект с атрибутами
    return {'data' : data, 'sequence' : seq, 'name' : name , 'organism' : organism, 'protein' : protein, 'gene' : gene}

def getUniprotFromEntrez(gi, entrez_email):
    Entrez.email = entrez_email
    data = Entrez.efetch(db="protein", id=str(gi), rettype="gb")
    data = data.read().split('\n')
    uni = [i.split()[1] for i in data if 'ACCESSION' in i][0]
    return uni
#getUniFromEntrez('124267')

def getUniprotFromLocalPDB(pdb, chain=None):
    file = pdb_file_path(pdb)
    record = [i for i in SeqIO.parse(file, 'pdb-seqres')]

    dbxrefs = None
    if chain == None:
        dbxrefs = [i.dbxrefs for i in record]
    else:
        if chain in [i.annotations['chain'] for i in record]:
            dbxrefs = [i.dbxrefs for i in record if i.annotations['chain'] == chain]

    refs_list = None
    if dbxrefs != None:
        refs_list = [j.replace('UNP:', '') for i in dbxrefs for j in i if 'UNP:' in j]
        if refs_list == []:
            refs_list = [j for i in dbxrefs for j in i]
        refs_list = [refs_list[i] for i in range(0, len(refs_list), 2)]
    return refs_list
#print(getUniFromLocalPDB('1A0N'))
#print(getUniFromLocalPDB('1A0N', 'B'))

def getUniprotForPdb(pdb, chain=None):
    pdb = pdb.upper()
    url = 'https://rest.uniprot.org/uniprotkb/search?query=(xref:pdb-'+pdb+')'
    data = parseUrl(url)
    data = data[0].replace(':false', ':False').replace(':true', ':True').replace(':null', ':None')
    results = ast.literal_eval(data)['results']
    if results == []:
        print(f'No Uniprot found for {pdb} :(')
        return None 
    else:
        if chain != None:
            chains = []
            for result in results:
                crossrefs = result['uniProtKBCrossReferences']
                pdb_crossref = [crossref for crossref in crossrefs if crossref['database'] == 'PDB' and crossref['id'] == pdb][0]
                chain_region = [i['value'] for i in pdb_crossref['properties'] if i['key'] == 'Chains'][0]
                chains.append(chain_region.split('=')[0])
            index = [i for i,c in enumerate(chains) if chain in c]
            if index != []:
                uni = results[index[0]]['primaryAccession']
            else:
                print(f'Chain {chain} not found in {pdb}. Available chains: {chains}. Returning Uniprots for all chains.')
                uni = [result['primaryAccession'] for result in results]
                if len(uni) == 1:
                    uni = uni[0]
        else:
            uni = [result['primaryAccession'] for result in results]
        return uni
    
    
#### notes #####
'''
по такой ссылке последовательность выводится с нестандартными остатками в скобках https://www.rcsb.org/fasta/chain/1brz.A/display
'''
