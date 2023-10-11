#!/usr/bin/env python

import os

PDB_PATH = '/home/m.pak/db/pdb/'
pdb_file_path = lambda pdb: os.path.join(PDB_PATH, 'pdb'+pdb.lower()+'.pdb') #PDB_PATH/pdb1hho.pdb

CIF_PATH = '/home/m.pak/db/pdb_mmcif/mmcif_files/'
cif_file_path = lambda pdb: os.path.join(CIF_PATH, pdb.lower()+'.cif') #PDB_PATH/1hho.cif
