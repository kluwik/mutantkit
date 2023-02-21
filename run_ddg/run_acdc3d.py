# -*- coding: utf-8 -*-
from absl import app
from absl import flags
from absl import logging
import os, sys

logging.set_verbosity(logging.INFO)
FLAGS = flags.FLAGS

flags.DEFINE_string("mutation", None, "Mutation")
flags.DEFINE_string("pdb_file", None, "PDB file")
flags.DEFINE_string("chain", None, "Chain in the PDB file")
flags.DEFINE_string("profile", None, "Directory with profiles")
flags.DEFINE_string("output_tag", None, "Output name")
flags.DEFINE_string("output_dir", os.getcwd(), "Output directory")

# python run_acdc3d.py --mutation=A35N --pdb_file=../pdb/2ami_Repair.pdb --chain=A --output_tag=2AMI_A35N --output_dir=../tmp --profile=../fasta/wt_seq/msa/2AMI.profile

def run_acdc(mutation, pdb_file, chain, profile, output_tag, output_dir): 

    pdb_file = os.path.abspath(pdb_file)
    profile = os.path.abspath(profile)
    
    if os.path.isfile(pdb_file) == False:
        print(f"Oops! PDB file {pdb_file} doesn't exist.")
        sys.exit()
    else:
        # if output_tag not specified make output_tag pdb_file_chain_mutation
        if output_tag == None:
            output_tag = '_'.join(pdb_file.replace('.pdb', ''), chain, mutation)
        
        output_dir = os.path.abspath(output_dir)
        out_file = os.path.join(output_dir, output_tag+'.acdc')
        
        # run 
        cmd = ' '.join(['acdc-nn struct', mutation, profile, pdb_file, chain, '>', out_file])
        logging.info('Computing ddG: '+output_tag)
        logging.info(cmd)
        os.system(cmd)

def main(argv):
    run_acdc(FLAGS.mutation, FLAGS.pdb_file, FLAGS.chain, FLAGS.profile, FLAGS.output_tag, FLAGS.output_dir)

if __name__ == '__main__':
    app.run(main)

