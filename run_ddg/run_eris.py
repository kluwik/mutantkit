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
flags.DEFINE_string("output_tag", None, "Output name")
flags.DEFINE_string("output_dir", os.getcwd(), "Output directory")

flags.mark_flag_as_required("mutation")
flags.mark_flag_as_required("pdb_file")
flags.mark_flag_as_required("chain")

eris='/home/m.pak/bin/eris_standalone-Linux-1.0/bin/eris'

def runEris(mutation, pdb_file, chain, output_tag, output_dir): 
    
    # python run_eris.py --mutation=A35N --pdb_file=../pdb/2ami_Repair.pdb --chain=A --output_tag=2AMI_A35N
    
    pdb_file = os.path.abspath(pdb_file)
    
    if os.path.isfile(pdb_file) == False:
        print(f"Oops! PDB file {pdb_file} doesn't exist.")
        sys.exit()
    else:
        # if output_tag not specified make output_tag: pdb_file_chain_mutation
        if output_tag == None:
            output_tag = '_'.join(pdb_file.replace('.pdb', ''), chain, mutation)
        
        # ДОБАВИТЬ ЗАПИСЬ В УКАЗАННУЮ output_dir
        
        # run 
        cmd = [eris, 'ddg', '-m', mutation, '-f', pdb_file, '-j', output_tag]
        cmd = ' '.join(cmd)
        logging.info('Computing ddG: '+output_tag)
        #print(cmd)
        os.system(cmd)

def main(argv):
    runEris(FLAGS.mutation, FLAGS.pdb_file, FLAGS.chain, FLAGS.output_tag, FLAGS.output_dir)

if __name__ == '__main__':
    app.run(main)

