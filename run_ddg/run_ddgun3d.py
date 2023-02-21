# -*- coding: utf-8 -*-

# python run_ddgun3d.py --mutation=A35N --pdb_file=/home/m.pak/_prj/mega/inputs/structures/wt/pdb/repaired/2ami_Repair.pdb --chain=A --output_tag=2AMI_A35N --output_dir=.

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
    
def mutList(path, output_tag, mutation):
    mut_file = os.path.join(path, output_tag+'.mut')
    with open(mut_file, 'w') as f:
        f.write(mutation + '\n')
    return mut_file

def run_ddgun(mutation, pdb_file, chain, output_tag, output_dir): 

    pdb_file = os.path.abspath(pdb_file)
    
    if os.path.isfile(pdb_file) == False:
        logging.info("Oops! PDB file doesn't exist: "+pdb_file)
        sys.exit()
    else:
        # if output_tag not specified make output_tag: pdb_file_chain_mutation
        if output_tag == None:
            output_tag = '_'.join(pdb_file.replace('.pdb', ''), chain, mutation)
        
        output_dir = os.path.abspath(output_dir)
        output_dir = os.path.join(output_dir, output_tag)
        os.system('mkdir -p '+output_dir)
        out_file = os.path.join(output_dir, output_tag+'.ddgun')
        mut_file = mutList(output_dir, output_tag, mutation)
        
        if os.path.isfile(out_file) and os.stat(out_file).st_size != 0:
            logging.info('Output file exists: '+output_tag)
            sys.exit()
            
        # run 
        cmd = ['~/bin/ddgun/ddgun_3d.py', pdb_file, chain, mut_file, '>', out_file]
        cmd = ' '.join(cmd)
        logging.info('Computing ddG: '+output_tag)
        logging.info(cmd)
        os.system(cmd)

def main(argv):
    run_ddgun(FLAGS.mutation, FLAGS.pdb_file, FLAGS.chain, FLAGS.output_tag, FLAGS.output_dir)

if __name__ == '__main__':
    app.run(main)

