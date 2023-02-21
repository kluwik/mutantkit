# -*- coding: utf-8 -*-
from absl import app
from absl import flags
from absl import logging
import os, sys

logging.set_verbosity(logging.INFO)
FLAGS = flags.FLAGS

flags.DEFINE_string("mutation", None, "Mutation")
flags.DEFINE_string("fasta_file", None, "Sequence file")
flags.DEFINE_string("output_tag", None, "Output name")
flags.DEFINE_string("output_dir", os.getcwd(), "Output directory")
#flags.DEFINE_boolean('reverse', False, 'Calculate ddG of reverse mutation.')

flags.mark_flag_as_required("mutation")
flags.mark_flag_as_required("fasta_file")

# python run_ddgun.py --mutation=A35N --fasta_file=../fasta/wt_seq/2AMI.fasta --output_tag=2AMI_A33N --output_dir=../tmp

def mutList(path, output_tag, mutation):
    mut_file = os.path.join(path, output_tag+'.mut')
    with open(mut_file, 'w') as f:
        f.write(mutation + '\n')
    return mut_file

def run_ddgun(mutation, fasta_file, output_tag, output_dir): 

    fasta_file = os.path.abspath(fasta_file)
    
    if os.path.isfile(fasta_file) == False:
        print(f"Oops! Sequence file {fasta_file} doesn't exist.")
        sys.exit()
    else:
        # if output_tag not specified make output_tag: pdb_file_chain_mutation
        if output_tag == None:
            output_tag = '_'.join(fasta_file.replace('.fasta', ''), mutation)
        
        output_dir = os.path.abspath(output_dir)
        out_file = os.path.join(output_dir, output_tag+'.ddgun')
        mut_file = mutList(output_dir, output_tag, mutation)
        
        # run 
        cmd = ['~/bin/ddgun/ddgun_seq.py', fasta_file, mut_file, '>', out_file]
        cmd = ' '.join(cmd)
        logging.info('Computing ddG: '+output_tag)
        logging.info(cmd)
        os.system(cmd)

def main(argv):
    run_ddgun(FLAGS.mutation, FLAGS.fasta_file, FLAGS.output_tag, FLAGS.output_dir)

if __name__ == '__main__':
    app.run(main)

