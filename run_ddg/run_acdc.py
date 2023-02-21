# -*- coding: utf-8 -*-
from absl import app
from absl import flags
from absl import logging
import os

logging.set_verbosity(logging.INFO)
FLAGS = flags.FLAGS

flags.DEFINE_string("mutation", None, "Mutation")
flags.DEFINE_string("profile", None, "Directory with profiles")
flags.DEFINE_string("output_tag", None, "Output file name")
flags.DEFINE_string("output_dir", os.getcwd(), "Output directory")

# python run_acdc.py --mutation=A35N --output_tag=2AMI_A35N --output_dir=../tmp --profile=../fasta/wt_seq/msa/2AMI.profile

def run_acdc(mutation, profile, output_tag, output_dir): 

    # if output_tag not specified make output_tag from profile file name
    if output_tag == None:
        output_tag = '_'.join(profile.replace('.profile', ''), mutation)

    output_dir = os.path.abspath(output_dir)
    out_file = os.path.join(output_dir, output_tag+'.acdc')

    # run 
    cmd = ' '.join(['acdc-nn seq', mutation, profile, '>', out_file])
    logging.info('Computing ddG: '+output_tag)
    logging.info(cmd)
    os.system(cmd)

def main(argv):
    run_acdc(FLAGS.mutation, FLAGS.profile, FLAGS.output_tag, FLAGS.output_dir)

if __name__ == '__main__':
    app.run(main)

