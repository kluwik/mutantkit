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

def writeIndList(tag, chain, wt_aa, pos, mut_aa, path):
    '''
    Creates individual list file with a single mutation in the specified path.
    To create mutation file for mutation S68G in PDB ID 1f21 chain A in the directory /wrk/il:
    writeIndList(68, 'S', 'G', '1f21', 'A', '/wrk/il')
    The mutation file will be named 'individual_list_1f21_S68G.txt' and contain a single line with 'SA68G;'.
    '''
    mutation = wt_aa + chain + str(pos) + mut_aa
    il_file = os.path.join(path, 'individual_list_'+tag+'.txt')
    with open(il_file, 'w') as f:
        f.write(mutation + ';' + '\n')
    return il_file

def runFoldX(mutation, pdb_file, chain, output_tag, output_dir): 
    #output_tag: pdbid_chain_mut; output_dir: pwd
    
    # python run_foldx.py --mutation=A35N --pdb_file=../pdb/2ami_Repair.pdb --chain=A --output_tag=2AMI_A35N --output_dir=../tmp
    
    pdb_file = os.path.abspath(pdb_file)
    output_dir = os.path.abspath(output_dir)
    
    if os.path.isfile(pdb_file) == False:
        print(f"Oops! PDB file {pdb_file} doesn't exist.")
        sys.exit()
    else:
        # if output_tag not specified make output_tag: pdb_file_chain_mutation
        if output_tag == None:
            output_tag = '_'.join(pdb_file.replace('.pdb', ''), chain, mutation)
        # create output directory
        output_dir = os.path.join(output_dir, output_tag)
        os.system('mkdir -p ' + output_dir)
        
        # create individual_list
        wt_aa, pos, mut_aa = mutation[0], int(mutation[1:-1]), mutation[-1] # сделать проверку формата
        il_file = writeIndList(output_tag, chain, wt_aa, pos, mut_aa, output_dir)
        
        # run 
        pdb_dir = os.path.dirname(os.path.abspath(pdb_file))
        cmd1 = "cd " + pdb_dir
        cmd2 = " ~/bin/foldx --command=BuildModel --pdb=" + os.path.basename(pdb_file) + " --mutant-file=" + il_file + " --output-dir=" + output_dir + " --output-file=" + output_tag
        cmd = ';'.join([cmd1,cmd2])
        logging.info('Computing ddG: '+output_tag)
        #print(cmd)
        os.system(cmd)
    
def main(argv):
    runFoldX(FLAGS.mutation, FLAGS.pdb_file, FLAGS.chain, FLAGS.output_tag, FLAGS.output_dir)

if __name__ == '__main__':
    app.run(main)

