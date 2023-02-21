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

rosetta='~/bin/rosetta/rosetta_src_2020.08.61146_bundle/main'
rosetta_ddg=os.path.join(rosetta, 'source/bin/ddg_monomer.default.linuxgccrelease')
rosetta_db=os.path.join(rosetta, 'database')

def resfile(path, output_tag, wt_aa, pos, mut_aa):
    resfile_file = os.path.join(path, output_tag+'.resfile')
    with open(resfile_file, 'w') as f:
        f.write('total 1\n')
        f.write('1\n')
        f.write(f'{wt_aa} {pos} {mut_aa}' + '\n') 
    return resfile_file
        
def runRosetta(mutation, pdb_file, chain, output_tag, output_dir): 
    
    # python run_rosetta.py --mutation=A35N --pdb_file=../pdb/2ami_Repair.pdb --chain=A --output_tag=2AMI_A35N --output_dir=../tmp
    
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
        out_file = os.path.join(output_dir, output_tag+'.rosetta')
        
        # create resfile
        wt_aa, pos, mut_aa = mutation[0], int(mutation[1:-1]), mutation[-1] # сделать проверку формата
        resfile_file = resfile(output_dir, output_tag, wt_aa, pos, mut_aa)
        
        # run 
        cmd1 = "cd " + output_dir
        cmd2 = [rosetta_ddg, '-in:file:s', pdb_file, '-ddg::mut_file', resfile_file, '-ddg::out ', out_file, '-in::file::fullatom -database', rosetta_db, '-fa_max_dis 9.0 -ignore_unrecognized_res -ddg:weight_file soft_rep_design -ddg::iterations 50 -ddg::dump_pdbs false -ddg::local_opt_only true -ddg::suppress_checkpointing true -ddg::mean true -ddg::min false -ddg::output_silent false']
        cmd2 = ' '.join(cmd2)
        cmd = ';'.join([cmd1,cmd2])
        logging.info('Computing ddG: '+output_tag)
        #print(cmd)
        os.system(cmd)

def main(argv):
    runRosetta(FLAGS.mutation, FLAGS.pdb_file, FLAGS.chain, FLAGS.output_tag, FLAGS.output_dir)

if __name__ == '__main__':
    app.run(main)

