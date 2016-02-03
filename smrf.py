#!/usr/bin/env python
#
# smrf
#   Coevolution analysis tool by using structure-based Markov random field
#
# Written by Chan-Seok Jeong <csjeong@kaist.ac.kr>, 2014-2016.
#
# If you use this software for any published work, please cite the accompanying paper.
#
#   Chan-Seok Jeong and Donsup Kim (2016) Structure-based Markov random field model for representing evolutionary constraints on functional sites. submitted.
#

import sys
import tempfile
import os.path
import ConfigParser

from edge import build_edge
from mrf import run_pmrf
from score import calc_pos_score
from score import calc_pair_score
from util import message
from util import write_file

VERSION = '0.2.1-beta.1'
CONFIG_FILE = '%s/smrf.cfg'%os.path.dirname(__file__)
DEFAULT_CONFIG = """\
[external]
pmrf_path =                ; absolute path for PMRF software
"""

def load_config():
    if not os.path.exists(CONFIG_FILE):
        write_file(DEFAULT_CONFIG, CONFIG_FILE)
        print 'Configuration file is created as %s'%(CONFIG_FILE)
        raise SystemExit
    cfg = ConfigParser.ConfigParser()
    cfg.read(CONFIG_FILE)
    return {'pmrf_path': os.path.expanduser(cfg.get('external', 'pmrf_path').split(';')[0])}

def smrf(afa_file, pdb_file, out_file1, out_file2, pmrf_path):
    edge_tf = tempfile.NamedTemporaryFile()
    edge_file = edge_tf.name
    mrf_tf = tempfile.NamedTemporaryFile()
    mrf_file = mrf_tf.name

    edge_list = build_edge(afa_file, pdb_file)
    write_file('\n'.join(['%s\t%s'%(i, j) for i, j in edge_list]), edge_file)
    message('MRF edge is determined.')

    run_pmrf(afa_file, edge_file, mrf_file, pmrf_path)
    message('MRF model is parameterized.')

    pos_score = calc_pos_score(mrf_file, afa_file)
    write_file('\n'.join(['%d\t%f'%(idx, val) for idx, val in pos_score]), out_file1)
    message('Positional coevolution scores are written in %s.'%(out_file1))

    pair_score = calc_pair_score(mrf_file)
    write_file('\n'.join(['%d\t%d\t%f'%(x.i, x.j, x.score) for x in pair_score]), out_file2)
    message('Pairwise coevolution scores are written in %s.'%(out_file2))

def print_help():
    print """\
Usage: %s [msa_file] [pdb_file] [out_file1] [out_file2]

Input:
  msa_file      Multiple sequence alignment (FASTA format)
  pdb_file      Protein structure (PDB format)

Output:
  out_file1     Positional coevolution score
                Column #1: residue position
                Column #2: score
  out_file2     Pairwise coevolution score
                Column #1: residue position #1
                Column #2: residue position #2
                Column #3: score
"""%sys.argv[0]

if __name__ == '__main__':
    cfg = load_config()
    if len(sys.argv[1:]) != 4:
        print_help()
        sys.exit(1)
    afa_file, pdb_file, out_file1, out_file2 = sys.argv[1:]
    smrf(afa_file, pdb_file, out_file1, out_file2, cfg['pmrf_path'])
