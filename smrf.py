#!/usr/bin/env python
#
# smrf
#   Coevolution analysis tool by using structure-based Markov random field
#
# Written by Chan-Seok Jeong <csjeong@kaist.ac.kr>, 2014-2016.
#
# If you use this software for any published work, please cite the accompanying paper.
#
#   Chan-Seok Jeong, and Dongsup Kim (2016) Structure-based Markov random field model for representing evolutionary constraints on functional sites. BMC Bioinformatics. 17, 99.
#

import sys
import tempfile
import os.path
import optparse
import ConfigParser
import subprocess
import shlex

from edge import build_edge
from util import message
from util import write_file

VERSION = '0.4.0'
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

def smrf(afa_file, pdb_file, options, pmrf_path):
    edge_file, mrf_file = options.edge_file, options.mrf_file
    pmrf_exec = '%s/pmrf'%pmrf_path
    if not os.path.exists(pmrf_exec):
        print 'Cannot find the PMRF executable in the directory %s.'%(pmrf_path)
        sys.exit(1)

    ## Determine graph structure
    edge_list = build_edge(afa_file, pdb_file)
    write_file('\n'.join(['%s\t%s'%(i, j) for i, j in edge_list]), edge_file)
    message('MRF edge is determined.')

    ## Build MRF model
    cmd = '%s build %s --edge %s -o %s'%(pmrf_exec, afa_file, edge_file, mrf_file)
    subprocess.check_call(cmd.split())
    message('MRF model is parameterized.')

    ## Estimate positional coevolution
    cmd = '%s stat %s --mode pos'%(pmrf_exec, mrf_file)
    fp = sys.stdout if options.score_file1 == "stdout" else open(options.score_file1, 'w')
    p = subprocess.Popen(shlex.split(cmd), stdout=fp)
    p.wait()
    message('Positional coevolution scores are estimated.')

    ## Estimate pairwise coevolution
    if options.score_file2:
        cmd = '%s stat %s --mode pair'%(pmrf_exec, mrf_file)
        p = subprocess.Popen(shlex.split(cmd), stdout=open(options.score_file2, 'w'))
        p.wait()
        message('Pairwise coevolution scores are estimated.')

def parse_options():
    parser = optparse.OptionParser(usage = "%prog <msa_file> <pdb_file> [options]",
                                   version="%prog " + VERSION)
    grp = optparse.OptionGroup(parser, "Output options")
    grp.add_option("--pos", dest="score_file1", default="stdout", help="write positional coevolution scores to <score_file> (default: stdout)", metavar="<score_file>")
    grp.add_option("--pair", dest="score_file2", default=None, help="write pairwise coevolution scores to <score_file>", metavar="<score_file>")
    grp.add_option("--edge", dest="edge_file", default=None, help="write MRF edge list to <edge_file>", metavar="<edge_file>")
    grp.add_option("--mrf", dest="mrf_file", default=None, help="write MRF model to <mrf_file>", metavar="<mrf_file>")
    parser.add_option_group(grp)
    (options, args) = parser.parse_args()
    if len(args) != 2: parser.error("Incorrect number of arguments. See '%s --help'."%os.path.basename(sys.argv[0]))
    for arg in args:
        if not os.path.exists(arg): parser.error("The file %s does not exist"%arg)
    (msa_file, pdb_file) = args
    return msa_file, pdb_file, options

if __name__ == '__main__':
    cfg = load_config()
    afa_file, pdb_file, options = parse_options()
    if not options.edge_file:
        edge_tf = tempfile.NamedTemporaryFile()
        options.edge_file = edge_tf.name
    if not options.mrf_file:
        mrf_tf = tempfile.NamedTemporaryFile()
        options.mrf_file = mrf_tf.name
    smrf(afa_file, pdb_file, options, cfg['pmrf_path'])
