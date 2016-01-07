import subprocess
import os.path
import sys

import numpy as np

def run_pmrf(afa_file, edge_file, mrf_file, pmrf_path):
    opts = '%s --edge %s -o %s'%(afa_file, edge_file, mrf_file)
    pmrf_exec = '%s/pmrf'%pmrf_path
    if not os.path.exists(pmrf_exec):
        print 'Cannot find the PMRF executable in the directory %s.'%(pmrf_path)
        sys.exit(1)
    cmd = '%s build %s'%(pmrf_exec, opts)
    subprocess.check_call(cmd.split())

class MRFParser(object):

    def __init__(self, mrf_file, idx_beg):
        self._mrf_file = mrf_file
        self._idx_beg = idx_beg

    def _open_mrf_file(self):
        return open(self._mrf_file)

    def _parse_idx(self, x):
        return int(x) - 1 + self._idx_beg

    def iternodes(self):
        conv = lambda x: 0. if x == '*' else float(x)
        fp = self._open_mrf_file()
        for line in fp:
            if line.strip() == '# NODE': break
        fp.next()
        for line in fp:
            if line.strip() == '# EDGE': break
            cols = line.strip().split('\t')
            idx = self._parse_idx(cols[0].split()[1])
            v = np.array(map(conv, cols[1:]))
            yield MRFNode(idx, v)
        fp.close()

    def iteredges(self):
        conv = lambda x: 0. if x == '*' else float(x)
        fp = self._open_mrf_file()
        for line in fp:
            if line.strip() == '# EDGE': break
        fp.next()
        for line in fp:
            if line.strip() == '//': break
            cols = line.strip().split('\t')
            idx1 = self._parse_idx(cols[0].split()[1])
            idx2 = self._parse_idx(cols[1].split()[1])
            v = np.array(map(conv, cols[2:]))
            n = int(np.sqrt(len(v)))
            yield MRFEdge(idx1, idx2, v.reshape((n, n), order='C'))
        fp.close()

class MRFNode(object):

    def __init__(self, idx, weight):
        self.idx = idx
        self.weight = weight

class MRFEdge(object):

    def __init__(self, idx1, idx2, weight):
        self.idx1 = idx1
        self.idx2 = idx2
        self.weight = weight
