from Bio import SeqIO
import numpy as np

from util import open_file
from mrf import MRFParser

def calc_pos_score(mrf_file, afa_file):
    frac, psfm = _calc_psfm(afa_file)
    mrf_parser = MRFParser(mrf_file, idx_beg=1)
    score_list = [_calc_pos_score(node.weight[:-1], f) for node, f in zip(mrf_parser.iternodes(), psfm)]
    score_list = _rescale_score(score_list)
    vname_list = [node.idx for node in mrf_parser.iternodes()]
    return zip(vname_list, score_list)

def _rescale_score(scr_list):
    scr = np.array([x for x in scr_list if x != None])
    mn, sd = scr.mean(), scr.std()
    if sd > 0: sc = lambda x: (x - mn) / sd if x != None else None
    else: sc = lambda x: x - mn if x != None else None
    return [sc(x) for x in scr_list]

def _kl_divergence(p, q):
    return sum([x * np.log(x / y) for x, y in zip(p, q) if x > 0 and y > 0])

def _calc_pos_score(w, f):
    p = np.exp(w) / np.exp(w).sum()
    return _kl_divergence(f, p)

def _calc_psfm(afa_file):
    AA2Index = {aa:idx for aa, idx in zip('ACDEFGHIKLMNPQRSTVWY', range(20))}
    seqs = [str(r.seq) for r in SeqIO.parse(open_file(afa_file), 'fasta')]
    md_seq = lambda seq: ''.join([x for x, y in zip(seq, seqs[0]) if y != '-'])
    seqs = map(md_seq, seqs)
    n = len(seqs[0])
    f = np.zeros((n, 20), dtype=float)
    frac = np.zeros(n, dtype=float)
    for seq in seqs:
        for i, x in enumerate(seq):
            if AA2Index.has_key(x):
                if x == 'B': f[i, (AA2Index['N'], AA2Index['D'])] += .5
                elif x == 'J': f[i, (AA2Index['I'], AA2Index['L'])] += .5
                elif x == 'Z': f[i, (AA2Index['Q'], AA2Index['E'])] += .5
                elif x == 'O': f[i, AA2Index['K']] += 1.
                elif x == 'U': f[i, AA2Index['C']] += 1.
                else: f[i, AA2Index[x]] += 1.
                frac[i] += 1.
    f = f / np.reshape(np.repeat(f.sum(1), 20), (n, 20))
    frac /= len(seqs)
    return frac, f

def calc_pair_score(mrf_file):
    mrf_parser = MRFParser(mrf_file, idx_beg=1)
    score_list = _calc_pair_score(mrf_parser)
    zscore_list = _rescale_score([x.score for x in score_list])
    return [EdgeScore(x.i, x.j, zscore) for x, zscore in zip(score_list, zscore_list)]

def _calc_pair_score(mrf_parser):
    edge_score = lambda edge: EdgeScore(edge.idx1, edge.idx2, np.linalg.norm(edge.weight[:-1, :-1]))
    return _apc_correct([edge_score(edge) for edge in mrf_parser.iteredges()])

def _apc_correct(score_list):
    mn = {}
    for x in score_list:
        mn.setdefault(x.i, []).append(x.score)
        mn.setdefault(x.j, []).append(x.score)
    mn = {k:np.mean(v) for k, v in mn.iteritems()}
    tot_mn = np.mean([x.score for x in score_list])
    apc_func = lambda scr: scr.score - mn[scr.i] * mn[scr.j] / tot_mn
    ret = [EdgeScore(x.i, x.j, apc_func(x)) for x in score_list]
    return ret

class EdgeScore(object):

    def __init__(self, i, j, score):
        self.i = i
        self.j = j
        self.score = score
