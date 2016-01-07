import itertools

from Bio import SeqIO
from Bio import pairwise2
import Bio.PDB.PDBParser

from util import open_file
from util import message

ONECODE = {'ALA':'A', 'CYS':'C', 'ASX':'B', 'GLU':'E', 'ASP':'D',
           'GLY':'G', 'PHE':'F', 'ILE':'I', 'HIS':'H', 'LYS':'K', 
           'XLE':'J', 'MET':'M', 'LEU':'L', 'PYL':'O', 'ASN':'N', 
           'GLN':'Q', 'PRO':'P', 'SER':'S', 'ARG':'R', 'SEL':'U', 
           'THR':'T', 'TRP':'W', 'VAL':'V', 'TYR':'Y', 'XAA':'X', 
           'GLX':'Z'}

class PDBParser(object):

    def parse_atom_seq(self, fp):
        residue_list = self.parse_residue_list(fp)
        return ''.join([ONECODE[x.get_resname().upper()] for x in residue_list])

    def parse_atom_cb(self, fp):
        residue_list = self.parse_residue_list(fp)
        ret = []
        for x in residue_list:
            if x.get_resname() == 'GLY': ret.append(self._get_atom(x, 'CA', None))
            elif x.has_id('CB'): ret.append(self._get_atom(x, 'CB', None))
            else: ret.append(self._get_atom(x, 'CA', None))
        return ret

    def parse_residue_list(self, fp):
        parser = Bio.PDB.PDBParser()
        structure = parser.get_structure('Read by PDBParser instance', fp)
        residue_list = structure.get_list()[0].get_list()[0].get_list()
        return [x for x in residue_list if x.get_id()[0].strip() == '']

    def _get_atom(self, res, k, d=None):
        return res[k] if res.has_id(k) else d

class IndexMapper(object):

    def __init__(self):
        self.MSK = 'X'

    def assign_index(self, seq, ref_seq, beg_idx=0):
        aln_seq, aln_ref_seq = pairwise2.align.globalms(seq.upper(), ref_seq.upper(), 1, -1, -.5, -.1)[0][:2]
        ref_idx_list = range(beg_idx, beg_idx + len(ref_seq))
        aln_ref_idx_list = self.build_aln_idx_list(ref_idx_list, aln_ref_seq)
        return self.assign_index_from_aln(aln_seq, aln_ref_seq, aln_ref_idx_list=aln_ref_idx_list)

    def build_aln_idx_list(self, idx_list, aln_seq):
        aln_idx_list, i = [], 0
        for y in aln_seq:
            if y != '-':
                aln_idx_list.append(idx_list[i])
                i += 1
            else: aln_idx_list.append(None)
        return aln_idx_list

    def assign_index_from_aln(self, aln_seq, aln_ref_seq, aln_ref_idx_list=[], msk_size=3, nonidentical=False):
        msk_aln_seq = self.mask_ambiguous_aln(aln_seq, msk_size)
        msk_aln_ref_seq = self.mask_ambiguous_aln(aln_ref_seq, msk_size)
        idx_list = []
        for x, y, curr in zip(msk_aln_seq, msk_aln_ref_seq, aln_ref_idx_list):
            if x != '-':
                if y != '-' and x != self.MSK and (nonidentical or x == y): idx_list.append(curr)
                else: idx_list.append(None)
        return idx_list

    def mask_ambiguous_aln(self, aln_seq, msk_size):
        frag_list = []
        for i, x in enumerate(aln_seq):
            if i == 0:
                frag_list.append(x)
            else:
                if prev_x == '-':
                    if x == '-': frag_list[-1] += x
                    else: frag_list.append(x)
                else:
                    if x == '-': frag_list.append(x)
                    else: frag_list[-1] += x
            prev_x = x
        frag_leng_list = []
        for frag in frag_list:
            if frag[0] == '-': frag_leng = 0
            else: frag_leng = len(frag)
            frag_leng_list += [frag_leng] * len(frag)
        ret = ''
        for x, n in zip(aln_seq, frag_leng_list):
            if x != '-' and n < msk_size: ret += self.MSK
            else: ret += x
        return ret

class Distance(object):

    def __init__(self, i, ri, j, rj, dist):
        self.i = i
        self.ri = ri
        self.j = j
        self.rj = rj
        self.dist = dist

_pdb_parser = PDBParser()
_idx_mapper = IndexMapper()

def build_edge(afa_file, pdb_file, cutoff=8.):
    seq = str(SeqIO.parse(open_file(afa_file), 'fasta').next().seq).upper().replace('-', '')
    message('Target sequence is retrieved from MSA.\n%s'%seq)
    atm_seq = _pdb_parser.parse_atom_seq(open_file(pdb_file)).upper()
    idx_list = _idx_mapper.assign_index(atm_seq, seq, beg_idx=1)
    if atm_seq != seq:
        message('Sequence in PDB is not identical to target sequence.\n%s'%atm_seq)
        message('Only the following residues will be considered for the analysis.\n' + ', '.join(['%s%d'%(res, i) for i, res in zip(idx_list, atm_seq) if i != None]))
    dist_list = _calc_dist(idx_list, atm_seq, pdb_file)
    edge_list = [(x.i, x.j) for x in dist_list if x.dist < cutoff]
    return edge_list

def _calc_dist(idx_list, atm_seq, pdb_file):
    dist_list = []
    cb_list = _pdb_parser.parse_atom_cb(open_file(pdb_file))
    for (i, res_i, cb_i), (j, res_j, cb_j) in itertools.combinations(zip(idx_list, atm_seq, cb_list), 2):
        if i == None or cb_i == None or j == None or cb_j == None: continue
        d = cb_i - cb_j
        if i > j: i, res_i, j, res_j = j, res_j, i, res_i
        dist_list.append(Distance(i, res_i, j, res_j, d))
    return dist_list
