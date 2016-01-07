# smrf - Structure-based Markov random field

This is an evolutionary analysis tool using structure-based Markov random field


## Prerequisites

The standard Python package and the following packages are required.

  - [Biopython](http://biopython.org/)
  - [Numpy](http://www.numpy.org/)

The PMRF software is required.

  - [PMRF]()


## Installation

This tool is not packaged for the general installation, yet. Instead, download it from the git repository.

    $ git clone https://github.com/jeongchans/smrf


## Usage

1. Create configuration file

    $ smrf.py

2. Configure smrf.cfg with the directory where the PMRF software is installed 

3. Run smrf.py as follows

    Usage: smrf.py [msa_file] [pdb_file] [out_file1] [out_file2]

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


## Reference

If you use this software for any published work, please cite the accompanying paper.

    Chan-Seok Jeong and Donsup Kim (2016) Structure-based Markov random field model for representing evolutionary constraints on functional sites. submitted.

