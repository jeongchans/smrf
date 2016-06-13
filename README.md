# SMRF - Structure-based Markov random field
This is an evolutionary analysis tool using structure-based Markov random field.


## Installation

### Prerequisites
The standard Python package and the following packages are required.

- [Biopython]
- [Numpy]

The PMRF software (currently, compatible with v0.4.x) is required.

- [PMRF]

### Download
Download the the latest release from the [SMRF repository].

### Installation
This tool is not packaged for the general installation, yet. Instead, uncompress the downloaded file to the directory where the SMRF software will be installed.

### Configuring the Initial Setting
Create and configure smrf.cfg with the directory where the PMRF software is installed.

```sh
$ smrf.py   # This command creates smrf.cfg, if it does not exist in the directory where SMRF is installed
$ vi [SMRF_DIR]/smrf.cfg
```


## Usage

### Getting Started
Calculate the positional coevolution scores by using example.afa (MSA file formatted as FASTA) and example.pdb (PDB file).

```sh
$ smrf.py example.afa example.pdb
```

Show the help page.

```sh
$ smrf.py -h
```

### Writing coevolution scores
Write the positional and pairwise coevolution scores with `--pos` and `--pair` options, respectively.

```sh
$ smrf.py example.afa example.pdb --pos pos_score.txt --pair pair_score.txt
```

### Writing MRF model
Write the MRF model with `--mrf` option.

```sh
$ smrf.py example.afa example.pdb --mrf model.mrf
```

The MRF model is compatible with [PMRF] software.

## Reference
If you use this software for any published work, please cite the accompanying paper.

- Chan-Seok Jeong and Dongsup Kim (2016) Structure-based Markov random field model for representing evolutionary constraints on functional sites. BMC Bioinformatics. 17, 99.


[SMRF repository]: https://github.com/jeongchans/smrf/releases
[PMRF]: https://github.com/jeongchans/pmrf
[Biopython]: http://biopython.org/
[Numpy]: http://www.numpy.org/