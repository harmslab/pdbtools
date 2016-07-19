#pdbtools

A set of tools for manipulating and doing calculations on wwPDB macromolecule structure files

##Introduction
pdbTools is a set of command line [python](http://www.python.org) scripts that manipulate [wwPDB](http://www.wwpdb.org/) protein and nucleic acid structure files.  There are many programs, both open source and proprietary, that perform similar tasks; however, most of these tools are buried within programs of larger functionality.  Thus, relatively simple calculations often involve learning a new program, compiling modules, and installing libraries. To fill a niche (and get the tasks done that I needed done), I started writing my own toolset.  This has evolved into the pdbTools suite.  The suite of programs is characterized by the following philosophy:

  * Each program should run as a stand-alone application with a standard, GNU/POSIX style command line interface.
  * Each program should be written in such a way to allow it to be used as a library of functions for more complex programs.
  * Programs should require a minimum of external dependencies.

Most of the scripts will run "out of the box" using a python interpreter.  The command line parser is designed to be flexible.  It will take an arbitrarily long list of pdb files, pdb ids, text files with pdb ids, or some mixture of all three.  If the pdb file or id is not in the working directory, scripts will attempt to download the pdb file from [RCSB](http://www.rcsb.org/).  Depending on the type of operation being done, a program will either write output files in the working directory or will print to stdout.  All structure outputs are written in standard pdb format.  All data outputs are in fixed-width column format.  They were designed to be read by the statistics package [R](http://cran.r-project.org/); however, they should be easily parsed by other graphing programs.

*Note:* These scripts are only compatible with Python version 2.4-2.7.

##Current functions

###Miscellaneous
  * download pdb files from the RCSB database: [download.py](https://github.com/harmslab/pdbtools/blob/master/pdbTools/download.py)

###Structure-based calculations
####Geometry
  * calculate protein center of mass: [centermass.py](https://github.com/harmslab/pdbtools/blob/master/pdbTools/centermass.py)
  * calculate distance distributions: [dist-filter.py](https://github.com/harmslab/pdbtools/blob/master/pdbTools/dist-filter.py),  [ion_dist.py](https://github.com/harmslab/pdbtools/blob/master/pdbTools/ion-dist.py)
  * calculate backbone torsion angles: [torsion.py](https://github.com/harmslab/pdbtools/blob/master/pdbTools/torsion.py)
  * calculate atom-by-atom solvent accessibility [sasa.py](https://github.com/harmslab/pdbtools/blob/master/pdbTools/sasa.py) [requires NACCESS](http://www.bioinf.manchester.ac.uk/naccess/ NACCESS)
  * find disulfide bonds based on distance: [disfulfide.py](https://github.com/harmslab/pdbtools/blob/master/pdbTools/disulfide.py)
  * find residues within some distance of each other: [contact.py](https://github.com/harmslab/pdbtools/blob/master/pdbTools/contact.py), [water-contact.py](https://github.com/harmslab/pdbtools/blob/master/pdbTools/water-contact.py), [close-contacts.py](https://github.com/harmslab/pdbtools/blob/master/pdbTools/close-contacts.py)
  * find number of atoms neighboring another: [neighbors.py](https://github.com/harmslab/pdbtools/blob/master/pdbTools/neighbors.py)
  * find ligands in structure file (ignoring boring ligands like water): [ligand.py](https://github.com/harmslab/pdbtools/blob/master/pdbTools/ligand.py)
  * figure out oligomerization state of macrmolecule: [oligomer.py](https://github.com/harmslab/pdbtools/blob/master/pdbTools/oligomer.py)

####Energy calculation
  * calculate coulomb energy: [coulomb.py](https://github.com/harmslab/pdbtools/blob/master/pdbTools/coulomb.py)
  * calculate the dipole moment of the protein: [moment.py](https://github.com/harmslab/pdbtools/blob/master/pdbTools/moment.py)
  * calculate pKa of ionizable groups using the Solvent-Accessibility-modified Tanford-Kirkwood method [satk.py](https://github.com/harmslab/pdbtools/blob/master/pdbTools/satk.py) (requires fortran compiler)
####Structure properties
  * extract structure experiment properties: [exper.py](https://github.com/harmslab/pdbtools/blob/master/pdbTools/exper.py)
  * extract protein sequence from structure: [seq.py](https://github.com/harmslab/pdbtools/blob/master/pdbTools/seq.py)
  * calculate theoretical pI, MW, fraction titratable residues, charge: [param.py](https://github.com/harmslab/pdbtools/blob/master/pdbTools/param.py)

###File/structure manipulation
  * add polar hydrogens: [addH.py](https://github.com/harmslab/pdbtools/blob/master/pdbTools/addH.py) [requires CHARMM](http://www.charmm.org/)
  * add missing heavy atoms, remove alternate conformations, etc.: [clean.py](https://github.com/harmslab/pdbtools/blob/master/pdbTools/clean.py) [requires CHARMM](http://www.charmm.org/)
  * mutate a residue: [mutator.py](https://github.com/harmslab/pdbtools/blob/master/pdbTools/mutator.py) [requires CHARMM](http://www.charmm.org/)
  * renumber atoms: [atom-renumber.py](https://github.com/harmslab/pdbtools/blob/master/pdbTools/atom_renumber.py)
  * renumber residues: [residue-renumber.py](https://github.com/harmslab/pdbtools/blob/master/pdbTools/residue-renumber.py)
  * offset all residues by a fixed amount: [offset.py](https://github.com/harmslab/pdbtools/blob/master/pdbTools/offset.py)
  * center protein in xyz space: [centermass.py](https://github.com/harmslab/pdbtools/blob/master/pdbTools/centermass.py)
  * places the asymmetric unit inside the unit cell: [centerasu.py](https://github.com/harmslab/pdbtools/blob/master/pdbTools/centerasu.py)
  * take subset of residues from file: [subset.py](https://github.com/harmslab/pdbtools/blob/master/pdbTools/subset.py)
  * split an NMR ensemble structure into individual files: [splitnmr.py](https://github.com/harmslab/pdbtools/blob/master/pdbTools/splitnmr.py)
  * take a set of pdb files and create an individual directory for each one: [pdb2dir.py](https://github.com/harmslab/pdbtools/blob/master/pdbTools/pdb2dir.py)
  * load data into the b-factor column: [bfactor.py](https://github.com/harmslab/pdbtools/blob/master/pdbTools/bfactor.py)

Some of the programs are written as interfaces to other programs: [CHARMM](http://www.charmm.org/),  [NACCESS](http://www.bioinf.manchester.ac.uk/naccess/ NACCESS), which must be downloaded and installed separately if their functions are desired.  To use satk.py, a set of fortran packages must be compiled.


##Usage

Almost all programs in the pdbTools suite have the same usage:

XXXX.py pdb_input optional_args > output

**pdb_input** can be one of the following (in any arbitrary combination):
  * pdb files
  * directories of pdb files
  * four-character pdb ids
  * text files containing whitespace delimited (i.e. space, tab, carriage return) lists of any combination of the other allowed types of arguments. If the list of arguments contains pdb files or ids that do not exist locally, the parser will attempt to download the files from the RCSB database.  

**optional_args**: Although the arguments to each program are identical, the options are quite different depending on the program requirements.  The best way to learn how to use a particular program is to type `XXXX.py` --help.  This will spit out a list of available options.  In most cases, the options are actually optional: the program will use a sane default if none is specified.  In some cases (notably `mutator.py`), options must be specified for the program to run.

**output**: Most scripts dump out a pdb file to standard out.  This can be captured using the ">" redirect.   Some write an output file that uses the name of the input pdb file as a suffix (e.g. `close-contacts.py 1stn.pdb` creates a file called 1stn.pdb.close_contacts).  

##Third Party Software
Some scripts require installation of third-party programs.  These should be installed according to the instructions given by the third-party, then placed into the $PATH variable.  To use the scripts that require CHARMM, the `$CHARMM` environment variable must be set to the directory containing the `charmm` binary and the `$CHARMM_LIB` environment variable to the directory containing the charmm parameter files.  


##Contributing
If you find a bug or have an idea for a program you'd like in this package, feel free to open an issue.  Even better: feel free to make a pull request!

##Project Owner
Mike Harms (https://github.com/harmsm, http://harmslab.uoregon.edu)
