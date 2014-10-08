pdbtools
========

A set of tools for manipulating and doing calculations on wwPDB macromolecule structure files

Introduction
====
pdbTools is a set of command line [http://www.python.org python] scripts that manipulate [http://www.wwpdb.org/ wwPDB] protein and nucleic acid structure files.  There are many programs, both open source and proprietary, that perform similar tasks; however, most of these tools are buried within programs of larger functionality.  Thus, relatively simple calculations often involve learning a new program, compiling modules, and installing libraries. To fill a niche (and get the tasks done that I needed done), I started writing my own toolset.  This has evolved into the pdbTools suite.  The suite of programs is characterized by the following philosophy:

  # Each program should run as a stand-alone application with a standard, GNU/POSIX style command line interface.
  # Each program should be written in such a way to allow it to be used as a library of functions for more complex programs.
  # Programs should require a minimum of external dependencies.

Most of the scripts will run "out of the box" using a python interpreter.  The command line parser is designed to be flexible.  It will take an arbitrarily long list of pdb files, pdb ids, text files with pdb ids, or some mixture of all three.  If the pdb file or id is not in the working directory, scripts will attempt to download the pdb file from [http://www.rcsb.org/ RCSB].  Depending on the type of operation being done, a program will either write output files in the working directory or will print to stdout.  All structure outputs are written in standard pdb format.  All data outputs are in fixed-width column format.  They were designed to be read by the statistics package [http://cran.r-project.org/ R]; however, they should be easily parsed by other graphing programs.

*Note:* These scripts are only compatible with Python version 2.4-2.7.

Current functions
===

Miscellaneous
==
  * download pdb files from the RCSB database: [http://code.google.com/p/pdb-tools/source/browse/trunk/pdbTools/pdb_download.py pdb_download.py]

Structure-based calculations
==
Geometry
==
  * calculate protein center of mass: [http://code.google.com/p/pdb-tools/source/browse/trunk/pdbTools/pdb_centermass.py pdb_centermass.py]
  * calculate distance distributions: [http://code.google.com/p/pdb-tools/source/browse/trunk/pdbTools/pdb_dist-filter.py pdb_dist-filter.py],  [http://code.google.com/p/pdb-tools/source/browse/trunk/pdbTools/pdb_ion-dist.py pdb_ion_dist.py]
  * calculate backbone torsion angles: [http://code.google.com/p/pdb-tools/source/browse/trunk/pdbTools/pdb_torsion.py pdb_torsion.py]
  * calculate atom-by-atom solvent accessibility [http://code.google.com/p/pdb-tools/source/browse/trunk/pdbTools/pdb_sasa.py pdb_sasa.py](requires [http://www.bioinf.manchester.ac.uk/naccess/ NACCESS])
  * find disulfide bonds based on distance: [http://code.google.com/p/pdb-tools/source/browse/trunk/pdbTools/pdb_disulfide.py pdb_disfulfide.py]
  * find residues within some distance of each other: [http://code.google.com/p/pdb-tools/source/browse/trunk/pdbTools/pdb_contact.py pdb_contact.py], [http://code.google.com/p/pdb-tools/source/browse/trunk/pdbTools/pdb_water-contact.py pdb_water-contact.py], [http://code.google.com/p/pdb-tools/source/browse/trunk/pdbTools/pdb_close-contacts.py pdb_close-contacts.py]
  * find number of atoms neighboring another: [http://code.google.com/p/pdb-tools/source/browse/trunk/pdbTools/pdb_neighbors.py pdb_neighbors.py]
  * find ligands in structure file (ignoring boring ligands like water): [http://code.google.com/p/pdb-tools/source/browse/trunk/pdbTools/pdb_ligand.py pdb_ligand.py]
  * figure out oligomerization state of macrmolecule: [http://code.google.com/p/pdb-tools/source/browse/trunk/pdbTools/pdb_oligomer.py pdb_oligomer.py]

Energy calculation
==
  * calculate coulomb energy: [http://code.google.com/p/pdb-tools/source/browse/trunk/pdbTools/pdb_coulomb.py pdb_coulomb.py]
  * calculate the dipole moment of the protein: [http://code.google.com/p/pdb-tools/source/browse/trunk/pdbTools/pdb_moment.py pdb_moment.py]
  * calculate pKa of ionizable groups using the Solvent-Accessibility-modified Tanford-Kirkwood method [http://code.google.com/p/pdb-tools/source/browse/trunk/pdbTools/pdb_satk.py pdb_satk.py] (requires fortran compiler) 
Structure properties
==
  * extract structure experiment properties: [http://code.google.com/p/pdb-tools/source/browse/trunk/pdbTools/pdb_exper.py pdb_exper.py]
  * extract protein sequence from structure: [http://code.google.com/p/pdb-tools/source/browse/trunk/pdbTools/pdb_seq.py pdb_seq.py]
  * calculate theoretical pI, MW, fraction titratable residues, charge: [http://code.google.com/p/pdb-tools/source/browse/trunk/pdbTools/pdb_param.py pdb_param.py]

File/structure manipulation
==
  * add polar hydrogens: [http://code.google.com/p/pdb-tools/source/browse/trunk/pdbTools/pdb_addH.py pdb_addH.py] (requires [http://www.charmm.org/ CHARMM])
  * add missing heavy atoms, remove alternate conformations, etc.: [http://code.google.com/p/pdb-tools/source/browse/trunk/pdbTools/pdb_clean.py pdb_clean.py] (requires [http://www.charmm.org/ CHARMM])
  * mutate a residue: [http://code.google.com/p/pdb-tools/source/browse/trunk/pdbTools/pdb_mutator.py pdb_mutator.py] (requires [http://www.charmm.org/ CHARMM])
  * renumber atoms: [http://code.google.com/p/pdb-tools/source/browse/trunk/pdbTools/pdb_atom_renumber.py pdb_atom-renumber.py]
  * renumber residues: [http://code.google.com/p/pdb-tools/source/browse/trunk/pdbTools/pdb_residue-renumber.py pdb_residue-renumber.py]
  * offset all residues by a fixed amount: [http://code.google.com/p/pdb-tools/source/browse/trunk/pdbTools/pdb_offset.py pdb_offset.py]
  * center protein in xyz space: [http://code.google.com/p/pdb-tools/source/browse/trunk/pdbTools/pdb_centermass.py pdb_centermass.py]
  * places the asymmetric unit inside the unit cell: [http://code.google.com/p/pdb-tools/source/browse/trunk/pdbTools/pdb_centerasu.py pdb_centerasu.py]
  * take subset of residues from file: [http://code.google.com/p/pdb-tools/source/browse/trunk/pdbTools/pdb_subset.py pdb_subset.py]
  * split an NMR ensemble structure into individual files: [http://code.google.com/p/pdb-tools/source/browse/trunk/pdbTools/pdb_splitnmr.py pdb_splitnmr.py]
  * take a set of pdb files and create an individual directory for each one: [http://code.google.com/p/pdb-tools/source/browse/trunk/pdbTools/pdb_pdb2dir.py pdb_pdb2dir.py]
  * load data into the b-factor column: [http://code.google.com/p/pdb-tools/source/browse/trunk/pdbTools/pdb_bfactor.py pdb_bfactor.py]

Some of the programs are written as interfaces to other programs: [http://www.charmm.org/ CHARMM],  [http://www.bioinf.manchester.ac.uk/naccess/ NACCESS], which must be downloaded and installed separately if their functions are desired.  To use pdb_satk.py, a set of fortran packages must be compiled.

Contributing
==
If you find a bug or have an idea for a program you'd like in this package, feel free to make a pull request!

Project Owner
==
  * Mike Harms ([http://harmslab.uoregon.edu])
