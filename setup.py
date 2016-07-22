#!/usr/bin/env python

# Copyright 2007, Michael J. Harms
# This program is distributed under General Public License v. 3.  See the file
# COPYING for a copy of the license.

try:
    from setuptools import setup
except:
    from distutils.core import setup

setup(name='pdbtools',
      version='0.1',
      description='A set of tools for manipulating and doing calculations on wwPDB macromolecule structure files',
      author='Michael J. Harms',
      author_email='harms@uoregon.edu',
      packages=['pdbtools'],
      scripts=['scripts/pdb_addH',
         'scripts/pdb_atom_renumber',
         'scripts/pdb_bfactor',
         'scripts/pdb_centerasu',
         'scripts/pdb_centermass',
         'scripts/pdb_clean',
         'scripts/pdb_closecontacts',
         'scripts/pdb_contact',
         'scripts/pdb_contactplot',
         'scripts/pdb_coulomb',
         'scripts/pdb_dist_filter',
         'scripts/pdb_disulfide',
         'scripts/pdb_download',
         'scripts/pdb_exper',
         'scripts/pdb_iondist',
         'scripts/pdb_ligand',
         'scripts/pdb_moment',
         'scripts/pdb_mutator',
         'scripts/pdb_neighbors',
         'scripts/pdb_offset',
         'scripts/pdb_oligomer',
         'scripts/pdb_param',
         'scripts/pdb_pdb2dir',
         'scripts/pdb_residue_renumber',
         'scripts/pdb_sasa',
         'scripts/pdb_satk',
         'scripts/pdb_seq',
         'scripts/pdb_splitnmr',
         'scripts/pdb_subset',
         'scripts/pdb_torsion',
         'scripts/pdb_watercontact'
      ],
      #install_requires=[],
      zip_safe=False)
