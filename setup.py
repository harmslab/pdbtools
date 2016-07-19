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
      install_requires=[
      ],
      zip_safe=False)
