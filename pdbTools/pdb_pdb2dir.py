#!/usr/bin/env python

# Copyright 2007, Michael J. Harms
# This program is distributed under General Public License v. 3.  See the file
# COPYING for a copy of the license.  

__description__ = \
"""
pdb_pdb2dir.py

Make an individual subdirectory for each pdb file in a set of pdb files.
"""
__author__ = "Michael J. Harms"
__date__ = "080204"

import os, shutil, sys

def main():
    """
    Creates directories in output_dir for each pdb file in input_dir and copies
    pdb files into respective directories.
    """

    from helper import cmdline

    cmdline.initializeParser(__description__,__date__)
    cmdline.addOption(short_flag="o",
                      long_flag="output",
                      action="store",
                      default=".",
                      help="output directory",
                      nargs=1,
                      type=str)
    file_list, options = cmdline.parseCommandLine()
    
    if not os.path.isdir(options.output):
        err = "Output directory \"%s\" does not exist!" % options.output
        raise cmdline.parser.error(err)

    for pdb_file in file_list:
        
        short_pdb = os.path.split(pdb_file)[-1][:-4]

        # Create the directory if it does not already exist
        dir = os.path.join(options.output,short_pdb)
        if not os.path.isdir(dir):
            os.mkdir(dir)

        # Copy the pdb into its directory
        shutil.copy(pdb_file,dir)


# If called from command line
if __name__ == "__main__":
    main()

