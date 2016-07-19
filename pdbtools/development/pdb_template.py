#!/usr/bin/env python

# Copyright 2008, Michael J. Harms
# This program is distributed under General Public License v. 3.  See the file
# COPYING for a copy of the license.  

__description__ = \
"""
pdb_template.py

"""
__author__ = "Michael J. Harms"
__date__ = ""

def pdbTemplate(pdb):
    """
    """

    pass

def main():
    """
    Function to call if run from command line.
    """

    from helper import cmdline

    
    cmdline.initializeParser(__description__,__date__)
    cmdline.addOption(short_flag="s",
                          long_flag="some_option",
                          action="store",
                          default=None,
                          help="Set some option",
                          nargs=1,
                          type=int)


    file_list, options = cmdline.parseCommandLine()

    for pdb_file in file_list:
        
        f = open(pdb_file,'r')
        pdb = f.readlines()
        f.close()

        pdbTemplate(pdb)


# If run from command line...
if __name__ == "__main__":
    main()
        

