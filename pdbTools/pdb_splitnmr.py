#!/usr/bin/env python

# Copyright 2007, Michael J. Harms
# This program is distributed under General Public License v. 3.  See the file
# COPYING for a copy of the license.  

__description__ = \
"""
pdb_splitnmr.py

Takes an nmr pdb file and splits each model into its own pdb file.
"""
__author__ = "Michael J. Harms"
__date__ = "080204"

import sys, os

def splitNMR(pdb):
    """
    Split each model in an NMR pdb file into its own pdb.
    """

    to_strip = ["ENDMDL","MASTER"]
    pdb = [l for l in pdb if l[0:6] not in to_strip]
    pdb_splitter = [(l[0:6],i) for i, l in enumerate(pdb)]
    pdb_splitter = [x[1] for x in pdb_splitter if x[0] == "MODEL "]
    model_numbers = [pdb[i].split()[1].strip() for i in pdb_splitter]
    pdb_splitter.append(len(pdb))

    all_models = []
    for i in range(1,len(pdb_splitter)):
        all_models.append((model_numbers[i-1],
                           pdb[pdb_splitter[i-1]:pdb_splitter[i]]))
       
    return all_models 

def main():
    """
    Function to call if run from commmand line.
    """

    from helper import cmdline
    cmdline.initializeParser(__description__,__date__)
    file_list, options = cmdline.parseCommandLine()
    
    for pdb_file in file_list:

        f = open(pdb_file,'r')
        pdb = f.readlines()
        f.close()

        models = splitNMR(pdb)
        
        short_pdb = os.path.split(pdb_file)[-1][:-4]
        for model in models:
            index = model[0]
            pdb_lines = model[1]
            g = open("%s_%s.pdb" % (short_pdb,index),"w")
            g.writelines(pdb_lines)
            g.close()

if __name__ == "__main__":
    main()

