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
    pdb_hash = [(l[0:6],i) for i, l in enumerate(pdb)]
    pdb_hash = [x[1] for x in pdb_hash if x[0] == "MODEL "]
    pdb_hash.append(-1)

    all_models = []
    for i in range(1,len(pdb_hash)):
        all_models.append(pdb[pdb_hash[i-1]:pdb_hash[i]])
       
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
        for index, model in enumerate(models):
            g = open("%s_%i.pdb" % (short_pdb,index),"w")
            g.writelines(model)
            g.close()

if __name__ == "__main__":
    main()

