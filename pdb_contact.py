#!/usr/bin/env python

# Copyright 2008, Michael J. Harms
# This program is distributed under General Public License v. 3.  See the file
# COPYING for a copy of the license.  

__description__ = \
"""
pdb_contacting-residues.py

Find all residues within some cutoff of a target residue in a pdb file.
"""
__author__ = "Michael J. Harms"
__date__ = "111116"

import os

def pdbContacting(pdb,target,cutoff,target_type="resname"):
    """
    """
  
    cutoff_sq = (cutoff)**2

    to_take = ["ATOM  ","HETATM"] 
    all_coord = [[l[12:26],[float(l[30+8*i:38+8*i]) for i in range(3)]]
                 for l in pdb if l[0:6] in to_take] 

    if target_type == "resname":
        target_list = [a for a in all_coord if a[0][5:8].strip() == target]
    else:
        target_list = [a for a in all_coord if int(a[0][10:14]) == target]

    out = []
    for t in target_list:
        contacts = []
        for a in all_coord:
            if sum([(a[1][i]-t[1][i])**2 for i in range(3)]) < cutoff:

                # ignore self
                if t[0] == a[0]:
                    continue
                contacts.append(a[0][5:].strip())
      
        # Grab only unique contacts 
        contacts = dict([(c,()) for c in contacts]).keys() 
        out.append("%s\t%s\n" % (t[0],("\t".join(contacts))))
   
    return out 

def main():
    """
    Function to call if run from command line.
    """

    from helper import cmdline
    import pdb_splitnmr

    
    cmdline.initializeParser(__description__,__date__)
    cmdline.addOption(short_flag="r",
                      long_flag="resname",
                      action="store",
                      default="TRP",
                      help="look for contacts near specified resname",
                      nargs=1)
    cmdline.addOption(short_flag="n",
                      long_flag="resnum",
                      action="store",
                      default=None,
                      help="look for contacts near specified residue number",
                      nargs=1)
    cmdline.addOption(short_flag="d",
                      long_flag="distance",
                      action="store",
                      default=3.5,
                      help="distance cutoff for calling contacts",
                      nargs=1,
                      type=float)

    file_list, options = cmdline.parseCommandLine()

    if options.resnum != None:
        target_type = "resnum"
        target = int(options.resnum)
    else:
        target_type = "resname"
        target = options.resname


    out = []
    for pdb_file in file_list:
        
        f = open(pdb_file,'r')
        pdb = f.readlines()
        f.close()

        # Only take the first NMR model
        models = pdb_splitnmr.splitNMR(pdb)
        if len(models) > 0:
            pdb = models[0]

        tmp_out = pdbContacting(pdb,target,options.distance,target_type)

        out.extend(["%s\t%s" % (pdb_file[:-4],t) for t in tmp_out])

    return out


# If run from command line...
if __name__ == "__main__":
    
    out = main()
    print "".join(out)
        

