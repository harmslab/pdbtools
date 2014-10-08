#!/usr/bin/env python

# Copyright 2007, Michael J. Harms
# This program is distributed under General Public License v. 3.  See the file
# COPYING for a copy of the license.  

__description__ = \
"""
pdb_disulfide.py

Returns a list of cys in disulfide bonds in a protein using a simple distance
constraint.  
"""

__author__ = "Michael J. Harms"
__date__ = "070726"

import sys
from pdb_clean import stripACS


def pdbDisulfide(pdb,disulfide_cutoff=4.0):
    """
    Finds all disulfides in pdb file based on disulfide cutoff distance (A).
    """

    pdb, skipped = stripACS(pdb)
    
    # Grab sg_atoms, coordinates, and residue numbers
    sg_atoms = [l for l in pdb if l[0:6] == "ATOM  " and l[13:16] == "SG "]
    coord = [[float(sg[30+8*i:38+8*i]) for i in range(3)] for sg in sg_atoms]
    cys_resid = [l[21:26] for l in sg_atoms]

    # Define all cys within disulfide_cutoff angstroms of each other as in
    # disulfide bonds
    cutoff_squared = disulfide_cutoff**2
    
    # Find disulfide
    free_cys = cys_resid[:]
    disulfide = []
    num_cys = len(cys_resid)
    for i in range(num_cys):
        for j in range(i+1,num_cys):
            r = sum([(coord[i][k] - coord[j][k])**2 for k in range(3)])
            if r < cutoff_squared:

                # This is a hack (in a sense).  If more than two cysteines are
                # close to one another (as in a Cys-Fe cluster), this will 
                # treat them all as disulfide bonded.  If I don't have the ugly
                # try...except statements, the code breaks because each cys has
                # multiple binding partners.
                try:
                    free_cys.remove(cys_resid[i])
                    disulfide.append(cys_resid[i])
                except ValueError: 
                    pass
                try:
                    free_cys.remove(cys_resid[j])
                    disulfide.append(cys_resid[j])
                except ValueError:
                    pass

    return free_cys, disulfide


def main():
    """
    Main function, if called from command line.
    """

    from helper import cmdline

    # Parse command line
    cmdline.initializeParser(__description__,__date__)
    cmdline.addOption(short_flag="c",
                      long_flag="cutoff",
                      action="store",
                      default=4,
                      help="cutoff for disulfide bonds",
                      nargs=1,
                      type=float)
    

    file_list, options = cmdline.parseCommandLine()

    for pdb_file in file_list:

        f = open(pdb_file,'r')
        pdb = f.readlines()
        f.close()

        free_cys, disulfide_bonds = pdbDisulfide(pdb,options.cutoff)

        out = ["<<%s>>\nFree Cys:        " % pdb_file]
        out.extend(["%s, " % c.strip() for c in free_cys])
        out.append("\nDisulfide Bonds: ")
        out.extend(["%s, " % d.strip() for d in disulfide_bonds])

        print "".join(out)

if __name__ == "__main__":
    main()
    

