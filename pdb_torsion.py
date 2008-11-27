#!/usr/bin/env python

# Copyright 2007, Michael J. Harms
# This program is distributed under General Public License v. 3.  See the file
# COPYING for a copy of the license.  

__description__ = \
"""
Determines the dihedral angles (phi,psi) for each residue in a protein.
"""

__author__ = "Michael J. Harms"
__date__ = "070601"

import os, sys
from helper import cmdline, geometry

def pdbTorsion(pdb):
    """
    Calculate the backbone torsion angles for a pdb file.
    """
    
    # Filter pdb file taking only N, CA, and C atoms
    atoms = [line for line in pdb if line[0:4] == 'ATOM' and
             line[13:16] in ['N  ','CA ','C  ']]

    # Define arrays that will contain CA, CO, and N coordinates
    num_resid = len(atoms)/3
    N = [[0.,0.,0.] for i in range(num_resid)]
    CA = [[0.,0.,0.] for i in range(num_resid)]
    CO = [[0.,0.,0.] for i in range(num_resid)]

    # Read the list of atoms into the coordinate arrays
    for i in range(num_resid):
        for j in range(3):
            N[i][j] = float(atoms[i*3][30+8*j:39+8*j])
            CA[i][j] = float(atoms[1+i*3][30+8*j:39+8*j])
            CO[i][j] = float(atoms[2+i*3][30+8*j:39+8*j])

    # Calculate phi and psi for each residue
    labels = []
    dihedrals = []
    for i in range(1,num_resid-1):
        try:
            labels.append((atoms[i][17:20],atoms[i][21:26]))
            dihedrals.append(geometry.calcDihedrals(CO[i-1],N[i],CA[i],CO[i],
                                                    N[i+1]))
        except ValueError:
            pass

    return dihedrals, labels


def main():
    """
    Call if this is called from the command line.
    """

    cmdline.initializeParser(__description__,__date__)

    file_list, options = cmdline.parseCommandLine()

    out = []
    for pdb_file in file_list:

        # Read in input file
        f = open(pdb_file,'r')
        pdb = f.readlines()
        f.close()

        # Calculate torsion angles and secondary structure
        dihedrals, labels = pdbTorsion(pdb)

        # Print out results in pretty fashion
        short_pdb =  os.path.split(pdb_file)[-1][:-4] 
        for i in range(len(dihedrals)):
            out.append("%30s%4s \"%5s\"%10.2F%10.2F\n" %\
                       (short_pdb,labels[i][0],labels[i][1],
                        dihedrals[i][0],dihedrals[i][1]))

    out = ["%10i%s" % (i,x) for i, x in enumerate(out)]
    
    header = "%10s%30s%4s%8s%10s%10s\n" % (" ","pdb","aa","res","phi","psi")
    out.insert(0,header)
    print "".join(out) 

    
    

if __name__ == "__main__":
    main()

