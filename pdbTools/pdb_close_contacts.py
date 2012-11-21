#!/usr/bin/env python

# Copyright 2007, Michael J. Harms
# Modified 2011, Deborah L Crittenden
# This program is distributed under General Public License v. 3.  See the file
# COPYING for a copy of the license.  

__description__ =\
"""
pdb_close_contacts.py

Calculates close contacts for given residue/s
"""

__author__ = "Michael J. Harms"
__date__ = "070521"


import os, sys
import string
from helper import geometry

def pdbCloseContacts(pdb,distance_cutoff=3.5,sequence_cutoff=1,residue=0):
    """
    Finds and returns close contacts for atoms in a pdb file
    """

    pdb = [l for l in pdb if l[0:6] == "ATOM  "]

    crd = []
    coords = []
    splitcoords = []
    atm = []
    atoms = []
    splitatoms = []
    res = pdb[0][25]
    for j in range(0,len(pdb)):

	    line = pdb[j]
        coords.append([float(line[31+i*8:38+i*8]) for i in range(3)])
	    atoms.append(line[8:26])

        # Split pdb file into residues
	    if line[25] == res:
            crd.append([float(line[31+i*8:38+i*8]) for i in range(3)])
	        atm.append(line[8:26])
        else:
	    res = line[25]
	    splitcoords.append(crd)
	    splitatoms.append(atm)
	    crd = []
	    atm = []
            crd.append([float(line[31+i*8:38+i*8]) for i in range(3)])
	    atm.append(line[8:26])

    # Create list of atoms for which to find number of close contacts
    if residue == 0:
        print 'Warning: calculating all close contacts. This may take a'
	    print 'while and produce a lot of output. You should at least'
	    print 'consider setting distance cutoff to 3.5 (ultra-short)'
	compare_coord = coords
	compare_atoms = atoms
    else:
	compare_coord = splitcoords[residue-1]
	compare_atoms = splitatoms[residue-1]

    # Calculate array of distances between each coordinate
    num_compare = len(compare_coord)
    num_total = len(coords)

    dist_array = [[0. for j in range(num_total)]
                  for i in range(num_compare)]
    for i in range(num_compare):
        for j in range(num_total):
            dist_array[i][j] = geometry.dist(compare_coord[i],coords[j])

    # Count number of neighbors within distance_cutoff angstroms
    # for each position
    close_contacts = []
    for i in range(num_compare):
        for j in range(num_total):

            # Skip neighbors that are simply close in sequence
            if compare_atoms[i][13] == atoms[j][13]:
                seq_dist = int(compare_atoms[i][14:18])-int(atoms[j][14:18])
                if abs(seq_dist) < sequence_cutoff:
                    continue

            # If the neighboring atom is close in space, add to counter
            if dist_array[i][j] <= distance_cutoff:
                close_contacts.append([compare_atoms[i],atoms[j],str(dist_array[i][j])])
                
    return close_contacts


def main():
    """
    Call if program executed from command line.
    """

    from helper import cmdline

#    print 'Usage: pdb_close_contacts.py pdbfiles [-d distance cutoff]'
#    print '       -s [sequence cutoff] -r [residue number]'

    # Parse command line
    cmdline.initializeParser(__description__,__date__)
    cmdline.addOption(short_flag="d",
                      long_flag="dist",
                      action="store",
                      default=3.5,
                      help="distance cutoff for close contacts",
                      nargs=1,
                      type=float)
    cmdline.addOption(short_flag="s",
                      long_flag="seq",
                      action="store",
                      default=1,
                      help="minimum sequence separation between contacts",
                      nargs=1,
                      type=int)
    cmdline.addOption(short_flag="r",
                      long_flag="residue",
                      action="store",
                      default=0,
                      help="only calculate neighbors for residue number (e.g. 14)",
                      nargs=1,
                      type=int)

    file_list, options = cmdline.parseCommandLine()

    out = []
    for pdb_file in file_list:
    
        f = open(pdb_file,'r')
        pdb = f.readlines()
        f.close()

        close_contacts = pdbCloseContacts(pdb,options.dist,options.seq,
                                          options.residue)

        short_pdb = os.path.split(pdb_file)[-1][:-4]
        for i in range(len(close_contacts)):
            out.append("    ".join(close_contacts[i]) + "\n")

	f = open(short_pdb+".contacts",'w')
	f.writelines(out)
	f.close()


if __name__ == "__main__":
    main()

