#!/usr/bin/env python

# Copyright 2007, Michael J. Harms
# This program is distributed under General Public License v. 3.  See the file
# COPYING for a copy of the license.  

__description__ =\
"""
pdb_neighbors.py

Calculates the number of CB atoms surrounding an atom in a pdb file.  A fake CB
atom is added for GLY residues.
"""

__author__ = "Michael J. Harms"
__date__ = "070521"


import os, sys
from helper import geometry

def pdbNeighbors(pdb,distance_cutoff=10,sequence_cutoff=4,residue=None,
                 atom=None,calc_everything=False):
    """
    Finds number of neighbors for each residue in a pdb and returns a list of 
    residues and a list of number of neighbors.
    """

    pdb = [l for l in pdb if l[0:6] == "ATOM  "]

    cb_atoms = []
    cb_coord = []
    for line_counter, line in enumerate(pdb):

        # Look for Gly CA and use to create fake CB
        if line[17:20] == "GLY" and line[13:16] == "CA ":
            cb_atoms.append("\"CB%s\"" % line[15:26])
            
            n = []; ca = []; co = []
            for i in range(3):
                n.append(float(pdb[line_counter-1][30+i*8:38+i*8]))
                ca.append(float(pdb[line_counter][30+i*8:38+i*8]))
                co.append(float(pdb[line_counter + 1][30+i*8:38+i*8]))
            #cb_coord.append(ca)
            cb_coord.append(geometry.calcGlyCbeta(n,ca,co))

        # For all other residues, append a CB
        #elif line[13:16] == "CA ":
        elif line[13:16] == "CB ":
            cb_atoms.append("\"%s\"" % line[13:26])
            cb_coord.append([float(line[31+i*8:38+i*8]) for i in range(3)])
        else:
            continue

    # Create list of atoms for which to find number of cb atom neighbors.
    if calc_everything:
        compare_atoms = ["\"%s\"" % l[13:26] for l in pdb]
        compare_coord = [[float(l[30+i*8:38+i*8]) for i in range(3)]
                         for l in pdb]
    elif residue == None and atom == None:
        compare_atoms = cb_atoms[:]
        compare_coord = cb_coord[:]
    else:
        if residue != None:
            compare_list = [l for l in pdb if l[17:20] == residue]
        else:
            compare_list = pdb[:]
        if atom != None:
            compare_list = [l for l in compare_list if l[13:16].strip() == atom]

        compare_atoms = ["\"%s\"" % l[13:26] for l in compare_list]
        compare_coord = [[float(l[30+i*8:38+i*8]) for i in range(3)]
                         for l in compare_list]

    # Calculate array of distances between each coordinate
    num_compare = len(compare_coord)
    num_cb = len(cb_coord)

    dist_array = [[0. for j in range(num_cb)]
                  for i in range(num_compare)]
    for i in range(num_compare):
        for j in range(num_cb):
            dist_array[i][j] = geometry.dist(compare_coord[i],cb_coord[j])

    # Count number of neighbors within distance_cutoff angstroms
    # for each position
    num_neighbors = [0 for x in compare_atoms]
    for i in range(num_compare):
        neighbors = 0
        for j in range(num_cb):

            # Skip neighbors that are simply close in sequence
            if compare_atoms[i][9] == cb_atoms[j][9]:
                seq_dist = int(compare_atoms[i][10:14])-int(cb_atoms[j][10:14])
                if abs(seq_dist) < sequence_cutoff:
                    continue

            # If the neighboring atom is close in space, add to counter
            if dist_array[i][j] <= distance_cutoff:
                num_neighbors[i] += 1
                
    return compare_atoms, num_neighbors


def main():
    """
    Call if program executed from command line.
    """

    from helper import cmdline

    # Parse command line
    cmdline.initializeParser(__description__,__date__)
    cmdline.addOption(short_flag="d",
                      long_flag="dist",
                      action="store",
                      default=10.0,
                      help="distance cutoff for neighbors",
                      nargs=1,
                      type=float)
    cmdline.addOption(short_flag="s",
                      long_flag="seq",
                      action="store",
                      default=4,
                      help="minimum sequence separation between neighbors",
                      nargs=1,
                      type=int)
    cmdline.addOption(short_flag="r",
                      long_flag="residue",
                      action="store",
                      default=None,
                      help="only calculate neighbors for residue (i.e. LYS)",
                      nargs=1,
                      type=str)
    cmdline.addOption(short_flag="a",
                      long_flag="atom",
                      action="store",
                      default=None,
                      help="only calculate neighbors to an atom (i.e. CA)",
                      nargs=1,
                      type=str)
    cmdline.addOption(short_flag="e",
                      long_flag="everything",
                      action="store_true",
                      default=False,
                      help="find number of CB neighbors for all atoms in pdb")
    

    file_list, options = cmdline.parseCommandLine()

    out = []
    for pdb_file in file_list:
    
        f = open(pdb_file,'r')
        pdb = f.readlines()
        f.close()

        residue_list, num_neighbors = pdbNeighbors(pdb,options.dist,options.seq,
                                                   options.residue,options.atom,
                                                   options.everything)

        short_pdb = os.path.split(pdb_file)[-1][:-4]
        for i in range(len(residue_list)):
            out.append("%30s%17s%10i\n" % (short_pdb,residue_list[i],
                                           num_neighbors[i]))

    out = ["%10i%s" % (i,x) for i, x in enumerate(out)]
    out.insert(0,"%10s%30s%17s%10s\n" % (" ","pdb","residue","cb_neigh"))
    out.insert(0,"# Dist cutoff: %.2F\n# Seq cutoff:  %i\n" %
               (options.dist,options.seq))

    print "".join(out)

if __name__ == "__main__":
    main()

