#!/usr/bin/env python

# Copyright 2007, Michael J. Harms
# This program is distributed under General Public License v. 3.  See the file
# COPYING for a copy of the license.  

__description__ = \
"""
pdb_mutator.py

Mutates a residue in a pdb file.
"""

__author__ = "Michael J. Harms"
__date__ = "070729"
__description__ = "Mutates a residue in a pdb file"

import sys, time, string, os
import pdb_atom_renumber, pdb_clean
from helper import container
from pdb_data.common import *

class MutatorError(Exception):
    """
    General exception to raise if there is a problem with this module.
    """

    pass

def mutateResidue(pdb,residue_number,mutation,chain=None):
    """
    Renames residue defined by residue number to residue_name; designed for use
    with CHARMM.  The residue is "mutated" but all atoms besides N CA C O and
    CB are simply removed so that they will be added and minimized in CHARMM.
    """
    
    keep_atoms = ["N  ","CA ","C  ","O  ","CB "]
    
    # Find residue to mutate
    residue = [l for l in pdb if int(l[22:26]) == residue_number]
    if chain != None:
        residue = [l for l in residue if l[21] == chain]
    
    original_aa = residue[0][17:20]
    
    # Do mutation
    index = pdb.index(residue[0])
    for i, r in enumerate(residue):
        residue[i] = "%s%-4s%s" % (r[:17],mutation,r[21:])    
        pdb[index + i] = residue[i]
    
    # Remove non-backbone/CB atoms
    for atom in residue:
        if atom[13:16] not in keep_atoms:
            pdb.remove(atom)
    
    return pdb, original_aa

def pdbMutator(pdb,residue,mutation,chain=None,run_charmm=True):
    """
    Mutate a residue in the pdb file, energy minimizing with CHARMM if
    requested.
    """

    # grab header
    header = [l for l in pdb if l[0:6] not in pdb_clean.COORD_RECORDS and
                                l[0:6] not in pdb_clean.DEPRECATED_RECORDS]

    # Grab coordinates 
    coord = [l for l in pdb if l[0:6] == "ATOM  "]
    if pdb_clean.pdbCheck(coord):
        err = "There are no ATOM entries in this pdb file!" 
        raise MutatorError(err)

    coord, original_aa = mutateResidue(coord,residue,mutation,chain)
    mutation_string = "%s%i%s" % (AA3_TO_AA1[original_aa],residue,
                                  AA3_TO_AA1[mutation])
    
    # Set up log
    log = ["REMARK  %s introduced by pdb_mutator (harmsm@jhu.edu)\n" % \
           mutation_string]
    log_fmt = "REMARK   - %s\n"
    log.append(log_fmt % ("Process time: %s" % time.asctime()))
    if chain == None:
        log.append(log_fmt % ("Mutation introduced on all chains"))
    else:
        log.append(log_fmt % ("Mutation introduced on chain %s" % chain))
    
    # Add missing atoms using CHARMM
    if run_charmm:
        print log_fmt % "Adding mutated side chain using CHARMM",
        seqres = [l for l in header if l[0:6] == "SEQRES"]
        coord = pdb_clean.addMissingAtoms(coord,seqres)
        log.append(log_fmt % "Mutated sidechain built with CHARMM")
    
    # Renumber atoms from 1
    coord = pdb_atom_renumber.pdbAtomRenumber(coord)
    log.append(log_fmt % "Renumbered atoms from 1")
    print log[-1],

    # Standardize atom-type on far right pdb column
    coord = ["%s           %s  \n" % (c[:66],c[13]) for c in coord]
    log.append(log_fmt % "Atom types were standardized.")
    print log[-1],
    
    # Final check
    if pdb_clean.pdbCheck(coord):
        err = "Unknown error occured and pdb has been mangled!"
        raise MutatorError(err)

    # Return processed pdb file.
    out_pdb = []
    out_pdb.extend(log)
    out_pdb.extend(header)
    out_pdb.extend(coord)

    return out_pdb, mutation_string


def main():
    """
    To be called if module run from command line.
    """

    from helper import cmdline
    
    # Parse command line
    cmdline.initializeParser(__description__,__date__)
    cmdline.addOption(short_flag="c",
                      long_flag="chain",
                      action="store",
                      default=None,
                      help="CHAIN to mutate",
                      nargs=1)
    cmdline.addOption(short_flag="r",
                      long_flag="residue",
                      action="store",
                      type="int",
                      default=None,
                      help="Residue to mutate (REQUIRED)",
                      nargs=1)
    cmdline.addOption(short_flag="m",
                      long_flag="mutation",
                      action="store",
                      default=None,
                      help="Three-letter name of mutation (REQUIRED)",
                      nargs=1)
    cmdline.addOption(short_flag="s",
                      long_flag="simple",
                      action="store_true",
                      default=False,
                      help="No atoms beyond CB added (i.e. no CHARMM)")
                      
    
 
    file_list, options = cmdline.parseCommandLine()    
    
    # Parse command line options
    
    if options.residue == None:
        err = "Residue (-r) argument is required!"
        raise cmdline.parser.error(err)
    else:
        residue = options.residue
        
    if options.mutation == None:
        err = "Mutation (-m) argument is required!"
        raise cmdline.parser.error(err)
    else:
        mutation = options.mutation
    
    chain = options.chain
    run_charmm = not options.simple   

    for file in file_list:

        f = open(file,'r')
        pdb = f.readlines()
        f.close()
        
        print "Loading %s" % file
        pdb_id = file[:-4]
        pdb, mutation_string = pdbMutator(pdb,residue,mutation,chain,
                                          run_charmm)
    
        out_file = "%s_%s.pdb" % (pdb_id,mutation_string)
        g = open(out_file,"w")
        g.writelines(pdb)
        g.close()

        print "Mutated pdb written to %s" % out_file

if __name__ == "__main__":
    main()

