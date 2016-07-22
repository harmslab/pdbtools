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
import .atom_renumber
import .clean
from .helper import container
from .data.common import *

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
