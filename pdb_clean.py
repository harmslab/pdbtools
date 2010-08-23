#!/usr/bin/env python

# Copyright 2007, Michael J. Harms
# This program is distributed under General Public License v. 3.  See the file
# COPYING for a copy of the license.  

__description__ = \
"""
pdb_clean.py

Standardizes a Brookhaven pdb file
"""
__author__ = "Michael J. Harms"
__date__ = "070727"

import sys, time, string, os, shutil
import pdb_atom_renumber, charmm.interface
from helper import container
from pdb_data.common import *

class PdbCleanError(Exception):
    """
    General exception to raise if there is a problem with this module.
    """

    pass

def pdbCheck(coord):
    """
    Make sure the pdb file still has something in it after processing.
    """

    pdb_check = len([l for l in coord if l[0:6] == "ATOM  "])
    
    if pdb_check > 0:
        return 0
    else:
        return 1

def convertModifiedAA(coord,header):
    """
    Convert modified amino acids to their normal counterparts.
    """

    # See if there are any non-standard amino acids in the pdb file.  If there
    # are not, return
    modres = [l for l in header if l[0:6] == "MODRES"]
    if len(modres) == 0:
        return coord, header, []

    # Create list of modified residues
    mod_dict = dict([(l[12:15],l[24:27]) for l in modres])

    # Convert to ATOM entries, skipping non-backbone atoms.  These will be built
    # with CHARMM.
    backbone_atoms = ["N  ","CA ","C  ","O  "]
    new_coord = []
    for line in coord:
        if line[17:20] in mod_dict.keys():
            new = mod_dict[line[17:20]]
            if line[13:16] in backbone_atoms:
                new_line = "ATOM  %s%s%s" % (line[6:17],new,line[20:])
                new_coord.append(new_line)
        else:
            new_coord.append(line)

    # Convert non-standard atoms in the SEQRES entries
    converted_list = []
    new_header = []
    for line in header:
        if line[0:6] == "SEQRES":
            old_seq = line[19:70].split()
            new_seq = []   
            for aa in old_seq:
                if aa in mod_dict.keys():
                    new_seq.append(mod_dict[aa])
                else:
                    new_seq.append(aa)

            new_seq = "".join(["%s " % aa for aa in new_seq])
            new_seq.strip()
            new_seq = "%-50s" % new_seq
       
            new_header.append("%s%-50s%s" % (line[:19],new_seq,line[71:]))
        else:
            new_header.append(line)

    # Create output remarks
    conv = ["REMARK       converted %s to %s\n" % (k,mod_dict[k])
            for k in mod_dict.keys()]

    return new_coord, new_header, conv


def stripACS(coord):
    """
    Removes alternate confromations.
    """

    def removeLetters(line):
        """
        Mini function that removes letters that denote ACS.
        """

        if line[16] in string.letters:
            line = "%s %s" % (line[:16],line[17:])
        if line[26] in string.letters:
            line = "%s %s" % (line[:26],line[27:])

        return line

    # If a particular residue already has an atom, it will be in known_atom_dict
    # The second occurence of that atom in the same residue is assumed to be an
    # alternate conformation and is skipped.
    known_atom_dict = {}
    coord_out = []
    skipped = []
    for c in coord:
        residue = c[21:26]

        # If the residue is not known, update known_atom_dict and append line
        # to coordinate file
        if residue not in known_atom_dict.keys():
            out = removeLetters(c)
            coord_out.append(out)
            known_atom_dict.update([(residue,[c[13:16]])])

        # If the residue is known, determine if the atom has been seen before.  
        # If it has, skip it.  Otherwise, append to coord_out and 
        # known_atom_dict
        else:
            atom = c[13:16]
            if atom in known_atom_dict[residue]:
                skipped.append("REMARK%s" % c[6:])
            else:
                out = removeLetters(c)
                coord_out.append(out)
                known_atom_dict[residue].append(atom)
     
    return coord_out, skipped


def backboneCheck(coord):
    """
    Checks for duplicate residues (fatal) and missing backbone atoms.  If a 
    backbone atom is missing, the entire containing residue is deleted.
    """

    residue_numbers = []
    for line in coord:
        if line[17:26] not in residue_numbers:
            residue_numbers.append(line[17:26]) 
        
    to_remove = []
    for resid in residue_numbers:
        resid_atoms = [l for l in coord if l[17:26] == resid]

        # All backbone atoms in the protein
        backbone_atoms = [[l for l in resid_atoms if l[13:16] == "N  "],
                          [l for l in resid_atoms if l[13:16] == "CA "],
                          [l for l in resid_atoms if l[13:16] == "C  "],
                          [l for l in resid_atoms if l[13:16] == "O  "]]

        # If this is a proline, add CD to required backbone atoms
        if resid[0:3] == "PRO":
            backbone_atoms.append([l for l in resid_atoms if l[13:16] == "CD "])

        # If more than one of a backbone atom is found for a residue, we have
        # some sort of duplication.  If a backbone atom is missing, delete the
        # residue.
        for b in backbone_atoms:
            if len(b) > 1:
                err = "\%s\" is duplicated!" % resid
                raise PdbCleanError(err)
            if len(b) == 0:
                to_remove.append(resid)
                
    coord = [l for l in coord if l[17:26] not in to_remove]
    removed = ["REMARK       removed %s\n" % r for r in to_remove]

    return coord, removed


def addMissingAtoms(coord,seqres,keep_temp=False,renumber_residues=False,
                    pdb_id="",fix_atoms=True,num_steps=500):
    
    # Grab the b-factor and occupancy columns
    bfact_occ = dict([(l[13:26],l[54:67]) for l in coord])
    
    # Load pdb into pdb object to renumber for CHARMM
    pdb_obj = container.Structure("tmp",seqres,coord)
    pdb_obj.renumberAtoms()
    pdb_obj.dumpNumberConversion("numbering_conversion.txt")
    structure_list = pdb_obj.dumpStructures()
    
    # Do a charmm run to add missing atoms.
    try:
        new_coord = charmm.interface.charmmWash(structure_list,
            keep_temp=keep_temp,fix_atoms=fix_atoms,num_steps=num_steps)
    except charmm.interface.CharmmInterfaceError, (strerror):
        err = "CharmmInterfaceError\n%s\n" % strerror
        raise PdbCleanError(err)

    # Remove hydrogens
    new_coord = [l for l in new_coord if l[12] != "H" and l[13] != "H"]
    
    # Place charmm coordinates into new pdb container, and load in old numbering
    new_pdb = container.Structure("tmp",[],new_coord)
    if renumber_residues:
        shutil.move("numbering_conversion.txt",
                    "%s_resid-conversion.txt" % pdb_id)
    else:
        new_pdb.loadNumberConversion("numbering_conversion.txt","fixed")
        new_pdb.renumberAtoms()
        os.remove("numbering_conversion.txt")
   
    # Add bfactors, occupancies, and TER entries back in
    out = []
    for chain in new_pdb.chains:
        chain_atoms = chain.atom_lines
        
        for l in chain_atoms:
            try:
                out.append(3*"%s" % (l[:54],bfact_occ[l[13:26]],l[67:]))
            except KeyError:
                out.append(3*"%s" % (l[:54],"  1.00  1.00",l[67:]))
        
        ter = chain_atoms[-1]
        ter = "%s%s%54s\n" % ("TER   ",ter[6:26]," ")
        out.append(ter)
    
    out.append("%-80s\n" % "END")
    
    return out

def pdbClean(pdb,pdb_id="temp",chains="all",renumber_residues=False,
             keep_temp=False,fix_atoms=True,num_steps=500):
    """
    Standardize a pdb file:
        - Remove waters, ligands, and other HETATMS
        - Convert modified residues (i.e. Se-Met) to the normal residue
        - Remove alternate conformations (taking first in pdb file)
        - Find and remove residues with missing backbone atoms
        - 
        - Take only the specified chain
        - Renumber residues from 1
    """

    # Set up log
    log = ["REMARK  PDB processed using pdb_clean.py (harmsm@jhu.edu)\n"]
    log_fmt = "REMARK   - %s\n"
    log.append(log_fmt % ("Process time: %s" % time.asctime()))

    # Check pdb files for Brookhaven-added error warnings (CAVEAT and OBSLTE)
    error = [l for l in pdb if l[0:6] in ERROR_RECORDS]
    if len(error) != 0:
        err = "PDB might have problem!\n" + "".join(error)
        raise PdbCleanError(err)

    # Grab pdb header, excluding coordinates and deprecated records.
    header = [l for l in pdb if l[0:6] not in COORD_RECORDS]
    
    # Convert non-standard amino acids to standard ones
    coord = [l for l in pdb if l[0:6] in COORD_RECORDS]
    coord, header, converted = convertModifiedAA(coord,header)
    if len(converted) != 0:
        log.append(log_fmt % "Modified amino acids converted.")
        print log[-1],
        log.extend(converted)
    if pdbCheck(coord):
        err = "Modified amino acid converter removed all atoms!  Mangled pdb!"
        raise PdbCleanError(err)

    # Strip all entries in COORD_RECORDS except ATOM
    coord = [l for l in coord if l[0:6] == "ATOM  "]
    if pdbCheck(coord):
        err = "There are no ATOM entries in this pdb file!" 
        raise PdbCleanError(err)
    else:
        log.append(log_fmt % "HETATM entries removed.")
        print log[-1],

    # Grab only the chain we want, if specified 
    if chains != "all":
        coord = [l for l in coord if l[21] in chains]
        log.append(log_fmt % ("Took only chain %r." % chains))
        print log[-1],
        if pdbCheck(coord):
            err = "Chain filter (%r) removed all atoms in pdb file!" % chains
            raise PdbCleanError(err)

    # Strip alternate conformations 
    coord, skipped = stripACS(coord)
    if len(skipped) != 0:
        log.append(log_fmt % "Alternate conformations were removed.")
        print log[-1],
        log.extend(skipped)
    if pdbCheck(coord):
        err = "ACS stripper removed all atoms!  Mangled pdb file."
        raise PdbCleanError(err)

    # Check for missing backbone atoms; these residues are deleted
    coord, removed = backboneCheck(coord)
    if len(removed) != 0:
        log.append(log_fmt % "Residues with missing backbone atoms removed.")
        print log[-1],
        log.extend(removed)
    if pdbCheck(coord):
        err = "Backbone checker removed all atoms!  Mangled pdb file."
        raise PdbCleanError(err)
    
    # Add missing atoms using CHARMM
    print log_fmt % "Adding heavy atoms using CHARMM.",
    seqres = [l for l in header if l[0:6] == "SEQRES"]
    coord = addMissingAtoms(coord,seqres,keep_temp,renumber_residues,pdb_id,
                            fix_atoms,num_steps)
    log.append(log_fmt % "Missing heavy atoms were added with CHARMM.")
    
    # Renumber residues if requested
    if renumber_residues:
        log.append(log_fmt % "Residues renumbered from one.")
        print log[-1],
        
    # Renumber atoms from 1
    coord = pdb_atom_renumber.pdbAtomRenumber(coord)
    log.append(log_fmt % "Renumbered atoms from 1")
    print log[-1],

    # Standardize atom-type on far right pdb column
    coord = ["%s           %s  \n" % (c[:66],c[13]) for c in coord]
    log.append(log_fmt % "Atom types were standardized.")
    print log[-1],
    
    # Final check
    if pdbCheck(coord):
        err = "Unknown error occured and pdb has been mangled!"
        raise PdbCleanError(err)

    log = ["%-79s\n" % (l.strip()) for l in log]
    try:
        remark_pos = [l[0:6] for l in header].index("REMARK")
    except ValueError:
        remark_pos = 0

    # Return processed pdb file, placing log after preliminary remarks.
    out_pdb = []
    out_pdb.extend(header)
    out_pdb.extend(log)
    out_pdb.extend(coord)

    return out_pdb


def main():
    """
    To be called if module run from command line.
    """

    from helper import cmdline
    
    # Parse command line
    cmdline.initializeParser(__description__,__date__)
    cmdline.addOption(short_flag="c",
                      long_flag="chains",
                      action="store",
                      default=None,
                      help="File containing chains to take",
                      nargs=1)
    cmdline.addOption(short_flag="o",
                      long_flag="out_suffix",
                      action="store",
                      default="clean",
                      help="suffix to append to output pdb",
                      nargs=1)
    cmdline.addOption(short_flag="r",
                      long_flag="renumber_residues",
                      action="store_true",
                      default=False,
                      help="Renumber residues from 1")
    cmdline.addOption(short_flag="k",
                      long_flag="keep_temp",
                      action="store_true",
                      default=False,
                      help="Keep temporary files")
    cmdline.addOption(short_flag="s",
                      long_flag="skip",
                      action="store_true",
                      default=True,
                      help="skip messed up pdb files")
    cmdline.addOption(short_flag="f",
                      long_flag="fix_atoms",
                      action="store_false",
                      default=True,
                      help="fix atoms in original file")   
    cmdline.addOption(short_flag="n",
                      long_flag="num_steps",
                      action="store",
                      default=500,
                      help="number of minimization steps",
                      nargs=1)

 
    file_list, options = cmdline.parseCommandLine()    
   
    # Parse command line options
    if options.chains == None:
        chains = "all"
    else:
        chains = cmdline.readFile(options.chains)
    
    suffix = options.out_suffix
    renumber_residues = options.renumber_residues
    keep_temp = options.keep_temp
    fix_atoms = options.fix_atoms 
    num_steps = options.num_steps   
 
    for pdb_file in file_list:

        f = open(pdb_file,'r')
        pdb = f.readlines()
        f.close()
        
        print "Loading %s" % pdb_file
        pdb_id = pdb_file[:-4]

        try:
            pdb = pdbClean(pdb,pdb_id,chains,renumber_residues,keep_temp,
                           fix_atoms,num_steps)
        except PdbCleanError, (strerror):
            err = "Error cleaning \"%s\"\n%s\n" % (pdb_file,strerror)
            print err,

            if options.skip:
                g = open("error.log","a")
                g.write(err)
                g.close()
                continue
            else:
                sys.exit()

        out_file = "%s_%s.pdb" % (pdb_id,suffix)
        g = open(out_file,"w")
        g.writelines(pdb)
        g.close()

        print "Cleaned pdb written to %s" % out_file

if __name__ == "__main__":
    main()

