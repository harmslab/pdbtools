#!/usr/bin/env python

# Copyright 2007, Michael J. Harms
# This program is distributed under General Public License v. 3.  See the file
# COPYING for a copy of the license.  

__description__ = "Adds polar hydrogens to a pdb file for a UHBD calculation."
__author__ = "Michael J. Harms"
__date__ = "070727"

from helper import container, cmdline, geometry
from math import sqrt
from pdb_atom_renumber import pdbAtomRenumber
from pdb_disulfide import pdbDisulfide
import time, os
import charmm.interface

class PdbAddHError(Exception):
    """
    General error class for this module.
    """

    pass
    

def renameResidue(pdb,residue_name,residue_number):
    """
    Renames residue defined by residue number to residue_name.
    """
    
    residue = [l for l in pdb if l[21:26] == residue_number]
    index = pdb.index(residue[0])
    for i, r in enumerate(residue):
        pdb[index + i] = "%s%-4s%s" % (r[:17],residue_name,r[21:])
    
    return pdb


def convertHis(pdb,his_types=None):
    """
    Rename HIS residues to HSD or HIS (1 or 2 tautomers).  If his_types is 
    specified, use numbers in there.  Otherwise, set all to 2 tautomer.
    """
    
    # Find histidines
    all_his = ["HIS ","HSD ","HSE ","HISA","HISB"]
    his_dict = {1:"HSD ",2:"HIS ",3:"HSC "}
    
    # Change his types
    HIS_lines = [l for l in pdb if l[17:21] in all_his and l[13:16] == "N  "]
    if his_types == None:
        his_types = [2 for l in HIS_lines]

    # Create final list of hist types
    his_types = [his_dict[k] for k in his_types]
    
    # Change his in file to proper his types
    for index, HIS_line in enumerate(HIS_lines):    
        residue_number = HIS_line[21:26]
        pdb = renameResidue(pdb,his_types[index],residue_number)

    return pdb


def convertResidues(pdb,atom_conv={},resid_conv={},atom_skip=[],resid_skip=[]):
    """
    Convert a pdb file to a CHARMM readable input (i.e. set HIS tautomeric
    states and give every group a charge.
    """

    atom_to_convert = atom_conv.keys()
    res_to_convert = resid_conv.keys()

    new_pdb = []
    for line in pdb:

        atom = line[12:16]
        if atom in atom_skip:
            continue
        elif atom in atom_to_convert:
            line = "%s%-4s%s" % (line[:12],atom_conv[atom],line[16:])
            
        res = line[17:21]
        if atom in atom_skip:
            continue
        elif res in res_to_convert:
            line = "%s%-4s%s" % (line[:17],resid_conv[res],line[21:])

        new_pdb.append(line)  
    
    return new_pdb

def processCysteines(pdb):
    """
    Find disulfide bonds.  Add hydrogens to free cysteines.  Change name of
    disulfide bonded residues to CSS.
    """

    # Find disulfide bonds
    free_cys, disulfide = pdbDisulfide(pdb)
    
    # Add hydrogens to free cysteines
    for cys in free_cys:
        
        residue = [l for l in pdb if l[21:26] == cys]
        CB_line = [l for l in residue if l[13:16] == "CB "][0]
        SG_line = [l for l in residue if l[13:16] == "SG "][0]
        
        CB_coord = [float(CB_line[30+8*i:38+8*i]) for i in range(3)]
        SG_coord = [float(SG_line[30+8*i:38+8*i]) for i in range(3)]
    
        HG_coord = geometry.calcHG(CB_coord,SG_coord)
            
        HG_line = "ATOM  %5i  %-3s %-4s%5s    %8.3F%8.3F%8.3F%23s%s \n" % \
            (1,"HG","CYS",cys,HG_coord[0],HG_coord[1],HG_coord[2]," ","H")
        
        index = pdb.index(SG_line)
        pdb.insert(index+1,HG_line)
    
    # Rename CYS in disulfide bonds to CSS
    for cys in disulfide:
        pdb = renameResidue(pdb,"CSS",cys)
    
    return pdb

def processTerm(pdb):
    """
    Process termini.  Change name of N-terminal residues to RESN, name of
    C-terminal residues to RESC, and add C-terminal carboxyl hydrogen.
    """

    # Change name of N-terminal residue
    HT1_lines = [l for l in pdb if l[13:16] == "HT1"]
    for HT1_line in HT1_lines:
 
        # Change name of N terminal residue to RESN (i.e. ALAN, ARGN, etc.)
        residue_name = HT1_line[16:20].strip()
        residue_name = "%4s" % (residue_name + "N")
        residue_number = HT1_line[21:26]
        pdb = renameResidue(pdb,residue_name,residue_number)

    # Add C-terminal hydrogens
    OXT_lines = [l for l in pdb if l[13:16] == "OXT"]
    for OXT_line in OXT_lines:
        
        # Grab residues in C-terminus
        last_residue = [l for l in pdb if l[21:26] == OXT_line[21:26]]
        C_line = [l for l in last_residue if l[13:16] == "C  "][0]
        O_line = [l for l in last_residue if l[13:16] == "O  "][0]
    
        # Calculate position of HXT
        OXT_coord = [float(OXT_line[30+8*i:38+8*i]) for i in range(3)]
        C_coord = [float(C_line[30+8*i:38+8*i]) for i in range(3)]
        O_coord = [float(O_line[30+8*i:38+8*i]) for i in range(3)]
        HXT_coord = geometry.calcHXT(C_coord,O_coord,OXT_coord)

        # Convert position of HXT to a pdb line
        HXT_line = "ATOM  %5i  %-3s %9s    %8.3F%8.3F%8.3F%23s%s \n" % \
            (1,"HXT",OXT_line[17:26],
             HXT_coord[0],HXT_coord[1],HXT_coord[2]," ","H")

        # Insert HXT
        index = pdb.index(OXT_line)
        pdb.insert(index+1,HXT_line)
        
        # Change name of C terminal residue to RESC (i.e. ALAC, ARGC, etc.)
        residue_name = OXT_line[16:20].strip()
        residue_name = "%4s" % (residue_name + "C")
        residue_number = OXT_line[21:26]
        
        pdb = renameResidue(pdb,residue_name,residue_number)
        

    return pdb

def flipAtoms(pdb):
    """
    Flip the OX1/OX2 labels of GLU and ASP residues.
    """
   
    flip = {"ASP":("OD1","OD2"),
            "GLU":("OE1","OE2")}
    flip_keys = flip.keys()

    for index, line in enumerate(pdb):
        res = line[17:30]
        if res in flip_keys and line[13:16] in flip[res]:
            new_atom = flip[res].index(line[13:16]) - 1
            pdb[index] = "%s%s%s" % (line[0:13],new_atom,line[16:])

    return pdb
 


def pdbAddH(pdb,pdb_id,uhbd_style=False,his_types=None,calc_type="single",
            keep_temp=False,hbond=False):
    """
    Add polar hydrogens to the structure using CHARMM for a UHBD calculation.
    """
    
    # Residues to alter and skip during processing
    if calc_type == "single":
        pdb2charmm_resid = {"LYS ":"LYSN","ARG ":"ARGN","GLU ":"GLUH",
                            "ASP ":"ASPH","LYSH":"LYSN","LSN ":"LYSN"}
        charmm2pdb_resid = {"LYSN":"LYS ","ARGN":"ARG ","GLUH":"GLU ",
                            "ASPH":"ASP ","HIS ":"HISA","HSD ":"HISB"}
    elif calc_type == "full":
        pdb2charmm_resid = {"GLU ":"GLUH","ASP ":"ASPH","LYSH":"LYS ",
                            "LSN ":"LYS "}
        charmm2pdb_resid = {"GLUH":"GLU ","ASPH":"ASP ","HIS ":"HISA",
                            "HSD ":"HISB","HSC ":"HISA"}
    else:
        err = "Calculation type \"%s\" not recognized!" % calc_type
        raise PdbAddHError(err)        
 
    charmm2pdb_atom_skip = [" HT3"]
    all_his = ["HIS ","HSD ","HSE ","HISA","HISB"]

    # Grab sequence and atoms from pdb file
    seq_lines = [l for l in pdb if l[0:6] == "SEQRES"]
    atom_lines = [l for l in pdb if l[0:6] == "ATOM  "]    
    
    # Create a pdb object that will find termini and renumber all atoms.  The
    # renumbering scheme is dumped to pdb_id_resid-conversion.txt
    pdb_obj = container.Structure(pdb_id,seq_lines,atom_lines)
    pdb_obj.renumberAtoms()
    pdb_obj.dumpNumberConversion("%s_resid-conversion.txt" % pdb_id)
    structure_list = pdb_obj.dumpStructures()

    # Convert residue names in structure
    for index, struct in enumerate(structure_list):
        tmp_struct = convertResidues(struct[0],resid_conv=pdb2charmm_resid)
        structure_list[index][0] = tmp_struct
    
    # Convert histidines to correct type (speificied in tautomer file).  If 
    # not tautomer file is specified, default HIS is passed to charmm
    if calc_type == "single":
        his_list = his_types
        for index, struct in enumerate(structure_list):
            if his_types != None:
                num_his = len([l for l in struct[0] if l[17:21] in all_his]) 
                try:
                    his = his_list[:num_his]
                    his_list = his_list[num_his:]
                except IndexError:
                    err = "Number of HIS in pdb and tautomer file do not match!"
                    raise cmdline.parser.error(err)
            else:
                his = None
            
            tmp_struct = convertHis(struct[0][:],his)
            structure_list[index][0] = tmp_struct

        # Make sure that all his where used
        if his_types != None and len(his_list) != 0:
            raise cmdline.parser.error(err)

    # For full calculation, convert all histidines to charged form (HSC)
    elif calc_type == "full":
        for index, struct in enumerate(structure_list):
            num_his = len([l for l in struct[0] if l[17:21] in all_his]) 
            his = [3 for i in range(num_his)]
            tmp_struct = convertHis(struct[0][:],his)
            structure_list[index][0] = tmp_struct
    
    # Flip carboxyl atoms 
    if calc_type == "full":
        for index, struct in enumerate(structure_list):
            structure_list[index][0] = flipAtoms(struct[0])
 
    # User CHARMM to add hydrogens
    try:
        out_pdb = charmm.interface.charmmWash(structure_list,calc_type,
                                              keep_temp,hbond)
    except charmm.interface.CharmmInterfaceError, (strerr):
        err = "Error in charmm!\n" 
        err += "%s\n" % strerr
        raise PdbAddHError(err)
    
    # Deal with addH specific changes in cysteines, termini, and residue names
    out_pdb = processCysteines(out_pdb)
    out_pdb = convertResidues(out_pdb,resid_conv=charmm2pdb_resid,
                              atom_skip=charmm2pdb_atom_skip)
    out_pdb = processTerm(out_pdb)

    new_pdb = container.Structure("tmp",[],out_pdb)
    new_pdb.loadNumberConversion("%s_resid-conversion.txt" % pdb_id,"fixed")
    new_pdb.renumberAtoms()

    out = []
    for chain in new_pdb.chains:
        out.extend(chain.atom_lines)
        
        ter = out[-1]
        ter = "%s%s%54s\n" % ("TER   ",ter[6:26]," ")
        out.append(ter)
    
    out.append("%-80s\n" % "END")

    if uhbd_style:
        out = [l for l in out if l[0:3] != "TER"] 

        # UHBD takes a non-standard pdb file; atom names must be left-justified.
        out = ["%s%-4s%s" % (l[:12],l[12:16].strip(),l[16:]) for l in out]
    
        # UHDB also cannot handle chain identifiers, remove them
        out = ["%s %s" % (l[0:21],l[22:]) for l in out]

    # Add header and END
    out.insert(0,"%-79s\n" % "REMARK  Polar hydrogens added by pdb_addH.py")
    out.insert(1,"%-79s\n" % ("REMARK  Time: %s" % time.asctime()))
    
    return out


def main():
    """
    Call if program called from command line.
    """
   
    # Parse command line
    cmdline.initializeParser(__description__,__date__)
    cmdline.addOption(short_flag="t",
                      long_flag="his_tautomers",
                      action="store",
                      default=None,
                      help="File containing his tautomers to use",
                      nargs=1)
    cmdline.addOption(short_flag="k",
                      long_flag="keep_temp",
                      action="store_true",
                      default=False,
                      help="Keep temporary files")
    cmdline.addOption(short_flag="f",
                      long_flag="full",
                      action="store_true",
                      default=False,
                      help="Add hydrogens for UHBD full calculation")
    cmdline.addOption(short_flag="b",
                      long_flag="hbond",
                      action="store",
                      default=None,
                      help="Write out hbonds to file",
                      nargs=1,
                      type=str)
    cmdline.addOption(short_flag="u",
                      long_flag="uhbd_style",
                      action="store_true",
                      default=False,
                      help="Write out in non-standard uhbd format")
    cmdline.addOption(short_flag="s",
                      long_flag="skip",
                      action="store_true",
                      default=True,
                      help="skip messed up pdb files")
    
 
    file_list, options = cmdline.parseCommandLine()    
    
    # Deal with his_tautomers file
    if options.his_tautomers != None:
    
        # Read in as a standard ascii file.
        his_types = cmdline.readFile(options.his_tautomers)
        try:
            his_types = [int(h) for h in his_types]
        
            # Make sure that entries are valid
            check = [h == 1 or h == 2 for h in his_types]
            if False in check:
                raise ValueError
        except ValueError:
            err = "His tautomer file can contain only 1's (HISB) and 2's (HISA)"
            cmdline.parser.error(err)
    else:
        his_types = None

    # Decide whether to keep temp files and how to format output
    keep_temp = options.keep_temp
    uhbd_style = options.uhbd_style
    hbond = options.hbond

    # Decide whether to add "full" hydrogens.
    if options.full:
        calc_type = "full"
    else:
        calc_type = "single"

    # Add hydrogens for every file in file_list
    file_list.sort()
    for file in file_list:
        pdb_id = os.path.split(file)[-1][:-4]
        out_file = "%sH.pdb" % pdb_id

        print "%s --> %s" % (file,out_file)

        # Load in file
        f = open(file,'r')
        pdb = f.readlines()
        f.close()
        
        # Add hydrogens
        try:
            pdb_out = pdbAddH(pdb,pdb_id,uhbd_style=uhbd_style,
                              his_types=his_types,
                              keep_temp=keep_temp,
                              calc_type=calc_type,
                              hbond=hbond)
        except PdbAddHError, (strerror):
            err = "Addition of hydrogens failed for %s\n" % file
            err += "Problem was with CHARMM\n.%s" % strerror
            print err

            if options.skip:
                print "CHARMM output written to error.log"
                g = open("error.log","a")
                g.write(err)
                g.write("charmm.out copied below\n%s\n" % (79*"-"))                

                h = open("charmm.out",'r')
                charmm_out = h.readlines()
                h.close()

                g.writelines(charmm_out)
                g.write("\n\n")
                g.close()

                continue
            else:
                sys.exit()

        # Dump to output file
        g = open(out_file,'w')
        g.writelines(pdb_out)
        g.close()
    
if __name__ == "__main__":
    main()

