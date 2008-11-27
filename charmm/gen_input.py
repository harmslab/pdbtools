# Copyright 2007, Michael J. Harms
# This program is distributed under General Public License v. 3.  See the file
# COPYING for a copy of the license.  

__description__ = \
"""
gen_input.py

A module containing functions to generate a hbuild.inp file for input into
CHARMM.
"""
__author__ = "Michael J. Harms"
__date__ = "080422"

import time, string, sys

global divider
divider = "! %s !\n" % (76*"*")

# Some hard-coded constants; hopefully these will never change
SINGLE_RTF = ["MASSES.RTF","AMINO.RTF"]
SINGLE_PRM =["PARM.PRM"]
FULL_RTF = ["MASSES_TITRA.RTF","AMINO_TITRA.RTF"]
FULL_PRM =["PARM_TITRA.PRM"]
FULL_MINIMIZATION_CONSTANTS = {"NPRI":50,"TOLG":0.001,"STEP":0.005,"TOLS":0.0,
                               "TOLENR":0.0,"INBFRQ":50,"CUTNB":20.0,
                               "CTOFNB":18.0,"CTONNB":15.0,"EPS":4.0}
WRITE_HBONDS_CONSTANTS = {"CTONHB":3.5,"CTOFHB":4.0,"CUTHB":4.5,"CTONHA":50.0,
                          "CTOFHA":70.0,"CUTHA":90.0}

def generateHeader():
    """
    Generates comment header for pyUHBD generated CHARMM input file.
    """

    global divider

    out = [divider]
    out.append("! This CHARMM input file was generated automatically using ")
    out.append("pyUHBD.\n")

    time_tuple = time.localtime(time.time())
    out.append("! %s\n" % time.asctime(time_tuple))
    out.append(divider)
    out.append("\n")

    return "".join(out)

def generateLib(library,type,is_first):
    """
    Generates a block of code to import rtf and param files.
    """

    # Generate comment description
    lib_out = [divider]
    lib_out.append("! Read in %s as %s.  Append = %s\n\n" %
                   (library,type,(not is_first)))

    # Decide whether to add "append" flag
    if is_first:
        append_me = ""
    else:
        append_me = "appe"

    # Import library
    lib_out.append("open read unit 3 card name -\n")
    lib_out.append("   \"$CHARMM_LIB/%s\"\n" % (library))
    lib_out.append("read %s unit 3 card %s\n" % (type,append_me))
    lib_out.append("close unit 3\n\n")

    return "".join(lib_out)


def readSeq(chain,filename,n_terminus,c_terminus,charmm_card):
    """
    Reads a sequence into chain from file "filename" with termini, using
    card number "charmm_card."
    """

    # Create proper termini commands
    first, last = "", ""
    if not n_terminus:
        first = " first none"
    if not c_terminus:
        last = "last none"
    termini_string = "%s %s" % (first,last)

    # Open card and read in sequence
    out = []
    out.append("open unit %i read card name %s\n" % (charmm_card,filename))
    out.append("read sequ pdb unit %i\n" % charmm_card)
    out.append("gener %s setup warn%s\n" % (chain,termini_string))
    out.append("rewind unit %i\n\n" % charmm_card)

    return "".join(out)

def readCoor(charmm_card):
    """
    Reads coordinates from charmm_card with offset, then closes card.
    """

    out = []
    out.append("read coor univ unit %i\n" % charmm_card)
    out.append("close unit %i\n\n" % charmm_card)

    return "".join(out)


def importFragment(index,pdb_file,n_terminus,c_terminus):
    """
    Takes a chain, possibly made of multiple fragments, and creates CHARMM
    code to read in their sequences and coordinates.
    """

    global divider

    chain = index   

    # Generate comment description
    out = [divider]
    out.append("! Import sequence and coordinates of chain %s.\n" % chain)
    out.append("!     %s, amino = %s, carboxyl = %s\n" % \
               (pdb_file,n_terminus,c_terminus))
    out.append("\n")

    # Set up charmm card
    charmm_card = 10 + index

    # Read in sequence
    out.append(readSeq(chain,pdb_file,n_terminus,c_terminus,charmm_card))

    # Read in coordinates and close card
    out.append(readCoor(charmm_card))

    return "".join(out)


def addHydrogens():
    """
    Builds missing atoms and adds hydrogens.
    """

    # Generate comment description
    out = [divider]
    out.append("! Build in missing atoms and add hydrogens.\n")

    # Generate missing coordinates
    out.append("! Build missing atoms\n")
    out.append("ic purge\nic param\nic fill preserve\nic build\n\n")

    # Build hydrogens on this chain
    out.append("! Build hydrogens\n")
    out.append("define test sele ( .not. type H* ) .and. ")
    out.append("( .not. init) show end\n")
    out.append("coor init sele (type ALA:VAL .and. hydrogens) end\n")
    out.append("hbuild sele type H* end\n")
    out.append("define test sele .not. init show end\n\n")

    # Define atoms that should be fixed during minimization
    out.append("! Fix all atoms that were not added by CHARMM\n")
    out.append("define fixed sele ( .not. type H* ) .and. ( init ) end\n\n")
    
    return "".join(out)

def writeHbonds(filename):
    """
    Write hbonds to an output file.
    """

    out = [divider]
    out.append("! Find and write all hbonds\n")

    m = WRITE_HBONDS_CONSTANTS

    cmd = "hbond acce all ctonhb %.1F ctofhb %.1F cuthb %.1F ctonha %.1F" % \
          (m["CTONHB"],m["CTOFHB"],m["CUTHB"],m["CTONHA"])
    cmd = "%s -\nctofha %.1F cutha %.1F\n" % (cmd,m["CTOFHA"],m["CUTHA"])
    out.append(cmd)

    out.append("open unit 17 form write name %s\n" % filename)
    out.append("write hbonds card unit 17\n")

    return "".join(out)


def minimizeSingle(steps):
    """
    Generates CHARMM code to minimize a structure with $steps number of steps.
    """

    # Generate comment description
    out = [divider]
    out.append("! Minimize non-fixed atoms with %s steps.\n\n" % steps)

    # Minimize
    out.append("cons fix sele fixed end\n")
    out.append("mini sd nstep %s\n" % steps)
    out.append("cons fix sele none end\n\n")

    return "".join(out)

def minimizeFull(steps):
    """
    Generates CHARMM code to minimize a structure; used for full-site
    calculation.
    """

    # Grab constants for full minimization
    m = FULL_MINIMIZATION_CONSTANTS

    # Put in comment identifier
    out = [divider]
    out.append("! Minimize non-fixed atoms with %s steps.\n\n" % steps)

    # Create output.
    #out.append("cons fix sele fixed end\n")
    out.append("cons fix sele .not. type H* end\n")
    out.append("mini sd -\n")
    out.append("nste %i npri %i tolg %.3F step %.3F tols % 0.3F -\n" %
               (steps,m["NPRI"],m["TOLG"],m["STEP"],m["TOLS"]))
    out.append("tolenr %.7F inbfrq %i cutnb %.2F ctofnb %.2F -\n" %
               (m["TOLENR"],m["INBFRQ"],m["CUTNB"],m["CTOFNB"]))
    out.append("ctonnb %.2F vswitch switch cdie eps %.1F\n" % (m["CTONNB"],
                                                               m["EPS"]))
    out.append("\ncons fix sele none end\n\n")

    return "".join(out)

def writeCoord(output_file):
    """
    Writes CHARMM coordinates to pdb file.
    """

    # Generate comment description
    out = [divider]
    out.append("! Write minimized coordinates out to %s\n\n" % output_file)

    # Write out coordinates
    out.append("open write unit 1 card name \"%s\"\n" % output_file)
    out.append("write coor unit 1 card\n")
    out.append("* pyUHBD addition of hydrogens using CHARMM minimization\n")
    out.append("*\n\n")

    return "".join(out)


def createCharmmFile(pdb_files,calc_type="single",hbond=None,num_steps=500,
                     coord_out="out.cor"):
    """
    Generate a charmm input file to add hyrogens to structure.
    """

    if calc_type == "single":
        rtf_files = SINGLE_RTF
        param_files = SINGLE_PRM
        minimize = minimizeSingle
    elif calc_type == "full":
        rtf_files = FULL_RTF
        param_files = FULL_PRM
        minimize = minimizeFull
    else:
        err = "\"%s\" is not a valid calc_type" % calc_type
        raise IOError(err)
        
    # Dump out header
    out = [generateHeader()]
    out.append("bomlev -2\n\n")

    # Add RTF imports
    for counter, r in enumerate(rtf_files):
        if counter == 0:
            is_first = True
        else:
            is_first = False
        out.append(generateLib(r,"rtf",is_first))

    # Add PARAM imports
    for counter, p in enumerate(param_files):
        if counter == 0:
            is_first = True
        else:
            is_first = False
        out.append(generateLib(p,"para",is_first))

    # Define pdb format we are reading
    out.append(divider)
    out.append("! Define type of pdb file to be read.\n\n")
    out.append("read univ\n* modified PDB setup\n*\nPDB\nw 61 6\nend\n\n")
    
    # Add commands for pdb files to out
    out.extend([importFragment(i,p[0],p[1],p[2])
                for i, p in enumerate(pdb_files)])

    # Add hydrogens
    out.append(addHydrogens())

    # Write out hbond information
    if hbond != None:
        out.append(writeHbonds(hbond))

    # Perform minization (single or full)
    out.append(minimize(num_steps))

    # Write to file
    out.append(writeCoord(coord_out))

    # Close CHARMM input file
    out.append("stop\n")

    return "".join(out)
    

