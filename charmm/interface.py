# Copyright 2007, Michael J. Harms
# This program is distributed under General Public License v. 3.  See the file
# COPYING for a copy of the license.  

__description__ = \
"""
interface.py

Functions that allow a pdb to be run through CHARMM.  The general method is
pdb2charmm -> runCharmm -> charmm2pdb.
"""
__author__ = "Michael J. Harms"
__date__ = "080418"

import math, os
import gen_input

# Set up charmm binary
try:
    charmm_path = os.environ["CHARMM"]
except KeyError:
    raise OSError("Environment variable CHARMM not defined!")

global charmm_bin
charmm_bin = os.path.join(charmm_path,"charmm")


class CharmmInterfaceError(Exception):
    """
    General error class for this module.
    """
    
    pass


def convertCharmmLine(atom_number,line):
    """
    Converts an individual charmm.cor line to a .pdb line.
    """
    
    fmt = "ATOM  %5i %-4s %-4s %4i    %8.3F%8.3F%8.3F%23s%s \n"

    residue_number = int(line[5:11])
    residue_name = line[11:15]
    atom_name = line[16:20].strip()
    if len(atom_name) != 4:
        atom_name = " %-3s" % atom_name.strip()
        atom_type = atom_name[1]
    else:
        atom_type = atom_name[0]
    x = float(line[20:30])
    y = float(line[30:40])
    z = float(line[40:50])
    
    return fmt % (atom_number,atom_name,residue_name,residue_number,x,y,z," ",
                  atom_type)
    

def runCharmm(input):
    """
    Run charmm < input.
    """

    global charmm_bin
    if not os.path.isfile(charmm_bin):
        err = "charmm binary \"%s\" does not exist!" % charmm_bin
        raise CharmmInterfaceError(err)

    print "Running: %s" % (charmm_bin)

    cin, cout = os.popen2(charmm_bin)
    cin.write(input)
    cin.close()
    out = cout.read()
    cout.close()

    return out


def charmm2pdb(charmm_output):
    """
    Converts a CHARMM-generated .cor file to a standard pdb file.
    """

    # Remove header
    output = charmm_output[4:]
    
    # Convert to pdb format
    pdb = [convertCharmmLine(i+1,l) for i, l in enumerate(output)]

    # Rename oddball residues
    rename_dict = {"OCT1":" O  ","OCT2":" OXT"}
    rename_lines = [l for l in pdb if l[12:16] in rename_dict.keys()]
    for line in rename_lines:
        index = pdb.index(line)
        new_atom = rename_dict[line[12:16]]
        pdb[index] = "%s%s%s" % (pdb[index][:12],new_atom,pdb[index][16:])

    return pdb


def charmmWash(input_structures,calc_type="single",keep_temp=False,
               hbond=None):
    """
    Wash a structure through CHARMM, adding polar hydrogens and missing heavy
    atoms.
    
    Input:
        pdb_files = [(file_contents,n_terminus,c_terminus),...]
            file_contents are the coordinates of a pdb file
            n_terminus and c_terminus are boolean
        calc_type = type of hydrogen minimization to do
        keep_temp = whether or not to keep temporary files.
    """
       
    # Convert all structures in input_structures to charmm-readable files
    struct_input = []    
    for index, structure in enumerate(input_structures):
    
        pdb = structure[0]
        n_terminus = structure[1]
        c_terminus = structure[2]
        
        tmp_file = "%i_tmp.pdb" % index
        g = open(tmp_file,'w')
        g.writelines(pdb)
        g.close()
        
        struct_input.append((tmp_file,n_terminus,c_terminus))
  
    # Make sure that there are not too many input files for CHARMM.  If too 
    # many files are passed, CHARMM doesn't raise an error; it just idles. 
    # To prevent the stall, I've made it a hard error.
    if len(struct_input) > 25:
        err = "There are too many input files (%i) " % len(struct_input)
        err += "for CHARMM.  Try splitting the input pdb into parts."
        raise CharmmInterfaceError(err)
 
    # Create CHARMM
    input = gen_input.createCharmmFile(struct_input,calc_type,hbond)
    
    # Run CHARMM
    charmm_out = runCharmm(input)
    
    # Try to read CHARMM output coordinates
    try:
        f = open("out.cor","r")
        coord_out = f.readlines()
        f.close()
    except IOError:
        
        f = open("charmm.inp","w")
        f.write(input)
        f.close()
        
        f = open("charmm.out","w")
        f.writelines(charmm_out)
        f.close()
        
        err = ["It appears that CHARMM has failed.\n"]
        err.append("   Input written to: charmm.inp\n")
        err.append("   Output writen to: charmm.out")
        raise CharmmInterfaceError("".join(err))
    
    # Convert CHARMM coordinates to a standard pdb format
    pdb_out = charmm2pdb(coord_out)
    
    # Delete temporary files
    if not keep_temp:
        tmp_files = [s[0] for s in struct_input]
        tmp_files.append("out.cor")
        for f in tmp_files:
            try:
                os.remove(f)
            except OSError:
                pass
    else:
        f = open("charmm.inp","w")
        f.write(input)
        f.close()
        
        f = open("charmm.out","w")
        f.writelines(charmm_out)
        f.close()
    
    return pdb_out

