#!/usr/bin/env python

# Copyright 2007, Michael J. Harms
# This program is distributed under General Public License v. 3.  See the file
# COPYING for a copy of the license.  

__description__ = \
"""
pdb_satk.py

Run solvent-accessibility modified Tanford-Kirkwood calculations using the
structure in a pdb file.  This is actually a simple interface to a set of 
fortran-77 programs (satkelni) written by Bertrand Garcia-Moreno E. in 1988.

For more information, see:
Matthew JB, Gurd FR, Garcia-Moreno B, Flanagan MA, March KL, Shire SJ. 
"pH-dependent processes in proteins." CRC Crit Rev Biochem. 1985;18(2):91-197.
"""
__author__ = "Michael J. Harms"
__date__ = "080313"

import os, sys, shutil
from math import sqrt

# Grab the path to the satk binaries
satk_path = os.path.join(sys.path[0],"satk")

class SATKError(Exception):
    """
    General error class for this module.
    """

    pass

def runBin(bin,arg_list):
    """
    Run some binary with arguments in arg_list.
    """

    # Make sure the binary exists
    if not os.path.isfile(bin):
        err = "Specified binary \"%s\" does not exist!\n" % bin
        err += "Have the binaries been compiled?"
        raise SATKError(err)

    # Run the command
    arg_list.insert(0,bin)
    status = os.spawnvp(os.P_WAIT,bin,arg_list)

    # Make sure the command was sucessful
    if status != 0:
        err = "%s failed (returning %i)" % (bin,status)
        raise SATKError(err)

def runMkwij(wij_name,wij_id='nothing',Dint=4,Dsolv=78.5,T=298.16,a=20,b=18,g=0,
             salts=[0.001,0.005,0.01,0.05,0.1,0.5,1.0]):
    """
    Wrapper for mkwij (make Wij).  mkwij creates a file (prefix.wij) that 
    contains the work to bring two charges from infinity to some separation at
    the interface of a sphere with dielectric constant Dint embedded in a
    medium of dielectric constant Dsolv.  The data in this file are then used 
    as a lookup table for calculations of dGij terms in a protein.  The 
    calculation takes the following parameters:

    Dint    protein dielectric constant
    Dsolv   solvent dielectric constant
    T       temperature (Kelvin)
    a       protein radius + ion exclusion layer (angstroms)
    b       protein radius (angstroms)
    g       depth of charge burial (angstroms)
    salts   ionic strength values to calculate (M)
    """

    # Create wij input file
    if os.path.isfile(wij_name):
        print wij_name,'deleted.'
        os.remove(wij_name)

    out = []
    out.append("%s%s" % (wij_name,'\n'))
    out.append("n\n")
    out.append("%s%s" % (wij_id,'\n'))
    out.append("%.1F   %.1F   %.2F\n" % (Dint,Dsolv,T))
    out.append("%.2F      %.2F      %.2F\n" % (a,b,g))
    out.append("%s\n" % len(salts))
    for salt in salts:
        out.append("%.6F00\n" % salt)

    f = open('mkwij.inp','w')
    f.write("".join(out))
    f.close()

    mkwij_bin = os.path.join(satk_path,'mkwij.out')
    mkwij_arg = ["mkwij.inp"]
    runBin(mkwij_bin,mkwij_arg)

def runSetup(prefix):
    """
    Wrapper for crdfrmt, prescan, and statacc.  These are a set of programs 
    that prep for the primary calculation.  statcc is of particular interest 
    because it calculates the solvent accessibility of every ionizable group 
    in the protein using a Richards method.  
    """

    # Verify that pdb file exists
    if not os.path.isfile("%s%s" % (prefix,'.pdb')):
        err = "The specified structure \"%s\" does not exist." % prefix
        raise SATKError(err)

    # Remove everything but the ATOM and HETATM entries, then add a dummy
    # header to the pdb file.
    f = open("%s.pdb" % prefix,"r")
    pdb = f.readlines()
    f.close()

    keep_lines = ["ATOM  ","HETATM"]
    pdb = [l for l in pdb if l[0:6] in keep_lines]
    pdb.insert(0,5*("%79s\n" % " "))

    f = open("%s_tmp.pdb" % prefix,"w")
    f.writelines(pdb)
    f.close() 

    # Delete existing files that would cause fortran crash
    ext_check = ['.cor','.wat','.chr','.saout','.elc']
    found_files = ["%s%s" % (prefix,ext) for ext in ext_check
                   if os.path.isfile("%s%s" % (prefix,ext))]
    if found_files != []:
        print 'The following files were deleted:'
        for output_file in found_files:
            print '\t', output_file
            os.remove(output_file)

    crdfrmt_bin = os.path.join(satk_path,'crdfrmt.out')
    crdfrmt_arg = ["%s_tmp.pdb" % prefix,"%s.cor" % prefix,"%s.wat" % prefix]
    runBin(crdfrmt_bin,crdfrmt_arg)

    prescan_bin = os.path.join(satk_path,'prescan.out')
    prescan_arg = ["%s.cor" % prefix,"%s.chr" % prefix]
    runBin(prescan_bin,prescan_arg)

    f = open('statacc.inp','w')
    f.write('y\nn')
    f.close()
    statacc_bin = os.path.join(satk_path,'statacc.out')
    statacc_arg = ["%s.chr" % prefix,"%s.cor" % prefix,
                   "%s.saout" % prefix, "%s.elc" % prefix,"statacc.inp"]
    runBin(statacc_bin,statacc_arg)

def createSin(prefix,cutoff=3.):
    """
    Takes a .elc file and generates a .sin file on the basis of distance
    constraints and solvent accessibility.
    """

    # Read elc file
    f = open("%s.elc" % prefix,'r')
    elc = f.readlines()[5:]
    f.close()

    residue = [int(l[0:5]) for l in elc]
    coord = [[float(l[20+9*i:29+9*i]) for i in [0,1,2]] for l in elc]
    inverse_sa = [float(l[59:64]) for l in elc]
    
    # Find distance to nearest atom for each atom
    num_atoms = len(coord)
    nearest = [100. for i in xrange(num_atoms)]
    for i in xrange(num_atoms):
        for j in xrange(num_atoms):
            if i == j or residue[i] == residue[j]:
                continue
            d = sqrt(sum([(coord[i][k]-coord[j][k])**2 for k in [0,1,2]]))
            if d < nearest[i]:
                nearest[i] = d

    # Populate sin file with residues, selecting a single atom from each 
    # residue based on two criteria:
    #   1. Distance to nearest other titratable group < cutoff
    #   2. Otherwise, take group with the highest inverse_sa (that is, the 
    #      highest solvent accessibility.  The .elc file has the 1-SA term in
    #      the left-most column).
    sin = [5*("%79s\n" % " ")]
    i = 0
    while i < num_atoms - 1:
        if residue[i] == residue[i+1]:
            if nearest[i] < cutoff:
                sin.append(elc[i])
            elif nearest[i+1] < cutoff:
                sin.append(elc[i+1])
            else:
                sort_list = zip(inverse_sa[i:i+2],[i,i+1])
                sort_list.sort()
                sin.append(elc[sort_list[0][1]])
            i += 2
        else:
            sin.append(elc[i])
            i += 1
    
    # Try to grab last residue (will only exist if previous residue had only one
    # atom) 
    try: 
        sin.append(elc[i])
    except IndexError:
        pass

    # Write to file               
    outputfile = "%s.sin" % prefix
    g = open(outputfile,'w')
    g.write("".join(sin))
    g.close()

def runSatkelni(prefix,wij_table,high_pH=20,low_pH=0,interval=0.25,
                salts=['all'],all_resids=True):
    """
    Wrapper for satkelni
    """

    # Make sure required files are present
    if not (os.path.isfile("%s.sin" % prefix) and os.path.isfile(wij_table)):
        err = "Not all required files present.  Check .wij and .sin files."
        raise SATKError(err)

    # Make sure there are < 100 pH steps (to avoid fortran array overrun)
    num_pH = (high_pH - low_pH)/interval
    if num_pH > 100:
        err = "Too many pH steps!  (You can only have 100 steps)."
        raise SATKError(err)

    # Check for old files that will kill fortran and delete them.
    ext_check = ['.pot','.plt','.ijd','.pov','.sinout']
    found_files = []
    for extension in ext_check:
        output_file = "%s%s" % (prefix,extension)
        if os.access(output_file,os.R_OK) == 1:
            found_files.append(output_file)
    if found_files != []:
        print 'The following files were deleted:'
        for output_file in found_files:
            print '\t', output_file
            os.remove(output_file)

    # Count the number of atoms in the .sin file
    f = open("%s%s" % (prefix,'.sin'),'r')
    sinfile = f.readlines()[7:]
    f.close()
    
    num_atoms = len(sinfile)

    # Determine the salts to place in the satkelni.inp file
    new_salts = []
    if salts[0] == 'all':
        g = open(wij_table,'r')
        wij = g.readlines()
        g.close()
        new_salts = [l.split()[0] for l in wij if l[0:2] == "  "]
    else:
        new_salts = ["%.3F" % s for s in salts]

    # Generate satkelni.inp
    out = []
    if all_resids:
        out.append('y\ny\ny\nn\n')
    else:
        out.append('n\ny\nn\n')
    out.append("%s %.1F %.1F %.2F\n" % (num_atoms,high_pH,low_pH,interval))
    out.append("1\n%s\n%i\n" % (wij_table,len(new_salts)))
    out.extend(["%s\n" % s for s in new_salts])
    out.append('y\n1\n0\n0\n0\n')    

    h = open('satkelni.inp','w')
    h.write("".join(out))
    h.close()

    # Run satkelni
    satkelni_bin = os.path.join(satk_path,"satkelni.out")
    satkelni_args = ["satkelni.inp","%s.sin" % prefix,"%s.pot" % prefix,
                     "%s.plt" % prefix,"%s.ijd" % prefix,"%s.sinout" % prefix,
                     "%s.pov" % prefix]
    runBin(satkelni_bin,satkelni_args)

    print 'satkelni run complete for', prefix

def createPrettyOutput(prefix):
    """
    Convert full .sinout file into pretty R-readable output.
    """

    # Read in sinout file
    f = open("%s.sinout" % prefix,'r')
    all_lines = f.readlines()
    f.close()

    # Grab header information
    param = all_lines[8:14]
    header = []
    header.append("# Dint: %10.3F (protein dielecric constant)\n" %\
                  float(param[0][6:12]))
    header.append("# Dext: %10.3F (solvent dielectric constant)\n" %\
                  float(param[1][6:12]))
    header.append("# T:    %10.3F Kelvin (temperature)\n" %\
                  float(param[2][6:12]))
    header.append("# a:    %10.3F angstroms (sphere radius)\n" %\
                  float( param[3][2:8]))
    header.append("# b:    %10.3F angstroms (sphere + ion radius)\n" %\
                  float(param[4][2:8]))
    header.append("# g:    %10.3F angstroms (depth of burial) \n" %\
                  float(param[5][2:8]))

    # Find various salt outputs 
    hash = [l[0:12] for l in all_lines]
    salt_indexes = [i for i, h in enumerate(hash) if h == "          AA"]

    pka_dict = {}
    titr_dict = {}
    for start in salt_indexes:

        # Find salt value
        pka_end = all_lines[start:].index("\n")
        salt_value = float(all_lines[start+pka_end+2][23:33])
       
        # Extract pKa vlues 
        pka_lines = all_lines[start+1:start+pka_end]
        pka_list = [((int(l[15:28]),l[10:15].strip()),float(l[28:46]))
                    for l in pka_lines]
        pka_dict[salt_value] = dict(pka_list)
      
        # Extract total pH titration 
        titr_index = start + pka_end + 6
        titr_end = all_lines[titr_index:].index("\n")
        titr_lines = all_lines[titr_index:titr_index+titr_end]
        titr = [[l[0:5],l[10:19],l[19:27],l[27:35],l[35:43],l[43:51],l[51:61],
                 l[61:70]] for l in titr_lines]
        titr = [tuple([float(v) for v in l]) for l in titr] 

        titr_dict[salt_value] = titr[:]
 
    salts = pka_dict.keys()
    salts.sort()
    
    # Write out pKa values 
    out = header[:]
    out.append("%10s%12s%10s%10s\n" % (" ","residue","ionic_str","pKa"))
    residues = pka_dict[salts[0]].keys()
    residues.sort()
    counter = 0 
    for r in residues:
        resid = "%s_%i" % (r[1],r[0])
        resid = ["%s_" % x for x in resid.split()]
        resid = "".join(resid)[:-1]           
        
        for s in salts:

            pka = pka_dict[s][r] 
 
            out.append("%10i%12s%10.3F%10.3F\n" % (counter,resid,s,pka))
            counter += 1

    g = open("%s.pka" % prefix,"w")
    g.writelines(out)
    g.close()

    # Write out pH titration
    out = header[:]
    out.append(10*"%10s" % (" ","ionic_str","pH","cations","anions","ions",
                            "protons","charge","dG_elec","dG_proton"))
    out.append("\n")
    counter = 0
    for s in salts:
        titr = titr_dict[s]
        titr = [8*"%10.3F" % x for x in titr]
        for value in titr:
            out.append("%10i%10.3F%s\n" % (counter,s,value))
            counter += 1
    
    g = open("%s.titr" % prefix,"w")
    g.writelines(out)
    g.close()

def cleanUp(prefix):
    """
    Delete all of the random fortran temporary files.
    """

    extension_list = ["chr","cor","elc","plt","plt","pot","pov","wat"]
    input_list = ["mkwij.inp","statacc.inp","satkelni.inp",
                  "%s_tmp.pdb" % prefix]
   
    remove_list = ["%s.%s" % (prefix,ext) for ext in extension_list]
    remove_list.extend(input_list)

    for f in remove_list:
        try:
            os.remove(f)
        except OSError:
            pass
 
def pdbSatk(pdb_file,wij_file=None,keep_temp=False):
    """
    Run an satk calculation using pdb_file.
    """  
 
    prefix = pdb_file[:-4]
   
    # Generate wij file with all ij potentials
    if wij_file != None:
        if not os.path.isfile(wij_file):
            err = "Specified wij file \"%s\" does not exist!" % wij_file
            raise  SATKError(err)
    else:
        wij_file = "%s.wij" % prefix
        runMkwij(wij_file)

    # Process pdb file
    runSetup(prefix)

    # Determine on which atoms to place the charge based on solvent
    # accessibility and ion pairs.
    createSin(prefix)
   
    # Do primary calculation
    runSatkelni(prefix,wij_file)

    # Delete temporary files
    if not keep_temp:
        cleanUp(prefix)

    createPrettyOutput(prefix)


def main():
    """
    Function to call if run from command line.
    """

    from helper import cmdline

    cmdline.initializeParser(__description__,__date__)
    cmdline.addOption(short_flag="w",
                      long_flag="wij_file",
                      action="store",
                      default=None,
                      help="Specify a pre-generated wij file",
                      nargs=1,
                      type=str)
    cmdline.addOption(short_flag="p",
                      long_flag="protein_dielec",
                      action="store",
                      default=4.,
                      help="Protein dielectric constant",
                      nargs=1,
                      type=float)
    cmdline.addOption(short_flag="s",
                      long_flag="solvent_dielec",
                      action="store",
                      default=78.5,
                      help="Solvent dielectric constant",
                      nargs=1,
                      type=float)
    cmdline.addOption(short_flag="b",
                      long_flag="ion_radius",
                      action="store",
                      default=2.,
                      help="Radius of ion (angstroms)",
                      nargs=1,
                      type=float)
    cmdline.addOption(short_flag="r",
                      long_flag="protein_radius",
                      action="store",
                      default=18.,
                      help="Radius of protein (angstroms)",
                      nargs=1,
                      type=float)
    cmdline.addOption(short_flag="d",
                      long_flag="depth",
                      action="store",
                      default=0.,
                      help="Depth of sidechain burial (angstroms)",
                      nargs=1,
                      type=float)
    cmdline.addOption(short_flag="T",
                      long_flag="temperature",
                      action="store",
                      default=298.16,
                      help="Temperature (Kelvin)",
                      nargs=1,
                      type=float)
    cmdline.addOption(short_flag="i",
                      long_flag="ionic_strengths",
                      action="store",
                      default=None,
                      help="File with set of ionic strengths (M)",
                      nargs=1,
                      type=str)
    cmdline.addOption(short_flag="k",
                      long_flag="keep_temp",
                      action="store_true",
                      default=False,
                      help="Keep temporary files")

    file_list, options = cmdline.parseCommandLine()

    
    # Generate wij file with all ij potentials
    if options.wij_file != None:
        wij_file = options.wij_file
        if not os.path.isfile(wij_file):
            err = "Specified wij file \"%s\" does not exist!" % wij_file
            raise  SATKError(err)
    else:

        # If a file with ionic strengths is specified, use it.
        if options.ionic_strengths != None:
            if not os.path.isile(ionic_strengths):
                err = "Ionic strengths file \"%s\" could not be found!" % \
                    options.ionic_strengths
            
            f = open(ionic_strengths,'r')
            lines = f.readlines()
            f.close()

            salt_list = []
            lines = [l for l in lines if l[0] != "#" and l.strip() != ""]
            for line in lines:
                salt_list.extend([float(entry) for entry in line])
            salt_list.sort()

        else:
            salt_list = [0.001,0.005,0.01,0.05,0.1,0.5,1.0]


        # Create wij file using the command line options specified
        wij_file = "wij.wij"
        runMkwij(wij_file,Dint=options.protein_dielec,
                          Dsolv=options.solvent_dielec,
                          T=options.temperature,
                          a=options.ion_radius+options.protein_radius,
                          b=options.protein_radius,
                          g=options.depth,
                          salts=salt_list)


    # Run satk on all files in file_list
    for pdb_file in file_list:

        pdb_tuple = os.path.split(pdb_file)
        if pdb_tuple[0] != "":
            shutil.copy(pdb_file,pdb_tuple[1])

        pdb_file = pdb_tuple[1]

        pdbSatk(pdb_file,wij_file=wij_file,keep_temp=options.keep_temp)


def GetIndivTitration(prefix,residtype,residnum):
    """
    Takes individual titration information from .sinout and generates 
    list of pka, dg, q, and # iterations as a function of pH for every
    ionic strength in sinout file.
    """

    residnum = str(residnum)

    inputfile = "%s%s" % (prefix,'.sinout')
    outputfile = "%s%s%s%s%s" % (prefix,'_',residtype,residnum,'.indiv')

    f = open(inputfile,'r')
    all_lines = f.readlines()
    f.close()

    from string import split
    from standard import fmt
    pka = []; ionic = []; ph = []; dg = []; charge = []; iterat = []

    for line_counter, line in enumerate(all_lines):
        if line[0:13] == 'IONIC STRN. =':
            header = line.split()
            ionic.append(float(header[3]))
            ph.append(float(header[10]))
            iterat.append(float(line[90:93]))

        if line[0:3] == residtype:
            column = line.split()
            if column[1] == residnum:
                pka.append(float(column[4]))
                dg.append(float(column[6]))
                charge.append(float(column[5]))

    newionic = [ionic[0]]; bigionic = []
    newpka = [pka[0]]; bigpka = []
    newph = [ph[0]]; bigph = []
    newiterat = [iterat[0]]; bigiterat = []
    newdg = [dg[0]]; bigdg = []
    newcharge = [charge[0]];  bigcharge = []

    for i in range(1,len(ionic)):
        if ionic[i] == ionic[i-1]:
            newionic.append(ionic[i])
            newpka.append(pka[i])
            newph.append(ph[i])
            newiterat.append(iterat[i])
            newdg.append(dg[i])
            newcharge.append(charge[i])
        else:
            bigionic.append(newionic)
            newionic = []
            bigpka.append(newpka)
            newpka = []
            bigph.append(newph)
            newph = []
            bigiterat.append(newiterat)
            newiterat = []
            bigdg.append(newdg)
            newdg = []
            bigcharge.append(newcharge)
            newcharge = []

    bigionic.append(newionic)
    bigpka.append(newpka)
    bigph.append(newph)
    bigiterat.append(newiterat)
    bigdg.append(newdg)
    bigcharge.append(newcharge)

    out = []
    out.append("%s%s%s%s%s" % ('# Titration of ',residtype,' ',residnum,
                               ' at ionic strengths:\n'))
    for i in range(len(bigionic)):
        out.append("%s%s%s" % ('#\t',bigionic[i][0],'\n'))


    out.append('\t')
    for i in range(len(bigionic)):
        out.append("%s%s%s" % ('I',str(i),'_ph\t'))
        out.append("%s%s%s" % ('I',str(i),'_pka\t'))
        out.append("%s%s%s" % ('I',str(i),'_dg\t'))
        out.append("%s%s%s" % ('I',str(i),'_q\t'))
        out.append("%s%s%s" % ('I',str(i),'_iterat\t'))
    out.append('\n')

    for i in range(len(bigionic[0])-1):
        out.append(fmt(i))
        for j in range(len(bigpka)):
            out.append("%s%s%s%s%s" %(fmt(bigph[j][i]),fmt(bigpka[j][i]),
                                      fmt(bigdg[j][i]),fmt(bigcharge[j][i]),
                                      fmt(bigiterat[j][i])))
        out.append('\n')

    g = open(outputfile,'w')
    g.write("".join(out))
    g.close()

    print "%s%s%s%s%s" % (prefix,'_',residtype,residnum,'.indiv'), 'written'


def GetPHTitration(prefix):

    from numarray import resize, array
    inputfile = "%s%s" % (prefix,'.sinout')
    outputfile = "%s%s" % (prefix,'.pHtitr')

    f = open(inputfile)
    all_lines = f.readlines()
    f.close()

    IS = []; pH = []
    for line_counter, line in enumerate(all_lines):
        if line[0:12] == 'IONIC. STRN.':
            i = 1
            current_line = all_lines[line_counter + i ]
            while current_line != '\n':
                IS.append("%5.3F" % float(current_line[2:13]))
                i += 1
                current_line = all_lines[line_counter + i]

        if line[0:19] == '        IONIC STRN.':
            number_IS = len(IS)
            j = 4; number_pH = 0
            current_line = all_lines[line_counter + j]
            while current_line != '\n':
                pH.append("%5.2F" % float(current_line[0:6]))
                j += 1
                current_line = all_lines[line_counter + j]
            number_pH = len(pH)
            all_data = resize(array(0.0),[number_IS,number_pH,4])
            break

    i = 0
    for line_counter, line in enumerate(all_lines):
        if line[0:19] == '        IONIC STRN.':

            j = 4
            current_line = all_lines[line_counter + j]
            while current_line != '\n':
                tmp_H = current_line[36:44]
                tmp_q = current_line[44:52]
                if tmp_H == ' ****** ': tmp_H = 0
                if tmp_q == ' ****** ': tmp_q = 0
                all_data[i,j-4,0] = float(tmp_H)
                all_data[i,j-4,1] = float(tmp_q)
                all_data[i,j-4,2] = float(current_line[54:62])
                all_data[i,j-4,3] = float(current_line[63:71])
                j += 1
                current_line = all_lines[line_counter + j]
            i += 1
            
    out = ["%s%s%s" % ('#pH titration of ',prefix,' at IS values:\n')]
    for i in range(number_IS):
        out.append("%s%s%s" % ('#\t',IS[i],'\n'))

    out.append('\tpH\t')
    for i in range(number_IS):
        out.append("%s%s%s%s" % ('IS',i,'_boundH','\t'))
        out.append("%s%s%s%s" % ('IS',i,'_charge','\t'))
        out.append("%s%s%s%s" % ('IS',i,'_dGtot','\t'))
        out.append("%s%s%s%s" % ('IS',i,'_dGH','\t'))
    out.append('\n')
    
    for i in range(number_pH):
        out.append("%s%s" % (str(i),'\t'))
        out.append("%s%s" % (pH[i],'\t'))
        for j in range(number_IS):
            for k in range(4):
                out.append("%5.3F%s" % (all_data[j,i,k],'\t'))
        out.append('\n')

    g = open(outputfile,'w')
    g.write("".join(out))
    g.close()

    print outputfile, 'written.' 

if __name__ == "__main__":
    main()

