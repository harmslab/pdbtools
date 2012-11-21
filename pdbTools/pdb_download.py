#!/usr/bin/env python

# Copyright 2007, Michael J. Harms
# This program is distributed under General Public License v. 3.  See the file
# COPYING for a copy of the license.  

"""
pdb_download.py

Download a pdb or set of pdb files and uncompress them.
"""

__author__ = "Michael J. Harms"
__date__ = "070521"
__usage__ = "pdb_download.py pdb_id or file w/ list of ids"

import os, sys, ftplib, shutil, gzip

HOSTNAME="ftp.wwpdb.org"
DIRECTORY="/pub/pdb/data/structures/all/pdb/"
PREFIX="pdb"
SUFFIX=".ent.gz"

def unZip(some_file,some_output):
    """
    Unzip some_file using the gzip library and write to some_output.
    """

    f = gzip.open(some_file,'r')
    g = open(some_output,'w')
    g.writelines(f.readlines())
    f.close()
    g.close()

    os.remove(some_file)

def pdbDownload(file_list,hostname=HOSTNAME,directory=DIRECTORY,prefix=PREFIX,
                suffix=SUFFIX):
    """
    Download all pdb files in file_list and unzip them.
    """

    success = True

    # Log into server
    print "Connecting..."
    ftp = ftplib.FTP()
    ftp.connect(hostname)
    ftp.login()

    # Remove .pdb extensions from file_list
    for file_index, file in enumerate(file_list):
        try:
            file_list[file_index] = file[:file.index(".pdb")]
        except ValueError:
            pass

    # Download all files in file_list
    to_get = ["%s/%s%s%s" % (directory,prefix,f,suffix) for f in file_list]
    to_write = ["%s%s" % (f,suffix) for f in file_list]
    for i in range(len(to_get)):
        try:
            ftp.retrbinary("RETR %s" % to_get[i],open(to_write[i],"wb").write)
            final_name = "%s.pdb" % to_write[i][:to_write[i].index(".")]
            unZip(to_write[i],final_name) 
            print "%s retrieved successfully." % final_name
        except ftplib.error_perm:
            os.remove(to_write[i])
            print "ERROR!  %s could not be retrieved!" % file_list[i]
            success = False

    # Log out
    ftp.quit()

    if success:
        return True
    else: 
        return False


def main():
    """
    If the function is called from the command line.
    """

    try:
        file_input = sys.argv[1:]
    except IndexError:
        print __usage__
        sys.exit()

    pdb_list = []
    for arg in file_input:
        if os.path.isfile(arg) and arg[-4:] != ".pdb":

            # Read in input file
            f = open(arg,"r")
            tmp_list = f.readlines()
            f.close()
        
            # Strip comments and blank lines, recombine file, then split on all 
            # white space
            tmp_list = [p for p in tmp_list if p[0] != "#" and p.strip() != ""]
            tmp_list = "".join(tmp_list)
            tmp_list = tmp_list.split()
            tmp_list = [p.lower() for p in tmp_list]
            pdb_list.extend(tmp_list)
       
        else:
            
            # Lower case, remove .pdb if appended
            pdb_id = arg.lower()
            if pdb_id[-4:] == ".pdb":
                pdb_id = pdb_id[:-4]

            pdb_list.append(pdb_id)

    # Download pdbs
    pdbDownload(pdb_list)
        

if __name__ == "__main__":
    main()


