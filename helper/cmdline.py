# Copyright 2007, Michael J. Harms
# This program is distributed under General Public License v. 3.  See the file
# COPYING for a copy of the license.  

__description__ = \
"""
pdb_cmdline.py

A standard argument parser for scripts in pdb_tools.  Returns a list of pdb
files (either single file or taken from directory) and a dictionary of option
values.
"""
__author__ = "Michael J. Harms"
__usage__ = "[options] pdb files and directories with pdb files" 
__date__ = "080422"

# Load builtin modules
import sys, os
import pdb_download

try:
    from optparse import OptionParser
    from optparse import Option, OptionValueError
except:
    error = """
    optparse module is missing (introduced in Python 2.3). Please verify that
    your Python version is up to date and that optparse is available. You are
    using Python %s """ % sys.version.split()[0]

def initializeParser(description,version):
    """
    Initialize an instance of parser with program description and version.
    """
    global parser

    usage = "%prog " + __usage__
    parser = OptionParser(usage=usage,version=version,description=description)


def addOption(short_flag,long_flag,action,default,help,**kwargs):
    """
    Add an option to the parser.
    """

    global parser

    parser.add_option("-%s" % short_flag,"--%s" % long_flag,action=action,
                      default=default,help="%s [default " % help + "%default]",
                      **kwargs)

def parseArgs(args):
    """
    A recursive function that allows discovery of pdb files 1) on command line,
    2) in a directory specified on the command line, 3) in a file w/out .pdb
    extension, or 4) from the rcsb.
    """

    file_list = []
    to_download = []
    
    for arg in args:

        # If it is already a file, append it
        if os.path.isfile(arg):
            if arg[-4:] == ".pdb":
                file_list.append(arg)
            else:
                g = open(arg,'r')
                inp_file = g.readlines()
                g.close()

                # Strip comments and blank lines 
                inp_file = [l for l in inp_file
                            if l[0] != "#" and l.strip() != ""]

                # Read arguments
                new_args = []
                for line in inp_file:
                    new_args.extend(line.split())

                # Strip blank arguments and errant punctuation
                punc_strip = ["\'","\""]
                new_args = [a for a in new_args if a.strip() != ""]
                new_args = ["".join([c for c in a if c not in punc_strip])
                            for a in new_args]


                tmp_file, tmp_download = parseArgs(new_args)

                file_list.extend(tmp_file)
                to_download.extend(tmp_download)
                    
        # If it is a directory, check for pdb files in that directory
        elif os.path.isdir(arg):
            tmp_list = os.listdir(arg)
            tmp_list = [f for f in tmp_list if f[-4:] == ".pdb"]
            tmp_list = [os.path.join(arg,f) for f in tmp_list]
            file_list.extend(tmp_list)
        
        # If it is not file or directory, see if it is a pdb id.  If the
        # file does not exist, try do download it from rcsb.
        else:
            if os.path.isfile("%s.pdb" % arg):
                file_list.append("%s.pdb" % arg)
            else:
                to_download.append(arg.lower())

    return file_list, to_download


def parseCommandLine():
    """
    Parse command line.
    """

    global parser

    # ---------- Parse arguments --------------------

    options, args = parser.parse_args()

    # Generate list of pdb files on which to perform calculations, doing some
    # error checking along the way
    if len(args) < 1:
        err = "You must specify at least one pdb file!"
        parser.error(err)
    else:
        file_list, to_download = parseArgs(args)
           
    # Download missing pdb files
    if len(to_download) > 0:
        to_download.sort()
        if pdb_download.pdbDownload(to_download):
            file_list.extend(["%s.pdb" % f for f in to_download])
        else:
            err = "pdb could not be found on rcsb!" 
            parser.error(err)

    # Remove duplicates from file_list by placing in dictionary
    file_dict = dict([(f,"") for f in file_list])
    file_list = file_dict.keys()
    file_list.sort()

    return file_list, options


def readFile(some_file):
    """
    Reads an ascii file if it exists, removes blank lines, whitespace, and
    comments (denoted by #).  It then splits values by whitespace, giving a list
    of all values in the file.
    
    Example:
        The file:
        1 2 3
        4
        5
        6 7
        Gives:
        ['1','2','3','4','5','6','7']
    
    """    

    # Make sure the file is truly a file
    if os.path.isfile(some_file):

        # Read in the file
        f = open(some_file)
        lines = f.readlines()
        f.close()

        # Strip out comments, extra whitespace, and blank lines
        lines = [l for l in lines if l[0] != "#"]
        lines = [l.strip() for l in lines if len(l.strip()) != 0]

        values = "".join(lines)
        values = values.split()

        return values
    
    else:
        raise IOError("%s does not exist" % some_file)

