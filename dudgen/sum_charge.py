#!/bin/env python
"""Print charge counts from ligands and decoys .charge files.

Michael Mysinger 200910 Created
"""

import os
import sys
import fnmatch
from optparse import OptionParser

def print_counts(outfile, counts):
    charges = [int(x) for x in counts]
    charges.sort()
    for c in charges:
        outfile.write("%3d %10d\n" % (c, counts[str(c)]))

def sum_charges(filename, counts):
    f = open(filename)
    splitf = (x.split() for x in f)
    for spl in splitf:
        if len(spl) == 9:
            key = spl[8]
            counts[key] = counts.get(key, 0) + 1
    f.close()
    return counts

def search_dirs(indir='.', outfile=None):
    if outfile is None:
        outf = sys.stdout
    else:
        outf = open(outfile, 'w')
    ligands = {}
    decoys = {}
    for fn in os.listdir('.'):
        if fnmatch.fnmatch(fn, "ligands.*.charge"):
            sum_charges(fn, ligands)
        elif fnmatch.fnmatch(fn, "decoys.*.charge"):
            sum_charges(fn, decoys)
    outf.write("Ligands\n")
    print_counts(outf, ligands)
    outf.write("\n")
    outf.write("Decoys\n")
    print_counts(outf, decoys)

def main(argv):
    """Parse arguments."""
    description = "Print charge counts from ligands and decoys .charge files."
    usage = "%prog [options]"
    version = "%prog: version 200910 - created by Michael Mysinger"
    parser = OptionParser(usage=usage, description=description,
                          version=version)
    parser.set_defaults(indir='.', outfile=None)
    parser.add_option("-i", "--indir",
           help="input directory (default: %default)")  
    parser.add_option("-o", "--outfile",
           help="output file (default: stdout)")
    options, args = parser.parse_args(args=argv[1:])
    options, args = parser.parse_args(args=argv[1:])
    if len(args):
        parser.error("program takes no positional arguments.\n" +
                     "  Use --help for more information.")
    search_dirs(indir=options.indir, outfile=options.outfile)
    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv))

