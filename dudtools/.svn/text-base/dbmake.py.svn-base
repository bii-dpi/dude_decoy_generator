#!/bin/env python
"""Generate dockable database in parallel.

Michael Mysinger 200810 Created
"""

import os
import sys
import time
from optparse import OptionParser
import duddb
import processpool

# Constants
DB_SIZE = 4000
NUM_PROCESSES = 2

def create_database(infile=None, outdir='.', name="dbmake", prot="mid"):
    if infile is None:
        inf = sys.stdin
    else:
        inf = open(infile, 'r')
    pp = processpool.ProcessPool(NUM_PROCESSES)
    outdir = os.path.abspath(outdir)
    gendir = os.path.join(outdir, duddb.GEN_DIR)
    dbname = os.path.join(outdir, name)
    duddb.make_database(pp, dbname, gendir, inf, prot=prot)
    for res in pp:
        res.close()

def main(argv):
    """Parse arguments."""
    description = "Generate dockable database in parallel."
    usage = "%prog [options]"
    version = "%prog: version 200810 - created by Michael Mysinger"
    parser = OptionParser(usage=usage, description=description,
                          version=version)
    parser.set_defaults(infile=None, outdir='.', name="dbmake", prot="mid")
    parser.add_option("-i", "--infile",
           help="input smiles file (default: stdin)")
    parser.add_option("-o", "--outdir",
           help="output directory (default: %default)")
    parser.add_option("-n", "--name",
           help="name prefix for output databases (default: %default)")
    parser.add_option("-p", "--protonation", dest="prot", 
           help=("protonation type for databases: " +
                 "options = ref, mid, lo, hi, all  (default: %default)"))
    options, args = parser.parse_args(args=argv[1:])
    if len(args):
        parser.error("program takes no positional arguments.\n" +
                     "  Use --help for more information.")
    start_time = time.time()
    create_database(infile=options.infile, outdir=options.outdir, 
                    name=options.name, prot=options.prot)
    print 'Program took %.1f seconds.' % (time.time() - start_time)
    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv))
