#!/bin/env python
"""Unique charged properties file to generate ligands.charge.

Michael Mysinger 201205 Created
"""

import os
import sys
import logging
import os.path as op
from optparse import OptionParser
from dudzinc import write_uniqs

PID_START = 10000001

class ScriptError(Exception):
    def __init__(self, msg, value=99):
        self.value = value
        Exception.__init__(self, msg)

def charged2ligands(in_f, outfile):
    """Module code starts here."""
    splits = [line.strip().split('\t') for line in in_f]
    pids = range(PID_START, PID_START+len(splits))
    pids.reverse()
    split_by_id = {}
    [split_by_id.setdefault(spl[1], []).append(spl) for spl in splits]
    uniqs = {}
    for k, splits in split_by_id.iteritems():
        uniqs[int(k)] = dict([tuple(x[3:]), (pids.pop(), x[0])] for x in splits)
    write_uniqs(outfile, uniqs)

def handleio(infile=None, outfile=None, **kwargs):
    """I/O handling for the script."""
    in_f = open(infile, "r")
    try:
        try:
            charged2ligands(in_f, outfile)
        except ScriptError, message:
            logging.error(message)
            return message.value
    finally:
        in_f.close()
    return 0

def main(argv):
    """Parse arguments."""
    logging.basicConfig(level=logging.INFO,
                        format="%(levelname)s: %(message)s")
    description = "Unique charged properties file to generate ligands.charge."
    usage = "%prog [options]"
    version = "%prog: version 201205 - created by Michael Mysinger"
    parser = OptionParser(usage=usage, description=description,
                          version=version)
    parser.set_defaults(infile="charged.ism", outfile="ligands.charge")
    parser.add_option("-i", "--infile", 
                      help="input file (default: %default)")
    parser.add_option("-o", "--outfile", 
                      help="output file (default: %default)")
    options, args = parser.parse_args(args=argv[1:])
    if len(args):
        parser.error("program takes no positional arguments.\n" +
                     "  Use --help for more information.")
    return handleio(infile=options.infile, outfile=options.outfile)

if __name__ == "__main__":
    sys.exit(main(sys.argv))

