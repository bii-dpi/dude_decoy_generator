#!/bin/env python
"""Compute daylight tanimoto matrix and the maximum tanimoto for the query molecules over the target molecules.

Michael Mysinger 201103 Created
"""

import os
import sys
import logging
from optparse import OptionParser
import dudzinc
import dudtools
import tanimoto

def dostuff(qf, tf, outf=None, 
            fp_server=dudtools.FP_SERVER_OPTIONS["ecfp4"]):
    """Module code starts here."""
    main_thread = dudtools.EmptyContainer()
    main_thread.num = 17
    dudtools.init_FP_Thread(main_thread, fp_server)
    qraw = (x.split() for x in qf)
    qsmi = [(x[1], x[0]) for x in qraw]
    qfp, qbc = dudzinc.load_smiles(main_thread, qsmi)
    traw = (x.split() for x in tf)
    tsmi = [(x[1], x[0]) for x in traw]
    tfp, tbc = dudzinc.load_smiles(main_thread, tsmi)
    tm = tanimoto.matrix(qfp, qbc, tfp, tbc)
    outf.write('%s\n' %'\t'.join(['QueryID', 'RefID', 'Tc', 'QuerySMI', 'RefSMI']))
    for i in xrange(len(tm)):
        max_tc = max(tm[i])
        idx = tm[i].index(max_tc)
        outf.write('%s\n' %'\t'.join([qsmi[i][0], tsmi[idx][0], str(max_tc), qsmi[i][1], tsmi[idx][1]]))
    return 0

def handleio(queryfile, targetfile, outfile=None, **kwargs):
    """I/O handling for the script."""
    qf = open(queryfile, "r")
    tf = open(targetfile, "r")
    if outfile is None:
        outf = sys.stdout
    else:
        outf = open(outfile, "w")
    result = dostuff(qf, tf, outf, **kwargs)
    qf.close()
    tf.close()
    outf.close()
    return result

def main(argv):
    """Parse arguments."""
    logging.basicConfig(level=logging.INFO,
                        format="%(levelname)s: %(message)s")
    description = "Compute daylight tanimoto matrix and the maximum tanimoto for the query molecules over the target molecules."
    usage = "%prog [options] <query_file> <target_file>"
    version = "%prog: version 201103 - created by Michael Mysinger"
    parser = OptionParser(usage=usage, description=description,
                          version=version)
    parser.set_defaults(outfile=None, fp_server="ecfp4")
    parser.add_option("-o", "--outfile", 
                      help="output file (default: stdout)")
    fp_types = ', '.join(dudtools.FP_SERVER_OPTIONS.keys())
    parser.add_option("-f", "--fingerprint-server", dest="fp_server", 
           help="fingerprint server: options = %s" % fp_types +
                     " (default: %default)")
    options, args = parser.parse_args(args=argv[1:])
    if len(args) != 2:
        parser.error("program takes two positional arguments.\n" +
                     "  Use --help for more information.")
    options.fp_server = dudtools.FP_SERVER_OPTIONS[options.fp_server]
    return handleio(args[0], args[1], outfile=options.outfile, 
                    fp_server=options.fp_server)

if __name__ == "__main__":
    sys.exit(main(sys.argv))

