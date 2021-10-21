#!/bin/env python
"""Test decoy generation code.

Michael Mysinger 201109 Created
"""

import os
import sys
import logging
from optparse import OptionParser
import dudtools
import dudzinc
import tanimoto

def test_decoys(ligand_fn, decoy_fn):
    """Decoy generation test code."""
    lfile = open(ligand_fn)
    uniqs = dudzinc.read_uniqs(lfile)
    lfile.close()
    main_thread = dudtools.EmptyContainer()
    main_thread.num = 17
    dudtools.init_MySQL_Thread(main_thread)
    dudtools.init_FP_Thread(main_thread, dudtools.FP_SERVER_OPTIONS["ecfp4"])
    lfps, lbcs = dudzinc.process_ligands(main_thread, "trash", uniqs,
                                         loglevel=0)
    decoys = {}
    diter = dudtools.read_decoys(decoy_fn)
    ligand = diter.next()
    decoys[ligand] = list(diter)
    smis = [[dpid, smi.replace('\\\\', '\\')] for (smi, dzid, dpid) 
                in decoys[ligand]]
    res = dudzinc.get_fingerprints(main_thread, smis)
    asciifps = [fps for (dpid, fps) in res]
    dfps, dbcs = tanimoto.loadascii(asciifps)
    indices = tanimoto.fractionfilterB(0.50, lfps, lbcs, dfps, dbcs)
    import pdb
    pdb.set_trace()
    return 0

def main(argv):
    """Parse arguments."""
    logging.basicConfig(level=logging.INFO,
                        format="%(levelname)s: %(message)s")
    description = "Test decoy generation code."
    usage = "%prog [options] ligands.charge decoys.queried"
    version = "%prog: version 201109 - created by Michael Mysinger"
    parser = OptionParser(usage=usage, description=description,
                          version=version)
    options, args = parser.parse_args(args=argv[1:])
    if len(args) != 2:
        parser.error("program takes two positional arguments.\n" +
                     "  Use --help for more information.")
    return test_decoys(args[0], args[1])

if __name__ == "__main__":
    sys.exit(main(sys.argv))

