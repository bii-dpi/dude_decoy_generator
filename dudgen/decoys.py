#!/bin/env python
"""Given a list of ligand ZINC ids, generate DUD style decoys and 
       optional DOCK databases.

Michael Mysinger 200804 Created
Michael Mysinger 200808 Modified for web tool
Michael Mysinger 201109 Share command-line, read ligands.charge
"""

import os
import sys
import time
from optparse import OptionParser
import disk_server
import duddb
import dudtools
import dudzinc

def collect_zinc_dbs(outdb, pids, dbsplitsize=duddb.DB_SIZE):
    f = open(disk_server.DBHEADER, 'r')
    header = f.readlines()
    f.close()
    dbnum = 1
    dbcount = 0
    odb = disk_server.start_db(dbnum, outdb, header)
    for pid in pids:
        spid = str(pid).zfill(8)
        indb = os.path.join(dudtools.ZINCDIR, spid[4:6], spid[6:8],
                                str(pid)+'.db.bz2')
        dbnum, dbcount, odb = disk_server.append_db(dbnum, dbcount, odb, 
                                  indb, outdb, header, dbsplitsize)

def quick_decoy_database(infile=None, outdir='.', loglevel=1, 
                         fraction=dudzinc.FRACTION_TO_CUT, 
                         fp_server=dudtools.FP_SERVER_OPTIONS["ecfp4"], 
                         numdecoys=duddb.DECOY_SIZE, skipdb=False,  
                         column=0, charge_start=False):
    decoys = dudzinc.create_potential_decoys(infile, outdir, 
                 loglevel=loglevel, fraction=fraction, column=column,
                 charge_start=charge_start, fp_server=fp_server)
    print "Making decoy database."
    picked = duddb.select_decoys(decoys, outdir, numdecoys=numdecoys)
    if not skipdb:
        dpids = [d[2] for dlist in picked.itervalues() for d in dlist]
        outdb = os.path.join(outdir, duddb.DB_PREFIX + 'dec')
        collect_zinc_dbs(outdb, dpids)

def main(argv):
    """Parse arguments."""
    description = ( "Given a list of ligand ZINC ids, generate DUD " + 
                    "style decoys and optional DOCK databases." )
    usage = "%prog [options]"
    version = "%prog: version 201109 - created by Michael Mysinger"
    parser = duddb.option_parser(usage, description, version)
    parser.set_defaults(log_level=1, column=1, charge_start=False)
    parser.add_option("--column", type="int", 
           help="input file column containing ids (default: %default)")
    options = duddb.parse(parser, argv)

    start_time = time.time()
    quick_decoy_database(infile=options.infile, outdir=options.outdir, 
                         loglevel=options.log_level, fraction=options.fraction, 
                         fp_server=options.fp_server, 
                         numdecoys=options.numdecoys, skipdb=options.skipdb, 
                         column=(options.column-1),
                         charge_start=options.charge_start)
    print 'Program took %.1f seconds.' % (time.time() - start_time)
    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv))
