#!/bin/env python
"""Generate databases for charge outliers.

Michael Mysinger 200709 Created
Michael Mysinger 200812 Rewritten to use duddb backend
"""

import os
import sys
import random
import time
from optparse import OptionParser
import dudtools
import dudzinc
import duddb
import mmmutils

# Parameters
NUM_OUTLIERS = 500
CHARGES = [-4, -3, -2, +2, +3, +4]
OL_DIR = 'outliers'

def query_outliers(c, charge):
    sql = 'select mwt, xlogP, rb, n_h_acceptors, n_h_donors, net_charge, smiles, sub_id_fk, prot_id from protomer where ph_mod_fk=0 and net_charge=%s'
    c.execute(sql, charge)
    result = list(c.fetchall())
    uniq = dict(([x[:6], x[6:]] for x in result))
    return uniq

def select_outliers(inlist, numoutliers=NUM_OUTLIERS):
    try:
        picked = random.sample(inlist, numoutliers)
    except ValueError, mesg:
        print "Value Error: %s for charge %d!" % (mesg, charge)
        picked = inlist
    return picked

def make_outlier_db(ftp, pp, sdir, dbdir, genol, name, outliers, 
                    numoutliers=NUM_OUTLIERS, loglevel=2):
    picked = select_outliers(outliers.values(), numoutliers=numoutliers)
    smiles = [[dsmi, "C%08d" % dzid, "P%08d" % dpid] for (dsmi, dzid, dpid) 
                                                       in picked]
    if loglevel > 1:
        # Output redundant smiles file to search directory if it exists
        fn = os.path.join(sdir, "outliers." +  name + ".picked")
        mmmutils.write_splits(fn, smiles)
    smilesfile = os.path.join(genol, name+".smi")
    mmmutils.write_splits(smilesfile, smiles)
    odir = os.path.join(genol, name)
    dbname = os.path.join(dbdir, 'out_' + name)
    insmi = ('\t'.join(str(y) for y in x)+'\n' for x in smiles)
    duddb.make_database(pp, dbname, odir, insmi, prot="ref", ftp=ftp)
    return picked

def setup_threads():
    main_thread = dudtools.EmptyContainer()
    main_thread.num = 17
    dudtools.init_MySQL_Thread(main_thread)
    return main_thread

def setup_dirs(outdir='.'):
    outgen = os.path.join(outdir, duddb.GEN_DIR)
    genol = os.path.join(outgen, OL_DIR)
    for x in [outdir, outgen, genol]:
        if not os.path.exists(x):
            os.mkdir(x)
    return outdir, outgen, genol

def generate_outliers(outdir='.', charges=CHARGES, numoutliers=NUM_OUTLIERS, 
                      loglevel=3):
    # change above after debugging to loglevel=1
    main_thread = setup_threads()
    ftp, pp = duddb.setup_threads()
    outdir, sdir = dudzinc.setup_dirs(outdir, loglevel=loglevel)
    outdir, outgen, genol = setup_dirs(outdir)
    try:
        zids = []
        for charge in charges:
            target = "charge%+1d" % charge
            queried = query_outliers(main_thread.c, charge)
            if loglevel > 2:
                smiles = [[dsmi, "C%08d" % dzid, "P%08d" % dpid] 
                              for (dsmi, dzid, dpid) in queried.values()]
                fn = os.path.join(sdir, "outliers." +  target + ".queried")
                mmmutils.write_splits(fn, smiles)
            picked = make_outlier_db(ftp, pp, sdir, outdir, genol, 
                                     target, queried, numoutliers=numoutliers, 
                                     loglevel=loglevel)
            zids.extend("C%08d" % int(x[1]) for x in picked)
            ftp.poll()
        for res in pp:
            res.close()
        ftp.wait()
    finally:
        print "Shutting Down."
        ftp.shutdown()
    uzids = dudtools.unique(zids)
    uzids.sort()
    f = open(os.path.join(outdir, 'outliers.name'), 'w')
    for x in uzids:
        f.write(x+'\n')
    f.close()

def main(argv):
    """Parse arguments."""
    description = "Generate databases for charge outliers."
    usage = "%prog [options]"
    version = "%prog: version 200812 - created by Michael Mysinger"
    parser = OptionParser(usage=usage, description=description,
                          version=version)
    options, args = parser.parse_args(args=argv[1:])
    if len(args):
        parser.error("program takes no positional arguments." +
                     "  Use --help for more information.")
    start_time = time.time()
    generate_outliers()
    print 'Program took %.1f seconds.' % (time.time() - start_time)
    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv))
