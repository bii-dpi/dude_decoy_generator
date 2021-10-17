#!/bin/env python
"""Generate databases given potential decoys.

Michael Mysinger 200804 Created
Michael Mysinger 201109 Share command line code, add more options
"""

import os
import sys
import fnmatch
import random
import subprocess
import time
from optparse import OptionParser
import disk_server
import dudtools
import dudzinc
import mmmutils
import processpool
import threadpool

# Scripts
DB_START_SCRIPT = os.path.join(os.environ['DOCK_BASE'], 'etc', 'dbstart.csh')

# Constants
DECOY_SIZE = 50
SMI_SIZE = 1000
DB_SIZE = 4000

# Concurrency Constants
NUM_PROCESSES = 4
NUM_DISK_THREADS = 2

# File Names
DB_PREFIX = 'dbgen_'

GEN_DIR = 'dbgen'
LIG_DIR = 'ligands'
DEC_DIR = 'decoys'

LIG_SMI = 'ligands.smi'
LIG_NAME = 'ligands.name'
DEC_SMI = 'decoys.smi'
DEC_NAME = 'decoys.name'

def read_decoys(indir, glob='decoys.*.filtered'):
    """Read decoys into memory from input directory."""
    decoys = {}
    # decoy format is {(lsmi, lzid, lpid) -> [(smis, dzid, dpid)]}
    names = os.listdir(indir)
    for fn in fnmatch.filter(names, glob):
        diter = dudtools.read_decoys(os.path.join(indir, fn))
        ligand = diter.next()
        decoys[ligand] = list(diter)
    return decoys

def select_decoys(decoys, sdir, numdecoys=DECOY_SIZE, luid_label="P", 
                    loglevel=2):
    """Randomly pick numdecoys decoys from possibles."""
    picked = {}
    for ltup, dlist in decoys.iteritems():
        ddict = dict((x[2], x) for x in dlist) 
        pids = ddict.keys()
        try:
            rpids = random.sample(pids, numdecoys)
        except ValueError:
            print ( "Too few decoys available for C%08d %s%08d!" %
                                           (ltup[1], luid_label, ltup[2]) )
            rpids = pids
        rpids.sort(key=int)
        picked[ltup] = [ddict[x] for x in rpids]

    if loglevel > 1:
        # write picked decoys to search directory
        for ltup, dlist in picked.iteritems():
            mid = "%s%08d" % (luid_label, ltup[2])
            fn = os.path.join(sdir, 'decoys.' + mid + '.picked')
            dudtools.write_decoys(fn, ltup, dlist)
    return picked

def fix_smiles(ismi):
    """Smiles file should be (smiles, id)."""
    num = 1
    for line in ismi:
        splits = line.split()
        if len(splits) > 2:
            splits = splits[0:2]
        elif len(splits) == 0:
            continue
        elif len(splits) < 2:
            tempid = "G" + str(num).zfill(8)
            splits.append(tempid)
            num += 1
        yield "\t".join(splits) + "\n"

def split_smiles(ismi, tdir, smisplitsize=SMI_SIZE):
    """Split smiles file into chunks in separate subdirectories."""
    if not os.path.exists(tdir):
        os.mkdir(tdir)
    fnum = 0
    osmi = None
    for i, l in enumerate(ismi):
        if i % smisplitsize == 0:
            if osmi is not None:
                osmi.close()
            fnum += 1
            dirname = os.path.join(tdir, str(fnum).zfill(2))
            if not os.path.exists(dirname):
                os.mkdir(dirname)
            filename = '%s.smi' % (str(fnum).zfill(2))
            osmi = open(os.path.join(dirname, filename), 'w')
        osmi.write(l)
    osmi.close()
    # wait for nfs to catch up (is this necessary?)
    time.sleep(2)

def start_dbgen(pp, tdir, prot="mid"):
    """Start database generation on SGE cluster."""
    for sdir in os.listdir(tdir):
        fdir = os.path.join(tdir, sdir)
        if not os.path.isdir(fdir):
            continue
        print 'Starting database generation in %s.' % fdir
        outf = open(os.path.join(fdir, 'dbgen.log'), 'w')
        filename = os.path.join(fdir, sdir+'.smi')
        # is the shell necessary here?
        sub = subprocess.Popen(' '.join([DB_START_SCRIPT, filename, prot]),
                               stdout=outf, shell=True)
        res = pp.run(sub, outf)
        if res is not None:
            res.close()

def collect_results(mythread, dbname, *args):
    """Collect database on remote XML-RPC disk server."""
    print 'Starting database collection for %s.' % dbname 
    mythread.client.collect_results(dbname, *args)
    print 'Finished database collection for %s.' % dbname

def make_database(pp, dbname, gendir, insmiles, prot="mid", ftp=None):
    """Start and collect a given database."""
    split_smiles(fix_smiles(insmiles), gendir)
    start_dbgen(pp, gendir, prot=prot)
    if ftp is None:
        disk_server.collect_results(dbname, gendir, DB_SIZE)
    else:
        ftp.put(threadpool.Request(collect_results, (dbname, gendir, DB_SIZE), 
                                 sendself=True))

def make_decoy_db(pp, dbdir, gendec, picked, dbprefix=DB_PREFIX, 
                  ftp=None):
    """Generate decoy database."""
    smiles = [[dsmi, "C%08d" % dzid, "P%08d" % dpid] 
                      for dlist in picked.values() 
                          for (dsmi, dzid, dpid) in dlist]
    smilesfile = os.path.join(gendec, DEC_SMI)
    mmmutils.write_splits(smilesfile, smiles)
    dbname = os.path.join(dbdir, dbprefix + 'dec')
    insmi = ('\t'.join(str(y) for y in x)+'\n' for x in smiles)
    make_database(pp, dbname, gendec, insmi, prot="ref", ftp=ftp)

def setup_threads():
    """Setup threads and remote proxies."""
    ftp = threadpool.RemotePool(NUM_DISK_THREADS, dudtools.DISK_HOST, 
                                dudtools.DISK_PORT, dudtools.DISK_SCRIPT)
    pp = processpool.ProcessPool(NUM_PROCESSES)
    time.sleep(5)
    return ftp, pp

def setup_dirs(outdir='.'):
    """Setup output directories."""
    outdir = os.path.abspath(outdir)
    outgen = os.path.join(outdir, GEN_DIR)
    genlig = os.path.join(outgen, LIG_DIR)
    gendec = os.path.join(outgen, DEC_DIR)
    for x in (outdir, outgen, genlig, gendec):
        if not os.path.exists(x):
            os.mkdir(x)
    return outdir, outgen, genlig, gendec

def generate_database(indir='.', outdir='.', 
                      loglevel=2, numdecoys=DECOY_SIZE, 
                      decoys=None, luid_label="P"):
    """Generate databases given potential decoys.."""
    outdir, sdir = dudzinc.setup_dirs(outdir, loglevel)
    outdir, outgen, genlig, gendec = setup_dirs(outdir)
    ftp, pp = setup_threads()
    try:
        if decoys is None:
            decoys = read_decoys(sdir)
        picked = select_decoys(decoys, sdir, numdecoys=numdecoys, 
                               luid_label=luid_label, loglevel=loglevel)
        make_decoy_db(pp, outdir, gendec, picked, ftp=ftp)
        for res in pp:
            res.close()
        ftp.wait()
    finally:
        print "Shutting Down."
        ftp.shutdown()

def option_parser(usage, description, version):
    parser = dudzinc.option_parser(usage, description, version)
    parser.set_defaults(numdecoys=DECOY_SIZE, skipdb=False)
    parser.add_option("-n", "--numdecoys", type="int",
           help="number of decoys per ligand (default: %default)")
    parser.add_option("-s", "--skip-database", action="store_true", 
           dest="skipdb", help="skip DOCK database generation for decoys")
    return parser

def parse(parser, argv):
    return dudzinc.parse(parser, argv)

def main(argv):
    """Parse arguments."""
    description = "Generate databases given potential decoys."
    usage = "%prog [options]"
    version = "%prog: version 201109 - created by Michael Mysinger"
    parser = option_parser(usage, description, version)
    parser.remove_option("--tanimoto-cutoff")
    parser.remove_option("--fingerprint-server")
    parser.remove_option("--skip-database")
    options = parse(parser, argv)

    start_time = time.time()
    generate_database(indir=options.indir, outdir=options.outdir, 
                      loglevel=options.log_level, numdecoys=options.numdecoys)
    print 'Program took %.1f seconds.' % (time.time() - start_time)
    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv))

