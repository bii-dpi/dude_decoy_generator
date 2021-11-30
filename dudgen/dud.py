#!/usr/arch/bin/python -u
"""Generate DUD decoys for every DUD target.

Michael Mysinger 200805 Rewritten as wrapper around dudgen.py 
Michael Mysinger 201109 Share command-line, default to ECFP4 & 75% cut
"""

import os
import sys
import fnmatch
import shutil
import subprocess
import time
from optparse import OptionParser
import duddb
import dudgen
import dudtools
import dudzinc
import mmmutils
import processpool
import threadpool

def get_smiles(c, zids):
    """Get smiles using zids by querying ZINC database."""
    smiles_sql = 'select smiles, sub_id from substance where sub_id in %s'
    nzids = dudtools.translate_retired_zids(c, zids)
    smis = dudtools.one_to_one_query(c, smiles_sql, nzids)
    err = (i for i, s in enumerate(smis) if s is None)
    if err:
        for i in err:
            if zids[i] != nzids[i]:
                print ( "Ligand id %d has new id %d, " % (zids[i], nzids[i]) + 
                          "but the new smiles was not found!!!")
            else:
                print "Ligand id %d not found!!!" % zids[i]
    smis = [s for s in smis if s is not None]
    return smis

def read_enzymes(dirname):
    """Read enzyme targets from database directory."""
    dir_files = os.listdir(dirname)
    name_files = fnmatch.filter(dir_files, '*.name')
    enzymes = [x.split('.')[0] for x in name_files]
    enzymes.sort()
    return enzymes

def read_uniqs(filename):
    uniqs = {}
    # smi, lzid, luid, mwt, logP, rb, hba, hbd, net
    splits = mmmutils.read_splits(filename)
    typed = ([x[0], int(x[1][1:9]), int(x[2][1:9]), float(x[3]), float(x[4]),
              int(x[5]), int(x[6]), int(x[7]), int(x[8])] for x in splits)
    for x in typed:
        uniqs.setdefault(x[1], {}).update([[tuple(x[3:]), (x[2], x[0])]])
    return uniqs

def write_dot_name(outdir, target, zids, legacy=True):
    if legacy:
        if target == 'na':
            target = 'neua'
        elif target == 'ppar_gamma':
            target = 'ppar'
        elif target == 'rxr_alpha':
            target = 'rxr'
    zids = dudtools.unique(zids)
    zids.sort()
    splits = (["C%08d" % int(i)] for i in zids)
    fn = os.path.join(outdir, target+'.name')
    mmmutils.write_splits(fn, splits)

def count_lines(filename):
    p1 = subprocess.Popen(["zcat", filename], stdout=subprocess.PIPE)
    p2 = subprocess.Popen(["wc", "-l"], stdin=p1.stdout,
                                        stdout=subprocess.PIPE)
    return int(p2.communicate()[0].strip())

def count_generated(ddir):
    count = 0
    if os.path.exists(ddir):
        for sdir in os.listdir(ddir):
            fdir = os.path.join(ddir, sdir)
            if os.path.isdir(fdir):
                dbfile = os.path.join(fdir, "db.db.gz")
                if not os.path.exists(dbfile):
                    return None
                count += count_lines(dbfile) - 8
    return count

def count_decoys(outdir, prefix):
    count = 0
    names = os.listdir(outdir)
    decfiles = fnmatch.filter(names, prefix + "*.db.gz")
    for dfile in decfiles:
        count += count_lines(dfile) - 8
    return count

def make_ligand_db(mythread, pp, ligdir, genlig, outdir, name, ftp=None):
    ldir = os.path.join(genlig, name)
    if not os.path.exists(ldir):
        os.mkdir(ldir)
    # for now parse .name, in the future switch to smiles input
    fn = os.path.join(ligdir, name + ".name")
    lzids = dudtools.read_zids(fn, col=0)
    lzids = dudtools.unique(lzids)
    # query smiles from ZINC
    smis = get_smiles(mythread.c, lzids)
    smiles = [[x[0], "C%08d" % x[1]] for x in smis]
    smilesfile = os.path.join(ldir, name + ".smi")
    mmmutils.write_splits(smilesfile, smiles)
    dbname = "dud_" + name + "_lig"
    dbname = os.path.join(outdir, dbname)
    insmi = ('\t'.join(str(y) for y in x)+'\n' for x in smiles)
    duddb.make_database(pp, dbname, ldir, insmi, ftp=ftp)
    return ldir

def generate_decoys(indir=dudtools.DUDDIR, outdir='.', loglevel=2, 
                    fraction=dudzinc.FRACTION_TO_CUT, 
                    fp_server=dudtools.FP_SERVER_OPTIONS["ecfp4"], 
                    numdecoys=duddb.DECOY_SIZE, skipdb=False, 
                    restart=False):

    # ensure loglevel of 2 for restart and .name generation
    loglevel = max(loglevel, 2)
    
    ligdir = os.path.join(indir, 'ligands')
    enzymes = read_enzymes(ligdir)

    outdir, sdir = dudzinc.setup_dirs(outdir, loglevel)
    outdir, outgen, genlig, gendec = duddb.setup_dirs(outdir)
    main_thread, mysqlp, fpp, cpup, filep = dudzinc.setup_threads(fp_server)
    pp = processpool.ProcessPool(duddb.NUM_PROCESSES)
    try:
        lfps = []
        lbcs = []
        uniqdict = {}
        if not restart:
            for name in enzymes:
                ldir = make_ligand_db(main_thread, pp,
                                      ligdir, genlig, outdir, name)
                req = threadpool.Request(dudgen.collect_ligands, (ldir,))
                req.name = name
                cpup.put(req)
            for i in enzymes:
                req, uniqs = cpup.get()
                name = req.name
                print "Collecting %s ligands."  % name
                fn = ( dudzinc.CHARGE_FILE_LIG + "." + name +
                           dudzinc.CHARGE_FILE_EXT )
                fn = os.path.join(outdir, fn)
                fps, bcs = dudzinc.process_ligands(main_thread, fn, uniqs,
                               loglevel=loglevel, luid_label="F")
                lfps.extend(fps)
                lbcs.extend(bcs)
                uniqdict[name] = uniqs
                ndir = os.path.join(outdir, duddb.LIG_DIR)
                if not os.path.exists(ndir):
                    os.mkdir(ndir)
                zids = uniqs.iterkeys()
                write_dot_name(ndir, name, zids)
        # alternative restart: write ascii fps, then only need local loadascii
        else:
            print "Restart requested. Analyzing previous results..."
            for name in enzymes:
                fn = ( dudzinc.CHARGE_FILE_LIG + "." + name +
                           dudzinc.CHARGE_FILE_EXT )
                fn = os.path.join(outdir, fn)
                uniqs = read_uniqs(fn)
                smis = [v for uniq in uniqs.itervalues()
                            for v in uniq.itervalues()]
                fps, bcs = dudzinc.load_smiles(main_thread, smis)
                lfps.extend(fps)
                lbcs.extend(bcs)
                uniqdict[name] = uniqs
            incomplete = []
            for name in enzymes:
                ddir = os.path.join(gendec, name)
                numgen = count_generated(ddir)
                dbprefix = "dud_" + name + "_dec"
                numdec = count_decoys(outdir, dbprefix)
                if numgen and numgen == numdec:
                    print "%s is complete." % name
                else:
                    print "%s is incomplete." % name
                    incomplete.append(name)
            for name in incomplete:
                print "Cleaning up any %s remnants." % name
                ddir = os.path.join(gendec, name)
                if os.path.exists(ddir):
                    shutil.rmtree(ddir)
                tdir = os.path.join(sdir, name)
                if os.path.exists(tdir):
                    shutil.rmtree(tdir)
                dbprefix = "dud_" + name + "_dec"
                names = os.listdir(outdir)
                decfiles = fnmatch.filter(names, dbprefix + "*.db.gz")
                for dfile in decfiles:
                    os.remove(dfile)
            enzymes = incomplete
        for name in enzymes:
            print "Querying %s decoys."  % name
            ddir = os.path.join(gendec, name)
            tdir = os.path.join(sdir, name)
            ndir = os.path.join(outdir, duddb.DEC_DIR)
            for x in (ddir, tdir, ndir):
                if not os.path.exists(x):
                    os.mkdir(x)
            decoys = dudzinc.query_decoys(mysqlp, fpp, cpup, filep, tdir, 
                uniqdict[name], lfps, lbcs, loglevel=loglevel, 
                fraction=fraction, luid_label="F")
            dbprefix = "dud_" + name + "_"
            picked = select_decoys(decoys, tdir, numdecoys=numdecoys, 
                                   luid_label='F', loglevel=loglevel)
            if not skipdb:
                duddb.make_decoy_db(pp, outdir, ddir, picked, 
                                    dbprefix=dbprefix)
            zids = (zid for dlist in picked.itervalues()
                           for (smi, zid, pid) in dlist)
            write_dot_name(ndir, name, zids)
            #ftp.poll()
        print "Waiting for all databases."
        for res in pp:
            res.close()
        #ftp.wait()
    finally:
        print "Shutting down."
        dudtools.cleanup_MySQL_Thread(main_thread)
        #ftp.shutdown()
        mysqlp.shutdown()
        fpp.shutdown()
        cpup.shutdown()
        filep.shutdown()
    
def main(argv):
    """Parse arguments."""
    description = "Generate DUD decoys for every DUD target."
    usage = "%prog [options]"
    version = "%prog: version 201109 - created by Michael Mysinger"
    parser = duddb.option_parser(usage, description, version)
    parser.set_defaults(restart=False)
    parser.add_option("-r", "--restart", action="store_true", 
           help="restart a failed run")
    options = duddb.parse(parser, argv)

    start_time = time.time()
    generate_decoys(indir=options.indir, outdir=options.outdir,
                    loglevel=options.log_level, fraction=options.fraction,
                    fp_server=options.fp_server, 
                    numdecoys=options.numdecoys, skipdb=options.skipdb, 
                    restart=options.restart)
    print 'Program took %.1f seconds.' % (time.time() - start_time)
    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv))
