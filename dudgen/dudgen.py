#!/bin/env python
"""Given a ligand list, generate DUD style decoys and DOCK databasess.

Michael Mysinger 200804 Created
Michael Mysinger 201109 Share command-line, default to ECFP4 & 75% cut
"""

import os
import sys
import time
from optparse import OptionParser
from os.path import join, exists, isdir, basename, splitext
import duddb
import dudtools
import dudzinc
import mmmutils
import processpool

def read_properties(filename):
    """Read first five DUD properties from mitools output."""
    f = open(filename, 'r')
    # only need first line
    spl = f.readline().split('\t')
    f.close()
    # grab smiles, [mwt, logP, rb, hba, hbd]
    indices = (5, 2, 9, 6, 7)
    casts = (float, float, int, int, int)
    if len(spl) < 10:
        return None, None
    else:
        smi = spl[0].replace('\\\\', '\\')
        return smi, [cast(spl[i]) for cast, i in zip(casts, indices)]

def read_charge(filename):
    """Read charge from AMSOL output."""
    f = open(filename, 'r')
    spl = f.readline().split()
    f.close()
    return int(round(float(spl[2])))

def read_dict(inpath):
    """Read ids from translation dict."""
    splits = mmmutils.read_splits(join(inpath, 'dict'))
    iddict = dict((x[0], x[1]) for x in splits)
    return iddict

def collect_properties(subdir, sleep_time=15):
    """Collect ligand protomers in subdir."""
    while not exists(join(subdir, 'started')):
        time.sleep(sleep_time)
    protomers = {}
    subid = basename(subdir)[-1]
    for i, prot in enumerate(('ref', 'mid', 'hi', 'lo')):
        pdir = join(subdir, prot)
        if isdir(pdir):
            protid = str(9 - i)
            dirlist = open(join(pdir, 'dirlist'), 'r')
            working = [basename(line.strip()) for line in dirlist]
            dirlist.close()
            while working:
                waiting = []
                for tempname in working:
                    tempdir  = join(pdir, tempname)
                    # Check for completion in log file
                    logfn = join(tempdir, 'sge.log')
                    if not exists(logfn):
                        waiting.append(tempname)
                        continue
                    logf = open(logfn, 'r')
                    log = logf.read()
                    logf.close()
                    if "Run complete!" not in log:
                        waiting.append(tempname)
                        continue
                    # Run complete, check for errors
                    if "ailed in AMSOL" in log:
                        print "Skip: AMSOL failed for %s" % tempdir
                        continue
                    elif "ailed in mitools" in log:
                        print "Skip: mitools failed for %s" % tempdir
                        continue
                    elif "ailed in corina" in log:
                        print "Skip: corina failed for %s" % tempdir
                        continue
                    # Run complete, wait for files to be copied back from /tmp
                    for dummy in xrange(30):
                        if ( exists(join(tempdir, 'db.db.bz2')) and
                             exists(join(tempdir, "mitools.ism")) and
                             exists(join(tempdir, "solv.log")) ):
                            break
                        time.sleep(10)
                    else:
                        print ( "Skip: Timed out waiting for files in %s"
                                                                % tempdir )
                        continue
                    prefix = splitext(tempname)[0][:12]
                    prot_index = str(len(tempname)-17).zfill(2)
                    upid = int(protid + subid + prefix[8:12] + prot_index)
                    fname = join(tempdir, prefix)
                    smis, props = read_properties(join(tempdir, "mitools.ism"))
                    if smis is not None:
                        plist = [upid, smis]
                        plist.extend(props)
                        plist.append(read_charge(join(tempdir, "solv.log")))
                        protomers.setdefault(prefix, []).append(plist)
                    else:
                        print "Skip: Empty mitools.ism collecting %s" % tempdir
                working = waiting
                time.sleep(sleep_time)
    return protomers

def collect_ligands(gendir, sleep_time=15):
    """Collect and process ligands after database generation."""
    uniqs = {}
    for sdir in os.listdir(gendir):
        fdir = join(gendir, sdir)
        if isdir(fdir):
            print "Working in %s" % fdir
            protomers = collect_properties(fdir, sleep_time)
            iddict = read_dict(fdir)
            for k, prots in protomers.iteritems():
                uniq = dict((tuple(x[2:]), (x[0], x[1])) for x in prots)
                start, end = dudtools.get_id_position(iddict[k])
                try:
                    lzid = int(iddict[k][start:end])
                except ValueError, e:
                    print "Error: %s" % e
                    print "\nLikely caused by my rudimentary handling of molecule ids."
                    print "They should look like TEMP12345678, P12345678, or 12345678"
                    print "Sorry..."
                    sys.exit(5)
                uniqs[lzid] = uniq
    return uniqs

def generate_decoy_database(infile=None, outdir='.', loglevel=2,
                            fraction=dudzinc.FRACTION_TO_CUT,
                            fp_server=dudtools.FP_SERVER_OPTIONS["ecfp4"],
                            numdecoys=duddb.DECOY_SIZE, skipdb=False,
                            ligand_prot="mid", charge_start=False,
                            remote_disk=False):
    # Possible future optimizations
    # 1) build protomers and their properties ASAP
    #      a) change build procedure to run mitools early for every ligand
    # 2) minimize hits to disk
    #      a) change database generation - pipes?

    if infile is None:
        inf = sys.stdin
    else:
        inf = open(infile, 'r')
    outdir, sdir = dudzinc.setup_dirs(outdir, loglevel)
    outdir, outgen, genlig, gendec = duddb.setup_dirs(outdir)
    main_thread, mysqlp, fpp, cpup, filep = dudzinc.setup_threads(fp_server)
    if remote_disk:
        ftp, pp = duddb.setup_threads()
    else:
        ftp = None
        pp = processpool.ProcessPool(duddb.NUM_PROCESSES)
    try:
        if charge_start:
            print "Restart, skipping ligand database generation."
            uniqs = dudzinc.read_uniqs(inf)
        else:
            dbname = duddb.DB_PREFIX + os.path.basename(genlig)[:3]
            dbname = join(outdir, dbname)
            duddb.make_database(pp, dbname, genlig, inf, prot=ligand_prot,
                                ftp=ftp)
            uniqs = collect_ligands(genlig)
        fn = dudzinc.CHARGE_FILE_LIG + dudzinc.CHARGE_FILE_EXT
        propfile = os.path.join(outdir, fn)
        lfps, lbcs = dudzinc.process_ligands(main_thread, propfile, uniqs,
                                             loglevel=loglevel)
        decoys = dudzinc.query_decoys(mysqlp, fpp, cpup, filep, sdir,
                     uniqs, lfps, lbcs, loglevel=loglevel, fraction=fraction,
                     luid_label="F")
        picked = duddb.select_decoys(decoys, sdir, numdecoys=numdecoys,
                                       luid_label='F', loglevel=loglevel)
        if not skipdb:
            duddb.make_decoy_db(pp, outdir, gendec, picked, ftp=ftp)
        for res in pp:
            res.close()
        if ftp:
            ftp.wait()
    finally:
        print "Shutting Down."
        dudtools.cleanup_MySQL_Thread(main_thread)
        mysqlp.shutdown()
        fpp.shutdown()
        cpup.shutdown()
        filep.shutdown()
        if ftp:
            ftp.shutdown()

def main(argv):
    """Parse arguments."""
    description = ( "Given a ligand list, generate DUD style decoys" +
                    " and DOCK flexibases. " +
               "ligand list is a file with a smiles and an ID on each line." +
               "ID must be formatted in C12345678 style, one letter, 8 digits")
    usage = "%prog [options]"
    version = "%prog: version 201109 - created by Michael Mysinger"
    parser = duddb.option_parser(usage, description, version)
    parser.set_defaults(prot="mid", charge_start=False, mysinger=False)
    parser.add_option("-p", "--protonation", dest="prot",
           help=("protonation type for ligand databases: " +
                 "options = ref, mid, lo, hi, all  (default: %default)"))
    parser.add_option("-m", "--mysinger-disk", action="store_true",
           help="use (fragile) remote disk servers on raid5")
    options = duddb.parse(parser, argv)

    start_time = time.time()
    generate_decoy_database(infile=options.infile, outdir=options.outdir,
        loglevel=options.log_level, fraction=options.fraction,
        fp_server=options.fp_server,
        numdecoys=options.numdecoys, skipdb=options.skipdb,
        ligand_prot=options.prot, charge_start=options.charge_start,
        remote_disk=options.mysinger_disk)
    print 'Program took %.1f seconds.' % (time.time() - start_time)
    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv))

