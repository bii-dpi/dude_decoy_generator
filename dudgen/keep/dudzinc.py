
# File Names
SEARCH_DIR = 'search'
CHARGE_FILE_LIG = 'ligands'
CHARGE_FILE_EXT = '.charge'

# Constants
FRACTION_TO_CUT = 0.75

# Tuning Constants for decoy list size
UPPER_LIMIT = 9000
LOWER_LIMIT = 3000

# Threading Constants
NUM_FP_THREADS = 3
NUM_SQL_THREADS = 4
NUM_CPU_THREADS = 6

# Property Ranges for SQL Queries
CHG_RANGES =  [  0,   0,   0,   0,   0,   1,   2]
NHD_RANGES =  [  0,   0,   1,   1,   2,   2,   3]
NHA_RANGES =  [  0,   1,   2,   2,   3,   3,   4]
RB_RANGES =   [  1,   2,   2,   3,   3,   4,   5]
MWT_RANGES =  [ 20,  35,  50,  65,  80, 100, 125]
LOGP_RANGES = [0.4, 0.8, 1.2, 1.8, 2.4, 3.0, 3.6]

# SQL statements
RETIRED_SQL = "select new_id from retired where old_id=%s"
LIGAND_SQL = ( "select prot_id, smiles, mwt, xlogP, rb, " + 
               "n_h_acceptors, n_h_donors, net_charge from protomer " + 
               "where sub_id_fk=%s and ph_mod_fk in (0, 1)" )
# decoyprot is a protomer copy containing only protomer ids that are
#   free, purchasable, and ph_mod_fk = 0
END_DECOY_SQL = ( " from decoyprot where net_charge in %s and " + 
                  "n_h_donors in %s and n_h_acceptors in %s and rb in %s " + 
                  "and mwt>=%%s and mwt<=%%s and xlogP>=%%s and xlogP<=%%s " + 
                  "limit %%s" )
TEST_DECOY_SQL = "select count(*)" + END_DECOY_SQL
FULL_DECOY_SQL = "select smiles, sub_id_fk, prot_id" + END_DECOY_SQL

def get_retired_zid(c, lzid):
    """Get new id from retired table if available."""
    rnum = c.execute(RETIRED_SQL, (lzid))
    if rnum == 0:
        return lzid
    elif rnum == 1:
        return c.fetchone()[0]
    else:
        raise dudtools.DudError("Ligand %d has ambiguous retired " + 
                                "table entries." % lzid)

def get_ligand_protomers(mythread, lzid):
    """Query ligand for protomers and unique properties."""
    c = mythread.c
    nlzid = get_retired_zid(c, lzid)
    if c.execute(LIGAND_SQL, (nlzid)) == 0:
        if nlzid != lzid:
            raise dudtools.DudError("Ligand id %d has new id %d, " % 
                      (lzid, nlzid) + "but the new id was not found!!!")
        else:
            raise dudtools.DudError("Ligand id %d not found!!!" % nlzid)
    prots = c.fetchall()
    # Select only ligand protomers with unique tuples of properties
    #   { prop tuple -> (luid, smi) }
    uniq = dict([x[2:], (x[0], x[1])] for x in prots)
    return nlzid, uniq

def _in_range(start, end):
    """Construct SQL (in) expression given start and end integers."""
    values = [str(x) for x in range(int(start), int(end)+1)]
    return '(%s)' % ', '.join(values)

def _execute_idx(c, idx, sql, props):
    """Execute decoy SQL query for a given width idx."""
    mwt, logp, rb, nha, nhd, chg = props
    chg_in = _in_range(chg-CHG_RANGES[idx], chg+CHG_RANGES[idx])
    nhd_in = _in_range(nhd-NHD_RANGES[idx], nhd+NHD_RANGES[idx])
    nha_in = _in_range(nha-NHA_RANGES[idx], nha+NHA_RANGES[idx])
    rb_in = _in_range(rb-RB_RANGES[idx], rb+RB_RANGES[idx])
    c.execute(sql % (chg_in, nhd_in, nha_in, rb_in), 
                    (mwt-MWT_RANGES[idx], mwt+MWT_RANGES[idx],
                    logp-LOGP_RANGES[idx], logp+LOGP_RANGES[idx],
                    UPPER_LIMIT))

def get_decoys(mythread, uniq):
    """Query potential decoys matching properties of a unique ligand."""
    c = mythread.c
    decoys = {}  # holds decoys for each ligand protomer id
    myidx = {}   # holds ligand list and smiles
    # For each ligand protomer query ZINC for decoys
    max_index = len(CHG_RANGES) - 1
    for props, (luid, lsmi) in uniq.iteritems():
        # Adjust the query to get the right number of decoys
        # Seven levels (0-6) of query from narrowest in properties to widest
        # Start with index 2 and adjust the level from there
        idx = 2
        lidx = 2
        while True:
            _execute_idx(c, idx, TEST_DECOY_SQL, props)
            check = c.fetchone()[0]
            if check >= UPPER_LIMIT:
                # Use lidx to make sure we don't oscillate
                if idx > 0 and idx <= lidx:
                    lidx = idx
                    idx -= 1
                    continue
            elif check < LOWER_LIMIT:
                if idx < max_index:
                    lidx = idx
                    idx += 1
                    continue
                else:
                    if check == 0:
                        decoys[(luid, lsmi)] = {}
                        myidx[(luid, lsmi)] = max_index
                        break
            _execute_idx(c, idx, FULL_DECOY_SQL, props)
            #    { (luid, lsmi) -> [(dsmi, dzid, dpid)] }
            decoys[(luid, lsmi)] = c.fetchall()
            myidx[(luid, lsmi)] = idx
            break
    return myidx, decoys

def get_fingerprints(mythread, smiles):
    """Request fingperints from fp server."""
    ffunc = getattr(mythread.client, mythread.config[2])
    rsmiles = [(smi, str(uid)) for uid, smi in smiles]
    rfps = ffunc(rsmiles)
    fps = [(uid, smi) for smi, uid in rfps]
    result = []
    for i, fp in enumerate(fps):
        if not fp[1]:
            print "Error: Failed to generate fingerprint for %s." % str(smiles[i])
        else:
            result.append(fp)
    return result

def clean_decoys(decoys):
    """Remove duplicate decoys while evenly distrubuting them."""

    # dpid_dicts provides fast lookup of a decoy tuple on its dpid 
    # { (lsmi, lzid, luid) -> {dpid -> (smi, dzid, dpid)} } 
    dpid_dicts = {}
    for ltup, dlist in decoys.iteritems():
        dpid_dicts[ltup] = dict((x[2], x) for x in dlist)

    # dids reverse maps each decoy pid to a list of ligand id tuples
    # { dpid -> [(lsmi, lzid, luid)] }
    dids = {}
    for ltup, dlist in decoys.iteritems():
        for (smi, dzid, dpid) in dlist:
            dids.setdefault(dpid, []).append(ltup)

    # idcounts sorts from least duplicated decoy ids to most duplicated
    # So we can assign the decoys with the least options first
    idcounts = [(len(llist), dpid) for dpid, llist in dids.iteritems()]
    idcounts.sort()

    # To remove duplicates we start a new decoy dictionary and fill it
    # with the current decoys one at a time, assigning each decoy to
    # the ligand with the fewest decoys up to that point

    # Initialize new duplicate free decoy dictionary
    ndecoys = {}
    for ltup in decoys:
        ndecoys[ltup] = []
    for (length, dpid) in idcounts:
        # for this decoy, find which ligand list that is least full
        options = [(len(ndecoys[ltup]), ltup) for ltup in dids[dpid]]
        options.sort()
        ltup = options[0][1]
        # copy this decoy tuple to that ligand list
        ndecoys[ltup].append(dpid_dicts[ltup][dpid])
    return ndecoys

# Decoy callback chain starts here

def process_decoys(req, res):
    """Callback that retreives decoys and queues ascii fingerprints."""
    req.g.indexes.update(res[0])
    for (luid, lsmi), dlist in res[1].iteritems():
        smis = [[dpid, smi.replace('\\\\', '\\')] for (smi, dzid, dpid) 
                                                       in dlist]
        ltup = (lsmi, req.lzid, luid)
        if not smis:
            req.g.decoys[ltup] = []
        else:
            nreq = threadpool.Request(get_fingerprints, (smis, ), 
                       callback=load_decoys, sendself=True, g=req.g)
            nreq.ltup = ltup
            nreq.dlist = dlist
            req.g.fpp.put(nreq)
        # write queried decoys
        if req.g.loglevel > 3:
            mid = "%s%08d" % (req.g.luid_label, luid)
            fn = os.path.join(req.g.sdir, 'decoys.' + mid + '.queried')
            req.g.filep.put(threadpool.Request(dudtools.write_decoys, 
                                           (fn, ltup, dlist)))

def load_decoys(req, res):
    """Callback that retreives ascii fps and queues raw fps and bitcounts."""
    asciifps = [fps for (dpid, fps) in res]
    nreq = threadpool.Request(tanimoto.loadascii, (asciifps, ),
               callback=tanimoto_filter, g=req.g)
    nreq.ltup = req.ltup
    nreq.dlist = req.dlist
    req.g.cpup.put(nreq)


def tanimoto_filter(req, res):
    """Callback that retreives bitcounts and queues tanimoto filter."""
    dfps, dbcs = res
    if req.g.cutoff is not None:
        nreq = threadpool.Request(tanimoto.simfilterB,
                   (req.g.cutoff, req.g.lfps, req.g.lbcs, dfps, dbcs),
                   callback=finish_decoys, g=req.g)
    else:
        nreq = threadpool.Request(tanimoto.fractionfilterB, 
                   (req.g.fraction, req.g.lfps, req.g.lbcs, dfps, dbcs), 
                   callback=finish_decoys, g=req.g)
    nreq.ltup = req.ltup
    nreq.dlist = req.dlist
    req.g.cpup.put(nreq)

def finish_decoys(req, res):
    """Callback that retreives tanimoto cutoffs and writes output file."""
    smiles = [req.dlist[i] for i in xrange(len(req.dlist)) if i not in res]
    # decoy format is {(lsmi, lzid, luid) -> [(smis, dzid, dpid)]}
    req.g.decoys[req.ltup] = smiles

# Decoy callback chain ends here

def read_ligands(initer, column=1):
    "Read ligand ZINC ids from input iterator."""
    lzids = dudtools.read_zids(initer, col=column, raw=True)
    lzids = dudtools.unique(lzids)
    return lzids

def query_ligands(main_thread, lzids):
    """Query ligands in preperation for decoy generation."""
    uniqs = {}
    for lzid in lzids:
        try:
            nlzid, uniq = get_ligand_protomers(main_thread, lzid)
            uniqs[nlzid] = uniq
        except dudtools.DudError, err:
            print "Error:", err
    return uniqs

def process_ligands(mythread, filename, uniqs, loglevel=2, luid_label="P"):
    write_uniqs(filename, uniqs, loglevel=loglevel,
                        luid_label=luid_label)
    smis = [v for uniq in uniqs.itervalues()
                for v in uniq.itervalues()]
    lfps, lbcs = load_smiles(mythread, smis)
    return lfps, lbcs

def write_uniqs(filename, uniqs, loglevel=2, luid_label='P'):
    """Write unique ligand protomers."""
    if loglevel > 1:
        props = [[v[1], "C%08d" % lzid, "%s%08d" % (luid_label, v[0])] 
                    + list(k) for lzid, uniq in uniqs.iteritems() 
                                for k, v in uniq.iteritems()]
        mmmutils.write_splits(filename, props)

def read_uniqs(infile):
    """Read unique ligand protomers from prior ligands.charge file."""
    splits = mmmutils.read_splits(infile, raw_file=True)
    uniqs = {}
    for spl in splits:
        lzid = int(spl[1][1:])
        # hack to get around the lack of ternary operator in python 2.4
        key = tuple(('.' in x and [float(x)] or [long(x)])[0] for x in spl[3:])
        value = (long(spl[2][1:]), spl[0])
        uniqs.setdefault(lzid, {})[key] = value
    return uniqs

def load_smiles(mythread, smiles):
    """Get fingerprints and bitcounts for a list of smiles."""
    fps = get_fingerprints(mythread, smiles)
    asciifps = [x[1] for x in fps]
    lfps, lbcs = tanimoto.loadascii(asciifps)
    return lfps, lbcs

# Stub for more flexible query_decoys function
#   You tell it which queues to watch so that additional filters can easily
#   be added into the callback chain
#def query_decoys(filep, queue_pools, sdir, uniqs, ligand_sea_ids,
#                 loglevel=2, luid_label="P", callback=process_decoys):
#
def query_decoys(mysqlp, fpp, cpup, filep, sdir, uniqs, lfps, lbcs,
                 loglevel=2, cutoff=None, fraction=FRACTION_TO_CUT, 
                 luid_label="P"):
    """Mine potential decoys from ZINC database."""
    g = dudtools.EmptyContainer()
    g.sdir = sdir
    g.loglevel = loglevel
    g.cutoff = cutoff
    g.fraction = fraction
    g.luid_label = luid_label
#    for pool in queue_pools:
#        setattr(g, 
    g.fpp = fpp
    g.cpup = cpup
    g.filep = filep
    g.lfps = lfps
    g.lbcs = lbcs
    g.indexes = {}
    g.decoys = {}
    # Start decoy callback chain
    for lzid in uniqs.keys():
        uniq = uniqs.pop(lzid)
        req = threadpool.Request(get_decoys, (uniq,), 
                  sendself=True, callback=process_decoys, g=g)
        req.lzid = lzid
        mysqlp.put(req)
        # Poll queues to keep them going
        if len(uniqs) % 5 == 0:
            fpp.poll()
            mysqlp.poll()
            filep.poll()
            cpup.poll()
    # Poll queues until decoys are done
    count = 1
    while mysqlp.qsize or fpp.qsize or cpup.qsize or filep.qsize:
        fpp.poll()
        mysqlp.poll()
        filep.poll()
        cpup.poll()
        time.sleep(0.5)
        count += 1
    if loglevel > 5:
        fn = os.path.join(sdir, 'query.index')
        splits = (("%s%08d" % (luid_label, k[0]), "%1d" % v) 
                   for k, v in g.indexes.iteritems())
        filep.put(threadpool.Request(mmmutils.write_splits, (fn, splits)))
    if loglevel > 5:
        for ltup, dlist in g.decoys.iteritems():
            mid = "%s%08d" % (luid_label, ltup[2])
            fn = os.path.join(sdir, 'decoys.' + mid + '.tanimoto')
            filep.put(threadpool.Request(dudtools.write_decoys, 
                                         (fn, ltup, dlist)))
    cleaned = clean_decoys(g.decoys)
    if loglevel > 1:
        for ltup, dlist in cleaned.iteritems():
            mid = "%s%08d" % (luid_label, ltup[2])
            fn = os.path.join(sdir, 'decoys.' + mid + '.filtered')
            filep.put(threadpool.Request(dudtools.write_decoys, 
                                         (fn, ltup, dlist)))
    return cleaned

def setup_threads(fp_server):
    """Setup threads and remote proxies."""
    main_thread = dudtools.EmptyContainer()
    main_thread.num = 17
    dudtools.init_MySQL_Thread(main_thread)
    dudtools.init_FP_Thread(main_thread, fp_server)
    sqlinit = threadpool.Request(dudtools.init_MySQL_Thread, sendself=True)
    sqlcleanup = threadpool.Request(dudtools.cleanup_MySQL_Thread, 
                                    sendself=True)
    mysqlp = threadpool.ThreadPool(NUM_SQL_THREADS,
                                       init=sqlinit, cleanup=sqlcleanup)
    fpinit = threadpool.Request(dudtools.init_FP_Thread, args=(fp_server,), 
                                sendself=True)
    fpp = threadpool.ThreadPool(NUM_FP_THREADS, init=fpinit)
    cpup = threadpool.ThreadPool(NUM_CPU_THREADS)
    filep = threadpool.ThreadPool(1)
    time.sleep(2)
    return main_thread, mysqlp, fpp, cpup, filep

def setup_dirs(outdir, loglevel=2):
    outdir = os.path.abspath(outdir)
    if loglevel > 1:
        sdir = os.path.join(outdir, SEARCH_DIR)
    else:
        sdir = outdir
    for x in (outdir, sdir):
        if not os.path.exists(x):
            os.mkdir(x)
    return outdir, sdir

def create_potential_decoys(infile=None, outdir='.', loglevel=2, 
                            fraction=FRACTION_TO_CUT, column=1,
                            charge_start=False, cutoff=None, 
                            fp_server=dudtools.FP_SERVER_OPTIONS["ecfp4"]):
    """Generate lists of potential decoys by querying ZINC."""
    if infile is None:
        inf = sys.stdin
    else:
        inf = open(infile, 'r')
    outdir, sdir = setup_dirs(outdir, loglevel)
    main_thread, mysqlp, fpp, cpup, filep = setup_threads(fp_server)

    try:
        fn = os.path.join(outdir, CHARGE_FILE_LIG + CHARGE_FILE_EXT)
        save_loglevel = loglevel
        if charge_start:
            print "Reading Ligand Protomers."
            uniqs = read_uniqs(inf)
            if os.path.abspath(fn) == os.path.abspath(infile):
                loglevel = 0
        else:
            print "Querying Ligand Protomers."
            lzids = read_ligands(inf, column=column)
            uniqs = query_ligands(main_thread, lzids)
        lfps, lbcs = process_ligands(main_thread, fn, uniqs, loglevel=loglevel)
        loglevel = save_loglevel
    finally:
        inf.close()
        dudtools.cleanup_MySQL_Thread(main_thread)

    try:
        print "Querying Decoys."
        decoys = query_decoys(mysqlp, fpp, cpup, filep, sdir, 
                              uniqs, lfps, lbcs, loglevel=loglevel, 
                              cutoff=cutoff, fraction=fraction)
    finally:
        print "Shutting Down."
        mysqlp.shutdown()
        fpp.shutdown()
        cpup.shutdown()
        filep.shutdown()
    return decoys

def option_parser(usage, description, version):
    parser = OptionParser(usage=usage, description=description,
                          version=version)
    parser.set_defaults(infile=None, outdir='.', log_level=2, cutoff=None, 
                        fraction=FRACTION_TO_CUT, fp_server="ecfp4",
                        charge_start=False)
    parser.add_option("-i", "--infile",
           help="input file (default: stdin)")  
    parser.add_option("-o", "--outdir",
           help="output directory (default: %default)")
    parser.add_option("-l", "--log-level", type="int", 
           help="log level verbosity (from 1-9) (default: %default)")
    parser.add_option("-d", "--decoy-fraction", type="float", dest="fraction", 
           help="decoy fraction to remove by tanimoto similarity" +
                     " (default: %default)")
    parser.add_option("-c", "--decoy-cutoff", type="float", dest="cutoff", 
           help="tanimoto cutoff, disables decoy fraction" +
                     " (default is to use decoy fraction instead)")
    parser.add_option("-r", "--restart-ligands.charge",
           dest="charge_start", action="store_true", 
           help="input file is ligands.charge from a prior run")
    fp_types = ', '.join(dudtools.FP_SERVER_OPTIONS.keys())
    parser.add_option("-f", "--fingerprint-server", dest="fp_server", 
           help="fingerprint server: options = %s" % fp_types +
                     " (default: %default)")
    return parser

def parse(parser, argv):
    options, args = parser.parse_args(args=argv[1:])
    if len(args):
        parser.error("program takes no positional arguments.\n" +
                     "  Use --help for more information.")
    options.log_level=max(1, options.log_level)
    options.log_level=min(9, options.log_level)
    if options.fp_server not in dudtools.FP_SERVER_OPTIONS:
        parser.error("invalid fingerprint server.\n" +
            "  Use --help to get a current list of fingerprint servers.")
    options.fp_server = dudtools.FP_SERVER_OPTIONS[options.fp_server]
    return options

def main(argv):
    """Parse arguments."""
    description = "Generate lists of potential decoys by querying ZINC."
    usage = "%prog [options]"
    version = "%prog: version 201109 - created by Michael Mysinger"
    parser = option_parser(usage, description, version)
    options = parse(parser, argv)

    start_time = time.time()
    create_potential_decoys(infile=options.infile, outdir=options.outdir, 
        loglevel=options.log_level, cutoff=options.cutoff, 
        fraction=options.fraction, fp_server=options.fp_server,
        charge_start=options.charge_start)
    print 'Program took %.1f seconds.' % (time.time() - start_time)
    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv))
