#!/bin/env python
"""Select subset of ligands and matched decoys from previous dudgen results.

Michael Mysinger 201110 Created
"""

import os
import sys
import logging
import glob
import gzip
import shutil
import os.path as op
from optparse import OptionParser

SEARCH_DIR = "search"
LIGAND_FILE = "ligands.charge"
DB_FILES = "*.db.gz"
DECOY_FILES = "decoys.*.picked"
LIGAND_DATABASES = "*_lig_*.db.gz"
DECOY_DATABASES = "*_dec_*.db.gz"
DBDELIM = 'Family'
CHUNK_SIZE = 1000000
DATABASE_SIZE = 3000

DATABASE_HEADER="""DOCK 5.2 ligand_atoms
positive                       (1)
negative                       (2)
acceptor                       (3)
donor                          (4)
ester_o                        (5)
amide_o                        (6)
neutral                        (7)
"""

class DUDError(Exception):
    def __init__(self, msg, value=99):
        self.value = value
        Exception.__init__(self, msg)

def get_paths(indir, outdir, force=False):
    """Check input and output paths before starting."""
    if op.abspath(indir) == op.abspath(outdir):
        raise DUDError("Output directory must be different than input " +
                      "to prevent data loss.", 6)

    search_dir = op.join(indir, SEARCH_DIR)
    if not op.exists(search_dir):
        raise DUDError("Unable to find input mapping directory %s!" % 
                       search_dir, 7)
    ligand_f = op.join(indir, LIGAND_FILE)
    if not op.exists(ligand_f):
        raise DUDError("Unable to find input ligand protomer file %s!" % 
                       ligand_f, 7)

    out_ligand_f = op.join(outdir, LIGAND_FILE)
    out_search_dir = op.join(outdir, SEARCH_DIR)
    out_db_files = op.join(outdir, DB_FILES)
    if force:
        if op.exists(out_ligand_f):
            logging.warn("Forced removal of %s" % out_ligand_f)
            os.remove(out_ligand_f)
        if op.exists(out_search_dir):
            logging.warn("Forced removal of %s" % out_search_dir)
            if op.isdir(out_search_dir):
                shutil.rmtree(out_search_dir)
            else:
                os.remove(out_search_dir)
        for x in glob.glob(out_db_files):
            logging.warn("Forced removal of %s" % x)
            os.remove(x)

    if op.exists(out_search_dir):
        raise DUDError("Output dir contains mapping directory %s!\n" % 
                       out_search_dir + "  Use -f to force overwrite", 8)
    if op.exists(out_ligand_f):
        raise DUDError("Output dir contains ligand protomer file %s!\n" % 
                      out_ligand_f + "  Use -f to force overwrite", 8)
    if glob.glob(out_db_files):
        raise DUDError("Output dir contains docking databases!\n" +  
                      "  Use -f to force overwrite", 8)
    os.makedirs(out_search_dir)
    return ligand_f, search_dir, out_ligand_f, out_search_dir

def read_subset(subset_file):
    """Read ligand ids to use for subset selection from input file."""
    # subset_ids: (ligand_file_id, db_file_id)
    subset_ids = set()
    for line in subset_file:
        spl = line.split()
        if len(spl) > 1:
            raise DUDError("Unknown subset id format on %s!\n" % line +
                "  Ligand ids must be either pure numbers or ZINC ids\n" +
                "  containing 1 uppercase letter followed by 8 digits.", 5)
        if spl:
            iden = spl[0]
            if iden.isdigit():
                id_tuple = ("C%08d" % int(iden), iden)
            elif ( len(iden) != 9 or not iden[1:].isdigit() or 
                   not iden[0].isalpha() or not iden[0].isupper() ):
                raise DUDError("'%s' identifier malformed!\n" % iden +
                    "  Ligand ids must be either pure numbers or ZINC ids\n" +
                    "  containing 1 uppercase letter followed by 8 digits.", 5)
            else:
                id_tuple = (iden, "%d" % int(iden[1:]))
            if id_tuple not in subset_ids:
                subset_ids.add(id_tuple)
            else:
                logging.warn("Skipping repeated id %s." % iden)
    ligand_ids, db_ids = zip(*subset_ids)
    return ligand_ids, db_ids

def read_picked(decoy_filename):
    """Get decoy ZINC ids from picked decoy file."""
    decoys = []
    dec_f = open(decoy_filename)
    first_spl = dec_f.next().split()
    if not len(first_spl) == 4 or not first_spl[0] == "ligand":
        raise DUDError("Unknown format for decoy mapping file %s" % 
                          dec_f, 16) 
    for line in dec_f:
        spl = line.split()
        decoys.append(spl[1])
    dec_f.close()
    return decoys

def test_id(family, ids):
    """Return database family only if it matches selected ids."""
    lines = family.split("\n", 2)
    if len(lines) < 3:
        return None
    myid = lines[1][47:56].strip() 
    if myid in ids:
        return family
    return None

def filter_ids(splits, ids):
    """Subselect from a list database families using ids."""
    nsplits = []
    for family in splits[:-1]:
        if family:
            nfamily = test_id(family, ids)
            if nfamily:
                nsplits.append(nfamily)
    nsplits.append(splits[-1])
    return nsplits

def start_db(dbnum, outdb, header=DATABASE_HEADER):
    """Open new database and write header.""" 
    filename = "%s_%s.db.gz" % (outdb, str(dbnum).zfill(4))
    odb = gzip.GzipFile(filename, 'w')
    for x in header:
        odb.write(x)
    return odb    

def append_db(dbnum, dbcount, odb, indb, outdb, ids,
                dbsplitsize=DATABASE_SIZE, chunk_size=CHUNK_SIZE):
    """Filter one input database into output databases, chunking as needed."""
    idb = gzip.GzipFile(indb)
    if idb:
        past_header = False
        splits = [""]
        chunk = idb.read(chunk_size)
        while chunk:
            buffer = splits[-1] + chunk
            splits = buffer.split(DBDELIM)
            num = len(splits)
            if not past_header and num > 1:
                splits = splits[1:]
                num = len(splits)
                past_header = True
            splits = filter_ids(splits, ids)
            num = len(splits)
            if num > 1:
                edge = dbcount + num - dbsplitsize
                if edge > 0:
                    odb.write(DBDELIM.join([""]+splits[:-edge]))
                    odb.close()
                    dbnum += 1
                    odb = start_db(dbnum, outdb)
                    odb.write(DBDELIM.join([""]+splits[-edge:-1]))
                    dbcount = edge - 1
                else:
                    odb.write(DBDELIM.join([""] + splits[:-1]))
                    dbcount += num - 1
            chunk = idb.read(chunk_size)
        last = test_id(splits[-1], ids)
        if last:
            if dbcount < dbsplitsize:
                odb.write(DBDELIM.join(["", last]))
                dbcount += 1
            else:
                odb.close()
                dbnum += 1
                odb = start_db(dbnum, outdb)
                odb.write(DBDELIM.join(["", last]))
                dbcount = 1
    idb.close()
    return dbnum, dbcount, odb

def process_ligands(ligand_f, out_ligand_f, ligand_ids):
    """Get ligand protomers of selected ligand ids."""
    protomers = {}
    ligf = open(ligand_f)
    out_ligf = open(out_ligand_f, "w")
    for line in ligf:
        spl = line.split("\t")
        if len(spl) != 9:
            raise DUDError("Unknown format for ligand protomer file %s" % 
                          ligand_f, 11)
        if spl[1] in ligand_ids:
            protomers[spl[2]] = spl[1]
            out_ligf.write(line)
    ligf.close()
    out_ligf.close()
    return protomers

def get_decoys(search_dir, out_search_dir, protomers):
    """Get mapping from ligand protomer to decoy ZINC ids."""
    all_decoys = {}
    for dec_fn in glob.glob(op.join(search_dir, DECOY_FILES)):
        prot = dec_fn.split('.')[1].replace("F", "P")
        if prot in protomers:
            all_decoys[prot] = read_picked(dec_fn)
            out_dec_fn = op.join(out_search_dir, op.basename(dec_fn))
            os.symlink(op.abspath(dec_fn), op.abspath(out_dec_fn))
    return all_decoys

def make_set(*args):
    """Make a single combined set from any number of input iterators."""
    out_set = set()
    [out_set.add(y) for x in args for y in x]
    return out_set

def subset_databases(files, ids, outdir):
    """Subselect only some ids from docking database files."""
    try:
        first = True
        for db_f in glob.glob(files):
            if first:
                dbnum = 1
                dbcount = 0
                outdb = op.join(outdir, op.basename(db_f.rsplit("_", 1)[0]))
                odb = start_db(dbnum, outdb)
                first = False
            dbnum, dbcount, odb = append_db(dbnum, dbcount, odb, db_f,
                                            outdb, ids)
    finally:
        odb.close()

def decoy_picker(indir, outdir, subset_file=None, force=False):
    """Get ready to pick decoys."""
    if subset_file is None:
        subf = sys.stdin
    else:
        subf = open(subset_file)
    try:
        ligand_ids, db_ids = read_subset(subf)
        file_paths = get_paths(indir, outdir, force=force)
        ligand_f, search_dir, out_ligand_f, out_search_dir = file_paths
        protomers = process_ligands(ligand_f, out_ligand_f, ligand_ids)
        all_decoys = get_decoys(search_dir, out_search_dir, protomers)
        ligand_set = make_set(ligand_ids, db_ids)
        subset_databases(op.join(indir, LIGAND_DATABASES), ligand_set, outdir)
        decoy_set = make_set(*all_decoys.values())
        subset_databases(op.join(indir, DECOY_DATABASES), decoy_set, outdir)
    except DUDError, message:
        logging.error(message)
        return message.value
    return 0

def main(argv):
    """Parse arguments."""
    logging.basicConfig(level=logging.INFO,
                        format="%(levelname)s: %(message)s")
    description = "Select subset of ligands and matched decoys from previous dudgen results."
    usage = "%prog [options]"
    version = "%prog: version 201110 - created by Michael Mysinger"
    parser = OptionParser(usage=usage, description=description,
                          version=version)
    parser.set_defaults(subset_file=None, indir=None, outdir=None, force=False)
    parser.add_option("-i", "--indir", 
        help="input directory containing full dud dataset (required)")
    parser.add_option("-o", "--outdir", 
        help="output directory for new dud subset (required)")
    parser.add_option("-s", "--subset-file", 
        help="file listing selected subset ids (default: stdin)")
    parser.add_option("-f", "--force", action="store_true", 
        help="force overwrite of output directory")
    options, args = parser.parse_args(args=argv[1:])
    if len(args):
        parser.error("program takes no positional arguments.\n" +
                     "  Use --help for more information.")
    if not options.indir:
        parser.error("You must specify input directory containing\n" +
                     "  full dud dataset.\n" +
                     "  Use --help for more information.")
    if not options.outdir:
        parser.error("You must specify output directory for\n" +
                     "  new dud subset.\n" +
                     "  Use --help for more information.")
    return decoy_picker(options.indir, options.outdir, force=options.force, 
                        subset_file=options.subset_file)

if __name__ == "__main__":
    sys.exit(main(sys.argv))
