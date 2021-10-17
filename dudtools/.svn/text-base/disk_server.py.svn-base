#!/bin/env python
"""XML-RPC server to collect database generation results.

Often started on a remote disk server to allow local disk access.

Michael Mysinger 200709 Created
Michael Mysinger 201109 Collection entirely in python and non-sequential 
"""

import os
import sys
import bz2
import gzip
import subprocess
import time
import os.path as osp
from SimpleXMLRPCServer import SimpleXMLRPCServer
from SocketServer import ForkingMixIn
import mmmutils

# Constants
DBDELIM = 'Family'
DBHEADER = osp.join(osp.dirname(osp.abspath(__file__)), 'header.db')
QUIT = False
PROT_DIRS = ["ref", "mid", "hi", "lo"]
STARTED_FILE = "started"
FINISHED_FILE = "finished"
DICT_FILE = "dict"
TEMP_DIR_BZIP = "db.db.bz2"

def read_dict(initer):
    """Read and format DOCK database generation dict file."""
    splits = (x.strip().split(None, 4) for x in initer)
    res = []
    for x in splits:
        temp_id = x[0].replace('TEMP', 'P')[:9]
        if len(x) == 5:
            orig_id = x[1].replace('ZINC', 'C')
            molname = x[4].replace(' ', '_')
        else:
            # Handle the case where there was no original id provided
            nx = '\t'.join(x).split(None, 2)
            orig_id = temp_id
            molname = x[2].replace(' ', '_')            
        res.append([temp_id, [orig_id, molname]])
    return dict(res)

def read_subdirs(ndir):
    subdirs = set()
    for prot in PROT_DIRS:
        pdir = osp.join(ndir, prot)
        if osp.exists(pdir):
            subdirs.update(mmmutils.read_dirlist(pdir, skip_isdir=True))
    return subdirs

def translate_family(family, transdict):
    #translate name and id on first line after delimeter
    lines = family.split("\n", 2)
    if len(lines) < 3:
        return None
    oid = lines[1][47:56]
    nid, nname = transdict.get(oid, (oid, "Translation_error!"))
    nline = "%s  %s%s" % (nname[:45].ljust(45), nid[:9].ljust(9), 
                              lines[1][56:])
    return "\n".join([lines[0], nline, lines[2]])    

def translate_temp_ids(splits, transdict):
    nsplits = []
    for family in splits[:-1]:
        if family:
            nfamily = translate_family(family, transdict)
            if nfamily:
                nsplits.append(nfamily)
    nsplits.append(splits[-1])
    return nsplits

def start_db(dbnum, outdb, header):
    """Open new database and write header.""" 
    filename = "%s_%s.db.gz" % (outdb, str(dbnum).zfill(4))
    odb = gzip.GzipFile(filename, 'w')
    for x in header:
        odb.write(x)
    return odb    

def append_db(dbnum, dbcount, odb, indb, outdb, header, dbsplitsize, 
              transdict=None, chunk_size=1000000):
    try:
        idb = mmmutils.flex_open(indb)
    except mmmutils.MMMError:
        print "Skip: could not find file %s" % indb
        idb = None
    if idb:
        past_header = False
        splits = [""]
        try:
            chunk = idb.read(chunk_size)
        except EOFError:
            print "Truncated file %s, trying one more time." % indb
            time.sleep(10)
            idb.seek(0)
            chunk = idb.read(chunk_size)
        while chunk:
            buffer = splits[-1] + chunk
            splits = buffer.split(DBDELIM)
            num = len(splits)
            if not past_header and num > 1:
                splits = splits[1:]
                num = len(splits)
                past_header = True
            if transdict:
                splits = translate_temp_ids(splits, transdict)
                num = len(splits)
            if num > 1:
                edge = dbcount + num - dbsplitsize
                if edge > 0:
                    odb.write(DBDELIM.join([""]+splits[:-edge]))
                    odb.close()
                    dbnum += 1
                    odb = start_db(dbnum, outdb, header)
                    odb.write(DBDELIM.join([""]+splits[-edge:-1]))
                    dbcount = edge - 1
                else:
                    odb.write(DBDELIM.join([""] + splits[:-1]))
                    dbcount += num - 1
            chunk = idb.read(chunk_size)
        last = splits[-1]
        if transdict:
            last = translate_family(last, transdict)
        if last:
            if dbcount < dbsplitsize:
                odb.write(DBDELIM.join(["", last]))
                dbcount += 1
            else:
                odb.close()
                dbnum += 1
                odb = start_db(dbnum, outdb, header)
                odb.write(DBDELIM.join(["", last]))
                dbcount = 1
    return dbnum, dbcount, odb

def collect_results(dbname, gendir, dbsplitsize):
    """Collect generation results and accumulate into final database files."""
    print "Starting collection on %s." % gendir
    f = open(DBHEADER, 'r')
    header = f.readlines()
    f.close()
    ndirs = (osp.join(gendir, x) for x in os.listdir(gendir))
    # dirdict: ndir -> (started flag, translation dict, dict of temp subdirs)
    dirdict = dict((ndir, [False, None, None]) for ndir in ndirs
                                               if osp.isdir(ndir))
    dbnum = 1
    dbcount = 0
    outdb = dbname
    odb = start_db(dbnum, outdb, header)
    try:
        while dirdict:
            for ndir in dirdict.keys():
                if not dirdict[ndir][0]:
                    if osp.exists(osp.join(ndir, STARTED_FILE)):
                        dirdict[ndir][0] = True
                    else:
                        continue
                if dirdict[ndir][1] is None:
                    trans_dict_f = open(osp.join(ndir, DICT_FILE))
                    dirdict[ndir][1] = read_dict(trans_dict_f)
                    trans_dict_f.close()
                if dirdict[ndir][2] is None:
                    subdirs = read_subdirs(ndir)
                    dirdict[ndir][2] = subdirs

                for tdir in dirdict[ndir][2].copy():
                    db_fn = osp.join(tdir, TEMP_DIR_BZIP)
                    if osp.exists(db_fn):
                        dbnum, dbcount, odb = append_db(dbnum, 
                            dbcount, odb, db_fn, outdb, header, 
                            dbsplitsize, dirdict[ndir][1])
                        dirdict[ndir][2].remove(tdir)
                if not dirdict[ndir][2]:
                    open(osp.join(ndir, FINISHED_FILE), 'w').close()
                    dirdict.pop(ndir)
            print ".",
            time.sleep(60)
    finally:
        odb.close()
    return True

class ForkingXMLRPCServer(ForkingMixIn, SimpleXMLRPCServer):
    allow_reuse_address = True

# Main
if __name__ == '__main__':

    if len(sys.argv) == 3:
        host = sys.argv[1]
        port = int(sys.argv[2])
    else:
        host = raw_input('host: ')
        port = int(raw_input('port number: '))

    print 'Starting file server on %s:%d.' % (host, port)
    server = ForkingXMLRPCServer((host, port))
    server.register_function(collect_results)

    try:
        server.serve_forever()
    finally:
        server.server_close()
