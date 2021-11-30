"""dudtools.py holds common code for dud generation and analysis.

Michael Mysinger 200703 Created
"""

import os
import xmlrpclib
import MySQLdb

# Directory constants
DUDDIR = '/raid0/people/mysinger/dud3'
ZINCDIR='/raid2/db/zinc0510'
CODEDIR = os.path.dirname(os.path.abspath(__file__))

# MySQL server
MYSQL_HOST = 'zincdb6'
MYSQL_USER = 'lab'
MYSQL_DATABASE = 'zinc8'

# Fingerprint server configuration
FP_SCRIPT = os.path.join(CODEDIR, 'fp_server', 'fp_server.bash')

# Daylight fingerprint server (MMM)
DAY_HOST, DAY_PORT, DAY_FUNCTION = 'korn', 8121, "get_fingerprints"

# Daylight fingerprint server (MJK)
DAYK_HOST, DAYK_PORT, DAYK_FUNCTION = 'korn', 8181, "smi2daylight"

# OEchem fingerprint server
OE_HOST, OE_PORT, OE_FUNCTION = 'fingerprint', 8181, "smi2oechem"

# ChemAxon path-based fingerprint server
CAP_HOST, CAP_PORT, CAP_FUNCTION = 'fingerprint', 8182, "smi2axonpath"

# ChemAxon ECFP fingerprint server
CAE_HOST, CAE_PORT, CAE_FUNCTION = 'fingerprint', 8183, "smi2axonecfp"

# Pipeline pilot ECFP4 fingerprint server  
PP_HOST, PP_PORT, PP_FUNCTION = 'ppilot', 8181, "smi2ecfp"

# Fingerprint server dict
FP_SERVER_OPTIONS = {'daylight': (DAYK_HOST, DAYK_PORT, DAYK_FUNCTION),
                     'ecfp4':  (PP_HOST, PP_PORT, PP_FUNCTION),
                     'oechem':  (OE_HOST, OE_PORT, OE_FUNCTION),
                     'axonpath':  (CAP_HOST, CAP_PORT, CAP_FUNCTION),
                     'axonecfp':  (CAE_HOST, CAE_PORT, CAE_FUNCTION)}

# Disk server
DISK_HOST = 'nfshead3'  # disk server run local to raid5
DISK_PORT = 8141
DISK_SCRIPT = os.path.join(CODEDIR, 'disk_server.py')

class DudError(Exception):
    pass    

class EmptyContainer(object):
    pass

def unique(inlist):
    """Generate list of unique elements from inlist."""
    return list(set(inlist))

def notnone(inlist):
    return [x for x in inlist if x is not None]

def get_id_position(idstr):
    if idstr[:4].isalpha():
        start, end = 4, 12
    elif idstr[:1].isalpha():
        start, end = 1, 9
    else:
        start, end = 0, 8
    return start, end

def read_dirs(indir):
    """Read, sort, and return subdirectories of indir.""" 
    outdirs = [x for x in os.listdir(indir) if 
                  os.path.isdir(os.path.join(indir, x))]
    outdirs.sort()
    return outdirs

def read_zids(filename, col=1, raw=False):
    """Read ZINC ids from file."""
    if raw:
        f = filename
    else:
        f = open(filename, 'r')
    first = f.next().split()[col]
    start, end = get_id_position(first)
    yield int(first[start:end])    
    for x in f:
        yield int(x.split()[col][start:end])
    f.close()

def read_decoys(filename, raw=False):
    """Read decoy smiles file."""
    if raw:
        f = filename
    else:
        f = open(filename, 'r')
    # read ligand
    first = f.next().split()
    zstart, zend = get_id_position(first[2])
    pstart, pend = get_id_position(first[3])
    yield first[1], int(first[2][zstart:zend]), int(first[3][pstart:pend])
    # read decoy
    second = f.next().split()
    zstart, zend = get_id_position(second[1])
    pstart, pend = get_id_position(second[2])
    yield second[0], int(second[1][zstart:zend]), int(second[2][pstart:pend])
    # use first decoy format for the rest
    for x in f:
        spl = x.split()
        yield spl[0], int(spl[1][zstart:zend]), int(spl[2][pstart:pend])
    f.close()

def read_pids(filename):
    """Read protomer ids from second column of smiles type files.""" 
    f = open(filename, 'r')
    splits = (x.split() for x in f)
    pids = [s[1][1:].zfill(8) for s in splits if len(s) > 1]
    f.close()
    return pids

def read_smis(filename):
    """Read smiles file and return smiles and ids."""
    f = open(filename, 'r')
    splits = (x.split() for x in f)
    smis, ids = zip(*splits)
    f.close()
    return smis, ids

def write_zids(filename, zids):
    f = open(filename, 'w')
    [f.write('ZINC'+str(x).zfill(8)+'\n') for x in zids]
    f.close()

def write_dict(filename, outdict):
    """Write out dictionary of lists with one key and all associated values per line.""" 
    f = open(filename, 'w')
    for k, v in outdict.iteritems():
        f.write(str(k)+' ')
        [f.write(str(x)+' ') for x in v]
        f.write('\n')
    f.close()

def write_decoys(filename, ligand, decoys, zid_label="C", 
                          luid_label="P", dpid_label="P"):
    """Write decoy smiles file."""
    f = open(filename, 'w')
    # write ligand
    f.write("ligand\t%s\t%s%08d\t%s%08d\n" % (ligand[0], zid_label, 
                int(ligand[1]), luid_label, int(ligand[2])))
    # write decoys
    for s in decoys:
        f.write("%s\t%s%08d\t%s%08d\n" % (s[0], zid_label, int(s[1]), 
                                          dpid_label, int(s[2])))
    f.close()

# Threadpool functions

def init_FP_Thread(self, fp_server):
    """Initialize fingerprint xmlrpc connections."""
    self.client = xmlrpclib.ServerProxy('http://%s:%d' % (fp_server[0], 
                            fp_server[1]), allow_none=True)
    self.config = fp_server

def init_MySQL_Thread(self):
    """Initialize MySQL database connections."""
    self.db = MySQLdb.connect(host=MYSQL_HOST, user=MYSQL_USER, 
                              db=MYSQL_DATABASE)
    self.c = self.db.cursor()

def cleanup_MySQL_Thread(self):
    """Cleanup MySQL database connections."""
    self.c.close()
    self.db.close()

# Database functions below require MySQLdb cursor c

def sequential_query(c, insql, inlist):
    """Return one database row or None for each item in inlist."""
    query = []
    for initem in inlist:
        inexpr = "(%s)" % initem
        if c.execute(insql % inexpr) > 0:
            query.append(c.fetchone())
        else:
            query.append(None)
    return query

def multiple_query(c, insql, inlist, chunk_size=10000):
    """Perform SQL query on entire chunks of a list at once.

    Warning: in general, there is no one-to-one correspondence between 
    the output items and input items, because multiple database records 
    could share the same input value or a particular input value could 
    have no records. To get a one-to-one correspondence, use a key 
    field (or distinct query) and check for the correct length.
    """
    numloops = len(inlist) // chunk_size + 1
    query = []
    for i in xrange(numloops):
        strlist = (str(x) for x in inlist[i*chunk_size:(i+1)*chunk_size])
        inexpr = "(%s)" % ', '.join(strlist)
        c.execute(insql % inexpr)
        query.extend(c.fetchall())
    return query

def one_to_one_query(c, insql, inlist, chunk_size=10000):
    """Perform SQL query for every element in a list.

    Warning: Use a key field in SQL for one-to-one correspondence.
    """
    # try fast query
    query = multiple_query(c, insql, inlist, chunk_size=chunk_size)
    # if elements are missing, use slow query
    if len(query) != len(inlist):
        query = sequential_query(c, insql, inlist)
    return query    

def translate_pids_to_zids(c, pids):
    """Translate pids into zids by querying database."""
    zid_sql = 'select sub_id_fk from protomer where prot_id in %s'
    # Rewrote this for Python 2.4
    #zids = [x[0] if x is not None else x for x in 
    #             one_to_one_query(c, zid_sql, pids)]
    zids = [x is not None and x[0] or x for x in 
                 one_to_one_query(c, zid_sql, pids)]
    return zids

def translate_retired_zids(c, zids):
    """Translate any retired zids into new ids using database."""
    retired_sql = 'select new_id from retired where old_id in %s'
    query = sequential_query(c, retired_sql, zids)
    # Rewrote this for Python 2.4
    #nzids = [z if q is None else q[0] for z, q in zip(zids, query)]
    nzids = [q is None and z or q[0] for z, q in zip(zids, query)]
    return nzids
