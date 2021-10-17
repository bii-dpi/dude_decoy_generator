#!/bin/env python
"""Translate temp ids back to user ids.

Michael Mysinger 201205 Created
"""

import os
import sys
import logging
import os.path as op
from optparse import OptionParser

DICT_FILE = "input.ism"
DICT_IN_ID_COL = 2
DICT_OUT_ID_COL = 3
LIG_ID_COL = 2
DEC_ID_COL = 3

class ScriptError(Exception):
    def __init__(self, msg, value=99):
        self.value = value
        Exception.__init__(self, msg)

def read_translation(dict_f, icol=DICT_IN_ID_COL, ocol=DICT_OUT_ID_COL):
    trans_dict = {}
    for line in dict_f:
        splits = line.split()
        in_id = "C%08d" % int(splits[icol-1])
        out_id = splits[ocol-1]
        trans_dict[in_id] = out_id 
    return trans_dict

def translate_ids(indir, trans_dict, ligcol=LIG_ID_COL, deccol=DEC_ID_COL):
    """Translate ids in ligands.charge then the decoys*picked files."""
    lig_fn = os.path.join(indir, "ligands.charge")
    temp_fn = os.path.join(indir, "trans_temp.charge")
    os.rename(lig_fn, temp_fn)
    in_f = open(temp_fn)
    out_f = open(lig_fn, "w")
    for line in in_f:
        in_id = line.split()[ligcol-1]
        out_f.write(line.replace(in_id, trans_dict[in_id]))
    in_f.close()
    out_f.close()
    os.remove(temp_fn)
    dec_dir = os.path.join(indir, "decoys")
    for filename in os.listdir(dec_dir):
        dec_fn = os.path.join(dec_dir, filename)
        temp_fn = os.path.join(dec_dir, "trans_temp.picked")
        os.rename(dec_fn, temp_fn)
        in_f = open(temp_fn)
        out_f = open(dec_fn, "w")
        line = in_f.next()
        in_id = line.split()[deccol-1]
        out_f.write(line.replace(in_id, trans_dict[in_id]))
        for line in in_f:
            out_f.write(line)
        in_f.close()
        out_f.close()
        os.remove(temp_fn)

def handleio(indir='.', dictfile=DICT_FILE, **kwargs):
    """I/O handling for the script."""
    dict_f = open(dictfile)
    trans_dict = read_translation(dict_f)
    dict_f.close()
    try:
        translate_ids(indir, trans_dict)
    except ScriptError, message:
        logging.error(message)
        return message.value
    return 0

def main(argv):
    """Parse arguments."""
    logging.basicConfig(level=logging.INFO,
                        format="%(levelname)s: %(message)s")
    description = "Translate temp ids back to user ids. Translates both ligands.charge and all decoys.*.picked ligand lines."
    usage = "%prog [options]"
    version = "%prog: version 201205 - created by Michael Mysinger"
    parser = OptionParser(usage=usage, description=description,
                          version=version)
    parser.set_defaults(indir='.', dictfile=DICT_FILE)
    parser.add_option("-i", "--indir", 
        help="input directory (default: %default)")
    parser.add_option("-d", "--dictfile", 
                      help="dictionary file (default: %default)")
    options, args = parser.parse_args(args=argv[1:])
    if len(args):
        parser.error("program takes no positional arguments.\n" +
                     "  Use --help for more information.")
    return handleio(indir=options.indir, dictfile=options.dictfile)

if __name__ == "__main__":
    sys.exit(main(sys.argv))
