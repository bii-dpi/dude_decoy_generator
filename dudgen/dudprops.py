#!/usr/arch/bin/python -u
"""Generate property graphs and html pages for every dud target.

Michael Mysinger 200804 Rewritten as wrapper around dudprop.py
"""

import os
import sys
import shutil
import time
from optparse import OptionParser
import duddb
import dudhtml
import dudprop
import dudtools
import mmmutils

def do_html(outdir, targets):
    """Generate html pages with thumbnails and plots."""
    shutil.copy(dudhtml.STYLESHEET, outdir)
    stylesheet = os.path.basename(dudhtml.STYLESHEET)
    props = ['mw', 'logp', 'hba', 'hbd', 'rotb', 'net']
    menuitems = [[x, x] for x in targets]
    for target in targets:
        html = open(os.path.join(outdir, target+'.html'), 'w')
        title = target.upper()
        images = [os.path.join(target, x) for x in props]
        thtml = dudhtml.ThumbHTML(html, stylesheet, title, menuitems, images)
        thtml.page()
        html.close()
    index = os.path.join(outdir, 'index.html')
    if os.path.exists(index):
        os.remove(index)
    os.symlink(os.path.join(outdir, targets[0]+'.html'), index)

def gen_properties(indir='.', outdir='.', legacy=False, loglevel=2):
    """Calculate properties and display ligand versus decoy plots in html."""

    outdir, outgen, genlig, gendec = duddb.setup_dirs(outdir)
    targets = dudtools.read_dirs(genlig)
    
    pdir = os.path.join(outdir, dudprop.PROP_DIR)
    if not os.path.exists(pdir):
        os.mkdir(pdir)
    for i, target in enumerate(targets):
        if legacy:
            if target == 'na':
                target = 'neua'
                targets[i] = target
            elif target == 'ppar_gamma':
                target = 'ppar'
                targets[i] = target
            elif target == 'rxr_alpha':
                target = 'rxr'
                targets[i] = target
        print 'Processing output directory %s.' % target
        lproplist = dudprop.collect_dir(os.path.join(genlig, target))
        dproplist = dudprop.collect_dir(os.path.join(gendec, target))
        if loglevel > 1:
            fn = dudprop.CHARGE_FILE_DEC + "." + target + dudprop.CHARGE_FILE_EXT
            fn = os.path.join(outdir, fn)
            mmmutils.write_splits(fn, dproplist)
            fn = "ligands." + target + dudprop.CHARGE_FILE_EXT
            fn = os.path.join(outdir, fn)
            if not os.path.exists(fn):
                mmmutils.write_splits(fn, lproplist)
        tdir = os.path.join(pdir, target)
        if not os.path.exists(tdir):
            os.mkdir(tdir)
        dudprop.plot_props(tdir, lproplist, dproplist)
    do_html(pdir, targets)

def main(argv):
    """Parse arguments."""
    description = "Generate property graphs and html for every dud target."
    usage = "%prog [options]"
    version = "%prog: version 200804 - created by Michael Mysinger"
    parser = OptionParser(usage=usage, description=description,
                          version=version)
    parser.set_defaults(indir='.', outdir='.')
    parser.add_option("-i", "--indir",
           help="input directory (default: %default)")  
    parser.add_option("-o", "--outdir",
           help="output directory (default: %default)")
    options, args = parser.parse_args(args=argv[1:])
    if len(args):
        parser.error("program takes no positional arguments.\n" +
                     "  Use --help for more information.")
    start_time = time.time()
    gen_properties(indir=options.indir, outdir=options.outdir)
    print 'Program took %.1f seconds.' % (time.time() - start_time)
    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv))
