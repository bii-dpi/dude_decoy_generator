#!/bin/env python
"""Generate property graphs and html.

Michael Mysinger 200702 Created
Michael Mysinger 200703 Use dudtools
Michael Mysinger 200704 Use dudhtml
Michael Mysinger 200704 Use processqueue and fnmatch
Michael Mysinger 200711 Adapted to use dudtest.py generated database
Michael Mysinger 200804 Rewrote for single directory using dudgen
"""

import os
import sys
import shutil
import time
from optparse import OptionParser
import pylab
import scipy
import duddb
import dudgen
import dudhtml
import mmmutils

# Constants
PROP_DIR = 'properties'
CHARGE_FILE_DEC = 'decoys'
CHARGE_FILE_EXT = '.charge'

def collect_dir(gendir, uid_label="P"):
    """Collect all properties after database generation."""
    uniqs = dudgen.collect_ligands(gendir)
    props = [[v[1], "C%08d" % lzid, "%s%08d" % (uid_label, v[0])] 
                     + list(k) for lzid, uniq in uniqs.iteritems() 
                                 for k, v in uniq.iteritems()]
    # proplist is [smi, zid, pid, mwt, logP, rb, hba, hbd, net]
    return props

def parse_props(proplist):
    """Parse properties from propery splits."""
    trans = zip(*proplist)
    # return [mwt], [logP], [rb], [hba], [hbd], [net]
    return tuple(trans[3:9])

def histit(values, nbins=10, bounds=None, normalized=False):
    """Compute histogram counts and bin centers."""
    counts, bins = scipy.histogram(values, nbins, bounds, normalized, new=True)
    centers = (bins[0:-1] + bins[1:])/2
    return counts, centers

def basic_plot(ligands, decoys, params):
    """Create ligand and decoy plots given parameters."""
    pylab.hold(False)
    ly, lx = histit(ligands, params[0][0], params[1], True)
    lp = pylab.plot(lx, ly)
    pylab.hold(True)
    dy, dx = histit(decoys, params[0][1], params[1], True)
    dp = pylab.plot(dx, dy)
    return lp, dp

def gen_plot(pdir, ligands, decoys, params, labels):
    """Generate both thumbnail and large plot images."""
    pylab.gcf().set_size_inches(8.0/5, 6.0/5)
    lp, dp = basic_plot(ligands, decoys, params)
    pylab.subplots_adjust(top=0.87)
    gca = pylab.gca()
    pylab.setp(gca.get_xticklabels(), fontsize=6)
    pylab.setp(gca.get_yticklabels(), fontsize=5)
    pylab.title(labels[0], fontsize=8)
    pylab.savefig(os.path.join(pdir, labels[1]+'_thumbnail.png'))
    pylab.gcf().set_size_inches(8.0, 6.0)
    lp, dp = basic_plot(ligands, decoys, params)
    pylab.subplots_adjust(top=0.90)
    pylab.xlabel(labels[0])
    pylab.ylabel('Fraction')
    pylab.legend((lp[0], dp[0]), ('Ligands', 'Decoys'))
    pylab.savefig(os.path.join(pdir, labels[1]+'.png'))

def plot_props(pdir, lproplist, dproplist):
    """Parse properties and generate all property plots."""
    lmw, llogp, lrotb, lhba, lhbd, lnet = parse_props(lproplist)
    dmw, dlogp, drotb, dhba, dhbd, dnet = parse_props(dproplist)
    gen_plot(pdir, llogp, dlogp, [[8, 8], [-7., 9.]],
               ['LogP', 'logp'])
    gen_plot(pdir, lmw, dmw, [[11, 11], [125., 675.]],
               ['Molecular Weight', 'mw'])
    gen_plot(pdir, lhba, dhba, [[9, 9], [-1., 17.]],
               ['Hydrogen Bond Acceptors', 'hba'])
    gen_plot(pdir, lhbd, dhbd, [[10, 10], [-0.5, 9.5]],
               ['Hydrogen Bond Donors', 'hbd'])
    gen_plot(pdir, lrotb, drotb, [[8, 8], [-1., 15.]],
               ['# Rotatable Bonds', 'rotb'])
    gen_plot(pdir, lnet, dnet, [[10, 10], [-4.5, 5.5]],
               ['Net Charge', 'net'])

def do_html(outdir):
    """Generate html pages with thumbnails and plots."""
    shutil.copy(dudhtml.STYLESHEET, outdir)
    stylesheet = os.path.basename(dudhtml.STYLESHEET)
    props = ['mw', 'logp', 'hba', 'hbd', 'rotb', 'net']
    html = open(os.path.join(outdir, "properties.html"), 'w')
    title = "Properties"
    thtml = dudhtml.SingleThumbHTML(html, stylesheet, title, props)
    thtml.page()
    html.close()

def dud_props(indir='.', outdir='.', legacy=True, loglevel=2):
    """Calculate properties and display ligand versus decoy plots in html."""
    pdir = os.path.join(outdir, PROP_DIR)
    if not os.path.exists(pdir):
        os.mkdir(pdir)
    outdir, outgen, genlig, gendec = duddb.setup_dirs(outdir)
    lproplist = collect_dir(genlig)
    dproplist = collect_dir(gendec)
    if loglevel > 1:
        fn = os.path.join(outdir, CHARGE_FILE_DEC + CHARGE_FILE_EXT)
        mmmutils.write_splits(fn, dproplist)
    plot_props(pdir, lproplist, dproplist)
    do_html(pdir)

def main(argv):
    """Parse arguments."""
    description = "Generate property graphs and html."
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
    dud_props(indir=options.indir, outdir=options.outdir)
    print 'Program took %.1f seconds.' % (time.time() - start_time)
    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv))
