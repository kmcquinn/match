"""Class for reading the output.cmd file from calcsfh"""
from __future__ import print_function
import argparse
import os
import sys

import matplotlib.pyplot as plt
import numpy as np

from .config import EXT
from .fileio import read_match_cmd
from .graphics import match_plot
from .utils import parse_pipeline

__all__ = ['CMD']

class CMD(object):
    """
    A quikly made object to read the MATCH CMD file and hold paramters to
    automatically make plots with the same color scale as other MATCH CMD files.
    """
    def __init__(self, filename):
        self.base, self.name = os.path.split(os.path.abspath(filename))
        self.cmd, self.fit, self.colors, self.yfilter = read_match_cmd(filename)
        self.load_match_cmd()

    def set_labels(self):
        """Set up list of labels for pgpro"""
        strfmt = r'${{\rm {:s}}}$'
        labels = [strfmt.format(i) for i in ['data', 'model', 'd-m']]
        try:
            target, _ = parse_pipeline(self.name)
            labels[0] = strfmt.format(target)
        except:
            pass
        labels.extend(r'$-2\ln P = {:g}$'.format(self.fit))
        return labels

    def load_match_cmd(self):
        """
        pgcmd needs hesses and extent. Optional are max_* which set the vmins
        and vmaxs.
        """
        self.nmagbin = len(np.unique(self.cmd['mag']))
        self.ncolbin = len(np.unique(self.cmd['color']))
        self.data = self.cmd['Nobs'].reshape(self.nmagbin, self.ncolbin)
        self.model = self.cmd['Nsim'].reshape(self.nmagbin, self.ncolbin)
        self.diff = self.cmd['diff'].reshape(self.nmagbin, self.ncolbin)
        self.sig = self.cmd['sig'].reshape(self.nmagbin, self.ncolbin)
        self.hesses = [self.data, self.model, self.diff, self.sig]

        self.extent = [np.min(self.cmd['color']), np.max(self.cmd['color']),
                       np.max(self.cmd['mag']), np.min(self.cmd['mag'])]
        #self.max_counts = np.nanmax(np.concatenate([self.data, self.model]))
        #self.max_diff = np.nanmax(np.abs(self.diff))
        #self.max_sig = np.nanmax(np.abs(self.sig))


    def pgcmd(self, labels=None, outdir=None, logcounts=True, figname=None):
        '''produce the image that pgcmd.pro makes'''
        labels = labels or self.set_labels()

        if figname is None:
            base = outdir or self.base
            assert os.path.isdir(base), \
                '{} directory not found'.format(base)
            figname = os.path.join(base, os.path.split(self.name)[1] + EXT)

        hesses = self.hesses
        if logcounts:
            hesses[0] = np.log10(hesses[0])
            hesses[1] = np.log10(hesses[1])

        xlabel = r'${}$'.format(self.colors)
        ylabel = r'${}$'.format(self.yfilter)

        grid = match_plot(hesses, self.extent, labels=labels, ylabel=ylabel,
                          xlabel=xlabel)

        gates = self.cmd['gate']
        ugates = np.unique(gates)
        if len(ugates) > 1:
            dinds = np.digitize(gates, bins=np.unique(gates, right=True))
            _ = [grid.axes_all[0].plot(self.cmd['color'][dinds == i],
                                       self.cmd['mag'][dinds == i],
                                       '.', alpha=0.3)
                 for i in range(len(dinds)) if i > 0]

        plt.savefig(figname)
        plt.close()
        print('wrote {}'.format(figname))
        return grid


def call_pgcmd_byfit(cmdfns, nmax=5, outdir=None):
    """
    Call pgcmd with filename order them by increasing best fit value.
    cmdfns : string or list
        .cmd filename or list of filenames.

    nmax : int
        make best nmax plots.
    """
    if not isinstance(cmdfns, list):
        cmd = CMD(cmdfns)
        cmd.pgcmd()
        return

    cmds = np.array([CMD(cmdfn) for cmdfn in cmdfns])
    icmd = np.argsort(np.concatenate([cmd.fit for cmd in cmds]))
    for j, cmd in enumerate(cmds[icmd]):
        if j > nmax:
            break
        jstr = ('{}'.format(j)).zfill(4)
        cmd.pgcmd(figname=cmd.figname.replace(EXT, ''.join(jstr, EXT)),
                  outdir=outdir)
    return

def main(argv):
    """main function for cmd"""
    parser = argparse.ArgumentParser(description="plot cmd file")

    parser.add_argument('-o', '--outdir', type=str,
                        help='directory to place images')

    parser.add_argument('-f', '--byfit', action='store_true',
                        help='number filenames by best fit')

    parser.add_argument('cmdfiles', type=str, nargs='*',
                        help='.cmd files to plot')

    args = parser.parse_args(argv)

    if args.outdir is not None:
        assert os.path.isdir(args.outdir), \
            'directory {} not found'.format(args.outdir)

    if args.byfit:
        call_pgcmd_byfit(args.cmdfiles, nmax=len(args.cmdfiles),
                         outdir=args.outdir)
    else:
        for cmdfile in args.cmdfiles:
            cmd = CMD(cmdfile)
            cmd.pgcmd(outdir=args.outdir)

if __name__ == "__main__":
    main(sys.argv[1:])
