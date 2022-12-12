#!/usr/bin/env python

# PIGEON

# Module
from __future__ import division
import numpy as np
import pandas as pd
from subprocess import call
from itertools import product
import time, sys, traceback, argparse
import ldsc_mod.parse as ps
import ldsc_mod.sumstats as sumstats
import ldsc_mod.regressions as reg
from functools import reduce

TopHEAD = "*********************************************************************\n"
TopHEAD += "* Polygenic gene-environment interaction (PIGEON) \n"
TopHEAD += "* Version 1.0.0 \n"
TopHEAD += "* (C) Jiacheng Miao and Yixuan Wu \n"
TopHEAD += "* University of Wisconsin-Madison \n"
TopHEAD += "* https://github.com/qlu-lab/PIGEON \n"
TopHEAD += "* GNU General Public License v3\n"
TopHEAD += "*********************************************************************\n"
pd.set_option('display.max_rows', 1000)
pd.set_option('display.max_columns', 1000)
pd.set_option('display.width', 1000)
pd.set_option('precision', 4)
pd.set_option('max_colwidth',1000)
np.set_printoptions(linewidth=1000)
np.set_printoptions(precision=4)

class Logger(object):
    '''
    Lightweight logging.
    TODO: replace with logging module
    '''

    def __init__(self, fh):
        self.log_fh = open(fh, 'w')

    def log(self, msg):
        '''
        Print to log file and stdout with a single command.
        '''
        print(msg, file=self.log_fh)
        print(msg)

def sec_to_str(t):
    '''Convert seconds to days:hours:minutes:seconds'''
    [d, h, m, s, n] = reduce(lambda ll, b : divmod(ll[0], b) + ll[1:], [(t, 1), 60, 60, 24])
    f = ''
    if d > 0:
        f += '{D}d:'.format(D=d)
    if h > 0:
        f += '{H}h:'.format(H=h)
    if m > 0:
        f += '{M}m:'.format(M=m)

    f += '{S}s'.format(S=s)
    return f

######### Command line arguments #########
# Basic PIGEON analysis
parser = argparse.ArgumentParser()
parser.add_argument('--scan', default=False, action='store_true',
    help='This flag tells PIGEON to scan for PGSxE and covariant GxE')
parser.add_argument('--cond', default=False, action='store_true',
    help='This flag tells PIGEON to conduct conditional PGSxE analysis')
parser.add_argument('--gxe-var', default=False, action='store_true',
    help='This flag tells PIGEON to estimate GxE variance component')

# Required flag in PIGEON analysis
parser.add_argument('--out',type=str,help='Prefix of the output filename')
parser.add_argument('--gxe-sumstats',type=str,help='Path to GWIS (SNPxE) summary statistics')

# Optional flag in PIGEON analysis
parser.add_argument('--gwas-sumstats',type=str,help='Path to GWAS summary statistics')
parser.add_argument('--e-sumstats',type=str,help='Path to GWAS summary statistics of environment')
parser.add_argument('--covariates-R2',default=0.0, type=float, help='R-squared for the regression between the phenotype and environments/covaraites')
parser.add_argument('--e-sd', default=1.0, type=float,help='Standard deviation of the environment')
parser.add_argument('--no-intercept', action='store_true',
    help = 'Constrain the intercept for oraclePGSxE analysis to be zero.')
# Optional parameter in PIGEON analysis
parser.add_argument('--two-step', default=None, type=float,
    help='Test statistic bound for use with the two-step estimator. Not compatible with --no-intercept')
parser.add_argument('--chisq-max', default=None, type=float,
    help='Max chi^2.')
parser.add_argument('--M', default=None, type=str,
    help='# of SNPs (if you don\'t want to use the .l2.M files that came with your .l2.PIGEONore.gz files)')
parser.add_argument('--n-blocks', default=200, type=int,
    help='Number of block jackknife blocks.')
parser.add_argument('--not-M-5-50', default=False, action='store_true',
    help='This flag tells PIGEON to use the .l2.M file instead of the .l2.M_5_50 file.')

# Read the LD score
parser.add_argument('--ref-ld-chr',type=str,help='Index of SNP in genofile at which to finish computing test stats')
parser.add_argument('--w-ld-chr',type=str,help='Location of the phenotype file')
parser.add_argument('--ref-ld',type=str,help='Index of SNP in genofile at which to finish computing test stats', default=None)
parser.add_argument('--w-ld',type=str,help='Location of the phenotype file', default=None)

# Liability scale
parser.add_argument('--samp-prev',default=None,
    help='Sample prevalence of binary phenotype (for conversion to liability scale).')
parser.add_argument('--pop-prev',default=None,
    help='Population prevalence of binary phenotype (for conversion to liability scale).')

# Tmp
parser.add_argument('--no-check-alleles', default=False, action='store_true',
    help='For covariant GxE estimation, skip checking whether the alleles match. This check is '
    'redundant for pairs of chisq files generated using munge_sumstats.py and the '
    'same argument to the --merge-alleles flag.')
parser.add_argument('--intercept-h2', action='store', default=None,
    help = 'Intercepts for constrained-intercept single-trait.')
parser.add_argument('--intercept-gencov', action='store', default=None,
    help = 'Intercepts for constrained-intercept cross-trait.'
    ' Must have same length as --rg. The first entry is ignored.')
parser.add_argument('--stratified', default=False, action='store_true',
    help='This flag tells PIGEON to conduct')
parser.add_argument('--overlap-annot', default=False, action='store_true',
    help='This flag informs PIGEON that the partitioned LD Scores were generates using an '
    'annot matrix with overlapping categories (i.e., not all row sums equal 1), '
    'and prevents PIGEON from displaying output that is meaningless with overlapping categories.')
parser.add_argument('--annot', default=None, type=str,
    help='Filename prefix for annotation file for partitioned LD Score estimation. '
    'PIGEON will automatically append .annot or .annot.gz to the filename prefix. '
    'See docs/file_formats_ld for a definition of the .annot format.')
parser.add_argument('--thin-annot', action='store_true', default=False,
    help='This flag says your annot files have only annotations, with no SNP, CM, CHR, BP columns.')
parser.add_argument('--cts-bin', default=None, type=str,
    help='This flag tells PIGEON to compute partitioned LD Scores, where the partition '
    'is defined by cutting one or several continuous variable[s] into bins. '
    'The argument to this flag should be the name of a single file or a comma-separated '
    'list of files. The file format is two columns, with SNP IDs in the first column '
    'and the continuous variable in the second column. ')
parser.add_argument('--cts-breaks', default=None, type=str,
    help='Use this flag to specify names for the continuous variables cut into bins '
    'with --cts-bin. For each continuous variable, specify breaks as a comma-separated '
    'list of breakpoints, and separate the breakpoints for each variable with an x. '
    'For example, if binning on MAF and distance to gene (in kb), '
    'you might set --cts-breaks 0.1,0.25,0.4x10,100,1000 ')
parser.add_argument('--cts-names', default=None, type=str,
    help='Use this flag to specify names for the continuous variables cut into bins '
    'with --cts-bin. The argument to this flag should be a comma-separated list of '
    'names. For example, if binning on DAF and distance to gene, you might set '
    '--cts-bin DAF,DIST_TO_GENE ')
parser.add_argument('--print-coefficients',default=False,action='store_true',
    help='when categories are overlapping, print coefficients as well as heritabilities.')
parser.add_argument('--frqfile', type=str,
    help='For use with --overlap-annot. Provides allele frequencies to prune to common '
    'snps if --not-M-5-50 is not set.')
parser.add_argument('--frqfile-chr', type=str,
    help='Prefix for --frqfile files split over chromosome.')
parser.add_argument('--ref-ld-chr-cts', default=None, type=str,
    help='Name of a file that has a list of file name prefixes for cell-type-specific analysis.')
parser.add_argument('--print-all-cts', action='store_true', default=False)
parser.add_argument('--print-cov', default=False, action='store_true',
    help='For use with --h2/--rg. This flag tells PIGEON to print the '
    'covaraince matrix of the estimates.')
parser.add_argument('--print-delete-vals', default=False, action='store_true',
    help='If this flag is set, PIGEON will print the block jackknife delete-values ('
    'i.e., the regression coefficeints estimated from the data with a block removed). '
    'The delete-values are formatted as a matrix with (# of jackknife blocks) rows and '
    '(# of LD Scores) columns.')

if __name__ == '__main__':
    args = parser.parse_args()
    if args.out is None:
        raise ValueError('--out is required.')

    # Define the file output name
    out_name = args.out
    if args.scan: out_name += '_Scan'
    if args.cond: out_name += '_Cond'
    if args.gxe_var: out_name += '_GxE_Var'
    if args.e_sumstats is not None: out_name += '-rGE'
    
    log = Logger(out_name+'.log')
    defaults = vars(parser.parse_args(''))
    opts = vars(args)
    non_defaults = [x for x in opts.keys() if opts[x] != defaults[x]]
    header = TopHEAD    
    header += '\n### Options in effect ### \n'
    header += 'pigeon.py \\\n'
    options = ['--'+x.replace('_','-')+' '+str(opts[x])+' \\' for x in non_defaults]
    header += '\n'.join(options).replace('True','').replace('False','')
    header = header[0:-1]+'\n'
    log.log(header)
    log.log('### Analysis began at {T} ###'.format(T=time.ctime()))
    log.log('\n### Parse the input data ###')
    start_time = time.time()
    try:
        # PIGEON Oraclae PGSxE scan
        if args.scan: 
            sumstats.estimate_oPGSxE_scan(args, log) if args.e_sumstats is None else sumstats.estimate_oPGSxe_scan(args, log)
        # PIGEON cond analysis
        elif args.cond:
            sumstats.estimate_oPGSxE_cond(args, log)
        # PIGEON Estimate GxE Heritability
        elif args.gxe_var:
            sumstats.estimate_GxE_var(args, log) if args.e_sumstats is None else sumstats.estimate_GxE_var_e(args, log)
    except Exception:
        ex_type, ex, tb = sys.exc_info()
        log.log(traceback.format_exception_only(Exception, ex))
        raise
    finally:
        log.log('--- Analysis finished at {T}'.format(T=time.ctime()))
        time_elapsed = round(time.time()-start_time, 2)
        log.log('--- Total time elapsed: {T}'.format(T=sec_to_str(time_elapsed)))
