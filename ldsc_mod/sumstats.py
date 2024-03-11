import numpy as np
import pandas as pd
from scipy import stats
import itertools as it
from . import parse as ps
from . import regressions as reg
import sys
import traceback
import copy
import os
import re

_N_CHR = 22
# complementary bases
COMPLEMENT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
# bases
BASES = COMPLEMENT.keys()
# true iff strand ambiguous
STRAND_AMBIGUOUS = {''.join(x): x[0] == COMPLEMENT[x[1]]
                    for x in it.product(BASES, BASES)
                    if x[0] != x[1]}
# SNPS we want to keep (pairs of alleles)
VALID_SNPS = {x for x in map(lambda y: ''.join(y), it.product(BASES, BASES))
              if x[0] != x[1] and not STRAND_AMBIGUOUS[x]}
# T iff SNP 1 has the same alleles as SNP 2 (allowing for strand or ref allele flip).
MATCH_ALLELES = {x for x in map(lambda y: ''.join(y), it.product(VALID_SNPS, VALID_SNPS))
                 # strand and ref match
                 if ((x[0] == x[2]) and (x[1] == x[3])) or
                 # ref match, strand flip
                 ((x[0] == COMPLEMENT[x[2]]) and (x[1] == COMPLEMENT[x[3]])) or
                 # ref flip, strand match
                 ((x[0] == x[3]) and (x[1] == x[2])) or
                 ((x[0] == COMPLEMENT[x[3]]) and (x[1] == COMPLEMENT[x[2]]))}  # strand and ref flip
# T iff SNP 1 has the same alleles as SNP 2 w/ ref allele flip.
FLIP_ALLELES = {''.join(x):
                ((x[0] == x[3]) and (x[1] == x[2])) or  # strand match
                # strand flip
                ((x[0] == COMPLEMENT[x[3]]) and (x[1] == COMPLEMENT[x[2]]))
                for x in MATCH_ALLELES}


def _splitp(fstr):
    flist = fstr.split(',')
    flist = [os.path.expanduser(os.path.expandvars(x)) for x in flist]
    return flist


def _select_and_log(x, ii, log, msg):
    '''Fiter down to rows that are True in ii. Log # of SNPs removed.'''
    new_len = ii.sum()
    if new_len == 0:
        raise ValueError(msg.format(N=0))
    else:
        x = x[ii]
        log.log(msg.format(N=new_len))
    return x


def smart_merge(x, y):
    '''Check if SNP columns are equal. If so, save time by using concat instead of merge.'''
    if len(x) == len(y) and (x.index == y.index).all() and (x.SNP == y.SNP).all():
        x = x.reset_index(drop=True)
        y = y.reset_index(drop=True).drop('SNP', 1)
        out = pd.concat([x, y], axis=1)
    else:
        out = pd.merge(x, y, how='inner', on='SNP')
    return out


def _read_ref_ld(args, log):
    '''Read reference LD Scores.'''
    ref_ld = _read_chr_split_files(args.ref_ld_chr, args.ref_ld, log,
                                   'reference panel LD Score', ps.ldscore_fromlist)
    log.log(
        '--- {N} SNPs found in the reference panel LD Scores.'.format(N=len(ref_ld)))
    return ref_ld


def _read_annot(args, log):
    '''Read annot matrix.'''
    try:
        if args.ref_ld is not None:
            overlap_matrix, M_tot = _read_chr_split_files(args.ref_ld_chr, args.ref_ld, log,
                                                          'annot matrix', ps.annot, frqfile=args.frqfile)
        elif args.ref_ld_chr is not None:
            overlap_matrix, M_tot = _read_chr_split_files(args.ref_ld_chr, args.ref_ld, log,
                                                      'annot matrix', ps.annot, frqfile=args.frqfile_chr)
    except Exception:
        log.log('Error parsing .annot file.')
        raise

    return overlap_matrix, M_tot


def _read_M(args, log, n_annot):
    '''Read M (--M, --M-file, etc).'''
    if args.M:
        try:
            M_annot = [float(x) for x in _splitp(args.M)]
        except ValueError as e:
            raise ValueError('Could not cast --M to float: ' + str(e.args))
    else:
        if args.ref_ld:
            M_annot = ps.M_fromlist(
                _splitp(args.ref_ld), common=(not args.not_M_5_50))
        elif args.ref_ld_chr:
            M_annot = ps.M_fromlist(
                _splitp(args.ref_ld_chr), _N_CHR, common=(not args.not_M_5_50))

    try:
        M_annot = np.array(M_annot).reshape((1, n_annot))
    except ValueError as e:
        raise ValueError(
            '# terms in --M must match # of LD Scores in --ref-ld.\n' + str(e.args))

    return M_annot


def _read_w_ld(args, log):
    '''Read regression SNP LD.'''
    if (args.w_ld and ',' in args.w_ld) or (args.w_ld_chr and ',' in args.w_ld_chr):
        raise ValueError(
            '--w-ld must point to a single fileset (no commas allowed).')
    w_ld = _read_chr_split_files(args.w_ld_chr, args.w_ld, log,
                                 'regression weight LD Score', ps.ldscore_fromlist)
    if len(w_ld.columns) != 2:
        raise ValueError('--w-ld may only have one LD Score column.')
    w_ld.columns = ['SNP', 'LD_weights']  # prevent colname conflicts w/ ref ld
    log.log(
        '--- {N} SNPs found in regression weight LD Scores.'.format(N=len(w_ld)))
    return w_ld


def _read_chr_split_files(chr_arg, not_chr_arg, log, noun, parsefunc, **kwargs):
    '''Read files split across 22 chromosomes (annot, ref_ld, w_ld).'''
    try:
        if not_chr_arg:
            log.log('--- Reading {N} from {F} ... ({p})'.format(N=noun, F=not_chr_arg, p=parsefunc.__name__))
            out = parsefunc(_splitp(not_chr_arg), **kwargs)
        elif chr_arg:
            f = ps.sub_chr(chr_arg, '[1-22]')
            log.log('--- Reading {N} from {F} ... ({p})'.format(N=noun, F=f, p=parsefunc.__name__))
            out = parsefunc(_splitp(chr_arg), _N_CHR, **kwargs)
    except ValueError as e:
        log.log('Error parsing {N}.'.format(N=noun))
        raise e

    return out


def _read_sumstats(args, log, fh, alleles=False, dropna=False):
    '''Parse summary statistics.'''
    log.log('--- Reading summary statistics from {S}'.format(S=fh))
    sumstats = ps.sumstats(fh, alleles=alleles, dropna=dropna)
    log_msg = '--- {N} SNPs found in the summary statistics'
    log.log(log_msg.format(N=len(sumstats)))
    m = len(sumstats)
    sumstats = sumstats.drop_duplicates(subset='SNP')
    if m > len(sumstats):
        log.log(
            'Dropped {M} SNPs with duplicated rs numbers.'.format(M=m - len(sumstats)))

    return sumstats


def _check_ld_condnum(args, log, ref_ld):
    '''Check condition number of LD Score matrix.'''
    if len(ref_ld.shape) >= 2:
        cond_num = int(np.linalg.cond(ref_ld))
        if cond_num > 100000:
            if args.invert_anyway:
                warn = "WARNING: LD Score matrix condition number is {C}. "
                warn += "Inverting anyway because the --invert-anyway flag is set."
                log.log(warn.format(C=cond_num))
            else:
                warn = "WARNING: LD Score matrix condition number is {C}. "
                warn += "Remove collinear LD Scores. "
                raise ValueError(warn.format(C=cond_num))


def _check_variance(log, M_annot, ref_ld):
    '''Remove zero-variance LD Scores.'''
    ii = ref_ld.iloc[:, 1:].var() == 0  # NB there is a SNP column here
    if ii.all():
        raise ValueError('All LD Scores have zero variance.')
    else:
        # log.log('--- Removing partitioned LD Scores with zero variance.')
        ii_snp = np.array([True] + list(~ii))
        ii_m = np.array(~ii)
        ref_ld = ref_ld.iloc[:, ii_snp]
        M_annot = M_annot[:, ii_m]

    return M_annot, ref_ld, ii


def _warn_length(log, sumstats):
    if len(sumstats) < 200000:
        log.log(
            'WARNING: number of SNPs less than 200k; this is almost always bad.')


def _print_cov(ldscore_reg, ofh, log):
    '''Prints covariance matrix of slopes.'''
    log.log(
        'Printing covariance matrix of the estimates to {F}.'.format(F=ofh))
    np.savetxt(ofh, ldscore_reg.coef_cov)


def _print_delete_values(ldscore_reg, ofh, log):
    '''Prints block jackknife delete-k values'''
    log.log('Printing block jackknife delete values to {F}.'.format(F=ofh))
    np.savetxt(ofh, ldscore_reg.tot_delete_values)

def _print_part_delete_values(ldscore_reg, ofh, log):
    '''Prints partitioned block jackknife delete-k values'''
    log.log('Printing partitioned block jackknife delete values to {F}.'.format(F=ofh))
    np.savetxt(ofh, ldscore_reg.part_delete_values)


def _merge_and_log(ld, sumstats, noun, log):
    '''Wrap smart merge with log messages about # of SNPs.'''
    sumstats = smart_merge(ld, sumstats)
    msg = '--- {N} SNPs remain after merging with {F}.'
    if len(sumstats) == 0:
        raise ValueError(msg.format(N=len(sumstats), F=noun))
    else:
        log.log(msg.format(N=len(sumstats), F=noun))

    return sumstats


def _read_ld_sumstats(args, log, fh, alleles=False, dropna=True):
    sumstats = _read_sumstats(args, log, fh, alleles=alleles, dropna=dropna)
    ref_ld = _read_ref_ld(args, log)
    n_annot = len(ref_ld.columns) - 1
    M_annot = _read_M(args, log, n_annot)
    M_annot, ref_ld, novar_cols = _check_variance(log, M_annot, ref_ld)
    w_ld = _read_w_ld(args, log)
    sumstats = _merge_and_log(ref_ld, sumstats, 'reference panel LD', log)
    sumstats = _merge_and_log(sumstats, w_ld, 'regression SNP LD', log)
    w_ld_cname = sumstats.columns[-1]
    ref_ld_cnames = ref_ld.columns[1:len(ref_ld.columns)]
    return M_annot, w_ld_cname, ref_ld_cnames, sumstats, novar_cols

def estimate_h2(args, log):
    '''Estimate h2 and partitioned h2.'''
    args = copy.deepcopy(args)
    if args.samp_prev is not None and args.pop_prev is not None:
        args.samp_prev, args.pop_prev =map(
            float, [args.samp_prev, args.pop_prev])
    if args.intercept_h2 is not None:
        args.intercept_h2 = float(args.intercept_h2)
    if args.no_intercept:
        args.intercept_h2 = 1
    M_annot, w_ld_cname, ref_ld_cnames, sumstats, novar_cols = _read_ld_sumstats(
        args, log, args.h2)
    ref_ld = np.array(sumstats[ref_ld_cnames])
    _check_ld_condnum(args, log, ref_ld_cnames)
    _warn_length(log, sumstats)
    n_snp = len(sumstats)
    n_blocks = min(n_snp, args.n_blocks)
    n_annot = len(ref_ld_cnames)
    chisq_max = args.chisq_max
    old_weights = False
    if n_annot == 1:
        if args.two_step is None and args.intercept_h2 is None:
            args.two_step = 30
    else:
        old_weights = True
        if args.chisq_max is None:
            chisq_max = max(0.001*sumstats.N.max(), 80)

    def s(x): return np.array(x).reshape((n_snp, 1))
    chisq = s(sumstats.Z**2)
    if chisq_max is not None:
        ii = np.ravel(chisq < chisq_max)
        sumstats = sumstats.iloc[ii, :]
        log.log('Removed {M} SNPs with chi^2 > {C} ({N} SNPs remain)'.format(
                C=chisq_max, N=np.sum(ii), M=n_snp-np.sum(ii)))
        n_snp = np.sum(ii)  # lambdas are late-binding, so this works
        ref_ld = np.array(sumstats[ref_ld_cnames])
        chisq = chisq[ii].reshape((n_snp, 1))

    if args.two_step is not None:
        log.log('--- Using two-step estimator with cutoff at {M}.'.format(M=args.two_step))

    hsqhat = reg.Hsq(chisq, ref_ld, s(sumstats[w_ld_cname]), s(sumstats.N),
                     M_annot, n_blocks=n_blocks, intercept=args.intercept_h2,
                     twostep=args.two_step, old_weights=old_weights)

    if args.print_cov:
        _print_cov(hsqhat, args.out + '.cov', log)
    if args.print_delete_vals:
        _print_delete_values(hsqhat, args.out + '.delete', log)
        _print_part_delete_values(hsqhat, args.out + '.part_delete', log)

    log.log(hsqhat.summary(ref_ld_cnames, P=args.samp_prev, K=args.pop_prev, overlap = args.overlap_annot))
    if args.overlap_annot:
        overlap_matrix, M_tot = _read_annot(args, log)

        # overlap_matrix = overlap_matrix[np.array(~novar_cols), np.array(~novar_cols)]#np.logical_not
        df_results = hsqhat._overlap_output(ref_ld_cnames, overlap_matrix, M_annot, M_tot, args.print_coefficients)
        df_results.to_csv(args.out+'.results', sep="\t", index=False)
        log.log('Results printed to '+args.out+'.results')

    return hsqhat

# oPGSxE_scan
def estimate_oPGSxE_scan(args, log):
    '''Estimate oPGSxE between trait 1 and a list of other traits.'''
    args = copy.deepcopy(args)
    rg_paths, rg_files = _parse_rg(args.gxe_sumstats, args.gwas_sumstats)
    n_pheno = len(rg_paths)
    f = lambda x: _split_or_none(x, n_pheno)
    args.intercept_h2, args.intercept_gencov, args.samp_prev, args.pop_prev = map(f,
        (args.intercept_h2, args.intercept_gencov, args.samp_prev, args.pop_prev))
    map(lambda x: _check_arg_len(x, n_pheno), ((args.intercept_h2, '--intercept-h2'),
                                               (args.intercept_gencov, '--intercept-gencov'),
                                               (args.samp_prev, '--samp-prev'),
                                               (args.pop_prev, '--pop-prev')))
    if args.no_intercept:
        args.intercept_gencov = [0 for _ in range(n_pheno)]
    p1 = rg_paths[0]
    out_prefix = args.out + rg_files[0]
    M_annot, w_ld_cname, ref_ld_cnames, sumstats, _ = _read_ld_sumstats(args, log, p1,
                                                                        alleles=True, dropna=True)


    oPGSxE = []
    n_annot = M_annot.shape[1]
    if n_annot == 1 and args.two_step is None and args.intercept_h2 is None:
        args.two_step = 30
    if args.two_step is not None:
        log.log('--- Using two-step estimator with cutoff at {M}.\n'.format(M=args.two_step))
    # added
    if args.stratified:
        df_gencov_list = []
        df_oPGSxE_scan_list = []

    for i, p2 in enumerate(rg_paths[1:n_pheno]):
        log.log(
            '\n### Computing Oracle PGSxE for phenotype {I}/{N} ###'.format(I=i + 1, N=len(rg_paths)-1))
        try:
            loop = _read_other_sumstats(args, log, p2, sumstats, ref_ld_cnames) # This will output the table
            rghat = _oPGSxE_scan(loop, args, log, M_annot, ref_ld_cnames, w_ld_cname, i)
            oPGSxE.append(rghat)
            _print_gencor(args, log, rghat, ref_ld_cnames, i, rg_paths, i == 0)
            out_prefix_loop = out_prefix + '_' + rg_files[i + 1]
            if args.print_cov:
                _print_rg_cov(rghat, out_prefix_loop, log)
            if args.print_delete_vals:
                _print_rg_delete_values(rghat, out_prefix_loop, log)
            # added
            if args.stratified:
                df_gencov_list.append(rghat.gencov.format_output(ref_ld_cnames, gxe=rg_paths[0], gwas=p2))
                df_oPGSxE_scan_list.append(rghat.format_output(ref_ld_cnames, gxe=rg_paths[0], gwas=p2))                 

        except Exception:  # keep going if phenotype 50/100 causes an error
            msg = 'ERROR computing Oracle PGSxE for phenotype {I}/{N}, from file {F}.'
            log.log(msg.format(I=i + 2, N=len(rg_paths), F=rg_paths[i + 1]))
            ex_type, ex, tb = sys.exc_info()
            log.log(traceback.format_exc(ex) + '\n')
            if len(oPGSxE) <= i:  # if exception raised before appending to oPGSxE
                oPGSxE.append(None)

    if args.stratified:
        pd.concat(df_gencov_list, axis=0).to_csv(args.out + '_Covariant_GxE.txt', index=False, sep='\t')
        pd.concat(df_oPGSxE_scan_list, axis=0).to_csv(args.out + '_Oracle_PGSxE.txt', index=False, sep='\t')

    # log.log('\n### Summary of Oracle PGSxE Results ###\n' +
    #         _get_oPGSxE_scan_table(rg_paths, oPGSxE, args))
    log.log('### Writing out the results ###')
    _get_oPGSxE_scan_table(rg_paths, oPGSxE, args)
    log.log('--- The oracle PGSxE results have been writen to {M}_Scan.txt'.format(M=args.out))

    
    
    return oPGSxE

# oPGSxE_cond
def estimate_oPGSxE_cond(args, log):
    '''Estimate oPGSxE between trait 1 and a list of other traits.'''
    args = copy.deepcopy(args)
    rg_paths, rg_files = _parse_rg(args.gxe_sumstats, args.gwas_sumstats)
    if args.e_sumstats is not None:
        e_paths =  _parse_e_sumstats(args.e_sumstats) # add
        p3 = e_paths # add
    n_pheno = len(rg_paths)
    f = lambda x: _split_or_none(x, n_pheno)
    args.intercept_h2, args.intercept_gencov, args.samp_prev, args.pop_prev = map(f,
        (args.intercept_h2, args.intercept_gencov, args.samp_prev, args.pop_prev))
    map(lambda x: _check_arg_len(x, n_pheno), ((args.intercept_h2, '--intercept-h2'),
                                               (args.intercept_gencov, '--intercept-gencov'),
                                               (args.samp_prev, '--samp-prev'),
                                               (args.pop_prev, '--pop-prev')))
    if args.no_intercept:
        args.intercept_gencov = [0 for _ in range(n_pheno)]
    p1 = rg_paths[0]
    out_prefix = args.out + rg_files[0]
    M_annot, w_ld_cname, ref_ld_cnames, sumstats, _ = _read_ld_sumstats(args, log, p1,
                                                                        alleles=True, dropna=True)


    oPGSxE = []
    n_annot = M_annot.shape[1]
    if n_annot == 1 and args.two_step is None and args.intercept_h2 is None:
        args.two_step = 30
    if args.two_step is not None:
        log.log('--- Using two-step estimator with cutoff at {M}.\n'.format(M=args.two_step))
    
    ## calculating delta / using cond analysis
    loop = _read_multiple_sumstats(args, log, rg_paths[1:3], sumstats, ref_ld_cnames)
    log.log('\n### Performing conditional Oracle PGSxE analysis ###\n')
    rghat_test = _oPGSxE_cond(loop, args, log, M_annot, ref_ld_cnames, w_ld_cname, 0)
    oPGSxE.append(rghat_test)
    out_prefix_loop = 'out_prefix_loop'
    if args.print_cov:
            _print_rg_cov_mul(rghat_test, out_prefix_loop, log)
    if args.print_delete_vals:
            _print_rg_delete_values_mul(rghat_test, out_prefix_loop, log)

    # log.log('\n### Summary of Conditional Oracle PGSxE Results ###\n' +
    #          _get_oPGSxE_cond_table(rg_paths, oPGSxE, args))
    log.log('### Writing out the results ###')
    _get_oPGSxE_cond_table(rg_paths, oPGSxE, args)
    log.log('--- The conditional oracle PGSxE results have been writen to {M}_Cond.txt'.format(M=args.out))
    
    return oPGSxE

# oPGSxe scan
def estimate_oPGSxe_scan(args, log):
    '''Estimate oPGSxE between trait 1 and a list of other traits.'''
    args = copy.deepcopy(args)
    rg_paths, rg_files = _parse_rg(args.gxe_sumstats, args.gwas_sumstats)
    e_paths =  _parse_e_sumstats(args.e_sumstats) # add
    n_pheno = len(rg_paths)
    f = lambda x: _split_or_none(x, n_pheno)
    args.intercept_h2, args.intercept_gencov, args.samp_prev, args.pop_prev =map(f,
        (args.intercept_h2, args.intercept_gencov, args.samp_prev, args.pop_prev))
    map(lambda x: _check_arg_len(x, n_pheno), ((args.intercept_h2, '--intercept-h2'),
                                               (args.intercept_gencov, '--intercept-gencov'),
                                               (args.samp_prev, '--samp-prev'),
                                               (args.pop_prev, '--pop-prev')))
    if args.no_intercept:
        args.intercept_gencov = [0 for _ in range(n_pheno)]
    p1 = rg_paths[0]
    out_prefix = args.out + rg_files[0]
    M_annot, w_ld_cname, ref_ld_cnames, sumstats, _ = _read_ld_sumstats(args, log, p1,
                                                                        alleles=True, dropna=True)

    oPGSxe = []
    n_annot = M_annot.shape[1]
    if n_annot == 1 and args.two_step is None and args.intercept_h2 is None:
        args.two_step = 30
    if args.two_step is not None:
        log.log('--- Using two-step estimator with cutoff at {M}.\n'.format(M=args.two_step))
    
    loop = _read_multiple_sumstats(args, log, [e_paths]+rg_paths[1:], sumstats, ref_ld_cnames)
    rghat_test = _oPGSxe_scan(loop, args, log, M_annot, ref_ld_cnames, w_ld_cname, 0)
    oPGSxe.append(rghat_test)
    out_prefix_loop = 'out_prefix_loop'
    if args.print_cov:
            _print_rg_cov_mul(rghat_test, out_prefix_loop, log)
    if args.print_delete_vals:
            _print_rg_delete_values_mul(rghat_test, out_prefix_loop, log)
    
    # log.log('\n### Summary of Oracle PGSxE Results ###\n' +
    #         _get_oPGSxe_scan_table(rg_paths, e_paths, oPGSxe, args))
    log.log('### Writing out the results ###')
    _get_oPGSxe_scan_table(rg_paths, e_paths, oPGSxe, args)
    log.log('--- The oracle PGSxE results have been writen to {M}_Scan-rGE.txt'.format(M=args.out))
    return oPGSxe

def estimate_GxE_var(args, log):
    '''Estimate GxE Variance Component between trait 1 and a list of other traits.'''
    args = copy.deepcopy(args)

    
    M_annot, w_ld_cname, ref_ld_cnames, sumstats, _ = _read_ld_sumstats(args, log, args.gxe_sumstats,
                                                                        alleles=True, dropna=True)

    n_annot = M_annot.shape[1]

    log.log('\n### Estimating the GxE variance component ###')
    if n_annot == 1 and args.two_step is None and args.intercept_h2 is None:
        args.two_step = 30
    # if args.two_step is not None:
    #     log.log('--- Using two-step estimator with cutoff at {M}.'.format(M=args.two_step))
    
    hsq = _GxE_var(sumstats, args, log, M_annot, ref_ld_cnames, w_ld_cname)
    _print_hsq(args, log, hsq, ref_ld_cnames)

    out_prefix = args.out + args.gxe_sumstats.split('/')[-1]
    if args.print_cov:
            _print_cov(hsq, out_prefix + '.hsq.cov', log)
    if args.print_delete_vals:
            _print_delete_values(hsq, out_prefix + '.hsq.delete', log)
    
    if args.stratified:
        hsq.format_output(ref_ld_cnames, P=args.samp_prev, K=args.pop_prev, gxe=args.gxe_sumstats, gwas=args.gwas_sumstats).to_csv(args.out + '_GxE_Heritability.txt', index=False, sep='\t')

    # log.log('\n### Summary of GxE Variance Component Results ###\n' +
    #          _get_GxE_var_table(hsq, args))
    log.log('### Writing out the results ###')
    _get_GxE_var_table(hsq, args)
    log.log('--- The GxE variance component results have been writen to {M}_GxE_Var.txt'.format(M=args.out))

    return hsq

def estimate_GxE_var_e(args, log):
    '''Estimate GxE Variance Component between trait 1 and a list of other traits.'''
    args = copy.deepcopy(args)
    
    M_annot, w_ld_cname, ref_ld_cnames, sumstats, _ = _read_ld_sumstats(args, log, args.gxe_sumstats,
                                                                        alleles=True, dropna=True)

    n_pheno = 1
    f = lambda x: _split_or_none(x, n_pheno)
    args.intercept_h2, args.intercept_gencov = map(f,(args.intercept_h2, args.intercept_gencov))
    map(lambda x: _check_arg_len(x, n_pheno), ((args.intercept_h2, '--intercept-h2'), (args.intercept_gencov, '--intercept-gencov')))
    if len(args.intercept_h2) == 1:
        args.intercept_h2 = [None, None]
    n_annot = M_annot.shape[1]
    e_paths =  _parse_e_sumstats(args.e_sumstats) # add
    loop = _read_other_sumstats(args, log, e_paths, sumstats, ref_ld_cnames)

    if n_annot == 1 and args.two_step is None and args.intercept_h2 is None:
        args.two_step = 30
    # if args.two_step is not None:
        # log.log('--- Using two-step estimator with cutoff at {M}.\n'.format(M=args.two_step))
    
    log.log('\n### Estimating the GxE variance component')
    hsq = _GxE_var_e(loop, args, log, M_annot, ref_ld_cnames, w_ld_cname)
    _print_hsq(args, log, hsq, ref_ld_cnames)

    out_prefix = args.out + args.gxe_sumstats.split('/')[-1]
    if args.print_cov:
            _print_cov(hsq, out_prefix + '.hsq.cov', log)
    if args.print_delete_vals:
            _print_delete_values(hsq, out_prefix + '.hsq.delete', log)
    
    if args.stratified:
        hsq.format_output(ref_ld_cnames, P=args.samp_prev, K=args.pop_prev, gxe=args.gxe_sumstats, gwas=args.gwas_sumstats).to_csv(args.out + '_GxE_Heritability.txt', index=False, sep='\t')

    # log.log('\n### Summary of GxE Variance Component Results ###\n' +
    #          _get_GxE_var_e_table(hsq, args))
    log.log('### Writing out the results ###')
    _get_GxE_var_e_table(hsq, args) 
    log.log('--- The GxE variance component Results have been writen to {M}_GxE_Var-rGE.txt'.format(M=args.out))

    return hsq

def _read_other_sumstats(args, log, p2, sumstats, ref_ld_cnames):
    loop = _read_sumstats(args, log, p2, alleles=True, dropna=False)
    # print(sumstats.head(5))
    # print(loop.head(5))
    loop = _merge_sumstats_sumstats(args, sumstats, loop, log)
    loop = loop.dropna(how='any')
    alleles = loop.A1 + loop.A2 + loop.A1x + loop.A2x
    if not args.no_check_alleles:
        li = _filter_alleles(alleles)
        loop = _select_and_log(loop, li, log,
                       '{N} SNPs with valid alleles.')
        alleles = alleles[li]
        loop['Z2'] = _align_alleles(loop.Z2, alleles)

    # loop = loop.drop(['A1', 'A1x', 'A2', 'A2x'], axis=1)
    loop = loop.drop(['A1x', 'A2x'], axis=1)
    _check_ld_condnum(args, log, loop[ref_ld_cnames])
    _warn_length(log, loop)
    return loop

# newly added function, plist is a list of path
def _read_multiple_sumstats(args, log, plist, sumstats, ref_ld_cnames):
    for i, p2 in enumerate(plist):
        loop = _read_sumstats(args, log, p2, alleles=True, dropna=False)
        loop = _merge_new_sumstats(args, sumstats, loop, log, i)
        loop = loop.dropna(how='any')
        alleles = loop.A1 + loop.A2 + loop.A1x + loop.A2x
        if not args.no_check_alleles:
            loop = _select_and_log(loop, _filter_alleles(alleles), log,
                                '--- {N} SNPs with valid alleles.')
            loop['Z'+str(i+2)] = _align_alleles(loop['Z'+str(i+2)], alleles)
        loop = loop.drop(['A1x', 'A2x'], axis=1)
        sumstats = loop
    _check_ld_condnum(args, log, loop[ref_ld_cnames])
    _warn_length(log, loop)
    return loop

def _get_oPGSxE_scan_table(rg_paths, oPGSxE, args):
    '''Print a table of genetic correlations.'''
    t = lambda attr: lambda obj: getattr(obj, attr, np.nan)
    x = pd.DataFrame()
    x['gxe_sumstats'] = [rg_paths[0] for i in range(1, len(rg_paths))]
    x['gwas_sumstats'] = rg_paths[1:len(rg_paths)]
    x['oPGSxE'] = list(map(t('rg_ratio'), oPGSxE))
    x['oPGSxE_se'] = list(map(t('rg_se'), oPGSxE))
    x['oPGSxE_p'] = list(map(t('p'), oPGSxE))
    if args.e_sumstats is not None:
        x['oPGSxE'] = list(map(t('rg_gxe_ratio'), oPGSxE))
        x['oPGSxE_se'] = list(map(t('rg_gxe_se'), oPGSxE))
        x['oPGSxE_p'] = list(map(t('p_gxe'), oPGSxE))
    x['covar_gxe'] = list(map(t('tot'), list(map(t('gencov'), oPGSxE))))
    x['covar_gxe_se'] = list(map(t('tot_se'), list(map(t('gencov'), oPGSxE))))
    x['covar_gxe_p'] =list(map(t('p'), list(map(t('gencov'), oPGSxE))))
    if args.samp_prev is not None and \
            args.pop_prev is not None and \
            all((i is not None for i in args.samp_prev)) and \
            all((i is not None for it in args.pop_prev)):

        c = list(map(lambda x, y: reg.h2_obs_to_liab(1, x, y), args.samp_prev[1:], args.pop_prev[1:]))
        x['h2_gxe_liab'] = list(map(lambda x, y: x * y, c, list(map(t('tot'), list(map(t('hsq1'), oPGSxE))))))
        x['h2_gxe_liab_se'] = list(map(lambda x, y: x * y, c, list(map(t('tot_se'), list(map(t('hsq1'), oPGSxE))))))
        x['h2_gxe_liab_p'] = list(map(lambda x, y: x * y, c, list(map(t('p'), list(map(t('hsq1'), oPGSxE))))))
        x['h2_add_liab'] = list(map(lambda x, y: x * y, c, list(map(t('tot'), list(map(t('hsq2'), oPGSxE))))))
        x['h2_add_liab_se'] = list(map(lambda x, y: x * y, c, list(map(t('tot_se'), list(map(t('hsq2'), oPGSxE))))))
        x['h2_add_liab_p'] = list(map(lambda x, y: x * y, c, list(map(t('p'), list(map(t('hsq2'), oPGSxE))))))
    else:
        x['h2_gxe_obs'] = list(map(t('tot'), list(map(t('hsq1'), oPGSxE))))
        x['h2_gxe_obs_se'] = list(map(t('tot_se'), list(map(t('hsq1'), oPGSxE))))
        x['h2_gxe_obs_p'] = list(map(t('p'), list(map(t('hsq1'), oPGSxE))))
        x['h2_add_obs'] = list(map(t('tot'), list(map(t('hsq2'), oPGSxE))))
        x['h2_add_obs_se'] = list(map(t('tot_se'), list(map(t('hsq2'), oPGSxE))))
        x['h2_add_obs_p'] = list(map(t('p'), list(map(t('hsq2'), oPGSxE))))
    #     x['h2_E_obs'] = list(map(t('tot'), list(map(t('hsq3'), oPGSxE))))
    #     x['h2_E_obs_se'] = list(map(t('tot_se'), list(map(t('hsq3'), oPGSxE))))
    #     x['h2_E_obs_p'] = list(map(t('p'), list(map(t('hsq3'), oPGSxE))))
    # x['covar_gxe_int'] = list(map(t('intercept'), list(map(t('gencov'), oPGSxE))))
    # x['covar_gxe_int_se'] = list(map(t('intercept_se'), list(map(t('gencov'), oPGSxE))))
    # x['covar_gxe_int_p'] = list(map(t('intercept_p'), list(map(t('gencov'), oPGSxE))))
    # x['covar_gxe_e'] = list(map(t('covar_gxe_e'), oPGSxE))
    # x['covar_g_e'] = list(map(t('tot'), list(map(t('gencov_g_e'), oPGSxE))))
    # x['covar_g_e_se'] = list(map(t('tot_se'), list(map(t('gencov_g_e'), oPGSxE))))
    # x['covar_g_e_p'] = list(map(t('p'), list(map(t('gencov_g_e'), oPGSxE))))

    # reorder x by oPGSxE_p
    x = x.sort_values(by=['oPGSxE_p'], na_position='last')

    # write out the file
    outfile = args.out + '_Scan.txt'
    tfile = open(outfile, 'w')
    tfile.write(x.to_string(header=True, index=False))
    tfile.close()

    return x.to_string(header=True, index=False) + '\n'

def _get_oPGSxe_scan_table(rg_paths, e_paths, oPGSxe, args):
    '''Print a table of genetic correlations.'''
    t = lambda attr: lambda obj: getattr(obj, attr, np.nan)
    x = pd.DataFrame()
    x['gxe_sumstats'] = [rg_paths[0]] * (len(rg_paths) - 1)
    x['e_sumstats'] = [e_paths] * (len(rg_paths) - 1)
    x['gwas_sumstats'] = rg_paths[1:len(rg_paths)]

    x['oPGSxE'] = np.array(list(map(t('rg_ratio'), oPGSxe))).flatten()
    x['oPGSxE_se'] = np.array(list(map(t('rg_se'), oPGSxe))).flatten()
    x['oPGSxE_p'] = np.array(list(map(t('p'), oPGSxe))).flatten()
    x['covar_gxe'] = np.array(list(map(t('covar_gxe'), oPGSxe))).flatten()
    x['covar_gxe_se'] = np.array(list(map(t('covar_se'), oPGSxe))).flatten()
    x['covar_gxe_p'] = np.array(list(map(t('covar_p'), oPGSxe))).flatten()
    if args.samp_prev is not None and \
            args.pop_prev is not None and \
            all((i is not None for i in args.samp_prev)) and \
            all((i is not None for it in args.pop_prev)):

        c = np.array(list(map(lambda x, y: reg.h2_obs_to_liab(1, x, y), args.samp_prev[1:], args.pop_prev[1:]))).flatten()
        x['h2_gxe_liab'] = np.array(list(map(lambda x, y: x * y, c, list(map(t('h2_gxe'), oPGSxe))))).flatten().tolist() * (len(rg_paths) - 1)
        x['h2_gxe_liab_se'] = np.array(list(map(lambda x, y: x * y, c, list(map(t('h2_se'), oPGSxe))))).flatten().tolist() * (len(rg_paths) - 1)
        x['h2_gxe_liab_p'] = np.array(list(map(lambda x, y: x * y, c, list(map(t('h2_p'), oPGSxe))))).flatten().tolist() * (len(rg_paths) - 1)
        x['h2_add_liab'] = np.array(list(map(lambda x, y: x * y, c, list(map(t('tot'), list(map(t('hsq_add'), oPGSxe))))))).flatten()
        x['h2_add_liab_se'] = np.array(list(map(lambda x, y: x * y, c, list(map(t('tot_se'), list(map(t('hsq_add'), oPGSxe))))))).flatten()
        x['h2_add_liab_p'] = np.array(list(map(lambda x, y: x * y, c, list(map(t('p'), list(map(t('hsq_add'), oPGSxe))))))).flatten()
    else:
        x['h2_gxe_obs'] = np.array(list(map(t('h2_gxe'), oPGSxe))).flatten().tolist() * (len(rg_paths) - 1)
        x['h2_gxe_obs_se'] = np.array(list(map(t('h2_se'), oPGSxe))).flatten().tolist() * (len(rg_paths) - 1)
        x['h2_gxe_obs_p'] = np.array(list(map(t('h2_p'), oPGSxe))).flatten().tolist() * (len(rg_paths) - 1)
        x['h2_add_obs'] = list(map(t('tot'), np.array(list(map(t('hsq_add'), oPGSxe))).flatten()))
        x['h2_add_obs_se'] = list(map(t('tot_se'), np.array(list(map(t('hsq_add'), oPGSxe))).flatten()))
        x['h2_add_obs_p'] = list(map(t('p'), np.array(list(map(t('hsq_add'), oPGSxe))).flatten()))
    #     x['h2_E_obs'] = np.array(list(map(t('tot'), list(map(t('hsq_E'), oPGSxe))))).flatten().tolist() * (len(rg_paths) - 1)
    #     x['h2_E_obs_se'] = np.array(list(map(t('tot_se'), list(map(t('hsq_E'), oPGSxe))))).flatten().tolist() * (len(rg_paths) - 1)
    #     x['h2_E_obs_p'] = np.array(list(map(t('p'), list(map(t('hsq_E'), oPGSxe))))).flatten().tolist() * (len(rg_paths) - 1)
    # x['covar_gxe_int'] = list(map(t('intercept'), np.array(list(map(t('gencov'), oPGSxe))).flatten()))
    # x['covar_gxe_int_se'] = list(map(t('intercept_se'), np.array(list(map(t('gencov'), oPGSxe))).flatten()))
    # x['covar_gxe_int_p'] = list(map(t('intercept_p'), np.array(list(map(t('gencov'), oPGSxe))).flatten()))
    # x['covar_gxe_e'] = np.array(list(map(t('covar_gxe_e'), oPGSxe))).flatten()
    # x['covar_g_e'] = list(map(t('tot'), np.array(list(map(t('gencov_g_e'), oPGSxe))).flatten()))
    # x['covar_g_e_se'] = list(map(t('tot_se'), np.array(list(map(t('gencov_g_e'), oPGSxe))).flatten()))
    # x['covar_g_e_p'] = list(map(t('p'), np.array(list(map(t('gencov_g_e'), oPGSxe))).flatten()))

    # reorder x by oPGSxE_p
    x = x.sort_values(by=['oPGSxE_p'], na_position='last')

    # write out the file
    outfile = args.out + '_Scan-rGE.txt'
    tfile = open(outfile, 'w')
    tfile.write(x.to_string(header=True, index=False))
    tfile.close()

    return x.to_string(header=True, index=False) + '\n'

def _get_oPGSxE_cond_table(rg_paths, oPGSxE, args):
    '''Print a table of genetic correlations.'''
    t = lambda attr: lambda obj: getattr(obj, attr, np.nan)
    x = pd.DataFrame()
    x['gxe_sumstats'] = [rg_paths[0] for i in range(1, len(rg_paths))]
    x['gwas_sumstats'] = rg_paths[1:len(rg_paths)]
    x['oPGSxE'] = list(list(map(t('coeff_delta'), oPGSxE)))[0]
    x['oPGSxE_se'] = list((list(map(t('coeff_delta_se'), oPGSxE))))[0]
    x['oPGSxE_p'] = list((list(map(t('coeff_delta_p'), oPGSxE))))[0]
    # write out the file
    outfile = args.out + '_Cond.txt'
    tfile = open(outfile, 'w')
    tfile.write(x.to_string(header=True, index=False))
    tfile.close()

    return x.to_string(header=True, index=False) + '\n'

def _get_GxE_var_table(hsq, args):
    '''Print a table of genetic correlations.'''
    t = lambda attr: lambda obj: getattr(obj, attr, np.nan)
    x = pd.DataFrame()
    x['gxe_sumstats'] = [args.gxe_sumstats]

    if args.samp_prev is not None and \
            args.pop_prev is not None and \
            all((i is not None for i in args.samp_prev)) and \
            all((i is not None for it in args.pop_prev)):

        c = list(map(lambda x, y: reg.h2_obs_to_liab(1, x, y), args.samp_prev[1:], args.pop_prev[1:]))
        x['h2_gxe_liab'] = list(map(lambda x, y: x * y, c, [hsq.tot]))
        x['h2_gxe_liab_se'] = list(map(lambda x, y: x * y, c, [hsq.tot_se]))
        x['h2_gxe_liab_p'] = list(map(lambda x, y: x * y, c, [hsq.p]))
        x['h2_gxe_liab_int'] = list(map(lambda x, y: x * y, c, [hsq.intercept]))
        x['h2_gxe_liab_int_se'] = list(map(lambda x, y: x * y, c, [hsq.intercept_se]))
        x['h2_gxe_liab_int_p'] = list(map(lambda x, y: x * y, c, [hsq.intercept_p]))
        x = x.sort_values(by=['h2_gxe_liab_p'], na_position='last')
    else:
        x['h2_gxe_obs'] = [hsq.tot]
        x['h2_gxe_obs_se'] = [hsq.tot_se]
        x['h2_gxe_obs_p'] = [hsq.p]
        # x['h2_gxe_obs_int'] = [hsq.intercept]
        # x['h2_gxe_obs_int_se'] = [hsq.intercept_se]
        # x['h2_gxe_obs_int_p'] = [hsq.intercept_p]
        x = x.sort_values(by=['h2_gxe_obs_p'], na_position='last')

    # write out the file
    outfile = args.out + '_GxE_Var.txt'
    tfile = open(outfile, 'w')
    tfile.write(x.to_string(header=True, index=False))
    tfile.close()

    return x.to_string(header=True, index=False) + '\n'

def _get_GxE_var_e_table(hsq, args):
    '''Print a table of genetic correlations.'''
    t = lambda attr: lambda obj: getattr(obj, attr, np.nan)
    x = pd.DataFrame()
    x['gxe_sumstats'] = [args.gxe_sumstats]
    x['e_sumstats'] = [args.e_sumstats]

    if args.samp_prev is not None and \
            args.pop_prev is not None and \
            all((i is not None for i in args.samp_prev)) and \
            all((i is not None for it in args.pop_prev)):

        c = list(map(lambda x, y: reg.h2_obs_to_liab(1, x, y), args.samp_prev[1:], args.pop_prev[1:]))
        x['h2_gxe_liab'] = list(map(lambda x, y: x * y, c, [hsq.h2_gxe]))
        x['h2_gxe_liab_se'] = list(map(lambda x, y: x * y, c, [hsq.h2_se]))
        x['h2_gxe_liab_p'] = list(map(lambda x, y: x * y, c, [hsq.h2_p]))
        # x['h2_gxe_liab_int'] = list(map(lambda x, y: x * y, c, [hsq.h2_intercept]))
        # x['h2_gxe_liab_int_se'] = list(map(lambda x, y: x * y, c, [hsq.h2_intercept_se]))
        # x['h2_gxe_liab_int_p'] = list(map(lambda x, y: x * y, c, [hsq.h2_intercept_p]))
        x = x.sort_values(by=['h2_gxe_liab_p'], na_position='last')
    else:
        x['h2_gxe_obs'] = [hsq.h2_gxe]
        x['h2_gxe_obs_se'] = [hsq.h2_se]
        x['h2_gxe_obs_p'] = [hsq.h2_p]
        # x['h2_gxe_obs_int'] = [hsq.h2_intercept]
        # x['h2_gxe_obs_int_se'] = [hsq.h2_intercept_se]
        # x['h2_gxe_obs_int_p'] = [hsq.h2_intercept_p]
        # x['h2_covar_gxe_e'] = [hsq.covar_gxe_e]
        x = x.sort_values(by=['h2_gxe_obs_p'], na_position='last')

    # write out the file
    outfile = args.out + '_GxE_Var-rGE.txt'
    tfile = open(outfile, 'w')
    tfile.write(x.to_string(header=True, index=False))
    tfile.close()

    return x.to_string(header=True, index=False) + '\n'

## need to change for more appropriate output
def _print_gencor(args, log, rghat, ref_ld_cnames, i, rg_paths, print_hsq1):
    l = lambda x: x + ''.join(['-' for i in range(len(x.replace('\n', '')))])
    P = [args.samp_prev[0], args.samp_prev[i + 1]]
    K = [args.pop_prev[0], args.pop_prev[i + 1]]
    if args.samp_prev is None and args.pop_prev is None:
        args.samp_prev = [None, None]
        args.pop_prev = [None, None]
    if print_hsq1:
        log.log(l('\n*** GxE Variance Component ***\n'))
        log.log(rghat.hsq1.summary(ref_ld_cnames, P=P[0], K=K[0]))

    log.log(
        l('\n*** Additive Heritability of phenotype {I}/{N} ***\n '.format(I=i + 1, N=len(rg_paths)-1)))
    log.log(rghat.hsq2.summary(ref_ld_cnames, P=P[1], K=K[1]))
    log.log(l('\n*** Covariant GxE *** \n'))
    log.log(rghat.gencov.summary(ref_ld_cnames, P=P, K=K))
    log.log(l('\n*** Oracle PGSxE *** \n'))
    log.log(rghat.summary(ref_ld_cnames) + '\n')

def _print_hsq(args, log, hsq, ref_ld_cnames):
    l = lambda x: x + ''.join(['-' for i in range(len(x.replace('\n', '')))])
    log.log(l('\n*** GxE Variance Component ***\n'))
    log.log(hsq.summary(ref_ld_cnames, P=args.samp_prev, K=args.pop_prev) + '\n')

def _merge_sumstats_sumstats(args, sumstats1, sumstats2, log):
    '''Merge two sets of summary statistics.'''
    sumstats1.rename(columns={'N': 'N1', 'Z': 'Z1'}, inplace=True)
    sumstats2.rename(
        columns={'A1': 'A1x', 'A2': 'A2x', 'N': 'N2', 'Z': 'Z2'}, inplace=True)
    x = _merge_and_log(sumstats1, sumstats2, 'summary statistics', log)
    return x

def _merge_new_sumstats(args, sumstats1, sumstats2, log, i):
    '''Merge two sets of summary statistics.'''
    if 'N' in sumstats1.columns and 'Z' in sumstats1.columns:
        sumstats1.rename(columns={'N': 'N1', 'Z': 'Z1'}, inplace=True)
    sumstats2.rename(
        columns={'A1': 'A1x', 'A2': 'A2x', 'N': 'N'+str(i+2), 'Z': 'Z'+str(i+2)}, inplace=True)
    x = _merge_and_log(sumstats1, sumstats2, 'summary statistics', log)
    return x

def _filter_alleles(alleles):
    '''Remove bad variants (mismatched alleles, non-SNPs, strand ambiguous).'''
    ii = alleles.apply(lambda y: y in MATCH_ALLELES)
    return ii

def _align_alleles(z, alleles):
    '''Align Z1 and Z2 to same choice of ref allele (allowing for strand flip).'''
    try:
        z *= (-1) ** alleles.apply(lambda y: FLIP_ALLELES[y])
    except KeyError as e:
        msg = 'Incompatible alleles in .sumstats files: %s. ' % e.args
        msg += 'Did you forget to use --merge-alleles with munge_sumstats.py?'
        raise KeyError(msg)
    return z

def _oPGSxE_scan(sumstats, args, log, M_annot, ref_ld_cnames, w_ld_cname, i):
    '''Run the regressions.'''
    n_snp = len(sumstats)
    def s(x): return np.array(x).reshape((n_snp, 1))
    if args.chisq_max is not None:
        ii = sumstats.Z1**2*sumstats.Z2**2 < args.chisq_max**2
        n_snp = np.sum(ii)  # lambdas are late binding, so this works
        sumstats = sumstats[ii]
    n_blocks = min(args.n_blocks, n_snp)
    ref_ld = sumstats[ref_ld_cnames].to_numpy()
    intercepts = [args.intercept_h2[0], args.intercept_h2[
        i + 1], args.intercept_gencov[i + 1]]
    rghat = reg.oPGSxE_scan(s(sumstats.Z1), s(sumstats.Z2),
                   ref_ld, s(sumstats[w_ld_cname]), s(sumstats.N1), s(sumstats.N2), M_annot,
                   intercept_hsq1=intercepts[0], intercept_hsq2=intercepts[1],
                   intercept_gencov=intercepts[2], n_blocks=n_blocks, twostep=args.two_step, e_sd=args.e_sd, covariates_R2=args.covariates_R2)

    return rghat

def _oPGSxe_scan(sumstats, args, log, M_annot, ref_ld_cnames, w_ld_cname, i):
    '''Run the regressions.'''
    n_snp = len(sumstats)
    def s(x): return np.array(x).reshape((n_snp, 1))
    Z_list = [s(sumstats[col]) for col in sumstats.columns if re.match(r'^Z[0-9]+$', col)]
    N_list = [s(sumstats[col]) for col in sumstats.columns if re.match(r'^N[0-9]+$', col)] # length of Z_list and N_list should be equal
    if args.chisq_max is not None:
        ii = sumstats.Z1**2*sumstats.Z2**2 < args.chisq_max**2 ## may need change
        n_snp = np.sum(ii)  # lambdas are late binding, so this works
        sumstats = sumstats[ii]
    n_blocks = min(args.n_blocks, n_snp)
    ref_ld = sumstats[ref_ld_cnames].to_numpy()
    intercepts_hsq = args.intercept_h2 ## need accordingly number of intercepts, or using None
    intercepts_gencov = args.intercept_gencov ## need accordingly number of intercepts, or using None
    ## we will have 1+2+...+len(Z_list) covs

    rghat = reg.oPGSxe_scan(Z_list, ref_ld, s(sumstats[w_ld_cname]), N_list, M_annot, intercepts_hsq=intercepts_hsq, intercepts_gencov=intercepts_gencov, n_blocks=n_blocks, twostep=args.two_step, e_sd=args.e_sd, covariates_R2=args.covariates_R2)

    return rghat

def _GxE_var(sumstats, args, log, M_annot, ref_ld_cnames, w_ld_cname):
    '''Run the regressions.'''
    n_snp = len(sumstats)
    def s(x): return np.array(x).reshape((n_snp, 1))
    n_blocks = min(args.n_blocks, n_snp)
    ref_ld = sumstats[ref_ld_cnames].to_numpy()
    rghat = reg.GxE_var(np.square(s(sumstats.Z)), ref_ld, s(sumstats[w_ld_cname]), s(sumstats.N), M_annot,
                   intercept=args.intercept_h2, n_blocks=n_blocks, twostep=args.two_step, covariates_R2=args.covariates_R2)

    return rghat

def _GxE_var_e(sumstats, args, log, M_annot, ref_ld_cnames, w_ld_cname):
    '''Run the regressions.'''
    n_snp = len(sumstats)
    def s(x): return np.array(x).reshape((n_snp, 1))
    n_blocks = min(args.n_blocks, n_snp)
    ref_ld = sumstats[ref_ld_cnames].to_numpy()
    intercepts_hsq = args.intercept_h2
    intercepts_gencov = args.intercept_gencov

    rghat = reg.GxE_var_e(s(sumstats.Z1), s(sumstats.Z2),
                   ref_ld, s(sumstats[w_ld_cname]), s(sumstats.N1), s(sumstats.N2), M_annot,
                   intercept_hsq1=intercepts_hsq[0], intercept_hsq2=intercepts_hsq[1],
                   intercept_gencov=intercepts_gencov[0], n_blocks=n_blocks, twostep=args.two_step, covariates_R2=args.covariates_R2)

    return rghat

def _oPGSxE_cond(sumstats, args, log, M_annot, ref_ld_cnames, w_ld_cname, i):
    '''Run the regressions.'''
    n_snp = len(sumstats)
    def s(x): return np.array(x).reshape((n_snp, 1))
    Z_list = [s(sumstats[col]) for col in sumstats.columns if re.match(r'^Z[0-9]+$', col)]
    N_list = [s(sumstats[col]) for col in sumstats.columns if re.match(r'^N[0-9]+$', col)] # length of Z_list and N_list should be equal
    if args.chisq_max is not None:
        ii = sumstats.Z1**2*sumstats.Z2**2 < args.chisq_max**2 ##  need change
        n_snp = np.sum(ii)  # lambdas are late binding, so this works
        sumstats = sumstats[ii]
    n_blocks = min(args.n_blocks, n_snp)
    ref_ld = sumstats[ref_ld_cnames].to_numpy()
    intercepts_hsq = args.intercept_h2 ## need accordingly number of intercepts, or using None
    intercepts_gencov = args.intercept_gencov ## need accordingly number of intercepts, or using None
    ## we will have 1+2+...+len(Z_list) covs

    rghat = reg.oPGSxE_cond(Z_list, ref_ld, s(sumstats[w_ld_cname]), N_list, M_annot, intercepts_hsq=intercepts_hsq, intercepts_gencov=intercepts_gencov, n_blocks=n_blocks, twostep=args.two_step, e_sd=args.e_sd, covariates_R2=args.covariates_R2)

    return rghat

def _parse_rg(gxe_sumstats, gwas_sumstats):
    '''Parse args.oPGSxE.'''
    rg_paths = _splitp(gxe_sumstats) + _splitp(gwas_sumstats)
    rg_files = [x.split('/')[-1] for x in rg_paths]
    if len(rg_paths) < 2:
        raise ValueError(
            'Must specify at least two phenotypes for oPGSxE estimation.')

    return rg_paths, rg_files

def _parse_e_sumstats(e_sumstats):
    '''Parse args.oPGSxE.'''
    e_paths = e_sumstats

    return e_paths

def _print_rg_delete_values(oPGSxE, fh, log):
    '''Print block jackknife delete values.'''
    _print_delete_values(oPGSxE.hsq1, fh + '.hsq1.delete', log)
    _print_delete_values(oPGSxE.hsq2, fh + '.hsq2.delete', log)
    _print_delete_values(oPGSxE.gencov, fh + '.gencov.delete', log)

def _print_rg_cov(rghat, fh, log):
    '''Print covariance matrix of estimates.'''
    _print_cov(rghat.hsq1, fh + '.hsq1.cov', log)
    _print_cov(rghat.hsq2, fh + '.hsq2.cov', log)
    _print_cov(rghat.gencov, fh + '.gencov.cov', log)

def _print_rg_delete_values_mul(oPGSxE, fh, log):
    '''Print block jackknife delete values.'''
    for i, hsq in enumerate(oPGSxE.hsq_list):
        _print_delete_values(hsq, fh + '.hsq' + str(i+1) +'.delete', log)
    for i, gencov in enumerate(oPGSxE.gencov_list):
        _print_delete_values(gencov, fh + '.gencov (' + str(oPGSxE.convorder[i][0]+1) + ', ' + str(oPGSxE.convorder[i][1]+1) + ').delete', log)

def _print_rg_cov_mul(oPGSxE, fh, log):
    '''Print covariance matrix of estimates.'''
    for i, hsq in enumerate(oPGSxE.hsq_list):
        _print_cov(hsq, fh + '.hsq' + str(i+1) +'.cov', log)
    for i, gencov in enumerate(oPGSxE.gencov_list):
        _print_cov(gencov, fh + '.gencov_' + str(oPGSxE.convorder[i][0]+1) + '_' + str(oPGSxE.convorder[i][1]+1) + '.cov', log)

def _split_or_none(x, n):
    if x is not None:
        y = map(float, x.replace('N', '-').split(','))
    else:
        y = [None for _ in range(n)]
    return y

def _check_arg_len(x, n):
    x, m = x
    if len(x) != n:
        raise ValueError(
            '{M} must have the same number of arguments as --oPGSxE/--h2.'.format(M=m))
