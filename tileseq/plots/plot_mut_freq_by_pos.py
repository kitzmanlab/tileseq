import sys
import argparse

import numpy as np

import pandas as pd

import seaborn.objects as so

import matplotlib.pyplot as plt

def altplot_mutrate_1samp_by_pos_strat_type_ed(
    title,
    cts,
    cvg, 
    min_cvg,
    codon_rng_zoom = None,
    figsize = (20,6),
    freq_clip_rng=(1e-6,0.025) ):
    
    cc=pd.merge(cts,cvg,how='inner',left_on='aa_num',right_on='codon_num')
    cc['freq_singlemut'] = (cc['singlemut_reads']/cc['total_pass_reads']).fillna(0)
    cc['freq_singlemut_plussyn'] = (cc['singlemut_plussyn_reads']/cc['total_pass_reads']).fillna(0)
    
    ccg=cc.groupby( ['codon_num','class','ed_dist'] ).agg(
    { 'singlemut_reads': sum,
      'singlemut_plussyn_reads': sum,
      'total_pass_reads': min,
      'freq_singlemut': np.mean,
     'freq_singlemut_plussyn': np.mean,
     # 'codon_ref':lambda l:list(l)[0],
    } 
    )
    ccg2=ccg.reset_index()

    ccg2fl = ccg2.copy()
    ccg2fl['freq_singlemut_fl'] = ccg2fl['freq_singlemut'].clip(*freq_clip_rng)
    ccg2fl['freq_singlemut_plussyn_fl'] = ccg2fl['freq_singlemut_plussyn'].clip(*freq_clip_rng)
    
    ccg2fl_mincvg = ccg2fl.loc[ 
        ccg2fl['total_pass_reads']>=min_cvg
    ]

    fig = plt.figure(constrained_layout=True, figsize=figsize )
    # subfigs = fig.subfigures(1)

    p=(
        so.Plot(ccg2fl, x="codon_num",  )
        .facet( row='ed_dist', col='class' )
        .label(
            x="Codon position", y="per-mut mean freq", )
        .add(so.Dots(), y="freq_singlemut_plussyn_fl")
        .scale(y="log")
        .layout(size=(20, 6), )
    )

    if codon_rng_zoom:
        p=p.limit( x=codon_rng_zoom )

    fig.suptitle(title)

    return p, fig


if __name__ == '__main__':                
    main()

    