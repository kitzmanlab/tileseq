import sys
import argparse

import numpy as np

import pandas as pd

import seaborn.objects as so

import matplotlib.pyplot as plt

def mutratebycodon_1samp_bypos_strat_type_ed(
    title,
    cts,
    cvg, 
    min_cvg,
    codon_rng_zoom = None,
    figsize = (20,6),
    freq_clip_rng=(1e-6,0.025) ):
    
    cc=pd.merge(cts,cvg,how='inner',left_on='aa_num',right_on='codon_num')
    cc['freq_singlemut'] = (cc['singlemut_reads']/cc['total_pass_reads']).fillna(0)
    cc['freq_singlemut_allowsyn'] = (cc['singlemut_allowsyn_reads']/cc['total_pass_reads']).fillna(0)
    
    ccg=cc.groupby( ['codon_num','class','ed_dist'] ).agg(
    { 'singlemut_reads': sum,
      'singlemut_allowsyn_reads': sum,
      'total_pass_reads': min,
      'freq_singlemut': np.mean,
     'freq_singlemut_allowsyn': np.mean,
     # 'codon_ref':lambda l:list(l)[0],
    } 
    )
    ccg2=ccg.reset_index()

    ccg2fl = ccg2.copy()
    ccg2fl['freq_singlemut_fl'] = ccg2fl['freq_singlemut'].clip(*freq_clip_rng)
    ccg2fl['freq_singlemut_allowsyn_fl'] = ccg2fl['freq_singlemut_allowsyn'].clip(*freq_clip_rng)
    
    ccg2fl_mincvg = ccg2fl.loc[ 
        ccg2fl['total_pass_reads']>=min_cvg
    ]

    fig = plt.figure(constrained_layout=True, figsize=figsize )
    subfigs = fig.subfigures(1)

    p=(
        so.Plot(ccg2fl, x="codon_num")
        .label (
            x="Codon position", y="per-mut mean freq",  )
        .facet( row='ed_dist', col='class' )        
        .add(so.Dots(), y="freq_singlemut_allowsyn_fl", color="class")
        .scale(y="log")
        .layout(size=(8, 4), engine='constrained'  )
    )

    p.on(subfigs).plot()

    if codon_rng_zoom:
        p=p.limit( x=codon_rng_zoom )

    fig.suptitle("Variant freq by codon position, "+title)

    return p, fig
    # return p



def main():
    opts = argparse.ArgumentParser( description='by position mut freq plots' )
    # sp = opts.add_subparsers(dest='cmd')
    # mutratebycodon = sp.add_parser('mutratebycodon')

    opts.add_argument('--in_varcts', required=True, dest='in_varcts')
    opts.add_argument('--in_varcvg', required=True, dest='in_varcvg')
    opts.add_argument('--desc', required=True, dest='desc')
    opts.add_argument('--codon_range', default=None, dest='codon_range')
    opts.add_argument('--mincvg', type=int, default=1000, dest='mincvg')

    opts.add_argument('--out', required=True, dest='out')

    o = opts.parse_args()

    tblcts = pd.read_table(o.in_varcts)
    tblcvg = pd.read_table(o.in_varcvg)

    codon_rng_zoom = None
    if o.codon_range:
        codon_rng_zoom = (int(o.codon_range.split(',')[0]),int(o.codon_range.split(',')[1]))

    p,fig = mutratebycodon_1samp_bypos_strat_type_ed(
        o.desc,
        tblcts,
        tblcvg, 
        o.mincvg,
        codon_rng_zoom = codon_rng_zoom,
        figsize = (20,6),
        freq_clip_rng=(1e-6,0.025) )
    
    fig.savefig( o.out, facecolor='white' )
        

if __name__ == '__main__':                
    main()

    