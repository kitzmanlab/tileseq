import numpy as np
import pandas as pd
import matplotlib.pylab as plt
from collections import defaultdict
import argparse

import seaborn as sns
import altair as alt


def varagg_ntscore_2samp(
    vtbl_sel,
    vtbl_base,
    aarng,
    pcount=0.5,
    lcol_use=['singlemut_reads','singlemut_allowsyn_reads'],
):
    """ 
    aggregate variant count tables and prepare a simple two-sample nucleotide log2 ratio score (ntLOF score)
    
    Args:
        vtbl_sel (_type_): _description_
        vtbl_base (_type_): _description_
        aarng (_type_): _description_
        pcount (float, optional): _description_. Defaults to 0.5.
        lcol_use (list, optional): _description_. Defaults to ['singlemut_reads','singlemut_allowsyn_reads'].

    Returns:
        _type_: _description_
    """
    vsel = vtbl_sel.loc[ (vtbl_sel['aa_num']>=aarng[0]) & (vtbl_sel['aa_num']<=aarng[1]) ].copy()
    vbase = vtbl_base.loc[ (vtbl_base['aa_num']>=aarng[0]) & (vtbl_base['aa_num']<=aarng[1]) ].copy()

    vj = pd.merge(
        vsel.set_index( ['aa_num','codon_ref','codon_mut','aa_ref','aa_mut','class','ed_dist'] )[ lcol_use ],
        vbase.set_index( ['aa_num','codon_ref','codon_mut','aa_ref','aa_mut','class','ed_dist'] )[ lcol_use ],
        how='outer',
        left_index=True,
        right_index=True,
        suffixes=['_sel','_base']
    )
    
    vj = vj.fillna(0)
    vj = vj.reset_index()
    
    # suppress any mut codons w/ Ns (bug from tileesq module)
    vj = vj.loc[
        ~(vj.codon_mut.str.contains('N'))
    ].copy()

    vj['cts_sel'] = 0
    for cn in lcol_use:
        vj['cts_sel'] += vj['%s_sel'%cn] 
        
    vj['cts_base'] = 0
    for cn in lcol_use:
        vj['cts_base'] += vj['%s_base'%cn] 
    
        
    vj['cts_sel'] += pcount
    vj['cts_base'] += pcount
    
    vj['cpm_sel'] = 1e6 * vj['cts_sel']/vj['cts_sel'].sum()
    vj['cpm_base'] = 1e6 * vj['cts_base']/vj['cts_base'].sum()    
        
    vj['ntscore'] = np.log2( vj['cpm_base'] / vj['cpm_base'] )
        
    return vj

# sum(readcounts for mut codons condition 1)/sum(readcounts mut codons condition 2)
def calc_aascore_2samp_simpleagg(
    ntscores,
    mincpm_base = 100,
    cpm_sel = 'cpm_sel',
    cpm_baseline = 'cpm_base',
):
    """ very simplistic aa-level LOF score, summing across the counts for each codon that encodes that aa, and then taking log2ratio.  no filters whatsoever, just two-sample comparison.  to be used only for initial quick-n-dirty QC purposes. 

    Args:
        ntscores (_type_): _description_
        mincpm_2 (int, optional): _description_. Defaults to 100.

    Returns:
        _type_: _description_
    """

    aasc={ k:[] for k in ['aa_num','aa_ref','aa_mut','varclass','nv_e1','nv_e2','nv_e3',cpm_sel,cpm_baseline] }
    
    for aa, subtbl in ntscores.loc[
            ntscores[cpm_baseline]>=mincpm_base
        ].groupby(
            ['aa_num','aa_ref','aa_mut']
        ):
    
        aasc['aa_num'].append(aa[0])
        aasc['aa_ref'].append(aa[1])
        aasc['aa_mut'].append(aa[2])
        aasc['varclass'].append( 'NON' if aa[2]=='*' else 'SYN' if aa[2]==aa[1] else 'MIS' )
        aasc['nv_e1'].append( (subtbl['ed_dist'] == 1).sum() )
        aasc['nv_e2'].append( (subtbl['ed_dist'] == 2).sum() )
        aasc['nv_e3'].append( (subtbl['ed_dist'] == 3).sum() )
        
        aasc[cpm_sel].append( subtbl[cpm_sel].sum() )
        aasc[cpm_baseline].append( subtbl[cpm_baseline].sum() )
                
    aasc = pd.DataFrame(aasc)
    
    aasc['nv']=aasc['nv_e1']+aasc['nv_e2']+aasc['nv_e3']
    
    aasc['aascore']=np.log2( aasc[cpm_sel] / aasc[cpm_baseline] )
    
    return aasc


def varagg_ntscore_Nsamp(
    vartbl_rpt,
    lsampnames,
    lvtbl,
    loutsampnames,
    aarng,
    lcol_use_cts=['singlemut_reads','singlemut_allowsyn_reads'],
    lcol_use_summary=['summary_wt','summary_syn_only','summary_stopgain','summary_singlemis_only','summary_singlemis_inclsyn']  ,
    incl_raw_counts = False,
):
    """ 
    aggregate variant count tables and prepare a simple two-sample nucleotide log2 ratio score (ntLOF score)

    for each variant, take the sum of counts in lcol_use_cts columns 
    for each library, obtain freq by div by sum of counts in lcol_use_sumary col

    varcall_rpt is used to get count denominator

    """
    
    if vartbl_rpt.index.name != 'libname':
        vartbl_rpt = vartbl_rpt.set_index('libname')

    bysamp_ttl = {}
    for (sn,osn) in zip(lsampnames,loutsampnames):
        ttls = vartbl_rpt.loc[ sn, lcol_use_summary ]
        ttlsum = ttls.sum()
        bysamp_ttl[ osn ] = ttlsum

    lv_inrng = [ vt.loc[ (vt['aa_num']>=aarng[0]) & (vt['aa_num']<=aarng[1]) ].copy() 
                for vt in lvtbl ]

    lv_inrng = [ vt.set_index( ['aa_num','codon_ref','codon_mut','aa_ref','aa_mut','class','ed_dist'] )[ lcol_use_cts ]
                for vt in lv_inrng ]

    lv_inrng_sums = [ vt.sum(1) for vt in  lv_inrng ]

    vj = pd.concat( lv_inrng_sums, axis=1 )

    vj.columns = loutsampnames
    
    vj = vj.fillna(0)
    vj = vj.reset_index()
    
    # # suppress any mut codons w/ Ns (bug from tileesq module)
    vj = vj.loc[
        ~(vj.codon_mut.str.contains('N'))
    ].copy()

    for sn in loutsampnames:
        vj[f'cpm_{sn}'] = 1e6 * vj[sn] / bysamp_ttl[sn]

    if incl_raw_counts:    
        for sn in loutsampnames:
            vj[f'cts_{sn}'] = vj[sn]

    for sn in loutsampnames:
        del vj[sn]

    return vj
   
def apply_wt_correction(
    vt,
    lcolcpm_samp,
    colcpm_plas,
):
    # vt = nucleotide level tbl
    # returns a copy
    vtret=vt.copy()

    for colcpm_samp in lcolcpm_samp:
        X = np.array(vt[ colcpm_samp ] - vt[ colcpm_plas ])
        X = X.clip(0,1e9)
        vtret[colcpm_samp] = X
        
    return vtret


def mis_non_syn_plots( 
    aasc_tbl,
    pairname,
    sampnames,
    col_score='aascore',
    col_varclass='varclass',
    log2_score_clip = (-4,4),
 ):

    aasc_tbl = aasc_tbl.copy()
    aasc_tbl['aascore']=aasc_tbl['aascore'].clip(*log2_score_clip)
    
    f,ax=plt.subplots(2,1,figsize=(5,5),sharex=True,sharey=False)
    plt.suptitle(f'aa score distributions by var class\n{pairname}\n{sampnames[0]} / {sampnames[1]}' )

    plt.sca(ax[0])
    sns.ecdfplot(
        aasc_tbl,
        hue='varclass', x='aascore'
    )
    plt.ylabel('cumul frac ≥ x')
    plt.xlabel('aa score')
    plt.xlim(*log2_score_clip)

    plt.sca(ax[1])
    sns.histplot(
        aasc_tbl,
        common_norm=False,
        ec=None,
        stat='density',
        hue='varclass', x='aascore',
        bins=np.linspace(log2_score_clip[0],log2_score_clip[1],51)
    )
    plt.xlim(*log2_score_clip)
    
    return f


def plot_aasc_hmap_link_to_splot_indiv_vars(
    sctbl,
    ntsctbl,
    title,
    sn1,
    sn2,
    colcpm_num='cpm_sel',
    colcpm_denom='cpm_base',
    colntscore='ntscore',
    colplot='aascore',
    log2rng=[-5,5],
    dispinnerbrks=[-0.5,0.5],
    aasort=list('AVLIMFYWRHKDESTNQGCP*'),
    height=250,
    width=None
):        
    """ plot aa lofscore heatmap and link to individual variants' abundances on an interactive scatterplot 

    Args:
        sctbl (_type_): _description_
        ntsctbl (_type_): _description_
        title (_type_): _description_
        sn1 (_type_): _description_
        sn2 (_type_): _description_
        colplot (str, optional): _description_. Defaults to 'aascore'.
        log2rng (list, optional): _description_. Defaults to [-5,5].
        dispinnerbrks (list, optional): _description_. Defaults to [-0.5,0.5].
        aasort (_type_, optional): _description_. Defaults to list('AVLIMFYWRHKDESTNQGCP*').
        height (int, optional): _description_. Defaults to 250.
        width (_type_, optional): _description_. Defaults to None.

    Returns:
        _type_: _description_
    """

    sctbl = sctbl.copy()
    
    sctbl[colplot]=sctbl[colplot].clip( log2rng[0], log2rng[1] )
    
    if width==None:
        aarng=(sctbl.aa_num.min(),sctbl.aa_num.max())
        width=int(1000 * ((aarng[1]-aarng[0])/100))
                
    lcol=list(sctbl.columns)

    pl2base=alt.Chart(
        sctbl,
        height=height,
        width=width,
        title=title,
    )
    
    pl2r=pl2base.mark_rect(
    ).encode(
        alt.X('aa_num:O'),
        alt.Y('aa_mut:O',sort=aasort),
        alt.Color('%s:Q'%colplot, scale=alt.Scale( 
            clamp=True, 
            domain=[log2rng[0]]+dispinnerbrks+[log2rng[-1]], 
            range=['blue','white','white','darkred'],
            interpolate='hsl')),
        alt.Tooltip(lcol)
    )    
    
    pl2b=pl2base.mark_rect(
        color='#cdcdcd'
    ).encode(
    )

    selr = alt.selection_point( fields=['aa_mut','aa_num'] )
    pl2r = pl2r.add_params( selr )

    cpmLim=(0.1,10000)
    
    line = pd.DataFrame({
        colcpm_num: cpmLim,
        colcpm_denom: cpmLim,
    })    

   
    scplot = alt.Chart(
        ntsctbl,
        title=f'{sn1} vs {sn2}'
    ).transform_filter(
        selr
    ).mark_circle(
    ).encode(
        alt.X(colcpm_denom,scale=alt.Scale( type='log', domain=cpmLim,clamp=True)),
        alt.Y(colcpm_num,scale=alt.Scale( type='log', domain=cpmLim, clamp=True)),
        alt.Color('codon_mut'),
        tooltip=['aa_num','aa_ref','codon_ref','aa_mut','codon_mut','aa_ref',colcpm_num,colcpm_denom,colntscore]
    )
      
    scoverlay = alt.Chart(line).mark_line(color= 'black').encode(
        alt.X(colcpm_num).axis(title=f'cpm {sn2}'),
        alt.Y(colcpm_num).axis(title=f'cpm {sn1}'),
    )
          
    vhist = alt.Chart(
        sctbl,
        height=50
    ).mark_bar(
    ).encode(
        alt.X('%s:Q'%colplot,bin=alt.Bin(extent=[log2rng[0],log2rng[1]],maxbins=40)), #,bin=alt.Bin(extent=[-4, 4], step=0.5)),
        alt.Y('count()', axis=alt.Axis(title=None)),
        alt.Row('varclass'),
        alt.Color('nv:O', scale=alt.Scale(scheme='set2')),
    ).resolve_scale(
        y='independent',
    )
       
    p = (
        (pl2b + pl2r) &
        (( scplot+scoverlay ) | vhist ).resolve_scale(color='independent')
    ).resolve_scale(
        color='independent'
    )
    
    return p


def main():

    def intRange(s):
        try:
            lower,upper=list(map(int, s.split(',')))
            assert lower<=upper
            return lower,upper
        except:
            raise argparse.ArgumentTypeError('range must be x,y, x<=y')

    opts = argparse.ArgumentParser( description='quick-and-dirty qc plots for individual replicate pairs' )
    
    opts.add_argument('--tbl_varcall_rpt', required=True, dest='tbl_varcall_rpt',
                      help='input file with variant calling report from tileseq pipeline, usually named varcall_rpt.txt')
    
    opts.add_argument('--varcall_path_prefix', required=True, dest='varcall_path_prefix',
                      help='path to prepend to variant counts table paths referred to by the [tbl_varcall_rpt] table')

    opts.add_argument('--tbl_rep_info', required=True, dest='tbl_rep_info',
                      help='table with list of paired samples, required columns: pair_name,libname_sel,libname_base')

    opts.add_argument('--out_base', required=True, dest='out_base',
                      help='output base path')

    opts.add_argument('--aarange', required=True, type=intRange, dest='aarange') 

    opts.add_argument('--pcount', type=float, default=0.1, dest='pcount') 
    opts.add_argument('--min_base_cpm', type=float, default=50., dest='min_base_cpm') 

    opts.add_argument('--aasc_clip_range', default='-4,4', dest='aasc_clip_range') 
    opts.add_argument('--cols_counts_include',default='singlemut_reads,singlemut_allowsyn_reads',dest='cols_counts_include')

    o = opts.parse_args()

    tbl_varcall_rpt = pd.read_table(o.tbl_varcall_rpt)
    tbl_rep_info = pd.read_table(o.tbl_rep_info)

    aasc_clip_range = ( float(o.aasc_clip_range.split(',')[0]), 
                       float(o.aasc_clip_range.split(',')[1]) )
    
    cols_counts_include = o.cols_counts_include.split(',')
    
    tbl_varcall_rpt = tbl_varcall_rpt.set_index('libname')

    for _, reprow in tbl_rep_info.iterrows():

        pairname = reprow['pair_name']

        print(f'working on {pairname}')

        samprow_base = tbl_varcall_rpt.loc[ reprow['libname_base'] ]
        samprow_sel = tbl_varcall_rpt.loc[ reprow['libname_sel'] ]

        vartbl_base = pd.read_table(o.varcall_path_prefix + '/' +samprow_base['vartbl'])
        vartbl_sel = pd.read_table(o.varcall_path_prefix + '/' +samprow_sel['vartbl'])

        ntsc_sel_base = varagg_ntscore_2samp(
            vartbl_sel,
            vartbl_base,
            aarng=o.aarange,
            pcount=o.pcount,
            lcol_use=cols_counts_include,
        )
        
        ntsc_sel_base = ntsc_sel_base.loc[ ntsc_sel_base.aa_num.between( o.aarange[0], o.aarange[1] , inclusive='both') ]

        aasc_sel_base = calc_aascore_2samp_simpleagg(
            ntsc_sel_base,
            mincpm_base=o.min_base_cpm,
        )

        fnout_ntsc = o.out_base + f'{pairname}.ntscore.tsv'
        ntsc_sel_base.to_csv( fnout_ntsc, sep='\t', index=False )

        fnout_aasc = o.out_base + f'{pairname}.aascore.tsv'
        aasc_sel_base.to_csv( fnout_aasc, sep='\t', index=False )
        
        fig_score_dists = mis_non_syn_plots( aasc_sel_base,
            pairname,
            (reprow['libname_sel'],reprow['libname_base']),
            col_score='aascore',
            col_varclass='varclass',
            log2_score_clip = aasc_clip_range,
        )  

        fnout_dist_plot = o.out_base + f'{pairname}.aascore_hists.pdf'
        fig_score_dists.savefig(fnout_dist_plot)
       

        fig_heatmap = plot_aasc_hmap_link_to_splot_indiv_vars(
            aasc_sel_base,
            ntsc_sel_base,
            f"{pairname} : {reprow['libname_sel']} / {reprow['libname_base']} ",
            reprow['libname_sel'],
            reprow['libname_base'],
            colcpm_num='cpm_sel',
            colcpm_denom='cpm_base',
            colntscore='ntscore',
            colplot='aascore',
            log2rng=aasc_clip_range,
            dispinnerbrks=[-0.5,0.5],
            aasort=list('AVLIMFYWRHKDESTNQGCP*'),
            height=250,
            width=None
        )

        fnout_heatmap = o.out_base + f'{pairname}.heatmap.html'
        fig_heatmap.save(fnout_heatmap)







if __name__ == '__main__':                
    main()

    