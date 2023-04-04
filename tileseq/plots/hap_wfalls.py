import sys
import argparse
import numpy as np
import pandas as pd
import seaborn.objects as so
import matplotlib.pyplot as plt
import altair as alt
from statsmodels.distributions.empirical_distribution import ECDF

def plot_varfreq_wfall_allcodonposs_byed(
    varcounts,
    title,
    count_col='singlemut_allowsyn_reads',
    psuedocount=0.1,
    xrel=True,
    logy=False,
):
    f,ax=plt.subplots(2,1,figsize=(8,5),sharey=False,sharex=True)
    
    vt = varcounts.copy()
    
    if logy:
        vt[count_col]+=psuedocount
            
    grand_tot = vt[count_col].sum()
        
    if xrel:
        rescalex=lambda x:x/x.shape[0]
    else:
        rescalex=lambda x:x
    
    for ied,ed in enumerate(range(1,4)):
        vt_ed=vt.loc[vt['ed_dist']==ed]
        cts=np.array(vt_ed[count_col])
        ax[0].scatter( x=rescalex(np.arange(cts.shape[0])), y=cts[np.argsort(-cts)], label='%d edits'%ed, s=1 )
               
        cscts=np.cumsum( cts[np.argsort(-cts)] )
        ax[1].scatter( x=rescalex(np.arange(cts.shape[0])), y=cscts/cscts[-1], label='%d edits'%ed, s=1 )
    
    cts=np.array(vt[count_col])
    ax[0].scatter( x=rescalex(np.arange(cts.shape[0])), y=cts[np.argsort(-cts)], label='overall', s=1)
    cscts=np.cumsum( cts[np.argsort(-cts)] )
    ax[1].scatter( x=rescalex(np.arange(cts.shape[0])), y=cscts/cscts[-1], label='overall', s=1)
                
    ax[0].legend()
    
    ax[0].set_ylabel('reads')
    ax[1].set_ylabel('cumul. frac of reads')
    
    if xrel:
        ax[1].set_xlabel('per-variant fractional rank, more to less abundant')
    else:
        ax[1].set_xlabel('per-variant rank, more to less abundant')
    
    if logy:
        ax[0].set_yscale('log')
            
    plt.suptitle('{}\nOn-tile missense mutations, #reads by rank'.format(title))

    return f, ax


def plot_wfall_hapcts(
    hapcts,
    title,
    remove_wt=False,
    count_col='total_counts',
    xrel=False,
    logx=False,
    logy=False,
):
    f,ax=plt.subplots(2,1,figsize=(8,5),sharey=False,sharex=True)

    if remove_wt:
        hc = hapcts.loc[ hapcts.status!='wildtype' ]
    else:
        hc = hapcts.copy()
        
    
    grand_tot = hc[count_col].sum()
        
    if xrel:
        rescalex=lambda x:x/x.shape[0]
    else:
        rescalex=lambda x:x
        
    cts=np.array(hc[count_col])
    ax[0].scatter( x=rescalex(1+np.arange(cts.shape[0])), y=cts[np.argsort(-cts)], label='overall', s=1)
    cscts=np.cumsum( cts[np.argsort(-cts)] )
    ax[1].scatter( x=rescalex(1+np.arange(cts.shape[0])), y=cscts/cscts[-1], label='overall', s=1)

    cts=np.array(hc.loc[ hc['counted_codonwise'], count_col] )
    ax[0].scatter( x=rescalex(1+np.arange(cts.shape[0])), y=cts[np.argsort(-cts)], label='usable (1mis)', s=1)
    cscts=np.cumsum( cts[np.argsort(-cts)] )
    ax[1].scatter( x=rescalex(1+np.arange(cts.shape[0])), y=cscts/cscts[-1], label='usable (1mis)', s=1)
                
    ax[0].legend()
    
    ax[0].set_ylabel('reads')
    ax[1].set_ylabel('cumul. frac of reads')
    
    if xrel:
        ax[1].set_xlabel('per-haplotype fractional rank, more to less abundant')
    else:
        ax[1].set_xlabel('per-haplotype rank, more to less abundant')
    
    if logy:
        ax[0].set_yscale('log')
        
    if logx:
        ax[0].set_xscale('log')
        
    plt.suptitle('{}\nHaplotypes, #reads vs rank'.format(title))

    return f

def cross_samp_cumul_haps(
    samptbl,
    lqile,
    title,
    remove_wt=True,
    col_hapcts='haptbl',
    sortby=['target_name','libname'],
    count_col='total_counts',
    logy=False,
):
    # lcolqiles = ['nhaps_qile%02d'%q for q in lqile]
    # tblout = { k:[] for k in ['libname']+lcolqiles }

    lsamps_plot_order = list(samptbl.sort_values(by=sortby)['libname'])

    samptbl = samptbl.set_index('libname')

    tbl_hc_long = { k:[] for k in ['libname','hap_status','cumulthresh','hap_count'] }

    for libname,r in samptbl.iterrows():
        hapcts = pd.read_table( r[col_hapcts] )

        if remove_wt:
            hapcts = hapcts.loc[ hapcts['status']!='wildtype' ]

        cts_all = np.array(hapcts[count_col],dtype=int)
        li_usable = np.array(hapcts['counted_codonwise'],dtype=bool)
        li_sort = (-cts_all).argsort()
        cts_all,li_usable = cts_all[li_sort], li_usable[li_sort]
        
        ccts = cts_all.cumsum()
        ccts = ccts/ccts[-1]

        for qile in lqile:
            imin = np.where(ccts>qile)[0].min()
            
            n_haps_usable = li_usable[:imin].sum()
            n_haps_filtered = (~li_usable)[:imin].sum()

            tbl_hc_long['libname'].append(libname)
            tbl_hc_long['hap_status'].append('usable')
            tbl_hc_long['cumulthresh'].append(qile)
            tbl_hc_long['hap_count'].append(n_haps_usable)

            tbl_hc_long['libname'].append(libname)
            tbl_hc_long['hap_status'].append('filtered')
            tbl_hc_long['cumulthresh'].append(qile)
            tbl_hc_long['hap_count'].append(n_haps_filtered)
           
    tbl_hc_long = pd.DataFrame(tbl_hc_long)
    
    ch=alt.Chart(
        tbl_hc_long,
        title='Haplotype counts w/in cumul top % '+title,
        height=100
    ).mark_bar(
    ).encode(
        alt.X('libname',sort=list(tbl_hc_long['libname'])),
        alt.Y('hap_count'),
        alt.Color('hap_status', scale=alt.Scale(scheme='category20') ),
        alt.Row('cumulthresh'),
        alt.Tooltip( ['libname','hap_status','hap_count'] ) 
    )

    tbl_out = tbl_hc_long.pivot( index='libname', columns=['cumulthresh','hap_status'], values=['hap_count'] )

    return tbl_out, ch


def plot_countofhaps_by_rdthresh_of_poss_vars(
     vtbl, 
    haptbl,
    aarng,
    desc,
    eds_incl=[2,3],
    exclude_syn=True,
    readcount_thresh=[1,5,10],
    histo_hapcount_max=25,
    hmap_hapcount_rng=(1,50),
):
    """plot, out of all variants in a tile (subject to class & edit distance criteria), the number
    of corresponding haplotypes, at various min read count thresholds. 

    Args:
        vtbl (_type_): _description_
        haptbl (_type_): _description_
        aarng (_type_): _description_
        desc (_type_): _description_
        top_panel_cumul (bool, optional): _description_. Defaults to True.
        eds_incl (list, optional): _description_. Defaults to [2,3].
        exclude_syn (bool, optional): _description_. Defaults to True.
        readcount_thresh (list, optional): _description_. Defaults to [1,5,10].
        histo_hapcount_max (int, optional): _description_. Defaults to 25.
        hmap_hapcount_rng (tuple, optional): _description_. Defaults to (1,50).

    Returns:
        _type_: _description_
    """


    vtbl_ontgt = vtbl.loc[
        vtbl.aa_num.between(aarng[0],aarng[1],inclusive=True)
    ]
    
    vtbl_ontgt_ed = vtbl_ontgt.loc[
        vtbl_ontgt.ed_dist.isin(eds_incl)
    ]
    
    if exclude_syn:
        vtbl_ontgt_ed = vtbl_ontgt_ed.loc[ vtbl_ontgt_ed['class']!='SYN' ]
        
    
    Nvar = vtbl_ontgt_ed.shape[0]
    
    hapct_by_thresh = []
    
    pl_hcmtx = []
    
    for minrc in readcount_thresh:
        hcs = haptbl.query('counted_codonwise and total_counts>=%d'%(minrc)).groupby( ['aa_num','codon_mut'] )['hap'].agg(len)
        hcs = pd.DataFrame(  { 'Nhaps_at_thresh': hcs } )
    
        vtbl_ontgt_ed_cp = vtbl_ontgt_ed.copy()
        vtbl_ontgt_ed_cp['minrd_thresh']=minrc 
        
        vtbl_ontgt_ed_cp['aa_codon_mut'] = ['%s_%s'%( r.aa_mut,r.codon_mut) for _,r in vtbl_ontgt_ed_cp.iterrows() ]
               
        hapct_by_thresh.append( 
            pd.merge(vtbl_ontgt_ed_cp, hcs, left_on=['aa_num','codon_mut'], right_index=True,how='left')
        )
                        
        hapct_by_thresh[-1]['Nhaps_at_thresh']=hapct_by_thresh[-1]['Nhaps_at_thresh'].fillna(0).astype(int)        
                        
        mtxpl_hm = alt.Chart(
            hapct_by_thresh[-1],
            width=400,
            height=200,
            title='# haps with ≥%d reads, by position, codon'%(minrc)
        ).mark_rect(
        ).encode(
            alt.X('aa_num:O'),
            alt.Y('aa_codon_mut:O'),
            alt.Color('Nhaps_at_thresh', scale=alt.Scale( 
                clamp=True, 
                type='log',
                domain=hmap_hapcount_rng,
                scheme='viridis')),
            alt.Tooltip( ['aa_num','aa_codon_mut','minrd_thresh','Nhaps_at_thresh'])
        )
            
        mtxpl_missing_superimp = alt.Chart(
            hapct_by_thresh[-1].query('Nhaps_at_thresh==0'),
            width=400,
            height=200
        ).mark_rect(
            color='red'
        ).encode(
            alt.X('aa_num:O'),
            alt.Y('aa_codon_mut:O'),
            alt.Tooltip( ['aa_num','aa_codon_mut','minrd_thresh','Nhaps_at_thresh'])
        )
                
        pl_hcmtx.append(mtxpl_hm+mtxpl_missing_superimp)

        
    hapct_by_thresh = pd.concat( hapct_by_thresh )

    _hapct_by_thresh = hapct_by_thresh.copy()

    lout_hist_cumul = []
    # prep histogram tbl and ecdf tbl for output
    for minrd_thresh, subtbl in _hapct_by_thresh.groupby('minrd_thresh'):
    
        ec = ECDF( subtbl['Nhaps_at_thresh'], side='left' )

        vc= subtbl['Nhaps_at_thresh'].value_counts()

        cumul = 1 - ec( list(vc.index) )

        assert( np.isclose( 1-ec(2), (subtbl['Nhaps_at_thresh']>=2).mean() ) )

        oc = pd.DataFrame(
            {'minrd_thresh':minrd_thresh,
             'nvariants':list(vc),
             'nhaps_gte_thresh':list(vc.index),
             'cumul_frac_haps_gte_thresh':cumul,
            }
        )
        oc=oc.sort_values(by='nhaps_gte_thresh', ascending=True)
        
        # print(ec(0),ec(1),ec(2),(subtbl['Nhaps_at_thresh']==0).mean(),(subtbl['Nhaps_at_thresh']==1).mean(),(subtbl['Nhaps_at_thresh']==2).mean(),(subtbl['Nhaps_at_thresh']>=2).mean())

        lout_hist_cumul.append(oc)

    lout_hist_cumul=pd.concat(lout_hist_cumul)

    hapct_by_thresh.loc[
        hapct_by_thresh['Nhaps_at_thresh']<1,
        'Nhaps_at_thresh'
    ] = 0.1

    # phist=alt.Chart(
    #     hapct_by_thresh,
    #     title='# distinct haps per variant (ed%s %s vars aa %d-%d, N=%d) %s'%(
    #         ','.join([str(x) for x in eds_incl]),
    #         'mis+non+syn' if not exclude_syn else 'mis+non',
    #         aarng[0],aarng[1],vtbl_ontgt_ed.shape[0],
    #         desc ),
    #     height=150
    # ).mark_line(
    #     point=True,
    # ).encode(
    #     alt.X('Nhaps_at_thresh', scale=alt.Scale(type='log', domain=(0.1,histo_hapcount_max), clamp=True ),title='Count of distinct haplotypes'),
    #     alt.Y('count()',stack=False,title='Count of variants'),
    #     alt.Color('minrd_thresh:O'),
    # )

    # pcumul=alt.Chart(
    #     hapct_by_thresh,
    #     title='# distinct haps per variant (ed%s %s vars aa %d-%d, N=%d) %s'%(
    #         ','.join([str(x) for x in eds_incl]),
    #         'mis+non+syn' if not exclude_syn else 'mis+non',
    #         aarng[0],aarng[1],vtbl_ontgt_ed.shape[0],
    #         desc ),
    #         height=150
    # ).transform_window(
    #     ecdf="cume_dist()",
    #     groupby=['minrd_thresh'],
    #     sort=[{"field": "Nhaps_at_thresh"}],
    # ).transform_calculate(
    #     iecdf='1.0 - datum.ecdf'
    # ).mark_line(
    #     interpolate="step-after"
    # ).mark_line(
    #     point=True,
    # ).encode(
    #     alt.X('Nhaps_at_thresh', scale=alt.Scale(type='log', domain=(0.1,histo_hapcount_max), clamp=True ), title='Min count of distinct haplotypes'),
    #     alt.Y('iecdf:Q', title='Cumulative fraction of variants'),
    #     alt.Color('minrd_thresh:O'),
    # )

    toplot = lout_hist_cumul.copy()
    toplot = toplot.loc[ toplot.nhaps_gte_thresh>0 ]
    loc = list(lout_hist_cumul.columns)
    pcumul=alt.Chart(
        toplot,
        title='# distinct haps per variant (ed%s %s vars aa %d-%d, N=%d) %s'%(
            ','.join([str(x) for x in eds_incl]),
            'mis+non+syn' if not exclude_syn else 'mis+non',
            aarng[0],aarng[1],vtbl_ontgt_ed.shape[0],
            desc ),
            height=150
    ).mark_line(
        interpolate="step-after",
        point=True,
    ).encode(
        alt.X('nhaps_gte_thresh:Q', scale=alt.Scale(type='log', domain=(1,histo_hapcount_max), clamp=True ), title='Min count of distinct haplotypes'),
        alt.Y('cumul_frac_haps_gte_thresh:Q', title=['Cumulative fraction of variants','with ≥ X distinct haplotypes']),
        alt.Color('minrd_thresh:O'),
        alt.Tooltip(loc),
    )    

    pj = pcumul & alt.vconcat(*pl_hcmtx)

    return hapct_by_thresh,lout_hist_cumul,pj


def main():
    
    def intRange(s):
        try:
            lower,upper=list(map(int, s.split(',')))
            assert lower<=upper
            return lower,upper
        except:
            raise argparse.ArgumentTypeError('range must be x,y, x<=y')


    opts = argparse.ArgumentParser( description='haplotype waterfall plots' )
    
    sp = opts.add_subparsers(dest='cmd')
    
    cmd_hapwfall = sp.add_parser('hapwfall')

    cmd_hapwfall.add_argument('--in_hapcts', required=True, dest='in_hapcts')
    cmd_hapwfall.add_argument('--desc', required=True, dest='desc')
    cmd_hapwfall.add_argument('--remove_wt', default=False, action='store_true', dest='remove_wt')
    cmd_hapwfall.add_argument('--logx', default=False, action='store_true', dest='logx')
    cmd_hapwfall.add_argument('--logy', default=False, action='store_true', dest='logy')
    cmd_hapwfall.add_argument('--out', required=True, dest='out')

    cmd_varwfall = sp.add_parser('varwfall')

    cmd_varwfall.add_argument('--in_varcts', required=True, dest='in_varcts')
    cmd_varwfall.add_argument('--incl_nonsense', default=False, action='store_true', dest='incl_nonsense')
    cmd_varwfall.add_argument('--incl_syn', default=False, action='store_true', dest='incl_syn')
    cmd_varwfall.add_argument('--desc', required=True, dest='desc')
    cmd_varwfall.add_argument('--aarng', required=True, type=intRange, dest='aarng')
    # cmd_varwfall.add_argument('--xrel', default=False, action='store_true', dest='xrel')
    # cmd_varwfall.add_argument('--logy', default=False, action='store_true', dest='logy')
    cmd_varwfall.add_argument('--out', required=True, dest='out')

    cmd_comparesamps = sp.add_parser('comparesamps')

    cmd_comparesamps.add_argument('--in_samptbl', required=True, dest='in_samptbl')
    cmd_comparesamps.add_argument('--out_counttbl', required=True, dest='out_counttbl')
    cmd_comparesamps.add_argument('--include_wt', default=True, action='store_false', dest='remove_wt')
    cmd_comparesamps.add_argument('--desc', required=True, dest='desc')
    cmd_comparesamps.add_argument('--out_plot', required=True, dest='out_plot')
    cmd_comparesamps.add_argument('--col_haptbl', default='haptbl', dest='col_haptbl')
    cmd_comparesamps.add_argument('--cumul_cutoffs', default='75,90,95', dest='cumul_cutoffs')

    cmd_hapcvg = sp.add_parser('varcoverageplots')

    cmd_hapcvg.add_argument('--in_hapcts', required=True, dest='in_hapcts')
    cmd_hapcvg.add_argument('--in_varcts', required=True, dest='in_varcts')
    cmd_hapcvg.add_argument('--aa_range', required=True, type=intRange, dest='aa_range')
    cmd_hapcvg.add_argument('--desc', required=True, dest='desc')
    cmd_hapcvg.add_argument('--include_ed1_haps', default=False,action='store_true',dest='include_ed1_haps')
    cmd_hapcvg.add_argument('--out_plot', required=True, dest='out_plot')
    cmd_hapcvg.add_argument('--out_tbl', required=True, dest='out_tbl')

    o = opts.parse_args()

    if o.cmd == 'hapwfall':
        tblhapcts = pd.read_table(o.in_hapcts)

        f = plot_wfall_hapcts(
            tblhapcts,
            o.desc,
            remove_wt=o.remove_wt,
            logy=o.logy,
            logx=o.logx
        )
        
        f.savefig( o.out, facecolor='white', dpi=200)
    
    elif o.cmd == 'varwfall':
        tblvarcts = pd.read_table(o.in_varcts)

        tblvarcts = tblvarcts.loc[
            tblvarcts['aa_num'].between( o.aarng[0], o.aarng[1] )
        ]
        
        li_incl = tblvarcts['class']=='MIS'
        if o.incl_nonsense:
            li_incl |= tblvarcts['class']=='NON'
        if o.incl_syn:
            li_incl |= tblvarcts['class']=='SYN'

        tblvarcts = tblvarcts.loc[ li_incl ]

        f, ax = plot_varfreq_wfall_allcodonposs_byed(
            tblvarcts,
            o.desc,
            count_col = 'singlemut_allowsyn_reads',
            psuedocount=0.1,
            logy=True,
            xrel=True,
        )
        
        f.savefig( o.out, facecolor='white', dpi=200)
    
    elif o.cmd == 'comparesamps':
        samptbl = pd.read_table(o.in_samptbl)

        lqile = [ int(x)/100 for x in o.cumul_cutoffs.split(',') ]

        tbl_hapcts_byqile, ch = cross_samp_cumul_haps(
            samptbl,
            lqile,
            o.desc,
            remove_wt=o.remove_wt,
            col_hapcts='haptbl',
            count_col='total_counts',
            logy=False,
        )

        ch.save(o.out_plot)


        print(tbl_hapcts_byqile.columns)

        tbl_hapcts_byqile.columns = [
            '%s_%02d_%s'%(c[0],int(100*c[1]),c[2]) for c in tbl_hapcts_byqile.columns
        ]

        tbl_hapcts_byqile.to_csv( o.out_counttbl, index=True, sep='\t' )

    elif o.cmd == 'varcoverageplots':
        tblhapcts = pd.read_table(o.in_hapcts)
        tblvarcts = pd.read_table(o.in_varcts)

        _,hapcthists,p = plot_countofhaps_by_rdthresh_of_poss_vars(
            tblvarcts, 
            tblhapcts,
            o.aa_range,
            o.desc,
            eds_incl=[2,3] if not o.include_ed1_haps else [1,2,3],
            exclude_syn=True,
            readcount_thresh=[1,5,10],
            histo_hapcount_max=100,
            hmap_hapcount_rng=(1,50),
        )

        p.save(o.out_plot)

        hapcthists.to_csv(o.out_tbl, index=False, sep='\t')

    

if __name__ == '__main__':                
    main()

    