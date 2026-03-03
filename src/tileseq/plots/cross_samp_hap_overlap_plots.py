import sys
import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def join_haptbl_2samp(
    haptbl_sel,
    haptbl_base,
    aarng,
    pcount=0
):
    """ join haplotype count tables for two samples, rescale each to cpm

    Args:
        haptbl_sel (_type_): _description_
        haptbl_base (_type_): _description_
        aarng (_type_): _description_
        pcount (int, optional): _description_. Defaults to 0.

    Returns:
        _type_: _description_
    """
    h1 = haptbl_sel.loc[ (haptbl_sel['aa_num']>=aarng[0]) & (haptbl_sel['aa_num']<=aarng[1]) ].copy()
    h2 = haptbl_base.loc[ (haptbl_base['aa_num']>=aarng[0]) & (haptbl_base['aa_num']<=aarng[1]) ].copy()

    lcolidx = [cn for cn in h1.columns if cn!='total_counts' ]

    hj = pd.merge(
        h1.set_index(lcolidx)['total_counts'],
        h2.set_index(lcolidx)['total_counts'],
        how='outer',
        left_index=True,
        right_index=True,
        suffixes=['_1','_2']
    )
    
    hj = hj.fillna(0)
    hj = hj.reset_index()
    
    hj['total_counts_1'] += pcount
    hj['total_counts_2'] += pcount
    
    hj['cpm_1'] = 1e6 * hj['total_counts_1']/hj['total_counts_1'].sum()
    hj['cpm_2'] = 1e6 * hj['total_counts_2']/hj['total_counts_2'].sum()    
                
    return hj
        

def join_haptbl_Nsamp(
    lfnhaptbl,
    lsampname,
    aarng,
    pcount=0
):  
    h1 = pd.read_table(lfnhaptbl[0])
    h1 = h1.loc[ (h1['aa_num']>=aarng[0]) & (h1['aa_num']<=aarng[1]) ].copy()
    lcolidx = [cn for cn in h1.columns if cn!='total_counts' ]
    hj = h1.set_index(lcolidx)
    hj.columns=['total_counts_'+lsampname[0]]
    
    for isamp in range(1,len(lfnhaptbl)):
        hN = pd.read_table(lfnhaptbl[isamp])
        hN = hN.loc[ (hN['aa_num']>=aarng[0]) & (hN['aa_num']<=aarng[1]) ]
        assert list(hN.columns) == list(h1.columns)
        hN = hN.set_index( lcolidx )
        hN.columns=['total_counts_'+lsampname[isamp]]
        hj = pd.concat( [hj,hN], axis=1 )
        sys.stderr.write('%d haps...'%(hj.shape[0])); sys.stderr.flush()
    
    hj = hj.fillna(0)

    for sn in lsampname:
        hj['cpm_'+sn] = 1e6 * (hj['total_counts_'+sn]+pcount) / (hj['total_counts_'+sn]+pcount).sum()

    hj = hj.reset_index()

    return hj
        

def cross_sample_hap_splom(
    _hapsjoined, # should have col: 
    filt_to_status=None,
    mincpm_thresh=10,
    min_nsamps_at_thresh=1,
    cpm_cliprng=(-1,15),
):
    """ join haplotype tables across a list of samples
    can specify minimum cpm for each haplotype, and that it must appear in ≥ threshold nubmer of samples at that abundance

    Args:
        _hapsjoined (_type_): _description_
        mincpm_thresh (int, optional): _description_. Defaults to 10.
        min_nsamps_at_thresh (int, optional): _description_. Defaults to 1.
        cpm_cliprng (tuple, optional): _description_. Defaults to (-1,15).

    Returns:
        _type_: _description_
    """
    lcpmcol=[c for c in _hapsjoined.columns if c.startswith('cpm_')]
    hapsjoined = _hapsjoined[ ['status']+lcpmcol ]
    
    # remove haps that are very low abundance in all samples
    lihap_filt = (hapsjoined[lcpmcol]>=mincpm_thresh).sum(1)>=min_nsamps_at_thresh
    sys.stderr.write( 'filtered %d -> %d haplotypes\n'%(hapsjoined.shape[0],lihap_filt.sum() )); sys.stderr.flush()
    hapsjoined_filt = hapsjoined.loc[ lihap_filt ]

    hapsjoined2 = hapsjoined_filt[ ['status']+lcpmcol ]

    hapsjoined2_lg2sc = hapsjoined2.copy()
    for c in hapsjoined2.columns[1:]:
        hapsjoined2_lg2sc[c]=np.log2(hapsjoined2_lg2sc[c])
        
    if filt_to_status:
        hapsjoined2_lg2sc = hapsjoined2_lg2sc.loc[ hapsjoined2_lg2sc['status'].isin(filt_to_status) ]
        sys.stderr.write( '%d haplotypes w/ that status\n'%(hapsjoined2_lg2sc.shape[0])); sys.stderr.flush()
                
    toplot = hapsjoined2_lg2sc.copy()
    for cn in lcpmcol:
        toplot[cn]=toplot[cn].clip( cpm_cliprng[0], cpm_cliprng[1] )

    g=sns.PairGrid(toplot,
                vars=lcpmcol,
                hue='status')    

    g.map_offdiag(sns.scatterplot,s=1,ec='none')
    # g.map_diag(sns.ecdfplot)

    # weight the ecdf by the number of reads, such that y = fraction of READS in barcodes w/ cpm >= x, not fraction of barcodes w/ cpm>=x
    def wecdfplot(*args,**kwargs):
        cpm = 2**args[0]
        return sns.ecdfplot(x=np.array(args[0]),weights=cpm,**kwargs)
    
    g.map_diag(wecdfplot)

    for iax in range(len(g.axes)):
        for jax in range(len(g.axes)):
            if iax!=jax:
                g.axes[iax][jax].plot( cpm_cliprng, cpm_cliprng, lw=0.5, c='k', ls='--' )
                g.axes[iax][jax].set_xlim( cpm_cliprng[0]-0.5, cpm_cliprng[1]+0.5 )
                g.axes[iax][jax].set_ylim( cpm_cliprng[0]-0.5, cpm_cliprng[1]+0.5 )
        
    g.add_legend()

    return g

def overlapBeneathCumulSum( jh, _lsamps, frac=0.90 ):
    lsamps=list(_lsamps)
    Nsamps=len(lsamps)
    mtx = np.array( jh[ ['cpm_%s'%x for x in lsamps]  ] )

    mtxout = np.zeros( (Nsamps,Nsamps), dtype=np.float32 )
    for isamp in range(Nsamps):
        ttlSampI = mtx[ :, isamp ].sum()
        ipermSortByI = np.argsort( -(mtx[:, isamp]) )

        # print ttlSampI
        # print mtx[lsamps[isamp]][ipermSortByI].cumsum()

        ieltmin = np.argmax(mtx[:,isamp][ipermSortByI].cumsum()/float(ttlSampI) >= frac )

        for jsamp in range(Nsamps):
            ttlSampJ = mtx[ :, jsamp ].sum()
            ttlBelowThreshJ = mtx[  ipermSortByI, jsamp ][ :ieltmin ]

            mtxout[ isamp, jsamp ] = ttlBelowThreshJ.sum() / float(ttlSampJ)

    tblout = pd.DataFrame(mtxout, columns=lsamps)
    tblout['sample']=lsamps
    tblout=tblout[ ['sample']+lsamps ]
    return tblout

def make_heatmap( 
    mtx,
    subsetRows = None,
    subsetCols = None,
    cluster = False,
    title = None,
    figsize=(12,8), 
    qileused=0.90,
    vmin=0.0,
    vmax=1.0 ):

    tbls2s = mtx
    tbls2s = tbls2s.set_index('sample')

    if subsetRows is not None:
        if subsetCols is not None:
            tbls2s = tbls2s.loc[ subsetRows, subsetCols ]
        else:
            tbls2s = tbls2s.loc[ subsetRows ]
    else:
        if subsetCols is not None:
            tbls2s = tbls2s[ subsetCols ]

    lbls = np.asarray(
        ['%.2f'%v for v in np.array(tbls2s).flatten()]
    ).reshape(
        *(tbls2s.shape)
    )

    if not cluster:    
        plt.figure(figsize=figsize)
        hm = sns.heatmap( tbls2s, vmin=vmin, vmax=vmax, annot=lbls, fmt='', cmap='magma'  )
        plt.suptitle(title)
        plt.tight_layout()
    else:
        clgrid = sns.clustermap( tbls2s, figsize=figsize, vmin=vmin, vmax=vmax, annot=lbls, fmt='', cmap='magma'   )
        for tl in clgrid.ax_heatmap.yaxis.get_ticklabels():
            tl.set_rotation(0.)
        plt.suptitle(title)
        plt.tight_layout()

    f = plt.gcf()

def main():
    opts = argparse.ArgumentParser( description='weighted haplotype-haplotype overlap among sample sets' )

    def intRange(s):
        try:
            lower,upper=list(map(int, s.split(',')))
            assert lower<=upper
            return lower,upper
        except:
            raise argparse.ArgumentTypeError('range must be x,y, x<=y')
        
    opts.add_argument('--in_samptbl', required=True, help='var calling sample table', dest='in_samptbl')
    opts.add_argument('--haptbl_prepend_path', default='', help='base path to prepend to haplotype file', dest='haptbl_prepend_path')
    opts.add_argument('--desc', required=True, help='experiment description', dest='desc')
    
    opts.add_argument('--displ_vmin', default=0.20, help='min display value on heatmap', dest='displ_vmin')
    opts.add_argument('--displ_vmax', default=0.50, help='cap display value on heatmap', dest='displ_vmax')
    opts.add_argument('--qile', default=0.95, type=float, help='ref sample cumulative quantile (e.g., 0.95 = 95%ile)', dest='qile')

    opts.add_argument('--haps_filt_condition', default=None, help='query expression by which to filter haplotypes eg \'status=="singlemis_inclsyn"\' ', 
                      dest='haps_filt_condition')

    opts.add_argument('--aarng', required=True, type=intRange, help='restrict haplotypes to amino acid range x,y', dest='aarng')

    opts.add_argument('--cluster', action='store_true', default=False, dest='cluster')

    opts.add_argument('--out_base', required=True, dest='out_base')

    o = opts.parse_args()
        
    samptbl = pd.read_table( o.in_samptbl )
    samptbl['haptbl'] = [ o.haptbl_prepend_path + fn for fn in samptbl['haptbl'] ]

    hapsjoin = join_haptbl_Nsamp(
        samptbl['haptbl'],
        samptbl['libname'],
        o.aarng,
        pcount=0
    )
    
    if o.haps_filt_condition is not None:
        npre=hapsjoin.shape[0]
        hapsjoin = hapsjoin.query(o.haps_filt_condition)
        npost=hapsjoin.shape[0]
        print('filtered from %d -> %d haplotypes'%(npre,npost))

    ovlmtx = overlapBeneathCumulSum( hapsjoin, samptbl['libname'], frac=o.qile )
    ovlmtx.to_csv( o.out_base + 'samp2samp.txt', sep='\t',  )

    make_heatmap( ovlmtx, cluster=o.cluster,
                 title=o.desc+'\nSum in (column sample) of haps within (row sample) top %.1f%%ile'%(100*o.qile),
                 vmin = o.displ_vmin,
                 vmax = o.displ_vmax )

    plt.savefig( o.out_base + '.samp2samp.png' )

    # TODO make an altair version 

if __name__ == '__main__':                
    main()

    