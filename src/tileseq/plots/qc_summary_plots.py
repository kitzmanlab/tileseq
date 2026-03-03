import sys
import argparse

import numpy as np

import pandas as pd

import altair as alt

def tn5seq_overlap_report_plot(
    fn_summarytbl,
    exptDesc='',
    sort_by=None,
    rescale_to_pct=False,
    ):
    
    # if sortby==None then use default sort; otherwise sort by annotated reference fields given as a comma-separated string (eg "target_name,replicate,condition,passage" )

    tbl = pd.read_table(fn_summarytbl)
    if sort_by :
        sort_by = sort_by.split(',')
        tbl = tbl.sort_values(by=sort_by)

    tbl['nreads_fail_adtrim']=tbl['nreads_input']-tbl['nreads_postadtrim']  

    tbl2=tbl[
        ['libname','nreads_fail_adtrim']+
        'aligns_pass_ovlreads aligns_fail_multiparts_ovlreads aligns_fail_longdel_ovlreads aligns_fail_longins_ovlreads aligns_fail_longclip_ovlreads aligns_fail_unaligned_ovlreads aligns_pass_nonovlreads aligns_fail_notproppair_nonovlreads aligns_fail_longdel_nonovlreads aligns_fail_longins_nonovlreads aligns_fail_longclip_nonovlreads aligns_fail_unaligned_nonovlreads'.split(' ')
    ]

    if rescale_to_pct:
        tbl2pcts=tbl2.set_index('libname')
        for i,r in tbl2pcts.iterrows():
            tbl2pcts.loc[i] = tbl2pcts.loc[i]/tbl2pcts.loc[i].sum()
        tbl2pcts=tbl2pcts.reset_index()
        tbl2L=tbl2pcts.melt(id_vars='libname',var_name='status',value_name='readcount')

        ch=alt.Chart(
            tbl2L,
            title='Read status during overlap/align '+exptDesc
        ).mark_bar(
        ).encode(
            alt.X('libname',sort=list(tbl.index)),
            alt.Y('readcount',axis=alt.Axis(format='.0%')),
            alt.Color('status' ),
            alt.Tooltip( ['status','readcount'] ) 
        )       

    else:
        tbl2L=tbl2.melt(id_vars='libname',var_name='status',value_name='readcount')

        ch=alt.Chart(
            tbl2L,
            title='Read status during overlap/align '+exptDesc
        ).mark_bar(
        ).encode(
            alt.X('libname',sort=list(tbl.index)),
            alt.Y('readcount'),
            alt.Color('status' ),
            alt.Tooltip( ['libname','status','readcount'] ) 
        )
    
   
    return ch

def tileseq_overlap_report_plot(
    fn_summarytbl,
    exptDesc='',
    sort_by=None,
    rescale_to_pct=False,
    ):
    
    # if sortby==None then use default sort; otherwise sort by annotated reference fields given as a comma-separated string (eg "target_name,replicate,condition,passage" )

    tbl = pd.read_table(fn_summarytbl)
    if sort_by :
        sort_by = sort_by.split(',')
        tbl = tbl.sort_values(by=sort_by)

    tbl['nreads_fail_overlap']=tbl['nreads_input']-tbl['nreads_postoverlap']
    tbl['nreads_fail_stuffer']=tbl['nreads_postoverlap']-tbl['nreads_poststufferrmv']
    tbl['nreads_fail_endseqwrong']=tbl['nreads_poststufferrmv']-tbl['nreads_postprimerrmv']

    tbl2=tbl[
        ['libname','nreads_fail_overlap','nreads_fail_stuffer','nreads_fail_endseqwrong',
        'aligns_fail_multiparts','aligns_fail_longdel','aligns_fail_longins','aligns_fail_longclip','aligns_fail_unaligned',
        'aligns_pass']
    ]

    if rescale_to_pct:
        tbl2pcts=tbl2.set_index('libname')
        for i,r in tbl2pcts.iterrows():
            tbl2pcts.loc[i] = tbl2pcts.loc[i]/tbl2pcts.loc[i].sum()
        tbl2pcts=tbl2pcts.reset_index()
        tbl2L=tbl2pcts.melt(id_vars='libname',var_name='status',value_name='readcount')

        ch=alt.Chart(
            tbl2L,
            title='Read status during overlap/align '+exptDesc
        ).mark_bar(
        ).encode(
            alt.X('libname',sort=list(tbl.index)),
            alt.Y('readcount',axis=alt.Axis(format='.0%')),
            alt.Color('status' ),
            alt.Tooltip( ['status','readcount'] ) 
        )       

    else:
        tbl2L=tbl2.melt(id_vars='libname',var_name='status',value_name='readcount')

        ch=alt.Chart(
            tbl2L,
            title='Read status during overlap/align '+exptDesc
        ).mark_bar(
        ).encode(
            alt.X('libname',sort=list(tbl.index)),
            alt.Y('readcount'),
            alt.Color('status' ),
            alt.Tooltip( ['libname','status','readcount'] ) 
        )
    
    return ch


def tileseq_hapstatus_report_plot(
    fn_summarytbl,
    exptDesc='',
    sort_by=None,
    rescale_to_pct=False,
    ):
    
    # if sortby==None then use default sort; otherwise sort by annotated reference fields given as a comma-separated string (eg "target_name,replicate,condition,passage" )

    tbl = pd.read_table(fn_summarytbl)
    if sort_by :
        sort_by = sort_by.split(',')
        tbl = tbl.sort_values(by=sort_by)

    tbl2=tbl[
        ['libname']+
        'summary_skip_readfail summary_wt summary_indelfs summary_indelnonfs_only summary_syn_only summary_stopgain summary_singlemis_only summary_singlemis_inclsyn summary_multimis_only summary_multimis_inclsyn summary_other'.split(' ')
    ]

    if rescale_to_pct:
        tbl2pcts=tbl2.set_index('libname')
        for i,r in tbl2pcts.iterrows():
            tbl2pcts.loc[i] = tbl2pcts.loc[i]/tbl2pcts.loc[i].sum()
        tbl2pcts=tbl2pcts.reset_index()
        tbl2L=tbl2pcts.melt(id_vars='libname',var_name='status',value_name='readcount')

        ch=alt.Chart(
            tbl2L,
            title='Variant/haplotype status '+exptDesc
        ).mark_bar(
        ).encode(
            alt.X('libname',sort=list(tbl.index)),
            alt.Y('readcount',axis=alt.Axis(format='.0%')),
            alt.Color('status', scale=alt.Scale(scheme='category20') ),
            alt.Tooltip( ['libname','status','readcount'] ) 
        )       

    else:
        tbl2L=tbl2.melt(id_vars='libname',var_name='status',value_name='readcount')

        ch=alt.Chart(
            tbl2L,
            title='Variant/haplotype status '+exptDesc
        ).mark_bar(
        ).encode(
            alt.X('libname',sort=list(tbl.index), axis=alt.Axis(labelLimit=600)),
            alt.Y('readcount',axis=alt.Axis(labelLimit=600)),
            alt.Color('status', scale=alt.Scale(scheme='category20') ),
            alt.Tooltip( ['libname','status','readcount'] ) 
        )
    
    return ch


def main():
    opts = argparse.ArgumentParser( description='QC plots from tileseq overlap & align processing steps' )

    sp = opts.add_subparsers(dest='cmd')
    
    cmd_overlapqc = sp.add_parser('overlapqc')

    cmd_overlapqc.add_argument('--in_summary_tbl', required=True, dest='in_summary_tbl')
    cmd_overlapqc.add_argument('--desc', required=True, help='experiment description', dest='desc')
    cmd_overlapqc.add_argument('--sortby', default='target_name,libname', help='fields to sort by', dest='sortby')
    cmd_overlapqc.add_argument('--out_base', required=True, dest='out_base')

    cmd_overlapqc_tn5 = sp.add_parser('overlapqc_tn5')

    cmd_overlapqc_tn5.add_argument('--in_summary_tbl', required=True, dest='in_summary_tbl')
    cmd_overlapqc_tn5.add_argument('--desc', required=True, help='experiment description', dest='desc')
    cmd_overlapqc_tn5.add_argument('--sortby', default='target_name,libname', help='fields to sort by', dest='sortby')
    cmd_overlapqc_tn5.add_argument('--out_base', required=True, dest='out_base')

    cmd_varcallqc = sp.add_parser('varcallqc')

    cmd_varcallqc.add_argument('--in_summary_tbl', required=True, dest='in_summary_tbl')
    cmd_varcallqc.add_argument('--desc', required=True, help='experiment description', dest='desc')
    cmd_varcallqc.add_argument('--sortby', default='target_name,libname', help='fields to sort by', dest='sortby')
    cmd_varcallqc.add_argument('--out_base', required=True, dest='out_base')

    o = opts.parse_args()

    if o.cmd=='overlapqc':
        c=tileseq_overlap_report_plot(
            o.in_summary_tbl,
            exptDesc=o.desc,
            sort_by=o.sortby,
            rescale_to_pct=False,
            )
        c.save( o.out_base+'.read_cts.html' )

        c=tileseq_overlap_report_plot(
            o.in_summary_tbl,
            exptDesc=o.desc,
            sort_by=o.sortby,
            rescale_to_pct=True
            )
        c.save( o.out_base+'.read_pct.html' )

    elif o.cmd=='overlapqc_tn5':
        c=tn5seq_overlap_report_plot(
            o.in_summary_tbl,
            exptDesc=o.desc,
            sort_by=o.sortby,
            rescale_to_pct=False,
            )
        c.save( o.out_base+'.read_cts.html' )

        c=tn5seq_overlap_report_plot(
            o.in_summary_tbl,
            exptDesc=o.desc,
            sort_by=o.sortby,
            rescale_to_pct=True
            )
        c.save( o.out_base+'.read_pct.html' )

    elif o.cmd=='varcallqc':

        c=tileseq_hapstatus_report_plot(
            o.in_summary_tbl,
            exptDesc=o.desc,
            sort_by=o.sortby,
            rescale_to_pct=False,
            )
        c.save( o.out_base+'.var_cts.html' )

        c=tileseq_hapstatus_report_plot(
            o.in_summary_tbl,
            exptDesc=o.desc,
            sort_by=o.sortby,
            rescale_to_pct=True
            )
        c.save( o.out_base+'.var_pct.html' )


if __name__ == '__main__':                
    main()

    