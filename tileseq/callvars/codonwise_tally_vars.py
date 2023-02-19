import sys
import argparse

import pysam

import pandas as pd

import numpy as np
import numpy.random as rand

from collections import OrderedDict as odict
from collections import defaultdict

import Bio.Seq
import Bio.SeqUtils 

from tileseq.tileseq_shared import *


# main set:
# wt
# indelfs
# indelnonfs_only
# syn_only
# singlemis_only
# singlemis_inclsyn
# other 

# broken down in more detail
# singlemise1_only
# singlemise1_inclsyn
# singlemise2_only
# singlemise2_inclsyn
# singlemise3_only
# singlemise3_inclsyn



# options:
#   include_any_fail  TODO
#   include_fail_readbq  TODO
#   include_fail_varbq  TODO
#    
# (1) prepare summary
# (2) singlemut_tbl
#  	 - (mis + nonse + syn == 1)  --> record the single mut. 
# (3) singlemut_tbl_allowsyns
#    - (mis + nonse <= 1), syn==anything --> record the mis or non if it exists; if ONLY syns then dont count


def main():
    opts = argparse.ArgumentParser( description='Tally tileseq variant call table' )

    opts.add_argument('--ref_fasta', required=True, dest='ref_fasta')

    opts.add_argument('--per_read_tbl', required=True, dest='per_read_tbl')

    opts.add_argument('--orf_interval', required=True, type=chrStartStopArg_11incl_to_00incl, dest='orf_interval',
        help='interval for the ORF, 1-based, inclusive, e.g., chrom:1000-1299')

    opts.add_argument('--orf_strand', default='+', dest='orf_strand',
        help='ORF strand, should be + or -' )

    opts.add_argument('--out_count_tbl', required=True, dest='out_count_tbl')

    opts.add_argument('--out_cvg_tbl', required=True, dest='out_cvg_tbl')

    opts.add_argument('--multi_syn_resolve', default='dontcount', choices=['dontcount','random','first'], dest='multi_syn_resolve', help='when there are multiple syn variants, and only syn variants, how to resolve')


    # opts.add_argument('--tile_amplicon', required=True,  type=chrStartStopArg_11incl_to_00incl, dest='tile_amplicon',
    #     help='interval for teh amplicon itself, including constant primer regions, 1-based, inclusive')

    o = opts.parse_args()

    target_chrom = o.orf_interval[0]

    reffa = pysam.FastaFile( o.ref_fasta )

    def ed(s,t):
        res=0
        assert len(s)==len(t)
        for i in range(len(s)):
            if s[i].upper()!=t[i].upper():
                res+=1
        return res

    def splitvar(v):
        p0,p1=v.split(':')
        return ( int(p0), p1.split('>')[0], p1.split('>')[1] )

    def edvar(v):
        v=splitvar(v)
        return ed( v[1],v[2] )

    #initialize list of codons in the tile 
    m_codonnum_wtseq={}
    
    Ncodons = int( (o.orf_interval[2]-o.orf_interval[1]+1)/3 )

    for codon_num in range(1, Ncodons+1):
        corz_codon = ( o.orf_interval[1] + 3*(codon_num-1), 
                       o.orf_interval[1] + 3*(codon_num)-1 )

        seq_codon = reffa.fetch( target_chrom, corz_codon[0], corz_codon[1]+1 ).upper()
        res_aa = mTransTbl[seq_codon]

        m_codonnum_wtseq[ codon_num ] = seq_codon        


    summary_wt=0
    summary_indelfs=0
    summary_indelnonfs_only=0
    summary_syn_only=0
    summary_singlemis_only=0
    summary_stopgain=0
    summary_singlemis_inclsyn=0
    summary_other=0

    summary_total_counted=0
    summary_total_skipped=0

    tbl_muts_templ = { k:[] for k in 
        ['aa_num','codon_ref','codon_mut','cdna_coord','aa_ref','aa_mut','class','ed_dist'] }

    for codon_num in sorted(m_codonnum_wtseq.keys()):
        wtcodon=m_codonnum_wtseq[codon_num]
        for ocodon in mTransTbl:
            if wtcodon==ocodon: continue
            tbl_muts_templ['aa_num'].append(codon_num)
            tbl_muts_templ['codon_ref'].append(wtcodon)
            tbl_muts_templ['codon_mut'].append(ocodon)
            tbl_muts_templ['cdna_coord'].append( 1+3*(codon_num-1) )
            tbl_muts_templ['aa_ref'].append( mTransTbl[wtcodon] )
            tbl_muts_templ['aa_mut'].append( mTransTbl[ocodon] )
            if tbl_muts_templ['aa_mut'][-1] == '*': tbl_muts_templ['class'].append('NON')
            elif tbl_muts_templ['aa_mut'][-1] == tbl_muts_templ['aa_ref'][-1]: tbl_muts_templ['class'].append('SYN')
            else: tbl_muts_templ['class'].append('MIS')

            tbl_muts_templ['ed_dist'].append(ed(wtcodon,ocodon))

    # (codon_num,wt_codon,mut_codon)->count
    cts_singlemut = defaultdict(int)
    cts_singlemut_plussyn = defaultdict(int)

    codnumz_to_cvg = np.zeros(Ncodons, dtype='int')

    reader = pd.read_table(o.per_read_tbl, chunksize=10000)

    ctr=0
    for chunk in reader:

        #keep track of indices to add to each table.
        lacc_smtbl=[]
        lacc_as_smtbl=[]

        for _,l in chunk.iterrows():
            ctr+=1

            count_in_tally=False

            if ctr%10000 == 0:
                sys.stderr.write('%d..'%(ctr))
                sys.stderr.flush()

            if l['pass']!='pass': 
                summary_total_skipped+=1
                continue
            else:
                summary_total_counted+=1

                icstart=int(l['codon_start'])
                icend=int(l['codon_end'])
                codnumz_to_cvg[ icstart-1:icend ] += 1

            if l['is_wt']:
                summary_wt+=1
                continue

            if l['n_codon_indel_frameshift']>0:
                summary_indelfs+=1
                continue

            if l['n_codon_indel_nonfs']>0:
                if l['n_codon_indel_nonfs']==1 and\
                   (l['n_codon_syn']+l['n_codon_mis']+l['n_codon_stopgain'])==0:                    
                    
                    summary_indelnonfs_only+=1
                    continue
                else:
                    summary_other += 1
                    continue

            if l['n_codon_stopgain'] > 0:
                summary_stopgain += 1
            elif l['n_codon_mis'] == 1 :
                if l['n_codon_syn'] == 0:  # and stopgain==0
                    summary_singlemis_only  += 1                    
                else:
                    summary_singlemis_inclsyn  += 1
            elif l['n_codon_mis'] == 0 and l['n_codon_syn'] >= 1: # and stopgain==0
                summary_syn_only += 1
            else:
                summary_other += 1

            if (l['n_codon_stopgain']+l['n_codon_mis']+l['n_codon_syn'])==1:
                curvar = l['vars_cdna_eff'].split(',')
                assert len(curvar)==1
                curvar = curvar[0]
                scv = splitvar(curvar)
                scv = (1+int(scv[0]-1)/3,scv[1],scv[2])
                cts_singlemut[ scv ] += 1
            elif (l['n_codon_stopgain']+l['n_codon_mis'])==1 and l['n_codon_syn']>1:
                # one stopgain or mis
                for curvar,pro in zip( l['vars_cdna_eff'].split(','), 
                                  l['vars_pro_eff'].split(',') ):
                    #    - (mis + nonse <= 1), syn==anything --> record the mis or non if it exists; if ONLY syns then dont
                    # if its syn, or stop, or doesnt match missense spec
                    if '=' not in pro: # and '*' not in pro and '>' in pro:
                        cur_ed = edvar(curvar)
                        scv = splitvar(curvar) 
                        scv = (1+int(scv[0]-1)/3,scv[1],scv[2])
                        cts_singlemut_plussyn[ scv ] += 1 
            elif (l['n_codon_stopgain']+l['n_codon_mis'])==0 and l['n_codon_syn']>1:
                # no stopgains, no mis, only syns
                lcurvar,lpro = l['vars_cdna_eff'].split(','), l['vars_pro_eff'].split(',')

                if o.multi_syn_resolve == 'dontcount':
                    continue
                elif o.multi_syn_resolve == 'random':
                    isyn = rand.randint(0,len(lcurvar))
                    scv = splitvar(lcurvar[isyn])
                    scv = (1+int(scv[0]-1)/3,scv[1],scv[2])
                    cts_singlemut_plussyn[ scv ] += 1
                elif o.multi_syn_resolve == 'first':
                    isyn = 0
                    scv = splitvar(lcurvar[isyn])
                    scv = (1+int(scv[0]-1)/3,scv[1],scv[2])
                    cts_singlemut_plussyn[ scv ] += 1


    cts_singlemut = pd.Series( cts_singlemut )
    cts_singlemut_plussyn = pd.Series( cts_singlemut_plussyn )
    
    bycodon_tbl = pd.DataFrame( tbl_muts_templ ).set_index( ['aa_num','codon_ref','codon_mut'] )
    bycodon_tbl['singlemut_reads'] = cts_singlemut
    bycodon_tbl['singlemut_reads'] = bycodon_tbl['singlemut_reads'].fillna(0).astype(int)
    bycodon_tbl['singlemut_plussyn_reads'] = cts_singlemut_plussyn
    bycodon_tbl['singlemut_plussyn_reads'] = bycodon_tbl['singlemut_plussyn_reads'].fillna(0).astype(int)

    cvg_tbl=pd.DataFrame(
        { 'codon_num':np.arange(1,Ncodons+1),
          'total_pass_reads': codnumz_to_cvg }  )

    bycodon_tbl.to_csv( o.out_count_tbl, sep='\t' )

    if o.out_cvg_tbl is not None:
        cvg_tbl.to_csv( o.out_cvg_tbl, index=False, sep='\t' )

    sys.stderr.write('done\n')
    sys.stderr.flush()

    print('summary_wt',summary_wt)
    print('summary_indelfs',summary_indelfs)
    print('summary_indelnonfs_only',summary_indelnonfs_only)
    print('summary_syn_only',summary_syn_only)
    print('summary_singlemis_only',summary_singlemis_only)
    print('summary_stopgain',summary_stopgain)
    print('summary_singlemis_inclsyn',summary_singlemis_inclsyn)
    print('summary_other',summary_other)
    print('summary_total_counted',summary_total_counted)
    print('summary_total_skipped',summary_total_skipped)

if __name__ == '__main__':                
    main()

    # summary_wt=0
    # summary_indelfs=0
    # summary_indelnonfs_only=0
    # summary_syn_only=0
    # summary_singlemis_only=0
    # summary_singlemis_only_byed={1:0, 2:0, 3:0}
    # summary_singlemis_inclsyn=0
    # summary_singlemis_inclsyn_byed={1:0, 2:0, 3:0}
    # summary_other=0




#TODO need to handle aa ofs#
#TODO need to handle cdna ofs#



# (1) prepare summary
# (2) singlemut_tbl
#    - (mis + nonse + syn == 1)  --> record the single mut. 
# (3) singlemut_tbl_allowsyns
#    - (mis + nonse <= 1), syn==anything --> record the mis or non if it exists; if ONLY syns then dont count
