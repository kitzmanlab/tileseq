import sys
import argparse

import pysam

import pandas as pd

import numpy as np
import numpy.random as rand

rand.seed(0)

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

    opts.add_argument('--out_hap_tbl', required=True, dest='out_hap_tbl')

    opts.add_argument('--out_byread_status_tbl', required=True, dest='out_byread_status_tbl')

    opts.add_argument('--multi_syn_resolve', default='dontcount', choices=['dontcount','random','first','split'], dest='multi_syn_resolve', help='when there are multiple syn variants, and only syn variants, how to resolve')


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

    hdr_out_byreadstatus = ['readname','status']

    #initialize list of codons in the tile 
    m_codonnum_wtseq={}
    
    Ncodons = int( (o.orf_interval[2]-o.orf_interval[1]+1)/3 )

    codnumz_to_cvg = np.zeros(Ncodons, dtype='int')

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
    summary_stopgain=0
    summary_singlemis_only=0
    summary_singlemis_inclsyn=0
    summary_multimis_only=0
    summary_multimis_inclsyn=0
    summary_other=0

    #summary_total_counted=0
    summary_skip_readfail=0

    lcolout_mutstbl = ['aa_num','codon_ref','codon_mut','cdna_coord','aa_ref','aa_mut','class','ed_dist']

    # haplotype tracking
    m_hap_idx = {}
    tbl_haps_out = { k:[] for k in 
        ['hap','status','total_counts','counted_codonwise']+lcolout_mutstbl }

    mhapid_rand = {}

    tbl_muts_templ = { k:[] for k in 
        lcolout_mutstbl }

    for codon_num in sorted(m_codonnum_wtseq.keys()):
        wtcodon=m_codonnum_wtseq[codon_num]
        for ocodon in mTransTbl:
            if wtcodon==ocodon: continue
            elif 'N' in ocodon: continue
            elif 'N' in wtcodon:
                raise ValueError(f'wt codon {codon_num} has N, but degenerate bases not allowed in reference')
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

    tbl_muts_templ2 = pd.DataFrame(tbl_muts_templ).set_index( ['aa_num','codon_ref','codon_mut'] )

    # (codon_num,wt_codon,mut_codon)->count
    cts_singlemut = defaultdict(float)
    cts_singlemut_plussyn = defaultdict(float)

    reader = pd.read_table(o.per_read_tbl, chunksize=10000, compression='infer')

    first_chunk=True

    ctr=0
    for chunk in reader:

        #keep track of indices to add to each table.
        lacc_smtbl=[]
        lacc_as_smtbl=[]

        tbl_out_byreadstatus = { k:[] for k in  hdr_out_byreadstatus  }

        for _,l in chunk.iterrows():
            ctr+=1

            count_in_tally=False

            tbl_out_byreadstatus['readname'].append(l['readname'])

            if ctr%10000 == 0:
                sys.stderr.write('%d..'%(ctr))
                sys.stderr.flush()

            if l['pass']!='pass': 
                summary_skip_readfail+=1
                tbl_out_byreadstatus['status'].append('skipped')
                continue
            else:
                #summary_total_counted+=1

                icstart=int(l['codon_start'])
                icend=int(l['codon_end'])
                codnumz_to_cvg[ icstart-1:icend ] += 1

            hapid = l['vars_vcf_style']

            if l['is_wt']:
                summary_wt+=1

                tbl_out_byreadstatus['status'].append('wildtype')

                if hapid not in m_hap_idx:
                    ihap = len(tbl_haps_out['hap'])
                    m_hap_idx[hapid] = ihap
                    tbl_haps_out['hap'].append(hapid)
                    tbl_haps_out['status'].append('wildtype')
                    tbl_haps_out['total_counts'].append(0)
                    tbl_haps_out['counted_codonwise'].append(False)
                    for c in lcolout_mutstbl: 
                        tbl_haps_out[c].append(None)
                    tbl_haps_out['aa_num'][-1]=-1
                    tbl_haps_out['cdna_coord'][-1]=-1
                    tbl_haps_out['ed_dist'][-1]=0
                else:
                    ihap = m_hap_idx[hapid]

                tbl_haps_out['total_counts'][ihap]+=1

                continue

            if l['n_codon_indel_frameshift']>0:
                summary_indelfs+=1

                tbl_out_byreadstatus['status'].append('indelfs')

                if hapid not in m_hap_idx:
                    ihap = len(tbl_haps_out['hap'])
                    m_hap_idx[hapid] = ihap
                    tbl_haps_out['hap'].append(hapid)
                    tbl_haps_out['status'].append('indelfs')
                    tbl_haps_out['total_counts'].append(0)
                    tbl_haps_out['counted_codonwise'].append(False)
                    for c in lcolout_mutstbl: 
                        tbl_haps_out[c].append(None)
                    tbl_haps_out['aa_num'][-1]=-1
                    tbl_haps_out['cdna_coord'][-1]=-1
                    tbl_haps_out['ed_dist'][-1]=l['sum_ed_dist']
                else:
                    ihap = m_hap_idx[hapid]

                tbl_haps_out['total_counts'][ihap]+=1

                continue

            if l['n_codon_indel_nonfs']>0:
                if l['n_codon_indel_nonfs']==1 and\
                   (l['n_codon_syn']+l['n_codon_mis']+l['n_codon_stopgain'])==0:                    

                    tbl_out_byreadstatus['status'].append('indelnonfs')

                    if hapid not in m_hap_idx:
                        ihap = len(tbl_haps_out['hap'])
                        m_hap_idx[hapid] = ihap
                        tbl_haps_out['hap'].append(hapid)
                        tbl_haps_out['status'].append('indelnonfs')
                        tbl_haps_out['total_counts'].append(0)
                        tbl_haps_out['counted_codonwise'].append(False)
                        for c in lcolout_mutstbl: 
                            tbl_haps_out[c].append(None)
                        tbl_haps_out['aa_num'][-1]=-1
                        tbl_haps_out['cdna_coord'][-1]=-1
                        tbl_haps_out['ed_dist'][-1]=l['sum_ed_dist']
                    else:
                        ihap = m_hap_idx[hapid]

                    tbl_haps_out['total_counts'][ihap]+=1


                    summary_indelnonfs_only+=1
                    continue
                else:

                    tbl_out_byreadstatus['status'].append('other')

                    if hapid not in m_hap_idx:
                        ihap = len(tbl_haps_out['hap'])
                        m_hap_idx[hapid] = ihap
                        tbl_haps_out['hap'].append(hapid)
                        tbl_haps_out['status'].append('other')
                        tbl_haps_out['total_counts'].append(0)
                        tbl_haps_out['counted_codonwise'].append(False)
                        for c in lcolout_mutstbl: 
                            tbl_haps_out[c].append(None)
                        tbl_haps_out['aa_num'][-1]=-1
                        tbl_haps_out['cdna_coord'][-1]=-1
                        tbl_haps_out['ed_dist'][-1]=l['sum_ed_dist']
                    else:
                        ihap = m_hap_idx[hapid]

                    tbl_haps_out['total_counts'][ihap]+=1


                    summary_other += 1
                    continue

            if l['n_codon_stopgain'] > 0:
                curstatus='stopgain'
                summary_stopgain += 1
            elif l['n_codon_mis'] == 1 :
                if l['n_codon_syn'] == 0:  # and stopgain==0
                    curstatus='singlemis'
                    summary_singlemis_only  += 1              
                else:
                    curstatus = 'singlemis_inclsyn'
                    summary_singlemis_inclsyn  += 1
            elif l['n_codon_mis'] == 0 and l['n_codon_syn'] >= 1: # and stopgain==0
                curstatus = 'syn_only'
                summary_syn_only += 1
            elif l['n_codon_mis'] >=2 :
                if l['n_codon_syn'] >= 1:
                    curstatus = 'multimis_inclsyn'
                    summary_multimis_inclsyn += 1
                else:
                    curstatus = 'multimis'
                    summary_multimis_only += 1
            else:
                curstatus = 'other2'
                summary_other += 1

            scv=None

            if (l['n_codon_stopgain']+l['n_codon_mis']+l['n_codon_syn'])==1:
                curvar = l['vars_cdna_eff'].split(',')
                assert len(curvar)==1
                curvar = curvar[0]
                scv = splitvar(curvar)
                scv = (1+int((scv[0]-1)/3),scv[1],scv[2])
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
                        scv = (1+int((scv[0]-1)/3),scv[1],scv[2])
                        cts_singlemut_plussyn[ scv ] += 1 
            elif (l['n_codon_stopgain']+l['n_codon_mis'])==0 and l['n_codon_syn']>1:
                # no stopgains, no mis, only syns
                lcurvar,lpro = l['vars_cdna_eff'].split(','), l['vars_pro_eff'].split(',')

                if o.multi_syn_resolve == 'dontcount':
                    curstatus = curstatus + '_notcounted'
                    scv = None
                elif o.multi_syn_resolve == 'random':
                    # needs to resolve to the same randomly selected syn var each time this hap is encounterered
                    # if hapid not in mhapid_rand :
                    if l['vars_cdna_eff'] not in mhapid_rand:  # key on the cdna eff instead of hapid, b/c of edge case where sometimes shotgun read starts midway thru a codon and sometimes not, and if a var in that codon its sometimes counted other times not
                        isyn = rand.randint(0,len(lcurvar))
                        mhapid_rand[hapid]=isyn
                    else:
                        # isyn = mhapid_rand[hapid]
                        isyn = mhapid_rand[l['vars_cdna_eff']]

                    scv = splitvar(lcurvar[isyn])
                    scv = (1+int((scv[0]-1)/3),scv[1],scv[2])
                    cts_singlemut_plussyn[ scv ] += 1
                elif o.multi_syn_resolve == 'split':
                    for cv in lcurvar:
                        scv = splitvar(cv)
                        scv = (1+int((scv[0]-1)/3),scv[1],scv[2])
                        cts_singlemut_plussyn[ scv ] += 1./len(lcurvar)
                elif o.multi_syn_resolve == 'first':
                    isyn = 0
                    scv = splitvar(lcurvar[isyn])
                    scv = (1+int((scv[0]-1)/3),scv[1],scv[2])
                    cts_singlemut_plussyn[ scv ] += 1

            tbl_out_byreadstatus['status'].append(curstatus)

            if scv is not None:

                # as of now, some like double mis are not counted and fall through the above block with scv still None.  
                # 
                # need to modify this script to be able to track higher order combos of missense w/ syn haplotypes
                #
                if hapid not in m_hap_idx:
                    ihap = len(tbl_haps_out['hap'])
                    m_hap_idx[hapid] = ihap
                    tbl_haps_out['hap'].append(hapid)
                    tbl_haps_out['status'].append(curstatus)
                    tbl_haps_out['total_counts'].append(0)
                    tbl_haps_out['counted_codonwise'].append(True)
                    tbl_haps_out['aa_num'].append(scv[0])
                    tbl_haps_out['codon_ref'].append(scv[1])
                    tbl_haps_out['codon_mut'].append(scv[2])
                    rec_cod2aa = tbl_muts_templ2.loc[ scv ]
                    for c in ['cdna_coord','aa_ref','aa_mut','class']:
                        tbl_haps_out[c].append(rec_cod2aa[c])
                    tbl_haps_out['ed_dist'].append(l['sum_ed_dist'])
                else:
                    ihap = m_hap_idx[hapid]
            else:
                if hapid not in m_hap_idx:
                    ihap = len(tbl_haps_out['hap'])
                    m_hap_idx[hapid] = ihap
                    tbl_haps_out['hap'].append(hapid)
                    tbl_haps_out['status'].append(curstatus)
                    tbl_haps_out['total_counts'].append(0)
                    tbl_haps_out['counted_codonwise'].append(False)
                    tbl_haps_out['aa_num'].append(-1)
                    tbl_haps_out['codon_ref'].append('')
                    tbl_haps_out['codon_mut'].append('')
                    for c in ['cdna_coord','aa_ref','aa_mut','class']:
                        tbl_haps_out[c].append('')
                    tbl_haps_out['ed_dist'].append(l['sum_ed_dist'])
                else:
                    ihap = m_hap_idx[hapid]

            tbl_haps_out['total_counts'][ihap]+=1
        
        if len(tbl_out_byreadstatus['readname'])>0:
            tbl_out_byreadstatus = pd.DataFrame.from_dict( tbl_out_byreadstatus )
            tbl_out_byreadstatus.to_csv(o.out_byread_status_tbl, sep='\t', index=False, header=first_chunk, mode='w' if first_chunk else 'a' )
            first_chunk = False

    tbl_haps_out = pd.DataFrame(tbl_haps_out)
    tbl_haps_out['aa_num'] = tbl_haps_out['aa_num'].astype(int)
    tbl_haps_out = tbl_haps_out.sort_values(by=['aa_num','codon_mut','total_counts'])

    tbl_haps_out.to_csv( o.out_hap_tbl, sep='\t', index=False )

    cts_singlemut = pd.Series( cts_singlemut ).astype(int)
    cts_singlemut_plussyn = pd.Series( cts_singlemut_plussyn ).astype(int)

    bycodon_tbl = pd.DataFrame( tbl_muts_templ ).set_index( ['aa_num','codon_ref','codon_mut'] )
    bycodon_tbl['singlemut_reads'] = cts_singlemut
    bycodon_tbl['singlemut_reads'] = bycodon_tbl['singlemut_reads'].fillna(0).astype(int)
    bycodon_tbl['singlemut_allowsyn_reads'] = cts_singlemut_plussyn
    bycodon_tbl['singlemut_allowsyn_reads'] = bycodon_tbl['singlemut_allowsyn_reads'].fillna(0).astype(int)

    cvg_tbl=pd.DataFrame(
        { 'codon_num':np.arange(1,Ncodons+1),
          'total_pass_reads': codnumz_to_cvg }  )

    bycodon_tbl.to_csv( o.out_count_tbl, sep='\t' )

    if o.out_cvg_tbl is not None:
        cvg_tbl.to_csv( o.out_cvg_tbl, index=False, sep='\t' )

    sys.stderr.write('done\n')
    sys.stderr.flush()

    print('summary_skip_readfail',summary_skip_readfail)
    print('summary_wt',summary_wt)
    print('summary_indelfs',summary_indelfs)
    print('summary_indelnonfs_only',summary_indelnonfs_only)
    print('summary_syn_only',summary_syn_only)
    print('summary_stopgain',summary_stopgain)
    print('summary_singlemis_only',summary_singlemis_only)
    print('summary_singlemis_inclsyn',summary_singlemis_inclsyn)
    print('summary_multimis_only',summary_multimis_only)
    print('summary_multimis_inclsyn',summary_multimis_inclsyn)
    print('summary_other',summary_other)

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
