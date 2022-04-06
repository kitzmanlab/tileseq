import sys
import argparse

import pysam

import pandas as pd

from collections import OrderedDict as odict

from tileseq.tileseq_shared import *
from tileseq.callvars.varcall_common import *

# hdrOutRg += ['sumEdDist',
# 'nCodonSyn','nCodonNonsyn','nInframeCodonaliDel',
# 'nFSIndel','nFSdel','nFSins','nIndel',
# 'ptc',
# 'codingGenotype',
# 'synVars',
# 'nPalsEntriesMatched',
# 'nNonpgmdVars',  
# 'palsEntriesMatched',
# 'nCodonSynInclAfterFS','nCodonNonsynInclAfterFS',
# 'codingGenotype_withEdits',
# 'lAcceptedVars']

import Bio.Seq
import Bio.SeqUtils 


hdrOut = [

'readname',

'pass',
'pass_covers_target',
'pass_read_bq',
'pass_vars_bq',

'sum_ed_dist',
'sum_bp_ins',
'sum_bp_del',
'sum_bp_nti',
'sum_bp_ntv',

'is_wt',
'has_indel',

'n_codon_mm1',
'n_codon_mm2',
'n_codon_mm3',
'n_codon_syn',
'n_codon_mis',

'n_codon_stopgain',
'n_codon_indel_frameshift',
'n_codon_indel_nonfs',

'tile_mean_bq',
'tile_min_bq',
'var_min_meanbq',
'var_min_minbq',

'vars_vcf_style',
'vars_cdna_eff',
'vars_pro_eff',

'uncounted_vars',   # uncounted b/c codon is on edge of alignemnt, or b/c N in query codon
]

def rgOverlap(bounds1, bounds2):
    if (min(bounds1) >= min(bounds2) and min(bounds1) <= max(bounds2) ) or \
       (max(bounds1) >= min(bounds2) and max(bounds1) <= max(bounds2) ) or \
       (min(bounds2) >= min(bounds1) and min(bounds2) <= max(bounds1) ) or \
       (max(bounds2) >= min(bounds1) and max(bounds2) <= max(bounds1) ):
        return 1
    else:
        return 0

def rngIntersect(bounds1, bounds2): 
    if bounds1==None or bounds2==None: return None
    elif rgOverlap(bounds1,bounds2) : return tuple(sorted(list(bounds1)+list(bounds2))[1:3])
    else: return None

# assumes mapped to same chrom
# assumes target and amplicon are fully contained within ORF  
def process_cdna_reads( 
    bamline, 
    ref_fasta,
    orf_interval,
    orf_strand,
    tile_min_mean_bq,
    tile_min_min_bq,
    var_min_mean_bq,
    var_min_min_bq,
    codon_pos_offset=0,
    cdna_pos_offset=0 ):

    out = odict( [ (col,None) for col in hdrOut ] )

    out['readname'] = bamline.query_name

    # bamline.ref_end = 1 past last aligned base

    gal = GappedAlign( bamline, ref_fasta )

    corz_align_isect_orf = rngIntersect(
        (bamline.reference_start, bamline.reference_end-1),
        (orf_interval[1], orf_interval[2]) )

    # is alignment outside of the orf?
    if corz_align_isect_orf is None:
        return None

    # pass_and_ofsgap = gal.ref_to_gapped_interval( bamline.reference_start, bamline.reference_end-1 )
    # problem this is getting tripped up with variants outside orf 

    pass_and_ofsgap = gal.ref_to_gapped_interval( corz_align_isect_orf[0], corz_align_isect_orf[1] )
    ofsrng_gap = (pass_and_ofsgap[1],pass_and_ofsgap[2])

    # are edges clean? (ie no indels @ target edges)
    pass_edges_clean = pass_and_ofsgap[0]

    tile_quals = bamline.query_qualities[ ofsrng_gap[0]:ofsrng_gap[1]+1 ]

    if len(tile_quals)>0:
        tile_mean_bq = np.mean(tile_quals)
        tile_min_bq = np.min(tile_quals)
    else:
        tile_mean_bq=tile_min_bq=-1

    out['tile_mean_bq']=tile_mean_bq
    out['tile_min_bq']=tile_min_bq

    pass_read_bq = (tile_mean_bq >= tile_min_mean_bq) and (tile_min_bq >= tile_min_min_bq)
    
    varlist_noncodonsplit = gal.find_variants_aggbydist(
                # coz_target_start = bamline.reference_start,
                # coz_target_end = bamline.reference_end-1,
                coz_target_start = corz_align_isect_orf[0], 
                coz_target_end = corz_align_isect_orf[1],
                maxdist=1 )

    # output a VCF-style variant list. 
    out['vars_vcf_style'] = ','.join(
        [ '%d:%s:%s'%( var.coz_ref_start+1, var.gapped_ref_seq, var.gapped_query_seq ) 
        for var in varlist_noncodonsplit ] )

    out['vars_pro_eff']=[]
    out['vars_cdna_eff']=[]

    varlist = gal.find_variants_aggbydist_groupbycodon( 
                # coz_target_start=bamline.reference_start,
                # coz_target_end=bamline.reference_end-1,
                coz_target_start = corz_align_isect_orf[0], 
                coz_target_end = corz_align_isect_orf[1],
                coz_exon_start = orf_interval[1],
                coz_exon_end=orf_interval[2],
                strand=-1 if orf_strand=='-' else 1,
                frame_start=0,
                maxdist=2 )

    if varlist is not None and len(varlist)>0:
        lv_minbq = [ v.get_min_bq_incl_flanks() for v in varlist ]    
        lv_meanbq = [ v.get_mean_bq_incl_flanks() for v in varlist ]
        out['var_min_minbq'] = min(lv_minbq)
        out['var_min_meanbq'] = min(lv_meanbq)

        pass_vars_bq = min(lv_minbq)>= var_min_min_bq and min(lv_meanbq) >= var_min_mean_bq
    else:
        out['var_min_minbq'] = None
        out['var_min_meanbq'] = None

        out['is_wt']=True
        out['has_indel']=False

        pass_vars_bq = True


    for c in ['sum_ed_dist','sum_bp_ins','sum_bp_del','sum_bp_nti','sum_bp_ntv',
              'n_codon_mm1','n_codon_mm2','n_codon_mm3','n_codon_syn','n_codon_mis',
              'n_codon_stopgain', 'n_codon_indel_nonfs', 'n_codon_indel_frameshift' ]:
        out[c]=0

    out['uncounted_vars'] = False

    if varlist is not None:
        for ivar in range(len(varlist)):

            var = varlist[ivar]

            # accumulate ed#, #ins, #del, #ti, #tv

            ed_cur = (var.n_mm+var.n_ins+var.n_del)
            out['sum_ed_dist'] += ed_cur
            out['sum_bp_ins'] += var.n_ins
            out['sum_bp_del'] += var.n_del
            out['sum_bp_nti'] += var.n_mm_ti
            out['sum_bp_ntv'] += var.n_mm_tv

            if ed_cur == 0:
                out['is_wt']=True
            else:
                out['is_wt']=False

            # find codon # containing vairant
            icodon = int( (var.coz_ref_start - orf_interval[1])/3 )
            assert icodon == int(  (var.coz_ref_end - orf_interval[1])/3 )


            coz_codon = ((orf_interval[1] + icodon * 3),
                         (orf_interval[1] + (icodon+1) * 3 - 1) )

            # if codon is not entirely within aligned region, then skip it
            if not ((coz_codon[0] >= corz_align_isect_orf[0] ) and (coz_codon[1] <= corz_align_isect_orf[1] )):
                out['uncounted_vars'] = True
                continue

            frame_varstart = (var.coz_ref_start - orf_interval[1]) % 3
            frame_varend = (var.coz_ref_end + 1 - orf_interval[1]) % 3

            gofs_codon = gal.ref_to_gapped_interval( coz_codon[0], coz_codon[1] )

            ugofs_ref_codon = (gal.gali_ugofsRef[ gofs_codon[1] ], 
                               gal.gali_ugofsRef[ gofs_codon[2] ])

            ugofs_query_codon = (gal.gali_ugofsRead[ gofs_codon[1] ], 
                               gal.gali_ugofsRead[ gofs_codon[2] ])

            seq_ref_codon = gal.ref_seq[ ugofs_ref_codon[0]:ugofs_ref_codon[1]+1 ].upper()
            seq_query_codon = gal.query_seq[ ugofs_query_codon[0]:ugofs_query_codon[1]+1 ].upper()

            if len(seq_ref_codon)!=len(seq_query_codon):
                out['has_indel']=True

                if abs(len(seq_ref_codon)-len(seq_query_codon)) % 3 == 0:
                    out['n_codon_indel_nonfs'] += 1

                    out['vars_pro_eff'].append('%d:indel_nonfs'%( codon_pos_offset+icodon+1 ))
                    out['vars_cdna_eff'].append('%d:indel_nonfs'%( cdna_pos_offset + 3*(codon_pos_offset+icodon) + 1))
                else:
                    out['n_codon_indel_frameshift'] += 1

                    out['vars_pro_eff'].append('%d:indel_fs'%( codon_pos_offset+icodon+1 ))
                    out['vars_cdna_eff'].append('%d:indel_fs'%( cdna_pos_offset + 3*(codon_pos_offset+icodon) + 1))
            else:
                out['has_indel']=False

                assert len(seq_ref_codon) == len(seq_query_codon) == 3

                if var.n_mm == 1:
                    out['n_codon_mm1']+=1
                elif var.n_mm == 2:
                    out['n_codon_mm2']+=1
                elif var.n_mm == 3:
                    out['n_codon_mm3']+=1

                aa_ref_codon = mTransTbl[ seq_ref_codon ]

                if 'N' in seq_query_codon:
                    out['uncounted_vars'] = True
                else:
                    aa_query_codon = mTransTbl[ seq_query_codon ]

                    if aa_query_codon == '*' and aa_ref_codon != '*':
                        out['n_codon_stopgain'] += 1

                        out['vars_pro_eff'].append('%d:%s>*'%( codon_pos_offset+icodon+1, maa1to3[ aa_ref_codon ] ))
                        out['vars_cdna_eff'].append('%d:%s>%s'%( cdna_pos_offset + 3*(codon_pos_offset+icodon) + 1, seq_ref_codon, seq_query_codon))


                    elif aa_query_codon == aa_ref_codon :
                        out['n_codon_syn']+=1

                        out['vars_pro_eff'].append('%d:%s='%( codon_pos_offset+icodon+1, maa1to3[ aa_ref_codon ] ))
                        out['vars_cdna_eff'].append('%d:%s>%s'%( cdna_pos_offset + 3*(codon_pos_offset+icodon) + 1, seq_ref_codon, seq_query_codon))
                    else:
                        out['n_codon_mis'] += 1

                        out['vars_pro_eff'].append('%d:%s>%s'%( codon_pos_offset+icodon+1, maa1to3[ aa_ref_codon ], maa1to3[ aa_query_codon ] ))
                        out['vars_cdna_eff'].append('%d:%s>%s'%( cdna_pos_offset + 3*(codon_pos_offset+icodon) + 1, seq_ref_codon, seq_query_codon))


    out['vars_pro_eff']=','.join(out['vars_pro_eff'])
    out['vars_cdna_eff']=','.join(out['vars_cdna_eff'])

    out['pass']= 'pass' if (pass_edges_clean and pass_read_bq and pass_vars_bq) else 'FAIL'
    out['pass_edges_clean']= 'pass' if pass_edges_clean else 'FAIL'
    out['pass_read_bq']= 'pass' if pass_read_bq else 'FAIL'
    out['pass_vars_bq']= 'pass' if pass_vars_bq else 'FAIL'

    return out


def main():
    opts = argparse.ArgumentParser( description='Call variants from cDNA tileseq reads' )

    opts.add_argument('--ref_fasta', required=True, dest='ref_fasta')

    opts.add_argument('--bam', required=True, dest='bam')

    opts.add_argument('--per_read_tbl', required=True, dest='per_read_tbl')

    opts.add_argument('--orf_interval', required=True, type=chrStartStopArg_11incl_to_00incl, dest='orf_interval',
        help='interval for the ORF, 1-based, inclusive, e.g., chrom:1000-1299')

    opts.add_argument('--orf_strand', default='+', dest='orf_strand',
        help='ORF strand, should be + or -' )

    opts.add_argument('--tile_min_mean_bq', default=25, type=int, dest='tile_min_mean_bq')
    opts.add_argument('--tile_min_min_bq', default=25, type=int, dest='tile_min_min_bq')

    opts.add_argument('--var_min_mean_bq', default=25, type=int, dest='var_min_mean_bq')
    opts.add_argument('--var_min_min_bq', default=25, type=int, dest='var_min_min_bq')

    o = opts.parse_args()

    target_chrom = o.orf_interval[0]

    bam = pysam.AlignmentFile(o.bam,'rb')
    reffa = pysam.FastaFile( o.ref_fasta )

    ############################

    Nskipped, Nprocessed = 0, 0

    tbl_result = []
    first_chunk = True

    for al in bam:
        if al.reference_name != target_chrom:
            Nskipped += 1
        else:
            res = process_cdna_reads( al, 
                   reffa,
                   o.orf_interval,
                   o.orf_strand,
                   o.tile_min_mean_bq,
                   o.tile_min_min_bq,
                   o.var_min_mean_bq,
                   o.var_min_min_bq )

            if res is None:
                Nskipped+=1
            else:
                Nprocessed += 1
                tbl_result.append( res )

        if (Nskipped+Nprocessed)%10000 == 0:
            sys.stderr.write('%d (%d skipped)... '%(Nprocessed+Nskipped, Nskipped))
            sys.stderr.flush()

        if len(tbl_result)>=100000:
            tbl_result = pd.DataFrame.from_dict( tbl_result )

            tbl_result.to_csv(o.per_read_tbl, sep='\t', index=False, header=first_chunk, mode='w' if first_chunk else 'a' )

            if first_chunk:
                first_chunk=False



            tbl_result = []


    if len(tbl_result)>0:
        tbl_result = pd.DataFrame.from_dict( tbl_result )
        tbl_result.to_csv(o.per_read_tbl, sep='\t', index=False, header=first_chunk, mode='w' if first_chunk else 'a' )


if __name__ == '__main__':                
    main()

# t=20sec / 150k for this. that's pretty ok. 
