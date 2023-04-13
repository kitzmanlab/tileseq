import argparse
import sys
import pysam 
import pandas as pd
import itertools

def main():
    
    opts = argparse.ArgumentParser(description='filter out alignments with long indels or masked regions, or w/ supplementary alignments')

    opts.add_argument( '--bam_in', required=True, dest='bam_in' )
    opts.add_argument( '--bam_out', required=True, dest='bam_out' )
    opts.add_argument( '--tbl_out', required=True, dest='tbl_out' )
    opts.add_argument( '--libname', required=True, dest='libname' )

    opts.add_argument('--max_del_len',default=10,type=int,dest='max_del_len')
    opts.add_argument('--max_ins_len',default=10,type=int,dest='max_ins_len')
    opts.add_argument('--max_clip_len',default=0,type=int,dest='max_clip_len')

    o=opts.parse_args()

    maxdel,maxins,maxclip=o.max_del_len,o.max_ins_len,o.max_clip_len

    bamin = pysam.AlignmentFile(o.bam_in)
    bamout = pysam.AlignmentFile(o.bam_out,'wb',template=bamin)

    n_pass=0
    n_filt_longdel=0
    n_filt_longins=0
    n_filt_longclip=0
    n_filt_multiparts=0
    n_filt_unaligned=0    

    for qn, _aligns in itertools.groupby( bamin, key=lambda ln:ln.query_name ):
        aligns=list(_aligns)
        if len(aligns)>1:
            n_filt_multiparts+=1
        else:
            if aligns[0].is_unmapped:
                n_filt_unaligned+=1
                continue

            cigs = aligns[0].cigartuples
            curmaxdel = max( [ 0 ] + [ t[1] for t in cigs if t[0]==2 ]  ) #2:Del
            if curmaxdel>maxdel:
                n_filt_longdel+=1
                continue

            curmaxins = max( [ 0 ] + [ t[1] for t in cigs if t[0]==1 ]  ) #1:Ins
            if curmaxins>maxins:
                n_filt_longins+=1
                continue

            curmaxclip = max( [ 0 ] + [ t[1] for t in cigs if t[0]==4 or t[0]==5 ]  ) #4/5: Soft/Hard Clip
            if curmaxclip>maxclip:
                n_filt_longclip+=1
                continue
            
            bamout.write(aligns[0])
            n_pass+=1
    
    mout={
          'aligns_pass':n_pass,
          'aligns_fail_multiparts':n_filt_multiparts,
          'aligns_fail_longdel':n_filt_longdel,
          'aligns_fail_longins':n_filt_longins,
          'aligns_fail_longclip':n_filt_longclip,
          'aligns_fail_unaligned':n_filt_unaligned
    }

    pd.Series(mout, name=o.libname).to_csv(o.tbl_out)


if __name__=='__main__':

    main()
