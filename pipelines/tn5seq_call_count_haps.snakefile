######################################################
# Configuration values you must specify:
#
# > sample_table  
#   one row per sequencing library to be processed
#   required columns
#   - libname  (unique name for each library)
#   - bam
#   - target_name (which is the tile/target?
#
# ref_fasta - path to a fasta file containing the reference from which to get the tileseq primer sequences
#
# ref_seqname - name of the sequence record within (ref_fasta) to use 
#
# orf_interval - start-end , 1-based coordinates, inclusive, for the first and last base of the ORF
#
#######################################################
#  optional config values:
#
#   prefix - prefix to be added to the front of all filenames
#
#   keycol  - the name of the key column (e.g., "libname") from the input sample_table
#

#######################################################
# For example, to run:
#  
# snakemake -s /path/to/this/Snakefile \
#    --printshellcmds \
#    --cores 32 \
#    --config \

import os.path as op
import os
import pandas as pd
import pysam 
import Bio.Seq

########

assert 'sample_table' in config, 'must specify sample table'
assert 'ref_fasta' in config, 'must specify path to fasta of reference sequence'
assert 'ref_seqname' in config, 'must specify name of reference chromosome within [ref_fasta]'
assert 'orf_interval' in config, 'must specify orf coordinate range'

PREFIX = config['prefix'] if 'prefix' in config  else ''
OUT_DIR = config['outdir']

MIN_MIN_BQ = int(config['read_min_min_bq']) if 'read_min_min_bq' in config else 10
MIN_MEAN_BQ = int(config['read_min_mean_bq']) if 'read_min_mean_bq' in config else 30
VAR_MIN_BQ = int(config['var_min_min_bq']) if 'var_min_min_bq' in config else 10

########
# load and check sample table.

KEYCOL = 'libname' if 'keycol' not in config else config['keycol']

l_reqd_cols = [ KEYCOL, 'bam', 'target_name' ]
tblSamples = pd.read_table( config['sample_table'] )
assert all( [ col in tblSamples.columns for col in l_reqd_cols ] ), 'sample table must have columns: '+','.join(l_reqd_cols)
assert len(set(tblSamples[KEYCOL])) == tblSamples.shape[0], 'all libname entries must be unique'
    
lLibs = tblSamples[KEYCOL].unique()

tblSamples = tblSamples.set_index( KEYCOL,drop=False )

########
# expected output files

assert 'outdir' in config, 'must specify output directory'

lout_callvars =  expand('{}/callvars/{}{{libname}}.perread.txt'.format(OUT_DIR,PREFIX), libname=lLibs)
lout_haps =  expand('{}/counts/{}{{libname}}.byhap.txt'.format(OUT_DIR,PREFIX), libname=lLibs)
lout_counts =  expand('{}/counts/{}{{libname}}.byvar.txt'.format(OUT_DIR,PREFIX), libname=lLibs)

lout_plotvarbypos =  expand('{}/plots_varbypos/{}{{libname}}.varbypos.png'.format(OUT_DIR,PREFIX), libname=lLibs)
lout_plothapwfall =  expand('{}/plots_hapwfall/{}{{libname}}.hapwf.png'.format(OUT_DIR,PREFIX), libname=lLibs)

outRpt = OUT_DIR+'/'+PREFIX+'varcall_rpt.txt', 

lOutFiles = [outRpt] + lout_callvars + lout_haps + lout_counts + lout_plotvarbypos + lout_plothapwfall

########

rule all:
    input:
        lOutFiles

rule callvars:
    input:
        bam = lambda wc: tblSamples.loc[ wc.libname ][ 'bam' ],
    output:
        perreadtbl = op.join(OUT_DIR,'callvars/'+PREFIX+'{libname}.perread.txt'),
    params:
        ref_fasta = config['ref_fasta'],
        ref_seqname = config['ref_seqname'],
        orf_interval = config['orf_interval'],
        libname=lambda wc: wc.libname,
    threads: 1
    run:
        shell("""
            shotty_reads_call_vars \
                --ref_fasta {params.ref_fasta} \
                --bam {input.bam} \
                --orf_interval {params.ref_seqname}:{params.orf_interval} \
                --per_read_tbl {output.perreadtbl} \
                --tile_min_min_bq %d --tile_min_mean_bq %d --var_min_min_bq %d --var_min_mean_bq %d
        """%(
            MIN_MIN_BQ, MIN_MEAN_BQ, VAR_MIN_BQ, VAR_MIN_BQ
        )
        )


rule tallyvars:
    input:
        perreadtbl = rules.callvars.output.perreadtbl
    output:
        haptbl = op.join(OUT_DIR,'counts/'+PREFIX+'{libname}.byhap.txt'),
        vartbl = op.join(OUT_DIR,'counts/'+PREFIX+'{libname}.byvar.txt'),
        cvgtbl = op.join(OUT_DIR,'counts/'+PREFIX+'{libname}.cvg.txt'),
        byreadhapstatus = op.join(OUT_DIR,'counts/'+PREFIX+'{libname}.byreadstatus.txt'),
        tallylog = op.join(OUT_DIR,'counts/'+PREFIX+'{libname}.log'),
    params:
        ref_fasta = config['ref_fasta'],
        ref_seqname = config['ref_seqname'],
        orf_interval = config['orf_interval'],
        libname=lambda wc: wc.libname,
    threads: 1
    run:      
        shell("""
            codonwise_tally_vars_hapaware \
                --ref_fasta {params.ref_fasta} \
                --per_read_tbl {input.perreadtbl} \
                --orf_interval {params.ref_seqname}:{params.orf_interval} \
                --multi_syn_resolve random \
                --out_hap_tbl {output.haptbl} \
                --out_cvg_tbl {output.cvgtbl} \
                --out_byread_status_tbl {output.byreadhapstatus} \
                --out_count_tbl {output.vartbl} > {output.tallylog}
        """
        )

rule mkplots:
    input:
        vartbl = rules.tallyvars.output.vartbl,
        haptbl = rules.tallyvars.output.haptbl,
        cvgtbl = rules.tallyvars.output.cvgtbl,
    params:
        libname=lambda wc: wc.libname,
    output:
        varbypos = op.join(OUT_DIR,'plots_varbypos/'+PREFIX+'{libname}.varbypos.png'),
        # varwfall = op.join(OUT_DIR,'plots_varwfall/'+PREFIX+'{libname}.varwf.png'),
        hapwfall = op.join(OUT_DIR,'plots_hapwfall/'+PREFIX+'{libname}.hapwf.png'),
    run:
        # PADBP=120
        # corng_amp_pad = ( tblSamples.loc[params.libname, 'amp_start']-PADBP,
        #                   tblSamples.loc[params.libname, 'amp_end']+PADBP  )

        corng_orf = ( int(config['orf_interval'].split('-')[0]),int(config['orf_interval'].split('-')[1]) )

        codon_max = 1+int((1+corng_orf[1]-corng_orf[0]) / 3)

        # codonrng_amp_pad = ( int( 1+(corng_amp_pad[0]-corng_orf[0])/3 ), 
        #                     int( 1+(corng_amp_pad[1]-corng_orf[0])/3 ) )

        shell("""
            plot_mut_freq_by_pos --in_varcts {input.vartbl} --in_varcvg {input.cvgtbl} --desc "{params.libname}" --out {output.varbypos} --codon_range %d,%d
        """%(
            1,codon_max
        ))

        shell("""
            hap_wfalls hapwfall --in_hapcts {input.haptbl} --desc "{params.libname}" --remove_wt --logx --logy --out {output.hapwfall}
        """)        


rule outtbl:
    input:
        haptbl=expand(rules.tallyvars.output.haptbl,libname=lLibs),
        vartbl=expand(rules.tallyvars.output.vartbl,libname=lLibs),
        cvgtbl=expand(rules.tallyvars.output.cvgtbl,libname=lLibs),
        byreadhapstatus=expand(rules.tallyvars.output.byreadhapstatus,libname=lLibs),
        hapstatuscounts=expand(rules.tallyvars.output.tallylog,libname=lLibs),
    output:
        table_out=outRpt
    run:
        # keep a few items from the input talbe:
        # libname, target_name
        tbl_in = pd.read_table( config['sample_table'] )
        tbl_in = tbl_in[ ['libname','target_name'] ].set_index('libname',drop=False)
        
        # gather alignment stats
        hapstatus = [ pd.read_csv(fn, sep=' ',names=['hapstatus','count']) for fn in input.hapstatuscounts ]
        for i in range(len(hapstatus)):
            hapstatus[i]=hapstatus[i].set_index('hapstatus')
            hapstatus[i].columns = [ lLibs[i] ]
        hapstatus = pd.concat(hapstatus,1)

        tbl_out = pd.concat( [tbl_in, hapstatus.transpose()], 1 )       

        tbl_out['haptbl']=list(input.haptbl)
        tbl_out['vartbl']=list(input.vartbl)
        tbl_out['cvgtbl']=list(input.cvgtbl)
        tbl_out['byreadhapstatus']=list(input.byreadhapstatus)
        tbl_out['hapstatuscounts']=list(input.hapstatuscounts)
        
        tbl_out.to_csv(output['table_out'][0],sep='\t',index=False)