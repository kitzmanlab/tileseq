######################################################
# Configuration values you must specify:
#
# > sample_table  
#   one row per sequencing library to be processed
#   required columns
#   - libname  (unique name for each library)
#   - bam
#   - target_name (which is the tile/target?  This must match one entry in the tiledef_table)
#
# tiledef_table, e.g., the tileseq_primers.txt file output by dmschef.
#   one row per tile/target being mutagenized 
#   required columns
#   - target_name (name of this target, e.g., tile02)
#   - target_start,	target_end - start/end of the tile
#   - pri_left_start, pri_left_end - start/end (on the top strand) of the tileseq forward primer
#   - pri_right_start, pri_right_end - start/end (on the top strand) of the tileseq reverse primer 
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
assert 'tiledef_table' in config, 'must specify tile definition table'
assert 'ref_fasta' in config, 'must specify path to fasta of reference sequence'
assert 'ref_seqname' in config, 'must specify name of reference chromosome within [ref_fasta]'
assert 'orf_interval' in config, 'must specify orf coordinate range'

PREFIX = config['prefix'] if 'prefix' in config  else ''
OUT_DIR = config['outdir']

MIN_MIN_BQ = int(config['tile_min_min_bq']) if 'tile_min_min_bq' in config else 10
MIN_MEAN_BQ = int(config['tile_min_mean_bq']) if 'tile_min_mean_bq' in config else 30
VAR_MIN_BQ = int(config['var_min_min_bq']) if 'var_min_min_bq' in config else 10

########
# load and check sample table.

KEYCOL = 'libname' if 'keycol' not in config else config['keycol']

l_reqd_cols = [ KEYCOL, 'bam', 'target_name' ]
tblSamples = pd.read_table( config['sample_table'] )
assert all( [ col in tblSamples.columns for col in l_reqd_cols ] ), 'sample table must have columns: '+','.join(l_reqd_cols)
assert len(set(tblSamples[KEYCOL])) == tblSamples.shape[0], 'all libname entries must be unique'
    
l_reqd_cols_td = [ 'target_name', 'target_start','target_end','pri_left_start', 'pri_left_end', 'pri_right_start', 'pri_right_end' ]
tblTileDefs = pd.read_table( config['tiledef_table'] )
assert all( [ col in tblTileDefs.columns for col in l_reqd_cols_td ] ), 'tile definition table must have columns: '+','.join(l_reqd_cols_td)

assert ( len(set(tblSamples['target_name']) - set(tblTileDefs['target_name'])))==0, 'tiles %s are not in the tile defintion table, '+','.join(
    set(tblSamples['target_name']) - set(tblTileDefs['target_name'])) 

tblSamples = pd.merge(
    tblSamples,
    tblTileDefs,
    how='left',
    on='target_name'
)

lLibs = tblSamples[KEYCOL].unique()

tblSamples = tblSamples.set_index( KEYCOL,drop=False )

########
# expected output files

assert 'outdir' in config, 'must specify output directory'

lout_callvars =  expand('{}/callvars/{}{{libname}}.perread.txt.gz'.format(OUT_DIR,PREFIX), libname=lLibs)
lout_haps =  expand('{}/counts/{}{{libname}}.byhap.txt'.format(OUT_DIR,PREFIX), libname=lLibs)
lout_counts =  expand('{}/counts/{}{{libname}}.byvar.txt'.format(OUT_DIR,PREFIX), libname=lLibs)

lout_plotvarbypos =  expand('{}/plots_varbypos/{}{{libname}}.varbypos.png'.format(OUT_DIR,PREFIX), libname=lLibs)
lout_plothapwfall =  expand('{}/plots_hapwfall/{}{{libname}}.hapwf.png'.format(OUT_DIR,PREFIX), libname=lLibs)
lout_plotvarcoverage =  expand('{}/plots_varcoverage/{}{{libname}}.varcoverage.html'.format(OUT_DIR,PREFIX), libname=lLibs)

outRpt = OUT_DIR+'/'+PREFIX+'varcall_rpt.txt', 

lOutFiles = [outRpt] + lout_callvars + lout_haps + lout_counts + lout_plotvarbypos + lout_plothapwfall + lout_plotvarcoverage

########

rule all:
    input:
        lOutFiles

rule callvars:
    input:
        bam = lambda wc: tblSamples.loc[ wc.libname ][ 'bam' ],
    output:
        perreadtbl = op.join(OUT_DIR,'callvars/'+PREFIX+'{libname}.perread.txt.gz'),
    params:
        ref_fasta = config['ref_fasta'],
        ref_seqname = config['ref_seqname'],
        orf_interval = config['orf_interval'],
        stop_after_opt = f'--stop_after {config["stop_after"]}' if 'stop_after' in config else '',
        libname=lambda wc: wc.libname,
    threads: 1
    run:
        corng_amplicon = ( tblSamples.loc[params.libname, 'amp_start'],
                          tblSamples.loc[params.libname, 'amp_end']  )
        corng_target = ( tblSamples.loc[params.libname, 'target_start'],
                          tblSamples.loc[params.libname, 'target_end']  )
       
        shell("""
            tile_reads_call_vars \
                --ref_fasta {params.ref_fasta} \
                --bam {input.bam} \
                --orf_interval {params.ref_seqname}:{params.orf_interval} \
                --tile_amplicon {params.ref_seqname}:%d-%d \
                --tile_target {params.ref_seqname}:%d-%d  \
                --per_read_tbl {output.perreadtbl} \
                {params.stop_after_opt} \
                --ends_fudge 10 \
                --tile_min_min_bq %d --tile_min_mean_bq %d --var_min_min_bq %d --var_min_mean_bq %d
        """%(
            corng_amplicon[0], corng_amplicon[1],
            corng_target[0], corng_target[1],
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
        byreadhapstatus = op.join(OUT_DIR,'counts/'+PREFIX+'{libname}.byreadstatus.txt.gz'),
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
        hapwfall = op.join(OUT_DIR,'plots_hapwfall/'+PREFIX+'{libname}.hapwf.png'),
        varcoverage = op.join(OUT_DIR, 'plots_varcoverage/'+PREFIX+'{libname}.varcoverage.html' ),
        varcoverage_tbl = op.join(OUT_DIR, 'plots_varcoverage/'+PREFIX+'{libname}.varcoverage.txt' )
    run:
        PADBP=120

        corng_amp_pad = ( tblSamples.loc[params.libname, 'amp_start']-PADBP,
                          tblSamples.loc[params.libname, 'amp_end']+PADBP  )

        corng_orf = ( int(config['orf_interval'].split('-')[0]),int(config['orf_interval'].split('-')[1]) )

        codonrng_amp_pad = ( int( 1+(corng_amp_pad[0]-corng_orf[0])/3 ), 
                            int( 1+(corng_amp_pad[1]-corng_orf[0])/3 ) )

        corng_target = ( tblSamples.loc[params.libname,'target_start'], tblSamples.loc[params.libname,'target_end'] )

        # find the aa range for this target

        aarng_target = ( 1 + int((corng_target[0]-corng_orf[0])/3),
                         1 + int((corng_target[1]-2-corng_orf[0])/3) )

        shell("""
            plot_mut_freq_by_pos --in_varcts {input.vartbl} --in_varcvg {input.cvgtbl} --desc "{params.libname}" --out {output.varbypos} --codon_range %d,%d
        """%(
            corng_orf[0], corng_orf[1]
        ))

        shell("""
            hap_wfalls hapwfall --in_hapcts {input.haptbl} --desc "{params.libname}" --remove_wt --logx --logy --out {output.hapwfall}
        """)        

        shell("""
            hap_wfalls varcoverageplots --in_hapcts {input.haptbl} --in_varcts {input.vartbl} --aa_range %d,%d --desc "{params.libname}" --out_plot {output.varcoverage} --out_tbl {output.varcoverage_tbl}
        """%(  
            aarng_target[0],aarng_target[1]
        ))

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
        hapstatus = pd.concat(hapstatus,axis=1)

        tbl_out = pd.concat( [tbl_in, hapstatus.transpose()], axis=1 )       

        tbl_out['haptbl']=list(input.haptbl)
        tbl_out['vartbl']=list(input.vartbl)
        tbl_out['cvgtbl']=list(input.cvgtbl)
        tbl_out['byreadhapstatus']=list(input.byreadhapstatus)
        tbl_out['hapstatuscounts']=list(input.hapstatuscounts)
        
        tbl_out.to_csv(output['table_out'][0],sep='\t',index=False)