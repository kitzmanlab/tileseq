######################################################
# Configuration values you must specify:
#
#
# > sample_table  
#   one row per sequencing library to be processed
#   required columns
#   - libname  (unique name for each library)
#   - fastq_fwd  (path to forward read fastq file)
#   - fastq_rev  (path to reverse read fastq file)
#   - target_name (which is the tile/target?  This must match one entry in the tiledef_table)
#
# ref_fasta - path to a fasta file containing the reference from which to get the tileseq primer sequences
#
# ref_seqname - name of the sequence record within (ref_fasta) to use 
#
#
#######################################################
#  optional config values:
#
#   prefix - prefix to be added to the front of all filenames
#
#   keycol  - the name of the key column (e.g., "libname") from the input sample_table
#
#   discords_to_n  - set this to a value (any value) to direct PEAR to mark discordant overlapping bases w/ N
#
#   bwa_ref - path to bwa reference - if not listen then it is assumed that ref_fasta points to a bwa-indexed fasta file 
#
#   bwa_options - options to pass to bwa (surrounded by quotes); by default :  "-A2 -E1 -L50,50 -B2"
#
#   max_indel_len - max size of indels to pass thru into final bam

#######################################################
# For example, to run:
#  
# snakemake -s /nfs/kitzman1/jacob/dev/tileseq/pipelines/tileseq_overlap_map_pearbwa.snakefile \
#     --printshellcmds \
#     --cores 32 \
#     --config \
#         sample_table=testing0219.2.txt \
#         tiledef_table=msh6_tileseq_primers_wpilot_lenti492.txt \
#         outdir=tsprocess \
#         ref_fasta=jkp0492_plenti_MSH6_2Abl.fa \
#         ref_seqname=jkp0492_plenti_MSH6_2Abl \
#         keycol=libname \
#         ngmerge_options="-q 30 -t 30 -p 0.001 -v 15 -c 80" 


import os.path as op
import os
import pandas as pd
import pysam 
import Bio.Seq

########

assert 'sample_table' in config, 'must specify sample table'
assert 'ref_fasta' in config, 'must specify path to fasta of reference sequence'
assert 'ref_seqname' in config, 'must specify name of reference chromosome within [ref_fasta]'

BWA_REF = config['bwa_ref'] if 'bwa_ref' in config else config['ref_fasta'] 

# -B3 instead of -B2 for Tn5seq
BWA_OPTIONS = config['bwa_options'] if 'bwa_options' in config else '-A2 -E1 -L50,50 -B3 '

# TODOC
TRIM_ADAPTOR_OPTS = "-a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -A CTGTCTCTTATACACATCTGACGCTGCCGACGA --minimum-length 40:40"


MAX_INDEL_LEN = int(config['max_indel_len']) if 'max_indel_len' in config else 10

PREFIX = config['prefix'] if 'prefix' in config  else ''
OUT_DIR = config['outdir']
NGMERGE_OPTIONS = config['ngmerge_options'] if 'ngmerge_options' in config else '-m 20 -d '

########
# load and check sample table.

KEYCOL = 'libname' if 'keycol' not in config else config['keycol']

l_reqd_cols = [ KEYCOL, 'fastq_fwd', 'fastq_rev', 'target_name' ]
tblSamples = pd.read_table( config['sample_table'] )
assert all( [ col in tblSamples.columns for col in l_reqd_cols ] ), 'sample table must have columns: '+','.join(l_reqd_cols)
assert len(set(tblSamples[KEYCOL])) == tblSamples.shape[0], 'all libname entries must be unique'
    
lLibs = tblSamples[KEYCOL].unique()
tblSamples = tblSamples.set_index( KEYCOL,drop=False )

########
# expected output files

assert 'outdir' in config, 'must specify output directory'

lOutBam =  expand('{}/align/sorted/{}{{libname}}.s.bam'.format(OUT_DIR,PREFIX), libname=lLibs)

outRpt = OUT_DIR+'/'+PREFIX+'merge_report.txt'

lOutFiles = [outRpt] + lOutBam 

########

rule all:
    input:
        lOutFiles

rule trimadaptors:
    # just count the inputs
    input:
        fq_fwd = lambda wc: tblSamples.loc[ wc.libname ][ 'fastq_fwd' ],
        fq_rev = lambda wc: tblSamples.loc[ wc.libname ][ 'fastq_rev' ],
    output:
        fq_fwd_trim = temp(op.join(OUT_DIR,'fastq/'+PREFIX+'{libname}.trim.fwd.fq')),
        fq_rev_trim = temp(op.join(OUT_DIR,'fastq/'+PREFIX+'{libname}.trim.rev.fq')),
        trim_log = op.join(OUT_DIR,'fastq/'+PREFIX+'{libname}.destuff.log'),
    threads: 32
    run:
        shell("""
            cutadapt -j {threads} {TRIM_ADAPTOR_OPTS} -o {output.fq_fwd_trim} -p {output.fq_rev_trim} {input.fq_fwd} {input.fq_rev} > {output.trim_log}
        """)

rule overlap:
    # merge overlapping read pairs
    input:
        fq_fwd = rules.trimadaptors.output.fq_fwd_trim,
        fq_rev = rules.trimadaptors.output.fq_rev_trim
    output:
        counts_out = op.join(OUT_DIR,'fastq/'+PREFIX+'{libname}.overlapped.linects.txt'),
        temp_asm_fq=temp(op.join(OUT_DIR,'fastq/'+PREFIX+'{libname}.assembled.fastq')),
        temp_unasm_fwd_fq=temp(op.join(OUT_DIR,'fastq/'+PREFIX+'{libname}.unassembled_1.fastq')),
        temp_unasm_rev_fq=temp(op.join(OUT_DIR,'fastq/'+PREFIX+'{libname}.unassembled_2.fastq')),
        # temp_disc_fq=temp(op.join(OUT_DIR,'fastq/'+PREFIX+'{libname}.discarded.fastq'))
    params:
        fq_overlap_base = op.join(OUT_DIR,'fastq/'+PREFIX+'{libname}'),
        temp_unasm_base=temp(op.join(OUT_DIR,'fastq/'+PREFIX+'{libname}.unassembled')),
    threads: 30
    run:
        # shell("""   
        #     pear -j{threads} --forward-fastq {input.fq_fwd} --reverse-fastq {input.fq_rev} {PEAR_OPTIONS} -o {params.fq_overlap_base}
        #     wc -l {params.fq_overlap_base}.assembled.fastq > {output.counts_out}
        # """)
        shell("""
            NGmerge -n {threads} -1 {input.fq_fwd} -2 {input.fq_rev} -f {params.temp_unasm_base} {NGMERGE_OPTIONS} -o {output.temp_asm_fq} -v -y
            wc -l {output.temp_asm_fq} > {output.counts_out}
        """)

rule align:
    input:
        asm_fq = rules.overlap.output.temp_asm_fq,
        unasm_rev_fq = rules.overlap.output.temp_unasm_fwd_fq,
        unasm_fwd_fq = rules.overlap.output.temp_unasm_rev_fq
    output:
        bam_ovl = temp(op.join(OUT_DIR,'align/unsorted/'+PREFIX+'{libname}.ovl.raw.bam')),
        bam_nonovl = temp(op.join(OUT_DIR,'align/unsorted/'+PREFIX+'{libname}.nonovl.raw.bam')),
    params:
        options = BWA_REF + " " + BWA_OPTIONS
    threads: 30
    run:
        shell( r"""
            bwa mem  -t {threads} {params.options} {input.asm_fq} | samtools view -@ {threads} -bS - > {output.bam_ovl};
            bwa mem  -t {threads} {params.options} {input.unasm_rev_fq} {input.unasm_fwd_fq} | samtools view -@ {threads} -bS - > {output.bam_nonovl};
            """)

rule align_filt_ovl:
    input:
        bam = rules.align.output.bam_ovl,
    output:
        bam = temp(op.join(OUT_DIR,'align/unsorted/'+PREFIX+'{libname}.ovl.bam')),
        counts_out = op.join(OUT_DIR,'align/unsorted/'+PREFIX+'{libname}.ovl.counts.txt'),
    params:
        libname=lambda wc: wc.libname
    threads: 1
    run:
        shell("filt_indel_mask_aligns --bam_in {input.bam} --bam_out {output.bam} --libname {params.libname} --tbl_out {output.counts_out} --max_del_len {MAX_INDEL_LEN} --max_ins_len {MAX_INDEL_LEN} --max_clip_len 0")

rule align_filt_nonovl:
    input:
        bam = rules.align.output.bam_nonovl,
    output:
        bam = temp(op.join(OUT_DIR,'align/unsorted/'+PREFIX+'{libname}.nonovl.bam')),
        counts_out = op.join(OUT_DIR,'align/unsorted/'+PREFIX+'{libname}.nonovl.counts.txt'),
    params:
        libname=lambda wc: wc.libname
    threads: 1
    run:
        shell("filt_indel_mask_aligns_PE --bam_in {input.bam} --bam_out {output.bam} --libname {params.libname} --tbl_out {output.counts_out} --max_del_len {MAX_INDEL_LEN} --max_ins_len {MAX_INDEL_LEN} --max_clip_len 0;")

rule join_bam:
    input:
        bam_ovl = rules.align_filt_ovl.output.bam,
        bam_nonovl = rules.align_filt_nonovl.output.bam, 
    output:
        bam = op.join(OUT_DIR,'align/sorted/'+PREFIX+'{libname}.s.bam')
    threads: 8
    run:
        shell("samtools sort -@{threads} <(samtools merge -@ {threads} {input.bam_ovl} {input.bam_nonovl} -o {output.bam}; samtools sort -@{threads} {output.bam} )")

def findline_cutadapt_log(fn):
    ll=[l.rstrip() for l in open(fn,'r').readlines()]
    lnprocd = [ x.group(1).replace(',','') for x in [re.search('.*Total read pairs processed:.*?([0-9,]+)',l) for l in ll ] if x]
    lnwritten = [ x.group(1).replace(',','') for x in [re.search('.*Reads written \(passing filters\):.*?([0-9,]+)',l) for l in ll ] if x]
    assert len(lnwritten)==1, 'not exactly one line with num of reads written in cutadpat log '+fn
    return int(lnprocd[0]), int(lnwritten[0])

# gather up read counts and put in a summary table 
rule outtbl:
    input:
        trimadapt_log = expand(rules.trimadaptors.output.trim_log, libname=lLibs),
        count_ovl = expand(rules.overlap.output.counts_out,libname=lLibs),
        counts_ovl_bam_filt = expand(rules.align_filt_ovl.output.counts_out, libname=lLibs),
        counts_nonovl_bam_filt = expand(rules.align_filt_nonovl.output.counts_out, libname=lLibs),
        bam = expand( rules.join_bam.output.bam, libname=lLibs ),
    output:
        table_out=outRpt #op.join(OUT_DIR,PREFIX+'sample_table.txt')
    run:
        tbl_out = pd.read_table( config['sample_table'] )

        tbl_out['bam'] = list(input.bam)

        tbl_out['nreads_input'] = [ findline_cutadapt_log(fn)[0] for fn in input.trimadapt_log ]
        tbl_out['nreads_postadtrim'] = [ findline_cutadapt_log(fn)[1] for fn in input.trimadapt_log ]

        tbl_out['nreads_overlap'] = [int(open(fn,'r').readline().strip().split()[0]) for fn in input.count_ovl ]
        assert (all(tbl_out['nreads_overlap']%4)==0), '#lines in input one of the input fqs not div by 4!'
        tbl_out['nreads_overlap']=(tbl_out['nreads_overlap']/4).astype(int)

        # gather alignment stats
        align_stats_ovl = [ pd.read_csv(fn,index_col=0) for fn in input.counts_ovl_bam_filt ]
        align_stats_ovl = pd.concat( align_stats_ovl, 1 ).transpose()
        tbl_out = pd.concat( [tbl_out, align_stats_ovl], 1 )

        align_stats_nonovl = [ pd.read_csv(fn,index_col=0) for fn in input.counts_nonovl_bam_filt ]
        align_stats_nonovl = pd.concat( align_stats_nonovl, 1 ).transpose()
        tbl_out = pd.concat( [tbl_out, align_stats_nonovl], 1 )

        tbl_out.to_csv(output['table_out'],sep='\t',index=False)
        
        