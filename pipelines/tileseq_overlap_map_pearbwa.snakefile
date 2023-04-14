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
# tiledef_table, e.g., the tileseq_primers.txt file output by dmschef.
#   one row per tile/target being mutagenized 
#   required columns
#   - target_name (name of this target, e.g., tile02)
#   - pri_left_start, pri_left_end - start/end (on the top strand) of the tileseq forward primer
#   - pri_right_start, pri_right_end - start/end (on the top strand) of the tileseq reverse primer 
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
#   cutadapt_trimpri_opts - options passed to cutadapt for constant flanking tileseq primer trimming
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
#         pear_options="-q 30 -t 30 -p 0.001 -v 15 -c 80" 


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

BWA_REF = config['bwa_ref'] if 'bwa_ref' in config else config['ref_fasta'] 

BWA_OPTIONS = config['bwa_options'] if 'bwa_options' in config else '-A2 -E1 -L50,50 -B2'

MAX_INDEL_LEN = int(config['max_indel_len']) if 'max_indel_len' in config else 10

PREFIX = config['prefix'] if 'prefix' in config  else ''
OUT_DIR = config['outdir']
PEAR_OPTIONS = config['pear_options']
if 'discords_to_n' in config:
    if '-z' not in PEAR_OPTIONS: PEAR_OPTIONS+=' -z '

CUTADAPT_TRIMPRI_OPTS = config['cutadapt_trimpri_opts'] if 'cutadapt_trimpri_opts' in config else ' --error-rate 0.1 --overlap 24 '

STUFFERS = ['GGAGACGCTCAGTACGTCTAAAGCGGCCGCacgtagtgaCCTGCAGGTTTAATGCACGTAGTCGTCTCC',
            'ctta GCAGGTG CTCAGTACGTC ac GCGGCCGC acgtagtga CCTGCAGG at ATGCACGTAGT CACCTGC gcga'.lower().replace(' ','')]

########
# load and check sample table.

KEYCOL = 'libname' if 'keycol' not in config else config['keycol']

l_reqd_cols = [ KEYCOL, 'fastq_fwd', 'fastq_rev', 'target_name' ]
tblSamples = pd.read_table( config['sample_table'] )
assert all( [ col in tblSamples.columns for col in l_reqd_cols ] ), 'sample table must have columns: '+','.join(l_reqd_cols)
assert len(set(tblSamples[KEYCOL])) == tblSamples.shape[0], 'all libname entries must be unique'
    
l_reqd_cols_td = [ 'target_name', 'pri_left_start', 'pri_left_end', 'pri_right_start', 'pri_right_end' ]
tblTileDefs = pd.read_table( config['tiledef_table'] )
assert all( [ col in tblTileDefs.columns for col in l_reqd_cols_td ] ), 'tile definition table must have columns: '+','.join(l_reqd_cols_td)

assert ( len(set(tblSamples['target_name']) - set(tblTileDefs['target_name'])))==0, 'tiles %s are not in the tile defintion table, '+','.join(
    set(tblSamples['target_name']) - set(tblTileDefs['target_name'])) 

# retreive tile primer sequences
fa_ref = pysam.FastaFile(config['ref_fasta'])
tblTileDefs['seqtrim_left']=[ 
    fa_ref.fetch( config['ref_seqname'], r['pri_left_start']-1, r['pri_left_end'] )
    for _,r in tblTileDefs.iterrows() 
]
tblTileDefs['seqtrim_right']=[ 
    fa_ref.fetch( config['ref_seqname'], r['pri_right_start']-1, r['pri_right_end'] )
    for _,r in tblTileDefs.iterrows() 
]
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

# final output fastqs
# lOutFq = expand('{}/fastq/{}{{libname}}.fq.gz'.format(OUT_DIR,PREFIX), libname=lLibs)
lOutBam =  expand('{}/align/unsorted/{}{{libname}}.bam'.format(OUT_DIR,PREFIX), libname=lLibs)
lInCounts = expand('{}/fastq/{}{{libname}}.input.linects.txt'.format(OUT_DIR,PREFIX), libname=lLibs)

outRpt = OUT_DIR+'/'+PREFIX+'merge_report.txt'

lOutFiles = [outRpt] + lOutBam 

########

rule all:
    input:
        lOutFiles

rule countinputs:
    # just count the inputs
    input:
        fq_fwd = lambda wc: tblSamples.loc[ wc.libname ][ 'fastq_fwd' ],
    output:
        counts_input = op.join(OUT_DIR,'fastq/'+PREFIX+'{libname}.input.linects.txt')
    log: op.join(OUT_DIR,'fastq/'+PREFIX+'{libname}.countinputs.log')
    threads: 8
    resources:
        mem_per_cpu="4gb", 
        cpus="8", 
        time="0:30:00"
    run:
        shell("""
            pigz -d -c -p8 {input.fq_fwd} | wc -l > {output.counts_input}
        """)

rule overlap:
    # merge overlapping read pairs
    input:
        fq_fwd = lambda wc: tblSamples.loc[ wc.libname ][ 'fastq_fwd' ],
        fq_rev = lambda wc: tblSamples.loc[ wc.libname ][ 'fastq_rev' ]
    output:
        counts_out = op.join(OUT_DIR,'fastq/'+PREFIX+'{libname}.overlapped.linects.txt'),
        temp_asm_fq=temp(op.join(OUT_DIR,'fastq/'+PREFIX+'{libname}.assembled.fastq')),

        temp_unasm_fwd_fq=temp(op.join(OUT_DIR,'fastq/'+PREFIX+'{libname}.unassembled.forward.fastq')),
        temp_unasm_rev_fq=temp(op.join(OUT_DIR,'fastq/'+PREFIX+'{libname}.unassembled.reverse.fastq')),
        temp_disc_fq=temp(op.join(OUT_DIR,'fastq/'+PREFIX+'{libname}.discarded.fastq'))
    log: op.join(OUT_DIR,'fastq/'+PREFIX+'{libname}.overlap.log')
    params:
        fq_overlap_base = op.join(OUT_DIR,'fastq/'+PREFIX+'{libname}')
    threads: 24
    resources:
        mem_per_cpu="4gb", 
        cpus="24", 
        time="1:00:00"
    run:
        shell("""   
            pear -j{threads} --forward-fastq {input.fq_fwd} --reverse-fastq {input.fq_rev} {PEAR_OPTIONS} -o {params.fq_overlap_base}
            wc -l {params.fq_overlap_base}.assembled.fastq > {output.counts_out}
        """)

rule destuffer:
    input:
        fq = rules.overlap.output.temp_asm_fq
    output:
        counts_out = op.join(OUT_DIR,'fastq/'+PREFIX+'{libname}.destuff.log'),
        destuff_fq = temp(op.join(OUT_DIR,'fastq/'+PREFIX+'{libname}.destuff.fastq')),
    threads: 24
    log: op.join(OUT_DIR,'fastq/'+PREFIX+'{libname}.destuffer.log')
    resources:
        mem_per_cpu="4gb", 
        cpus="24", 
        time="1:00:00"
    run:
        stuffarg = ' '.join( ['-b %s '%stuff for stuff in STUFFERS] )

        shell("""
            cutadapt %s --error-rate 0.1 --overlap 24 --discard-trimmed -j {threads} -o {output.destuff_fq} {input.fq} > {output.counts_out}
        """%(stuffarg))
        # put log in sep file, get counts out

rule trim_tilepris:
    input:
        fq = rules.destuffer.output.destuff_fq
    output:
        counts_out = op.join(OUT_DIR,'fastq/'+PREFIX+'{libname}.trim.log'),
        trim_fq = op.join(OUT_DIR,'fastq/'+PREFIX+'{libname}.fq.gz'),
    params:
        libname=lambda wc: wc.libname
    log: op.join(OUT_DIR,'fastq/'+PREFIX+'{libname}.trim_tilepris.log')
    threads: 24
    resources:
        mem_per_cpu="4gb", 
        cpus="24", 
        time="1:00:00"
    run:
        seqtrim_left = tblSamples.loc[params.libname, 'seqtrim_left']
        seqtrim_right = tblSamples.loc[params.libname, 'seqtrim_right']
        # seqtrim_right_rc = str(Bio.Seq.Seq(seqtrim_right).reverse_complement())
        # don't reverse complement since we are on merged reads

        trimarg = '-a "^%s...%s"'%( seqtrim_left, seqtrim_right )
            
        shell("""
            cutadapt %s %s --discard-untrimmed -j {threads} -o {output.trim_fq} {input.fq} > {output.counts_out}
        """%(trimarg, CUTADAPT_TRIMPRI_OPTS ))
    

rule align:
    input:
        fq = rules.trim_tilepris.output.trim_fq
    output:
        bam = temp(op.join(OUT_DIR,'align/unsorted/'+PREFIX+'{libname}.raw.bam')),
    threads: 24
    log: op.join(OUT_DIR,'fastq/'+PREFIX+'{libname}.align.log')
    resources:
        mem_per_cpu="4gb", 
        cpus="24", 
        time="1:00:00"
    run:
        # if 'exclude_n_containing' in config:
        #     cmd=r"""
        #         zcat {input.fq_merged} | perl -ne 'push @a, $_; @a = @a[@a-4..$#a]; if ($. % 4 == 0) {{ if ($a[1] !~ /N/) {{ print "$a[0]$a[1]$a[2]$a[3]" }} }}' | \
        #         bwa mem  -t {threads} {params.ref} {params.options} /dev/stdin | samtools view -bS - > {output.bam_unsorted};
        #         """
        #     shell(cmd)
        # else: 
        shell( r"""
            bwa mem  -t {threads} {BWA_REF} {BWA_OPTIONS} {input.fq} | samtools view -bS - > {output.bam};
            """)

rule align_filt:
    input:
        bam = rules.align.output.bam
    output:
        bam = op.join(OUT_DIR,'align/unsorted/'+PREFIX+'{libname}.bam'),
        counts_out = op.join(OUT_DIR,'align/unsorted/'+PREFIX+'{libname}.bamfilt.counts.txt'),
    params:
        libname=lambda wc: wc.libname
    log: op.join(OUT_DIR,'fastq/'+PREFIX+'{libname}.align_filt.log')
    resources:
        mem_per_cpu="4gb", 
        cpus="2", 
        time="1:00:00"
    run:
        shell("filt_indel_mask_aligns --bam_in {input.bam} --bam_out {output.bam} --libname {params.libname} --tbl_out {output.counts_out} --max_del_len {MAX_INDEL_LEN} --max_ins_len {MAX_INDEL_LEN} --max_clip_len 0")


def findline_cutadapt_log(fn):
    ll=[l.rstrip() for l in open(fn,'r').readlines()]
    lnwritten = [ x.group(1).replace(',','') for x in [re.search('.*Reads written \(passing filters\):.*?([0-9,]+)',l) for l in ll ] if x]
    assert len(lnwritten)==1, 'not exactly one line with num of reads written in cutadpat log '+fn
    return int(lnwritten[0])

# gather up read counts and put in a summary table 
rule outtbl:
    input:
        fq_overlap=expand(rules.trim_tilepris.output.trim_fq, libname=lLibs),
        bam=expand(rules.align_filt.output.bam, libname=lLibs),
        counts_input = expand(rules.countinputs.output.counts_input, libname=lLibs),
        counts_postovl = expand(rules.overlap.output.counts_out, libname=lLibs),
        counts_postdestuff = expand(rules.destuffer.output.counts_out, libname=lLibs),
        counts_postpritrim = expand(rules.trim_tilepris.output.counts_out, libname=lLibs),
        counts_align = expand(rules.align_filt.output.counts_out, libname=lLibs),
    resources:
        mem_per_cpu="4gb", 
        cpus="2", 
        time="0:30:00"
    output:
        table_out=outRpt #op.join(OUT_DIR,PREFIX+'sample_table.txt')
    run:
        tbl_out = pd.read_table( config['sample_table'] )

        tbl_out['bam'] = list(input.bam)

        tbl_out['nreads_input'] = [int(open(fn,'r').readline().strip().split()[0]) for fn in input.counts_input ]
        assert (all(tbl_out['nreads_input']%4)==0), '#lines in input one of the input fqs not div by 4!'
        tbl_out['nreads_input']=(tbl_out['nreads_input']/4).astype(int)
        

        tbl_out['nreads_postoverlap'] = [int(open(fn,'r').readline().strip().split()[0]) for fn in input.counts_postovl ]
        assert (all(tbl_out['nreads_postoverlap']%4)==0), '#lines in input one of the input fqs not div by 4!'
        tbl_out['nreads_postoverlap']=(tbl_out['nreads_postoverlap']/4).astype(int)
        
        tbl_out['nreads_poststufferrmv'] = [ findline_cutadapt_log(fn) for fn in input.counts_postdestuff ]
        tbl_out['nreads_postprimerrmv'] = [ findline_cutadapt_log(fn) for fn in input.counts_postpritrim ]

        tbl_out = tbl_out.set_index(KEYCOL,drop=False)

        # gather alignment stats
        align_stats = [ pd.read_csv(fn,index_col=0) for fn in input.counts_align ]
        align_stats = pd.concat( align_stats, 1 ).transpose()
        tbl_out = pd.concat( [tbl_out, align_stats], 1 )

        #           'aligns_pass':n_pass,
        #   'aligns_fail_multiparts':n_filt_multiparts,
        #   'aligns_fail_longdel':n_filt_longdel,
        #   'aligns_fail_longins':n_filt_longins,
        #   'aligns_fail_longclip':n_filt_longclip,
        #   'aligns_fail_unaligned':n_filt_unaligned


        tbl_out['frac_fail_overlap'] = 1 - tbl_out['nreads_postoverlap'] / tbl_out['nreads_input']
        tbl_out['frac_fail_stuffer'] = 1 - tbl_out['nreads_poststufferrmv'] / tbl_out['nreads_postoverlap']
        tbl_out['frac_fail_primerrmv'] = 1 - tbl_out['nreads_postprimerrmv'] / tbl_out['nreads_poststufferrmv']

        tbl_out.to_csv(output['table_out'],sep='\t',index=False)
        
        