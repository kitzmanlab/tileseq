# eg
# from your data directory:
#
#  
# snakemake -s /path/to/this/Snakefile \
#    --printshellcmds \
#    --cores 24 \
#    --config \
#        sample_table=sample_key.txt \
#        prefix="" \
#        outdir=process \
#        pear_options="-q 30 -t 30 -p 0.001 -v 15 -c 80"


# config values you must specify:
#   sample_table - path to a file with at least the following columns: 
#                           - libname (must be unique) - or the value of the option keycol
#                           - fastq_fwd  (path to read1 fastq)
#                           - fastq_rev  (path to read2 fastq)
#   prefix - prefix added to the front of all filenames
#   outdir - base of output directory 
#   pear_options - string of options to pass to pear 
# 
#  optional ones:
#
#   keycol  - the name of the key column (e.g., "libname") from the input sample_table
#
#   count_perfects  - set this to a value (any value) to direct PEAR to mark discordant overlapping bases
#          with N; additionally this script will count the number of perfect overlapping reads. 

# makes:
#    {outdir}/pear_overlap/  
#    {outdir}/{prefix}.ovl_summary.txt

import os.path as op
import os
import pandas as pd

########

PREFIX = config['prefix']
OUT_DIR = config['outdir']
PEAR_OPTIONS = config['pear_options']
COUNT_PERFECTS = 'count_perfects' in config

if COUNT_PERFECTS:
    if '-z' not in PEAR_OPTIONS: PEAR_OPTIONS+=' -z '

########
# load and check sample table.

KEYCOL = 'libname' if 'keycol' not in config else config['keycol']

l_reqd_cols = [ KEYCOL, 'fastq_fwd', 'fastq_rev' ]

tblSamples = pd.read_csv( config['sample_table'], sep='\t' )

assert all( [ col in tblSamples.columns for col in l_reqd_cols ] ), 'sample table must have columns: '+','.join(l_reqd_cols)
assert len(set(tblSamples[KEYCOL])) == tblSamples.shape[0], 'all libname entries must be unique'

lLibs = tblSamples[KEYCOL].unique()

tblSamples = tblSamples.set_index( KEYCOL,drop=False )

########
# expected output files

assert 'prefix' in config, 'must specify a value for prefix (eg --config prefix="myexpt")'
assert 'outdir' in config, 'must specify output directory'

# final output fastqs
lOutFq = expand('{}/fastq/{}{{libname}}.fq.gz'.format(OUT_DIR,PREFIX), libname=lLibs)

lOutCounts = expand('{}/fastq/{}{{libname}}.counts.txt'.format(OUT_DIR,PREFIX), libname=lLibs)

outRpt = OUT_DIR+'/'+PREFIX+'merge_report.txt'

lOutFiles = lOutFq + lOutCounts + [outRpt]

########

rule all:
    input:
        lOutFiles

# call variants per sample
rule overlap:
    input:
        fq_fwd = lambda wc: tblSamples.loc[ wc.libname ][ 'fastq_fwd' ],
        fq_rev = lambda wc: tblSamples.loc[ wc.libname ][ 'fastq_rev' ]
    output:
        fq_overlap = op.join(OUT_DIR,'fastq/'+PREFIX+'{libname}.fq.gz'),
        out_counts = op.join(OUT_DIR,'fastq/'+PREFIX+'{libname}.counts.txt'),

        temp_asm_fq=temp(op.join(OUT_DIR,'fastq/'+PREFIX+'{libname}.assembled.fastq')),
        temp_unasm_fwd_fq=temp(op.join(OUT_DIR,'fastq/'+PREFIX+'{libname}.unassembled.forward.fastq')),
        temp_unasm_rev_fq=temp(op.join(OUT_DIR,'fastq/'+PREFIX+'{libname}.unassembled.reverse.fastq')),
        temp_disc_fq=temp(op.join(OUT_DIR,'fastq/'+PREFIX+'{libname}.discarded.fastq'))
    params:
        fq_overlap_base = op.join(OUT_DIR,'fastq/'+PREFIX+'{libname}')
    threads: 8
    run:
        if COUNT_PERFECTS:
            shell("""   
                pear -j{threads} --forward-fastq {input.fq_fwd} --reverse-fastq {input.fq_rev} {PEAR_OPTIONS} -o {params.fq_overlap_base};
                export NON_N=$(sed -ne "2~4p" {params.fq_overlap_base}.assembled.fastq | grep -v N |wc -l);
                export ALL=$(sed -ne "2~4p" {params.fq_overlap_base}.assembled.fastq | wc -l);
                cat {params.fq_overlap_base}.assembled.fastq | gzip > {output.fq_overlap}
                echo "{wildcards.libname} $ALL $NON_N" | perl -pi -e 's/ /\t/gi' > {output.out_counts}
            """)
        else:
            shell("""   
                pear -j{threads} --forward-fastq {input.fq_fwd} --reverse-fastq {input.fq_rev} {PEAR_OPTIONS} -o {params.fq_overlap_base};
                export ALL=$(sed -ne "2~4p" {params.fq_overlap_base}.assembled.fastq | wc -l);
                cat {params.fq_overlap_base}.assembled.fastq | gzip > {output.fq_overlap}
                echo "{wildcards.libname} $ALL" | perl -pi -e 's/ /\t/gi' > {output.out_counts}
            """)


# gather up read counts and put in a summary table 
rule outtbl:
    input:
        fq_overlap=expand(rules.overlap.output.fq_overlap, libname=lLibs)
    output:
        table_out=outRpt #op.join(OUT_DIR,PREFIX+'sample_table.txt')
    run:
        global tblSamples

        tbl_out = tblSamples.copy()
        tbl_out.set_index(KEYCOL,inplace=True)
        tbl_out['num_pair_merged']=0

        if COUNT_PERFECTS:
            tbl_out['num_pair_merged_perfect']=0
            tbl_out['frac_of_merged_perfect']=0

        for libname ,r in tbl_out.iterrows():
            path_count_file = '{}/fastq/{}{}.counts.txt'.format(OUT_DIR,PREFIX,libname )
            lcounts = open(path_count_file,'r').readline().rstrip().split('\t')
            tbl_out.loc[libname,'num_pair_merged']=int(lcounts[1])
            if COUNT_PERFECTS:
                tbl_out.loc[libname,'num_pair_merged_perfect'] = int(lcounts[2])
                tbl_out.loc[libname,'frac_of_merged_perfect'] = int(lcounts[2])/float(lcounts[1]) if int(lcounts[1])>0 else 0

        tbl_out.loc[ lLibs, 'fastq_merged' ] = [ op.abspath(fn) for fn in  input['fq_overlap'] ]

        tbl_out.reset_index(drop=False,inplace=True)

        tbl_out.to_csv(output['table_out'],sep='\t',index=False)
        