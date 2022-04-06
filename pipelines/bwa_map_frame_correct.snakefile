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
#        adpatorremoval_options=" --minalignmentlength 30 --adapter1 CTGTCTCTTATACACATCTCCGAGCCCACGAGAC --adapter2 CTGTCTCTTATACACATCTGACGCTGCCGACGA " \
#        bwa_options=" ref.fa " \
#        py2_environment="mypy2" \


# config values you must specify:
#   conda_env - name of the conda environment in which tools can be run
#   sample_table - must be a file w/ columns:
#       libname
#       fastq_merged
#   prefix - prefix added to the front of all filenames
#   outdir - base of output directory 
#   bwa_ref - path to bwa reference
#   bwa_options - options to pass to bwa (surrounded by quotes)

# optional config values
#   exclude_n_containing - set to true (or any value) to filter away N-containing reads before alignment
#   frameaware_realign - set to true (or any value) to filter away N-containing reads before alignment
#   refname - for frameaware-realignment, the name of the reference sequence for the CDS
#   corng_cds - for frameaware-realignment, the coordinate the reference sequence for the CDS, 1-based inclusive, e.g., 300,602

# makes:
#    {outdir}/align/unsorted/   
#    {outdir}/align/framecorrect/   
#    {output}/{PREFIX}sample_table.txt

import os.path as op
import os
import pandas as pd

########

PREFIX = config['prefix']
OUT_DIR = config['outdir']

########
# load and check sample table.

l_reqd_cols = [ 'libname', 'fastq_merged' ]

tblSamples = pd.read_csv( config['sample_table'], sep='\t' )

assert all( [ col in tblSamples.columns for col in l_reqd_cols ] ), 'sample table must have columns: '+','.join(l_reqd_cols)
assert len(set(tblSamples['libname'])) == tblSamples.shape[0], 'all libname entries must be unique'

lLibs = tblSamples['libname'].unique()

tblSamples = tblSamples.set_index( 'libname',drop=False )

########
# expected output files

assert 'prefix' in config, 'must specify a value for prefix (eg --config prefix="myexpt")'
assert 'outdir' in config, 'must specify output directory'

# final output 
lOutBams = expand('{}/align/unsorted/{}{{libname}}.bam'.format(OUT_DIR,PREFIX), libname=lLibs)
if 'frameaware_realign' in config:
    lOutBamsFar = expand('{}/align/framecorrect/{}{{libname}}.bam'.format(OUT_DIR,PREFIX), libname=lLibs)
else:
    lOutBamsFar = []

lOutOther = [op.join(OUT_DIR,PREFIX+'sample_table.txt')]

lOutFiles = lOutBams + lOutBamsFar + lOutOther

rule all:
    input:
        lOutFiles

# print(lOutFiles)

########

# map overlapped reads with bwa
rule bwa:  
    input:
        fq_merged=lambda wc:tblSamples.loc[ wc.libname ][ 'fastq_merged' ],
    output:
        bam_unsorted=op.join(OUT_DIR,'align/unsorted/'+PREFIX+'{libname}.bam'),
    params:
        ref = config['bwa_ref'],
        options = config['bwa_options']
    threads: 4
    run:
        if 'exclude_n_containing' in config:
            cmd=r"""
                zcat {input.fq_merged} | perl -ne 'push @a, $_; @a = @a[@a-4..$#a]; if ($. % 4 == 0) {{ if ($a[1] !~ /N/) {{ print "$a[0]$a[1]$a[2]$a[3]" }} }}' | \
                bwa mem  -t {threads} {params.ref} {params.options} /dev/stdin | samtools view -bS - > {output.bam_unsorted};
                """
            shell(cmd)
        else: 
            shell( r"""
                bwa mem  -t {threads} {params.ref} {params.options} {input.fq_merged} | samtools view -bS - > {output.bam_unsorted};
                """)

if 'frameaware_realign' in config:
    rule frameaware_realign:
        input: 
            bam_unsorted = rules.bwa.output.bam_unsorted
        output:
            bam_far = op.join(OUT_DIR,'align/framecorrect/'+PREFIX+'{libname}.bam')
        params:
            ref = config['bwa_ref'],
            refname = config['refname'],
            corng_cds = config['corng_cds']
        shell:
            """
            frameaware_realign --faRef {params.ref} --refName {params.refname} --corngCds {params.corng_cds} --inBam {input.bam_unsorted} --outBam {output.bam_far}
            """

    rule outtbl:
        input:
            bam=expand(rules.bwa.output.bam_unsorted, libname=lLibs),
            bamfar=expand(rules.frameaware_realign.output.bam_far, libname=lLibs)
        output:
            table_out=op.join(OUT_DIR,PREFIX+'sample_table.txt')
        run:
            global tblSamples

            tblSamples['bam']=''
            tblSamples.loc[lLibs,'bam']=[ op.abspath(fn) for fn in lOutBams ]
            tblSamples.loc[lLibs,'bam_out']=[ op.abspath(fn) for fn in lOutBamsFar]

            # write     out final table w/ output
            tblSamples.to_csv( output['table_out'], sep='\t', index=False )
else:

    rule outtbl:
        input:
            bam=expand(rules.bwa.output.bam_unsorted, libname=lLibs),
        output:
            table_out=op.join(OUT_DIR,PREFIX+'sample_table.txt')
        run:
            global tblSamples

            tblSamples['bam']=''
            tblSamples.loc[lLibs,'bam']=[ op.abspath(fn) for fn in lOutBams ]
            tblSamples.loc[lLibs,'bam_out']=[ op.abspath(fn) for fn in lOutBams ]

            # write     out final table w/ output
            tblSamples.to_csv( output['table_out'], sep='\t', index=False )
