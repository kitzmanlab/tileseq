from setuptools import setup
# from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext

import pysam
import numpy
import glob
import os.path as op

setup(
    name="tileseq",

    # cmdclass = {'build_ext':build_ext},

    packages=['tileseq','tileseq.align','tileseq.callvars'],

  #   ext_modules = [ 
  #     Extension('tile_callvar_shared',['tile_callvar_shared.pyx'], extra_objects=glob.glob(op.join( pysam.get_include()[0],'*.so')) ),
    # ],

# ext_modules =  cythonize(Extension("frameaware_realign_core", ["tileseq/align/frameaware_realign_core.pyx"]), annotate=True),

    ext_modules =  cythonize(["tileseq/callvars/varcall_common.pyx"], language_level="3"),

    include_dirs = [numpy.get_include()]+pysam.get_include(),
    define_macros = pysam.get_defines(),

    entry_points = {
        'console_scripts': [ 'filt_indel_mask_aligns = tileseq.align.filt_indel_mask_aligns:main',
                             'filt_indel_mask_aligns_PE = tileseq.align.filt_indel_mask_aligns_PE:main',
                             'tile_reads_call_vars = tileseq.callvars.tile_reads_call_vars:main',
                             'shotty_reads_call_vars = tileseq.callvars.shotty_reads_call_vars:main',
                             'codonwise_tally_vars_hapaware = tileseq.callvars.codonwise_tally_vars_hapaware:main',

                             'qc_summary_plots = tileseq.plots.qc_summary_plots:main',
                             'plot_mut_freq_by_pos = tileseq.plots.plot_mut_freq_by_pos:main',
                             'hap_wfalls = tileseq.plots.hap_wfalls:main'
                           ]
    }
)
