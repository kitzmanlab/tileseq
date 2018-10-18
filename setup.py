from setuptools import setup
#from distutils.core import setup
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

    packages=['tileseq'],#,'tileseq.align','tileseq.callvars'],

  #   ext_modules = [ 
  #     Extension('tile_callvar_shared',['tile_callvar_shared.pyx'], extra_objects=glob.glob(op.join( pysam.get_include()[0],'*.so')) ),
    # ],

# ext_modules =  cythonize(Extension("frameaware_realign_core", ["tileseq/align/frameaware_realign_core.pyx"]), annotate=True),

    ext_modules =  cythonize(["tileseq/align/frameaware_realign_core.pyx",
                              "tileseq/callvars/varcall_common.pyx"]),

    include_dirs = [numpy.get_include()]+pysam.get_include(),
    define_macros = pysam.get_defines(),

    entry_points = {
        'console_scripts': [ 'frameaware_realign = tileseq.align.frameaware_realign:main',
                           ]
    }
)
