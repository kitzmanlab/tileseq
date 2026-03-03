from setuptools import setup, Extension
import numpy
import pysam

setup(
    ext_modules=[
        Extension(
            "tileseq.callvars.varcall_common",
            ["src/tileseq/callvars/varcall_common.pyx"],
            include_dirs=[numpy.get_include()]+pysam.get_include(),
        )
    ]
)
# ext_modules =  cythonize(["tileseq/callvars/varcall_common.pyx"], language_level="3"),
# include_dirs = [numpy.get_include()]+pysam.get_include(),
#     define_macros = pysam.get_defines(),