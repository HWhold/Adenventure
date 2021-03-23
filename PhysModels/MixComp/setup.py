# Falls du irgendwelche Änderungen an der .pyx Datei vornimmst, musst du die 
# Datei neu kompilieren: python setup.py build_ext --inplace

from distutils.core import setup
from Cython.Build import cythonize
import numpy


setup(
      ext_modules = cythonize('cython_translation.pyx'),
      include_dirs=[numpy.get_include()]
)