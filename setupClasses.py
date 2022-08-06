from Cython.Build import cythonize
from distutils.core import setup

setup(
    ext_modules=cythonize("./GoGraph/classes/*.pyx", compiler_directives={'embedsignature': True, 'language_level': 3})
)
