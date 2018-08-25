from distutils.core import setup, Extension
import os

SurfaceCode = Extension('SurfaceCode',
                    include_dirs=['./blossom5-v2.05'],
                    # libraries=['libblossom.a'],
                    # library_dirs=['./blossom5-v2.05/lib'],
                    extra_objects=['./blossom5-v2.05/lib/libblossom.a'],
                    extra_compile_args=["-std=c++11"],
                    sources=['main.cpp'])

setup(name='SurfaceCode', ext_modules=[SurfaceCode])
