from setuptools import setup
from setuptools import Extension

nnls_module = Extension(name='libNNLS',  # module name
                           sources=['NNLS.cpp'],    # c++ source code
                           include_dirs=[r'.\eigen-3.3.8',
                                        r'.\env\Scripts',     # Header files of dependent third-party libraries
                                         r'.\env\Lib\site-packages\pybind11\include']
                           )

setup(ext_modules=[nnls_module])