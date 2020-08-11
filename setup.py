#!/usr/bin/env python
# from setuptools import setup

# setup(use_scm_version=True)

import platform
import os
import sys
import shutil
from setuptools import setup, Extension
from distutils.command.clean import clean as Clean
import numpy

# Version number
version = "0.0.1"


def readme():
    with open("README.md") as f:
        return f.read()


try:
    from Cython.Distutils import build_ext
except ImportError:
    use_cython = False
else:
    use_cython = True

# use_cython=False



# set up macro
if platform.system() == "Darwin":
    macros = [("__APPLE__", "1")]
    intel_root = os.path.join(os.path.dirname(__file__),"external/intel/linux")
    mp5lib = 'iomp5'
    extra_compile_args = ["-fopenmp"]  #!!cmk '-fpermissive'

elif "win" in platform.system().lower():
    macros = [("_WIN32", "1")]
    intel_root = os.path.join(os.path.dirname(__file__),"external/intel/windows")
    mp5lib = 'libiomp5md'
    extra_compile_args = ["/EHsc", "/openmp"]
else:
    macros = [("_UNIX", "1")]
    intel_root = os.path.join(os.path.dirname(__file__),"external/intel/linux")
    mp5lib = 'iomp5'
    extra_compile_args = ["-fopenmp"]  #!!cmk '-fpermissive'

library_list = [intel_root+"/compiler/lib/intel64"]
runtime_library_dirs = None if "win" in platform.system().lower() else library_list


# see http://stackoverflow.com/questions/4505747/how-should-i-structure-a-python-package-that-contains-cython-code
if use_cython:
    ext_modules = [
        Extension(
            name="sgkit_plink.wrap_plink_parser",
            language="c++",
            sources=[
                "sgkit_plink/wrap_plink_parser.pyx",
                "sgkit_plink/CPlinkBedFile.cpp",
            ],
            libraries = [mp5lib],
            library_dirs = library_list,
            runtime_library_dirs = runtime_library_dirs,
            include_dirs=library_list+[numpy.get_include()],
            extra_compile_args=extra_compile_args,
            define_macros=macros,
        ),
        Extension(
            name="sgkit_plink.wrap_matrix_subset",
            language="c++",
            sources=[
                "sgkit_plink/wrap_matrix_subset.pyx",
                "sgkit_plink/MatrixSubset.cpp",
            ],
            include_dirs=[numpy.get_include()],
            define_macros=macros,
        ),
    ]
    cmdclass = {"build_ext": build_ext, "clean": CleanCommand}
else:
    ext_modules = [
        Extension(
            name="sgkit_plink.wrap_plink_parser",
            language="c++",
            sources=[
                "sgkit_plink/wrap_plink_parser.cpp",
                "sgkit_plink/CPlinkBedFile.cpp",
            ],
            libraries = [mp5lib],
            library_dirs = library_list,
            runtime_library_dirs = runtime_library_dirs,
            include_dirs=library_list+[numpy.get_include()],
            extra_compile_args=extra_compile_args,
            define_macros=macros,
        ),
        Extension(
            name="sgkit_plink.wrap_matrix_subset",
            language="c++",
            sources=[
                "sgkit_plink/wrap_matrix_subset.cpp",
                "sgkit_plink/MatrixSubset.cpp",
            ],
            include_dirs=[numpy.get_include()],
            define_macros=macros,
        ),
    ]
    cmdclass = {}

install_requires = ["numpy>=1.11.3", "wheel>=0.34.2"]


class CleanCommand(Clean):
    description = "Remove build directories, and compiled files (including .pyc)"

    def run(self):
        Clean.run(self)
        if os.path.exists("build"):
            shutil.rmtree("build")
        for dirpath, _, filenames in os.walk("."):
            for filename in filenames:
                if (
                    filename.endswith(".so")
                    or filename.endswith(".pyd")
                    # or filename.find("wrap_plink_parser.cpp") != -1 # remove automatically generated source file
                    # or filename.find("wrap_matrix_subset.cpp") != -1 # remove automatically generated source file
                    or filename.endswith(".pyc")
                ):
                    tmp_fn = os.path.join(dirpath, filename)
                    print("removing", tmp_fn)
                    os.unlink(tmp_fn)




# !!!cmk see FIXUP's
setup(
    name="sgkit_plink",
    version=version,
    description="sgkit_plink",
    long_description=readme(),
    long_description_content_type="text/markdown",
    keywords="FIXUPgwas bioinformatics sets intervals ranges regions",
    url="https://fastlmm.github.io/",
    author="FIXUPFaST-LMM Team",
    author_email="FIXUPfastlmm-dev@python.org",
    license="FIXUPApache 2.0",
    classifiers=[
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python",
    ],
    packages=["sgkit_plink",],  # basically everything with a __init__.py
    install_requires=install_requires,
    # extensions
    cmdclass=cmdclass,
    ext_modules=ext_modules,
)

