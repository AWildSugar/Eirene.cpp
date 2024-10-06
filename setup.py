# Ripped from ripser.py
import sys
import os
import platform
# python setup.py build_ext --inplace

from setuptools import setup, find_packages
from setuptools.extension import Extension

try:
    from Cython.Build import cythonize
    from Cython.Distutils import build_ext
except:
    print("You don't seem to have Cython installed. Please get a")
    print("copy from www.cython.org or install it with `pip install Cython`")
    sys.exit(1)

class CustomBuildExtCommand(build_ext):
    """This extension command lets us not require numpy be installed before running pip install ripser
    build_ext command for use when numpy headers are needed.
    """

    def run(self):
        # Import numpy here, only when headers are needed
        import numpy

        # Add numpy headers to include_dirs
        self.include_dirs.append(numpy.get_include())
        # Call original build_ext command
        build_ext.run(self)


#extra_compile_args = ["-Ofast", "-I/opt/homebrew/include/"]
extra_compile_args = ["-Og", "-g", "-I/opt/homebrew/include/"]
extra_link_args = ["-lfmt", "-L/opt/homebrew/lib"]

if platform.system() == "Windows":
    extra_compile_args.extend(
        [
            # Supported by Visual C++ >=14.1
            "/std:c++20"
        ]
    )
elif platform.system() == "Darwin":
    extra_compile_args.extend(["-std=c++20", "-mmacosx-version-min=10.13"])
    extra_link_args.extend(["-stdlib=libc++", "-mmacosx-version-min=10.13"])
else:
    extra_compile_args.extend(["-std=c++20"])

ext_modules = Extension(
    "pyrene",
    sources=["src/pyrene.pyx", "src/eirene_impl.cpp", "src/util.cpp"],
    extra_compile_args=extra_compile_args,
    extra_link_args=extra_link_args,
    language="c++",
)


setup(
    name="pyrene",
    #version=verstr,
    description="A Python API for a C++ port of eirene.jl.",
    #long_description=long_description,
    long_description_content_type="text/markdown",
    author="Alex Sugar",
    #author_email="chris.tralie@gmail.com, nat@riverasaul.com",
    #url="https://ripser.scikit-tda.org",
    #license="MIT",
    packages=find_packages(),#["pyrene"],
    ext_modules=cythonize(ext_modules),
    install_requires=["Cython", "numpy", "scipy"],
    cmdclass={"build_ext": CustomBuildExtCommand},
    python_requires=">=3.6",
)
