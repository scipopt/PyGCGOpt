from setuptools import setup, Extension
import os, platform, sys, re

# look for environment variable that specifies path to SCIP and GCG
scipoptdir = os.environ.get('SCIPOPTDIR', '').strip('"')
gcgoptdir = os.environ.get('GCGOPTDIR', scipoptdir).strip('"')

extra_compile_args = []
extra_link_args = []

includedirs = []

for optdir in set([scipoptdir, gcgoptdir]):
    # determine include directory
    if os.path.exists(os.path.join(optdir, 'src')):
        # SCIP seems to be installed in place
        includedirs.append(os.path.abspath(os.path.join(optdir, 'src')))
    else:
        # assume that SCIP is installed on the system
        includedirs.append(os.path.abspath(os.path.join(optdir, 'include')))

if not gcgoptdir:
    if platform.system() == 'Linux':
        includedirs.append("/usr/include/gcg")

includedirs = list(set(includedirs))

print('Using include path <%s>.' % ", ".join(includedirs))


# determine scip library
if os.path.exists(os.path.join(scipoptdir, "lib", "shared", "libscip.so")):
    # SCIP seems to be created with make
    sciplibdir = os.path.abspath(os.path.join(scipoptdir, "lib", "shared"))
    sciplibname = "scip"
    extra_compile_args.append("-DNO_CONFIG_HEADER")
    # the following is a temporary hack to make it compile with SCIP/make:
    extra_compile_args.append("-DTPI_NONE")  # if other TPIs are used, please modify
else:
    # assume that SCIP is installed on the system
    sciplibdir = os.path.abspath(os.path.join(scipoptdir, "lib"))
    sciplibname = "libscip" if platform.system() in ["Windows"] else "scip"

# determine gcg library
if os.path.exists(os.path.join(gcgoptdir, "lib", "shared", "libgcg.so")):
    # SCIP seems to be created with make
    gcglibdir = os.path.abspath(os.path.join(gcgoptdir, "lib", "shared"))
    gcglibname = "gcg"
else:
    # assume that SCIP is installed on the system
    gcglibdir = os.path.abspath(os.path.join(gcgoptdir, "lib"))
    gcglibname = "libgcg" if platform.system() in ["Windows"] else "gcg"

print('Using SCIP library <%s> at <%s>.' % (sciplibname, sciplibdir))
print('Using GCG library <%s> at <%s>.' % (gcglibname, gcglibdir))

# set runtime libraries
if platform.system() in ['Linux', 'Darwin']:
    for libdir in set([sciplibdir, gcglibdir]):
        extra_link_args.append('-Wl,-rpath,'+libdir)

# enable debug mode if requested
if "--debug" in sys.argv:
    extra_compile_args.append('-UNDEBUG')
    sys.argv.remove("--debug")

use_cython = True

packagedir = os.path.join('src', 'pygcgopt')

with open(os.path.join(packagedir, '__init__.py'), 'r') as initfile:
    version = re.search(r'^__version__\s*=\s*[\'"]([^\'"]*)[\'"]',
                        initfile.read(), re.MULTILINE).group(1)

try:
    from Cython.Build import cythonize
except ImportError:
    if not os.path.exists(os.path.join(packagedir, 'gcg.cpp')):
        print('Cython is required')
        quit(1)
    use_cython = False

if not os.path.exists(os.path.join(packagedir, 'gcg.pyx')):
    use_cython = False

ext = '.pyx' if use_cython else '.cpp'

extra_compile_args.append("-std=c++11")

extensions = [Extension('pygcgopt.gcg', [os.path.join(packagedir, 'gcg'+ext)],
                          include_dirs=includedirs,
                          library_dirs=list(set([sciplibdir, gcglibdir])),
                          libraries=[sciplibname, gcglibname],
                          extra_compile_args = extra_compile_args,
                          extra_link_args=extra_link_args
                          )]

if use_cython:
    # Compiler directives needed for documentation, see https://stackoverflow.com/a/10060115/11362041
    extensions = cythonize(extensions, compiler_directives={'language_level': 3, 'embedsignature': True})

with open('README.md') as f:
    long_description = f.read()

setup(
    name='PyGCGOpt',
    version=version,
    description='Python interface and modeling environment for GCG',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/scipopt/PyGCGOpt',
    author='Lehrstuhl für Operations Research - RWTH Aachen University',
    author_email='gcg-bugs@or.rwth-aachen.de',
    license='MIT',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Education',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Cython',
        'Topic :: Scientific/Engineering :: Mathematics'],
    ext_modules=extensions,
    install_requires=[
        'wheel',
        'pyscipopt>=5.0.0'
    ],
    packages=['pygcgopt'],
    package_dir={'pygcgopt': packagedir},
    package_data = {'pygcgopt': ['gcg.pyx', 'gcg.pxd', '*.pxi']}
)
