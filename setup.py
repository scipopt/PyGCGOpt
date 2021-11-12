from setuptools import setup, Extension
import os, platform, sys, re

# look for environment variable that specifies path to SCIP and GCG
scipoptdir = os.environ.get('SCIPOPTDIR', '').strip('"')
gcgoptdir = os.environ.get('GCGOPTDIR', '').strip('"')

extra_compile_args = []
extra_link_args = []

includedirs = []

for optdir in [scipoptdir, gcgoptdir]:
    # determine include directory
    if os.path.exists(os.path.join(optdir, 'src')):
        # SCIP seems to be installed in place
        includedirs.append(os.path.abspath(os.path.join(optdir, 'src')))
    else:
        # assume that SCIP is installed on the system
        includedirs.append(os.path.abspath(os.path.join(optdir, 'include')))
    # Current release of GCG has broken folder structure when installed. All GCG sources *should* be in a `gcg` subfolder in the sources. Until this is fixed, we need to include the `gcg` folder so that imports work.
    if os.path.exists(os.path.join(optdir, 'include')) and not os.path.exists(os.path.join(optdir, 'include', 'graph')):
        includedirs.append(os.path.abspath(os.path.join(optdir, 'include', 'gcg')))

includedirs = list(set(includedirs))

print('Using include path <%s>.' % ", ".join(includedirs))


# determine library
if os.path.exists(os.path.join(scipoptdir, 'lib/shared/libscipsolver.so')):
    # SCIP seems to be created with make
    sciplibdir = os.path.abspath(os.path.join(scipoptdir, 'lib/shared'))
    sciplibname = 'scipsolver'
    extra_compile_args.append('-DNO_CONFIG_HEADER')
else:
    # assume that SCIP is installed on the system
    sciplibdir = os.path.abspath(os.path.join(scipoptdir, 'lib'))
    sciplibname = 'scip'
    if platform.system() in ['Windows']:
        sciplibname = 'libscip'

gcglibdir = os.path.abspath(os.path.join(gcgoptdir, "lib/shared"))
gcglibname = "gcg"

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

if platform.system() == 'Darwin':
    extra_compile_args.append("-std=c++11")

extensions = [Extension('pygcgopt.gcg', [os.path.join(packagedir, 'gcg'+ext)],
                          include_dirs=includedirs,
                          library_dirs=list(set([sciplibdir, gcglibdir])),
                          libraries=[sciplibname, gcglibname],
                          extra_compile_args = extra_compile_args,
                          extra_link_args=extra_link_args
                          )]

if use_cython:
    extensions = cythonize(extensions, compiler_directives={'language_level': 3})

with open('README.md') as f:
    long_description = f.read()

setup(
    name='PyGCGOpt',
    version=version,
    description='Python interface and modeling environment for GCG',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://git.or.rwth-aachen.de/gcg/PyGCGOpt',
    author='Lehrstuhl für Operations Research - RWTH Aachen University',
    author_email='',
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
        'pyscipopt>=3.4.0'
    ],
    packages=['pygcgopt'],
    package_dir={'pygcgopt': packagedir},
    package_data = {'pygcgopt': ['gcg.pyx', 'gcg.pxd', '*.pxi']}
)
