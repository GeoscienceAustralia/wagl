# The image processor

This repository was moved (with no history) from SVN on June 28 2013. The original SVN repository can be found at.

http://www.ga.gov.au:9080/svn/gemd/neo/landsat-nbar-processor/branches/ULA3

**THAT CODE SHOULD NEVER BE USED AGAIN AND IS ONLY KEPT FOR REFERENCE.*

## Getting started

Documentation is managed using [Sphinx](http://sphinx-doc.org). The configuration can be found in the directory 'docs/source'. That directory was created running the command `sphinx-quickstart` which created (see http://sphinx-doc.org/tutorial.html):

- Makefile,
- make_doc.bat,
- source
  - conf.py,
  - index.rst,
  - _templates,
  - _static 

docs/source/conf.py contains all the settings and alike. (see http://sphinx-doc.org/config.html and http://pythonhosted.org/an_example_pypi_project/sphinx.html for details).

Sphinx autodoc needs to be able to load every module in order to parse the docstrings. This requires that all modules used in ULA3 be available. On VAYU, this can be achieved with the following cst is (at the time of writing)

    MODULEPATH=/g/data/v10/opt/modules/modulefiles:$MODULEPATH
    module load intel-fc
    module load intel-cc
    module load intel-mkl/10.3.0
    module load geotiff/1.3.0
    module load python/2.7.2
    module load hdf4/4.2.6_2012
    module load gdal/1.9.0_HDF5
    module load proj/4.7.0
    module load pyephem/3.7.5.1
    module load numexpr/2.0.1
    module load geotiff/1.3.0
    module load sphinx/1.1.3

The only trick here is that gcc is also required for compiling Fuqin's code... actually, there is a bit of system specific stuff in the Makefile in the top level of the ULA3 package which you may need to take a look at.

To generate the API documentation for ULA3, one can then run:

make ula3 # required to build fortran based modules (using F2py - which you should get with hte python module)
make .docs

in the top level directory (i.e. this directory) which will put a bunch of *.rst files in the directory 'docs/source' (this is also a dependency of `make all`). Note that Sphinx needs to load each python file for which documentation will be generated as a module. Presently, this won't work for some ULA3 modules because they attempt to instantiate an instance of ULA3.common.processor_config.ProcessorConfig in their namespace, which will fail unless certain arguments (usually passed to process.py) are present on the command line. The fix is to not generate documentation for any packages/modules where this occurs by deleting the appropropriate *.rst files (and possibly modifying others). This is done by the call to make, which first calls docs/build-sphinx.sh to generate the *.rst files and delete the ones that cause errors. Once those errors are eventually removed, someone will need to modify docs/build-sphinx.sh to generate the rest of the documentation.

Documentation will be produced in docs/build/html (the root of which is, of course, index.html). This should be copied to an appropriate place on a web server or wherever else you feel like putting it.

Note that there are many other targets in docs/Makefile, but I am not sure what all of these do or if they have other dependences that need to be present on the system. The Sphinx home page - sphinx-doc.org - mentions some of the other formats which can be generated and more information is available at http://sphinx-doc.org/latest/builders.html (the list is pretty extensive and includes man pages, windows help files and pdf).

All further documentation is generated with Sphinx... so best go generate it now!

**P.S. PLEASE ONLY CHECK DOCUMENTATION (*.rst files) THAT ARE CREATED MANUALLY INTO THE REPOSITORY** (i.e. do not check in files that are auto generated!
