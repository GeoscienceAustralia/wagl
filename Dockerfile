# Note, this file excludes any AWS dependencies as it's designed to solely run
# up the most basic version of WAGL (unit tests, for instance - but could be used
# for development as well)
FROM ubuntu:focal
SHELL ["/bin/bash", "-c"]

ENV BUILD_DIR=/build
ENV PATH="${PATH}:${BUILD_DIR}/conda/bin"
ENV PYTHONPATH=${BUILD_DIR}/conda/lib/python3.8/site-packages/
ENV WAGL_DIR=${BUILD_DIR}/wagl

USER root

# Build deps
RUN apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y --fix-missing --no-install-recommends \
        git bzip2 ca-certificates gfortran-10 gcc-10 make software-properties-common libpq-dev

RUN ln -s $(which gfortran-10) $(which gfortran-10 | sed 's/\(.*\)\/\gfortran-10/\1\/gfortran/') \
    && ln -s $(which gcc-10) $(which gcc-10 | sed 's/\(.*\)\/\gcc-10/\1\/gcc/')

WORKDIR ${BUILD_DIR}

# Bump this when newer versions of python are required
ADD https://repo.continuum.io/miniconda/Miniconda3-py38_4.8.2-Linux-x86_64.sh /root/miniconda.sh
RUN chmod +x /root/miniconda.sh && /root/miniconda.sh -b -f -p conda

# GDAL 3.1 is being used because https://gdal.org/api/python.html#usage
RUN conda install -c conda-forge \
        gdal==3.1.4 \
        python-fmask==0.5.5

WORKDIR ${WAGL_DIR}
ADD . ./

# Install dependencies required for unittests
RUN pip install -r requirements.txt
RUN pip install git+https://github.com/sixy6e/idl-functions.git#egg=master

# include basic details for diagnostics
RUN which python; python --version;

# Warning: "setup.py test" is deprecated:
# https://setuptools.readthedocs.io/en/latest/userguide/commands.html#test-build-package-and-run-a-unittest-suite
# "setup.py test" creates a local build without installing it
#
# TODO: replace with a cleaner build system + pytest
RUN time python setup.py test
