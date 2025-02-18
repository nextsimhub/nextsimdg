#!/bin/env bash

set -eu

# Checkout (clone) current version of XIOS using Subversion (svn)
INSTALLDIR="/xios"
if [ ! -d "${INSTALLDIR}" ]; then
  cd /
  svn checkout http://forge.ipsl.fr/ioserver/svn/XIOS3/trunk ${INSTALLDIR}
fi
cd ${INSTALLDIR}

# Create a file for exporting paths to the include and lib directories for XIOS' dependencies
cat <<EOF >arch/arch-GCC_LINUX.env
export HDF5_INC_DIR=\$(pkg-config --variable=prefix hdf5)/include
export HDF5_LIB_DIR=\$(pkg-config --variable=prefix hdf5)/lib

export NETCDF_INC_DIR=\$(pkg-config --variable=prefix netcdf)/include
export NETCDF_LIB_DIR=\$(pkg-config --variable=prefix netcdf)/lib

export NETCDFFORT_INC_DIR=\$(pkg-config --variable=prefix netcdf-fortran)/include
export NETCDFFORT_LIB_DIR=\$(pkg-config --variable=prefix netcdf-fortran)/lib

export BOOST_INC_DIR=\$HOME/boost
export BOOST_LIB_DIR=\$HOME/boost
EOF

# Create a file for defining compiler and linking flags for building XIOS
cat <<EOF >arch/arch-GCC_LINUX.path
NETCDF_INCDIR="\$(pkg-config --cflags-only-I netcdf-fortran)"
NETCDF_LIBDIR="\$(pkg-config --libs-only-L netcdf) \$(pkg-config --libs-only-L netcdf-fortran)"
NETCDF_LIB="-lnetcdff -lnetcdf"

MPI_INCDIR=""
MPI_LIBDIR="\$(pkg-config --libs-only-L libcurl)"
MPI_LIB="-lcurl"

HDF5_INCDIR="-I \$HDF5_INC_DIR"
HDF5_LIBDIR="-L \$HDF5_LIB_DIR \$(pkg-config --libs-only-L zlib)"
HDF5_LIB="-lhdf5_hl -lhdf5 -lhdf5 -lz"

BOOST_INCDIR="-I \$BOOST_INC_DIR"
BOOST_LIBDIR="-L \$BOOST_LIB_DIR"
BOOST_LIB=""

OASIS_INCDIR="-I\$PWD/../../oasis3-mct/BLD/build/lib/psmile.MPI1"
OASIS_LIBDIR="-L\$PWD/../../oasis3-mct/BLD/lib"
OASIS_LIB="-lpsmile.MPI1 -lscrip -lmct -lmpeu"
EOF

# Create a file containing parameters to be passed to the fcm-make build system
cat <<EOF >arch/arch-GCC_LINUX.fcm
################################################################################
###################                Projet XIOS               ###################
################################################################################

%CCOMPILER      mpicc
%FCOMPILER      mpif90
%LINKER         mpif90

%BASE_CFLAGS    -w -std=c++11 -D__XIOS_EXCEPTION
%PROD_CFLAGS    -fPIC -O3 -DBOOST_DISABLE_ASSERTS
%DEV_CFLAGS     -g -O2
%DEBUG_CFLAGS   -fPIC -g

%BASE_FFLAGS    -D__NONE__ -ffree-line-length-312
%PROD_FFLAGS    -fPIC -O3
%DEV_FFLAGS     -g -O2
%DEBUG_FFLAGS   -fPIC -g

%BASE_INC       -D__NONE__
%BASE_LD        -lstdc++

%CPP            cpp
%FPP            cpp -P
%MAKE           gmake
EOF

./make_xios --arch GCC_LINUX --job 8 --full --debug
rm -r \
  ${INSTALL_DIR}/arch \
  ${INSTALL_DIR}/bin/generic_testcase.exe \
  ${INSTALL_DIR}/doc \
  ${INSTALL_DIR}/done \
  ${INSTALL_DIR}/flags \
  ${INSTALL_DIR}/generic_testcase \
  ${INSTALL_DIR}/inputs \
  ${INSTALL_DIR}/obj \
  ${INSTALL_DIR}/ppsrc \
  ${INSTALL_DIR}/tools \
  ${INSTALL_DIR}/xios_test_suite
