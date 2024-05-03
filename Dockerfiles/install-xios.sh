#!/usr/bin/env bash

# use svn to obtain current version of xios
cd /
installdir="xios"
svn checkout http://forge.ipsl.jussieu.fr/ioserver/svn/XIOS/trunk $installdir
cd $installdir

cat <<EOF > arch/arch-GCC_LINUX.path
NETCDF_INCDIR="-I \$NETCDF_INC_DIR -I \$NETCDFFORT_INC_DIR"
NETCDF_LIBDIR="-L \$NETCDF_LIB_DIR -L \$NETCDFFORT_LIB_DIR"
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

cat <<EOF > arch/arch-GCC_LINUX.env
export HDF5_INC_DIR=\$(pkg-config --variable=prefix hdf5)/include
export HDF5_LIB_DIR=\$(pkg-config --variable=prefix hdf5)/lib

export NETCDF_INC_DIR=\$(pkg-config --variable=prefix netcdf)/include
export NETCDF_LIB_DIR=\$(pkg-config --variable=prefix netcdf)/lib

export NETCDFFORT_INC_DIR=\$(pkg-config --variable=prefix netcdf-fortran)/include
export NETCDFFORT_LIB_DIR=\$(pkg-config --variable=prefix netcdf-fortran)/lib

export BOOST_INC_DIR=\$HOME/boost
export BOOST_LIB_DIR=\$HOME/boost
EOF

cat <<EOF > arch/arch-GCC_LINUX.fcm
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

%BASE_FFLAGS    -D__NONE__
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
rm -r /xios/obj /xios/bin/generic_testcase.exe /xios/src /xios/tools \
    /xios/inputs /xios/doc /xios/arch /xios/xios_test_suite /xios/flags \
    /xios/generic_testcase /xios/ppsrc /xios/done
