spack env activate -p nextsim

source ~/.envs/nextsim/bin/activate # TODO change this depending on where your nextsim python venv is located
XIOS_DIR="/home/tdm39/nextsimdg/Dockerfiles/XIOS" # TODO Change this to your XIOS install path

rm -rf build
mkdir -p build
cd build

cmake .. \
    -DCMAKE_BUILD_TYPE=Debug \
    -DCMAKE_C_COMPILER=$(which mpiicc) \
    -DCMAKE_CXX_COMPILER=$(which mpiicpc) \
    -DCMAKE_Fortran_COMPILER=$(which mpiifort) \
    -DENABLE_XIOS=ON \
    -Dxios_DIR="${XIOS_DIR}" \
    -DENABLE_MPI=ON

make -j 12
