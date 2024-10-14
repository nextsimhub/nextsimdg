# Build apptainer container for neXtSIM-DG
These can be used to compile and run neXtSIM on HPC computers

## Build the container
The input file is `nextsim_dg.def`. We use it to build an image file `nextsim_dg.sif`.
In the `apptainer` directory, type
```
sudo apptainer build nextsim_dg.sif nextsim_dg.def
```
You can save the image file anywhere

## Run the container
To run a command inside the container, type
```
apptainer exec --cleanenv nextsim_dg.sif COMMAND
```
To open an ubuntu shell inside the container, type
```
apptainer shell --cleanenv nextsim_dg.sif
```

## Compile the code
In the root directory of `nextsimdg`, do (can change the options passed to cmake)
```
apptainer shell --cleanenv apptainer/nextsim_dg.sif
mkdir -p build
cd build
cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DWITH_THREADS=ON
make -j [NUM_JOBS]
```

# Run the model
```
cd ../run
```
Link `./nextsim` from `build`
```
ln -s ../build/nextsim
```
Run with
```
./nextsim --config-file config_benchmark.cfg
```
Can also run a more realistic case, by downloading some initial conditions
```
wget ftp://ftp.nersc.no/nextsim/netCDF/init_25km_NH.nc
```
and running with the config file
```
./nextsim --config-file config_column.cfg
```
