# customLAMMPS
A collection of custom scripts to modify LAMMPS

## Contents
1. fix_activeBD
2. fix_active3D
3. pair_lj_cut_dipole2_long

## Adding code to LAMMPS
There are two ways to do this: `make` and `CMAKE`.
### Adding code to LAMMPS with `make`
- Download script and paste into the `/src/` directory of *your* copy of LAMMPS. For example:
```
cp fix_activeBD.* mylammps/src
```
 - For activeBd and active3D we need the DIPOLE package, pair_lj_cut_dipole2_long requires the KSPACE package. 
 To add these packages to your build of LAMMPS, type the following in the LAMMPS `/src/` directory.
 ```
 make yes-dipole
 make yes-kspace
 ```
 - Now we can build. For running on a single processor you can use serial, for parrallel you need mpi:
 ```
 make serial
 make mpi
 ```
 This will produce the LAMMPS executable `lmp_serial` and `lmp_mpi` in `lammps/src`.
 
 ### Adding code to LAMMPS with `CMAKE`
 - Download script and paste into the `/src/` directory of *your* copy of LAMMPS. For example:
```
cp fix_activeBD.* mylammps/src
```
- set up the build directory:
```
cd mylammps
mkdir build
cd build
```
 - For activeBd and active3D we need the DIPOLE package, pair_lj_cut_dipole2_long requires the KSPACE package.:
 ```
 cmake -D PKG_DIPOLE=on -D PKG_KSPACE=on ../cmake
 ```
 - Now build:
 ```
 cmake --build .
 ```
 - This will create the `lmp` executable in the `mylammps/build` directory. You can install this executable to your system with:
 ```
 make install
 ```
