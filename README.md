# ABP_LAMMPS
A collection of custom scripts to modify LAMMPS for the modelling of active Brownian particles

### Contents
1. fix_activeBD
2. fix_active3D


## The models
### Active Brownian Particles in 3D (fix_active3D)
We model active colloids as active Brownian particles, which propel with a constant velocity *v<sub>p</sub>*, along their individual direction vectors **e**, which in turn are subject to rotational diffusion. 
We implement this model through molecular dynamics simulations using a customised version of the open source LAMMPS package, which integrates the following
equations of motion:


<img src="https://github.com/FJMOORE/ABP_LAMMPS/blob/ece414769a2f352a0dc00552e9575ecdbf46367a/abp3D_eom1.png" width = "400">

<img src="https://github.com/FJMOORE/ABP_LAMMPS/blob/ece414769a2f352a0dc00552e9575ecdbf46367a/abp3D_eom2.png" width = "230">

Here <img src="https://render.githubusercontent.com/render/math?math=\dot{\mathbf{r}}"> is the particle velocity, <img src="https://render.githubusercontent.com/render/math?math=v_p">  is the magnitude of the constant applied active velocity, and **F** is the inter-particle force. The thermal fluctuations promoting translational diffusion are included in the Gaussian white-noise term
<img src="https://render.githubusercontent.com/render/math?math=\boldsymbol{\xi_T}">, where <img src="https://render.githubusercontent.com/render/math?math=\langle\boldsymbol{\xi_T}\rangle=0">, and D<sub>T</sub> is the translational diffusion coefficient. Thermal noise driving rotational diffusion of the direction vector **e** is represented by <img src="https://render.githubusercontent.com/render/math?math=\boldsymbol{\xi_R}">,  where <img src="https://render.githubusercontent.com/render/math?math=\langle\boldsymbol{\xi_R}\rangle=0">, and D<sub>R</sub> is the rotational diffusion coefficient. The two diffusion coefficients are related via  <img src="https://render.githubusercontent.com/render/math?math=D_{T} = D_{R}\sigma^{2}/3">.


### Induced-charge electrophoresis (ICEP) Janus particles as active Brownian particles (fix_activeBD)
We model active ICEP Janus colloids as active Brownian particles. Assuming the external field driving particle motion is in the *z* direction, the particles undergo active motion in the *xy* plane and are diffusive in the *z* plane. For a particle being driven at a velocity <img src="https://render.githubusercontent.com/render/math?math=v_p">, the particle velocity evolves according to:

<img src="https://github.com/FJMOORE/ABP_LAMMPS/blob/ece414769a2f352a0dc00552e9575ecdbf46367a/abpJanus_eom.png" width = "300">

and the orientation of the active propulsion is given by:

<img src="https://github.com/FJMOORE/ABP_LAMMPS/blob/6f6fcb5f23c6cbad653b70d618d29218ede3c36b/abpJanus_eom_theta.png" width = "150">

___

# Adding code to LAMMPS
There are two ways to do this: `make` and `CMAKE`.
### Adding code to LAMMPS with `make`
- Download script and paste into the `/src/` directory of *your* copy of LAMMPS. For example:
```
cp fix_activeBD.* mylammps/src
```
 - For activeBd and active3D we need the DIPOLE package. 
 To add these packages to your build of LAMMPS, type the following in the LAMMPS `/src/` directory.
 ```
 make yes-dipole
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
 - For activeBD and active3D we need the DIPOLE package:
 ```
 cmake -D PKG_DIPOLE=on ../cmake
 ```
 - Now build:
 ```
 cmake --build .
 ```
 - This will create the `lmp` executable in the `mylammps/build` directory. You can install this executable to your system with:
 ```
 make install
 ```
