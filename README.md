# OpenCAEPoroX

OpenCAEPoroX is an open-source, highly-efficient, cross-platform, and massively distributed parallel simulator for the multi-component and multi-phase fluid flows in porous media. The software employs a universal multi-component model equation and an abstract method framework design, implemented in C++. Its object-oriented programming characteristics provide excellent scalability and readability. As illustrated in the [Top-Level Architecture](#top-level-architecture-diagram) diagram, the simulator consists of five major modules:

1. **Preprocessing module**: Responsible for reading discrete grids and reservoir parameters, and completing the parallel partitioning of grids and distributed deployment of reservoir information. It handles the preprocessing work for the simulation. Supported input files include Eclipse-style input files and Gmsh grid files. The parallel partitioning part utilizes the open-source software package ParMetis.

2. **Reservoir core module**: Contains all the physical information of the reservoir, corresponding models, and solving modules, such as the PVT module, phase permeation module, rock module, flow module, and well module, etc. It is the core component of the simulator.

3. **Solver module**: Contains solvers for different problems, such as isothermal solver and non-isothermal solver. These solvers are responsible for providing and scheduling solving methods, such as semi-implicit methods, fully implicit methods, adaptive implicit methods, etc. Each method collaborates with various submodules in the reservoir core module to complete the simulation process.

4. **Solving control module**: Responsible for controlling the solving process, including the selection of time steps, solving nonlinear iterations, and convergence control.

5. **Postprocessing module**: Responsible for outputting the required results for subsequent visualization and data analysis. Currently, it can output VTK format files describing the dynamic changes of all reservoir grid physical information, Eclipse-style "SUMMARY" files reflecting reservoir overall state and well performance, "FastReview" files for runtime solving analysis, and "statistic" files reflecting parallel efficiency and load balancing.

#### Top-Level Architecture diagram
![simulation_framework](docs/img/simulation_framework.png)

Key features of the simulator include:

- Supported Eclipse-style input files and Gmsh grid files.
- Various discrete methods: semi-implicit methods, fully implicit methods, adaptive implicit methods.
- Scalar or coupled block systems.
- Modular structure for easy implementation of your own methods.
- Linux, Windows, and macOS support.

---

# Table of Contents

- [Quickstart](#quickstart)
  - [Requirements and Dependencies](#requirements-and-dependencies)
  - [Build Dependencies](#build-dependencies)
  - [Build OpenCAEPoroX](#build-opencaeporox)
  - [Running examples](#running-examples)
- [License](#license)

---

# Quickstart

Welcome to the OpenCAEPoroX installation guide. Here are the instructions on how to build OpenCAEPoroX and run a benchmark example SPE1 provided by the Society of Petroleum Engineers.

---

## Requirements and Dependencies

### 1. Requirements

Before starting the installation process, make sure you have the following system requirements:

- A supported operating system (e.g., Linux, Windows, and macOS).
- C++ compiler (e.g., [GCC](https://gcc.gnu.org/) 7.3.0).
- [CMake](https://cmake.org/) 3.17 or later.
- MPI compiler (e.g., [OpenMPI](https://www.open-mpi.org/) for Linux or [MPICH](https://www.mpich.org/downloads/) for Windows) or Intel compiler (e.g., [oneAPI](https://www.intel.com/content/www/us/en/developer/tools/oneapi/toolkits.html#base-kit)). The latter is recommended for optimal performance.

### 2. Dependencies

OpenCAEPoroX requires the following external open-source libraries:

- [lapack-3.11](https://netlib.org/lapack/lapack-3.11.0.html)
- [parmetis-4.0.3](http://glaros.dtc.umn.edu/gkhome/metis/parmetis/download)
- [hypre-2.28.0](https://github.com/hypre-space/hypre)
- [petsc-3.19.3](https://web.cels.anl.gov/projects/petsc/vault/petsc-3.19/docs/install/download.html)
- [petsc_solver]() (custom solver built on top of PETSc)

Ensure that these libraries are installed and accessible before proceeding with the installation of OpenCAEPoroX.

---

## Build Dependencies

Build the required dependencies in the following order:

### 1. lapack-3.11

```
   cd lapack-3.11
   cp make.inc.example make.inc
   make blaslib -j 16
   make cblaslib -j 16
   make lapacklib -j 16
   make lapackelib -j 16
```

### 2. parmetis-4.0.3

```
   cd parmetis-4.0.3
   make config cc=mpiicc prefix=ROOT_DIR/parmetis-4.0.3
   make -j 16
   make install
```

### 3. hypre-2.28.0

```
   cd hypre-2.28.0
   ./configure --prefix=ROOT_DIR/hypre-2.28.0 --with-MPI --enable-shared
   make -j 16
   make install   
```

### 4. petsc-3.19.3

```
   cd petsc-3.19.3

   export PETSC_DIR=ROOT_DIR/petsc-3.19.3
   export PETSC_ARCH=petsc_install

   ./configure CC=mpiicc CXX=mpiicpc \
      --with-fortran-bindings=0 \
      --with-hypre-dir=ROOT_DIR/hypre-2.28.0 \ 
      --with-debugging=0 \
      COPTFLAGS="-O3" \
      CXXOPTFLAGS="-O3" \

   make -j 20 PETSC_DIR=ROOT_DIR/petsc-3.19.3 PETSC_ARCH=petsc_install all
   make all check 
```

### 5. petsc_solver

```
   cd petsc_solver

   export CC=mpiicc
   export CXX=mpiicpc

   export CPATH=ROOT_DIR/lapack-3.11/CBLAS/include:ROOT_DIR/lapack-3.11/LAPACKE/include:$CPATH
   export LD_LIBRARY_PATH=ROOT_DIR/lapack-3.11:$LD_LIBRARY_PATH

   mkdir build
   cd build
   cmake ..
   make
```

Please note that this is a general outline, and you may need to adjust paths and commands based on your specific system configuration and directory structure. Make sure to replace placeholders such as `ROOT_DIR` with the actual root directory set by the user.

---

## Build OpenCAEPoroX

Once the dependencies are built, navigate to the OpenCAEPoroX directory and build library:

```
cd OpenCAEPoro

# users specific compilers
export CC=mpiicc
export CXX=mpiicpc

# users specific directory paths
export PARMETIS_DIR=ROOT_DIR/parmetis-4.0.3
export PARMETIS_BUILD_DIR=ROOT_DIR/parmetis-4.0.3/build/Linux-x86_64
export METIS_DIR=ROOT_DIR/parmetis-4.0.3/metis
export METIS_BUILD_DIR=ROOT_DIR/parmetis-4.0.3/build/Linux-x86_64
export PETSC_DIR=ROOT_DIR/petsc-3.19.3
export PETSC_ARCH=petsc_install
export PETSCSOLVER_DIR=ROOT_DIR/petsc_solver

export CPATH=ROOT_DIR/petsc-3.19.3/include/:$CPATH
export CPATH=ROOT_DIR/petsc-3.19.3/petsc_install/include/:ROOT_DIR/parmetis-4.0.3/metis/include:ROOT_DIR/parmetis-4.0.3/include:$CPATH
export CPATH=ROOT_DIR/lapack-3.11/CBLAS/include/:$CPATH

mkdir build
cd build

cmake -DUSE_PARDISO=ON -DUSE_PETSCSOLVER=ON -DUSE_PARMETIS=ON -DUSE_METIS=ON -DCMAKE_BUILD_TYPE=Release ..

make -j 16
make install
```

---

## Running examples

After installation, you can test the setup by running the following command in the terminal under the OpenCAEPoroX main directory:

```
mpirun -n p ./testOpenCAEPoro ./data/spe1a/spe1a.data
```

Replace `p` with the number of processes you want to use. Check the output on the screen and the newly generated files in `./data/spe1a/`. You should find files like `SUMMARY.out` and `FastReview.out`. If you run this command with more than one process, `statistics.out` will also be generated.

---

# License

This software is free software distributed under the Lesser General Public License or LGPL, version 3.0 or any later versions. This software distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with OpenCAEPoroX. If not, see [http://www.gnu.org/licenses/](http://www.gnu.org/licenses/).
