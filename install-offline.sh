# 1. Find a computer with internet access 
# install the ocp command line
# you will see a ocp folder in home
sh <(curl -s https://ocp-download.oss-cn-hongkong.aliyuncs.com/install.sh)

# Download the necessary packages
# the -d option will only download the package, because some fixes are needed for intel compiler 
# you will see a external folder in the ocp folder
ocp install -d

# 2. Now compress all everything in the ocp folder
tar -cJf ocp.tar.xz ~/ocp

# 3. Copy and decompress the file on the offline computer
tar -xJf ocp.tar.xz

# 4. Make sure it's still located in the home folder
# add the ocp cli to your path
PATH=$HOME/ocp/cli/latest:$PATH

# 5. finally, compile and install the dependencies
export OCP_CC=mpiicx
export OCP_CXX=mpiicpx
export OCP_FC=mpiifort

export OCP_COMPILER=intel
export OCP_HYPRE_DIR=~/ocp/external/hypre/2.28.0/install/intel/default
export OCP_PETSC_DIR=~/ocp/external/petsc/3.20.1/install/intel/int64
export OCP_PARMETIS_DIR=~/ocp/external/parmetis/4.0.3/install/intel/int64
export OCP_LAPACK_DIR=~/ocp/external/lapack/3.11.0/install/intel/default
export OCP_FASP_DIR=~/ocp/external/fasp/2.2.1/install/intel/default
export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:$OCP_PETSC_DIR/lib/pkgconfig:$OCP_LAPACK_DIR/lib/pkgconfig
# you will see external and toolkit folders in the ocp folder
ocp install

# 6. Now compile the OpenCAEPoroX code
export CC=$OCP_CC
export CXX=$OCP_CXX
export FC=$OCP_FC
cmake -S . -B build -DUSE_METIS=ON -DUSE_PARMETIS=ON \
    -DLAPACK_DIR=$OCP_LAPACK_DIR/lib/cmake -DMETIS_DIR=$OCP_PARMETIS_DIR  \
    -DPARMETIS_DIR=$OCP_PARMETIS_DIR -DFASP_DIR=$OCP_FASP_DIR -GNinja
cmake --build build 
