# WaComM
WaComM: A parallel Water quality Community Model for pollutant transport and dispersion operational predictions

# Compile

1. Create the wacomm directory and its subs

It is assumed that WaComM will be installed in the user's home directory.

cd $HOME
mkdir wacomm
cd wacomm
mkdir ext
mkdir bin
mkdir etc
mkdir dist

1. Clone the repository

git clone https://github.com/CCMMMA/wacomm.git

2. Create the profile script

In this example the compiler is set as PGI.

cat > etc/profile << EOF
export COMPILER_VERSION="pgi"
export WACOMM_ROOT="$HOME/wacomm"
export EXT_ROOT="$WACOMM_ROOT/ext/$COMPILER_VERSION"

export PATH=$EXT_ROOT/bin:$PATH
export LD_LIBRARY_PATH=$EXT_ROOT/lib:$LD_LIBRARY_PATH

export CPP=cpp
export CC=pgcc
export CFLAGS="-O -fPIC -m64 -tp=px"
export CXX=pgc++
export CXXFLAGS="-O -fPIC -m64 -tp=px"
export FC=pgfortran
export FCFLAGS="-O -fPIC -m64 -tp=px"
export F77=pgfortran
export FFLAGS="-O -fPIC -m64 -tp=px"
export CPPFLAGS="-DpgiFortran -I$EXT_ROOT/include"
export LDFLAGS="-L$EXT_ROOT/lib"

export HDF5_HOME=$EXT_ROOT
export HDF5=$EXT_ROOT
export NETCDF=$EXT_ROOT
export NETCDF_FORTRAN=$EXT_ROOT

export HDF5_DIR=$HDF5_HOME
export HDF5_INCDIR=$HDF5_HOME/include/
export HDF5_LIBDIR=$HDF5_HOME/lib/
EOF

3. Set the environment

source ext/profile

4. Download the dependencies

cd $WACOMM_HOME/dist
wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.2/src/hdf5-1.10.2.tar.gz
wget ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4.6.1.tar.gz
wget ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-fortran-4.4.4.tar.gz

5. Compile the dependencies

cd $WACOMM_HOME/dist
tar -xvzf hdf5-1.10.2.tar.gz
cd hdf5-1.10.2
./configure --prefix=$HDF5_HOME --enable-hl --enable-shared --enable-fortran --enable-build-mode=production --enable-unsupported --enable-cxx --with-zlib --with-szlib --enable-threadsafe --with-pthread
make
make install
cd ..
tar -xvzf netcdf-4.6.1.tar.gz
cd netcdf-4.6.1
export CPPFLAGS="-I$HDF5_HOME/include"
export LDFLAGS="-L$HDF5_HOME/lib"
./configure --prefix=$NETCDF --enable-netcdf-4 --enable-shared  --enable-dap --with-hdf5=$HDF5_HOME
make
make install
cd ..
tar -xvzf netcdf-fortran-4.4.4.tar.gz
cd netcdf-fortran-4.4.4
export CPPFLAGS="-I$NETCDF/include"
export LDFLAGS="-L$NETCDF/lib"
 ./configure --prefix=$NETCDF_FORTRAN
 
# 6. Compile WaComM

cd $WACOMM_HOME/wacomm
cd src
make COMPILER=pgi

# Cite WaComm

Montella, Raffaele, Alison Brizius, Diana Di Luccio, Cheryl Porter, Joshua Elliot, Ravi Madduri, David Kelly, Angelo Riccio, and Ian Foster. "Using the FACE-IT portal and workflow engine for operational food quality prediction and assessment: An application to mussel farms monitoring in the Bay of Napoli, Italy." Future Generation Computer Systems (2018).
https://www.sciencedirect.com/science/article/pii/S0167739X16308305

Galletti, Ardelio, Raffaele Montella, Livia Marcellino, Angelo Riccio, Diana Di Luccio, Alison Brizius, and Ian T. Foster. "Numerical and Implementation Issues in Food Quality Modeling for Human Diseases Prevention." In HEALTHINF, pp. 526-534. 2017.

Di Luccio, Diana, Ardelio Galletti, Livia Marcellino, Angelo Riccio, Raffaele Montella, and Alison Brizius. "Some remarks about a community open source Lagrangian pollutant transport and dispersion model." Procedia Computer Science 113 (2017): 490-495.
https://www.sciencedirect.com/science/article/pii/S1877050917317180

Montella, Raffaele, Diana Di Luccio, Pasquale Troiano, Angelo Riccio, Alison Brizius, and Ian Foster. "WaComM: A parallel Water quality Community Model for pollutant transport and dispersion operational predictions." In Signal-Image Technology & Internet-Based Systems (SITIS), 2016 12th International Conference on, pp. 717-724. IEEE, 2016.
https://ieeexplore.ieee.org/abstract/document/7907547/

Montella, Raffaele, Alison Brizius, Diana Di Luccio, Cheryl Porter, Joshua Elliot, Ravi Madduri, David Kelly, Angelo Riccio, and Ian Foster. "Applications of the FACE-IT portal and workflow engine for operational food quality prediction and assessment: Mussel farm monitoring in the Bay of Napoli, Italy." (2016).
http://ceur-ws.org/Vol-1800/short3.pdf
