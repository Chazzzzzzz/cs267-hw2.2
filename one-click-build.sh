rm -rf build  
module load cmake
module swap PrgEnv-intel PrgEnv-gnu
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
