#!/bin/sh
FORTRAN_ROUTING_DIR="/content/t-route/src/fortran_routing/mc_pylink_v00/MC_singleSeg_singleTS"
FORTRAN_RES_DIR="/content/t-route/src/fortran_routing/mc_pylink_v00/Reservoir_singleTS"
FAST_REACH_DIR="/content/t-route/src/python_routing_v02/fast_reach"

cd $FORTRAN_ROUTING_DIR
gfortran varPrecision.f90 -c -O0 -fPIC
gfortran -c -O0 -fPIC -o mc_single_seg.o MCsingleSegStime_f2py_NOLOOP.f90
gfortran pyMCsingleSegStime_NoLoop.f90 -c -o pymc_single_seg.o -O3 -fPIC
cp $FORTRAN_ROUTING_DIR/*.o $FAST_REACH_DIR 

cd $FORTRAN_RES_DIR
gfortran varPrecision.f90 -c -O0 -fPIC
gfortran -c -O2 -fPIC -o module_levelpool.o module_levelpool.f90
gfortran -c -O2 -fPIC -o pymodule_levelpool.o pymodule_levelpool.f90
cp $FORTRAN_RES_DIR/*.o $FAST_REACH_DIR 

numpy_I = "/usr/local/lib/python3.6/dist-packages/numpy/core/include"
py_I = "/usr/include/python3.6"
py_lib = "/usr/include/python3.6"

cd $FAST_REACH_DIR
rm *.so
cython -3 -v -p --gdb --line-directives -Wextra --cleanup 3 fortran_wrappers.pxd *.pyx
gcc -pthread -Wno-unused-result -Wsign-compare -DNDEBUG -fwrapv -O2 -Wall -Wstrict-prototypes -fPIC -fno-strict-aliasing -I$numpy_I -I$py_I -mtune=generic -c -o fortran_wrappers.o fortran_wrappers.c
gcc -pthread -Wno-unused-result -Wsign-compare -DNDEBUG -fwrapv -O2 -Wall -Wstrict-prototypes -fPIC -fno-strict-aliasing -I$numpy_I -I$py_I -mtune=generic -c -o mc_reach.o mc_reach.c
gcc -pthread -Wno-unused-result -Wsign-compare -DNDEBUG -fwrapv -O2 -Wall -Wstrict-prototypes -fPIC -fno-strict-aliasing -I$numpy_I -I$py_I -mtune=generic -c -o reservoir.o reservoir.c
gcc -pthread -Wno-unused-result -Wsign-compare -DNDEBUG -fwrapv -O2 -Wall -Wstrict-prototypes -fPIC -fno-strict-aliasing -I$numpy_I -I$py_I -mtune=generic -c -o reach.o reach.c
gcc -pthread -Wno-unused-result -Wsign-compare -DNDEBUG -fwrapv -O2 -Wall -Wstrict-prototypes -fPIC -fno-strict-aliasing -I$numpy_I -I$py_I -mtune=generic -c -o utils.o utils.c
gcc -pthread -shared -L $py_lib -lgfortran -o reach.cpython-36m-x86_64-linux-gnu.so mc_single_seg.o pymc_single_seg.o fortran_wrappers.o reach.o
gcc -pthread -shared -L $py_lib -lgfortran -o reservoir.cpython-36m-x86_64-linux-gnu.so module_levelpool.o pymodule_levelpool.o fortran_wrappers.o reservoir.o
gcc -pthread -shared -L $py_lib -o mc_reach.cpython-36m-x86_64-linux-gnu.so mc_reach.o
gcc -pthread -shared -L $py_lib -o utils.cpython-36m-x86_64-linux-gnu.so utils.o