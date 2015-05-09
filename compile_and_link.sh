export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
python2.7 setup.py build_ext --inplace
gfortran -O2 -c splev_bispeu.f90
g++ -O2 -c -I/usr/include/python2.7 demo.cpp
g++ -O2 -o demo splev_bispeu.o build/temp.linux-x86_64-2.7/scipy_tools.o demo.o -lpython2.7

