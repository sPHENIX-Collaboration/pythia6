# /bin/sh
gfortran -c -m64 -fPIC pythia*.f
gfortran -c -m64 -fPIC  -fno-second-underscore tpythia6_called_from_cc.F
echo 'void MAIN__() {}' >main.c
gcc -c -m64 -fPIC main.c
gcc -c -m64 -fPIC pythia6_common_address.c
gfortran -m64 -shared -Wl,-soname,libPythia6.so -o libPythia6.so \
 main.o  pythia*.o tpythia*.o
