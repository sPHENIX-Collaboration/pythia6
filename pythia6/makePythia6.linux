# /bin/sh
g77 -c -fPIC pythia*.f
g77 -c -fPIC  -fno-second-underscore tpythia6_called_from_cc.F
echo 'void MAIN__() {}' >main.c
gcc -c -fPIC main.c
gcc -c -fPIC pythia6_common_address.c
g77 -shared -Wl,-soname,libPythia6.so -o libPythia6.so main.o  pythia*.o tpythia*.o
