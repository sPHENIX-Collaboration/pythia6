# /bin/sh
f77 -c pythia*.f
f77 -c tpythia6_called_from_cc.F
echo 'void MAIN__() {}' >main.c
cc -c main.c
cc -c pythia6_common_address.c
ld -b -shared -o libPythia6.so \
 main.o tpythia6_called_from_cc.o pythia*.o \
 -L/usr/shlib  -lfor -lUfor  -lc -lm
