# /bin/sh
f77 -c pythia*.f
f77 -c tpythia6_called_from_cc.F
echo 'void MAIN__() {}' >main.c
CC -c main.c
CC -c pythia6_common_address.c
CC -Wl,-soname,libPythia.so -shared -g -o libPythia6.so \
 main.o tpythia6_called_from_cc.o pythia*.o -lftn
