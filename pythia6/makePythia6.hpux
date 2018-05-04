# /bin/sh
f77 -c  +Z +ppu pythia*.f
f77 -c  +Z +ppu tpythia6_called_from_cc.F
echo 'void MAIN__() {}' >main.c
CC -c main.c
CC -c pythia6_common_address.c
ld -b -z -o libPythia6.so \
 main.o tpythia6_called_from_cc.o pythia*.o
chmod 555 libPythia6.so  
