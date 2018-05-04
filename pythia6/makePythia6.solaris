# /bin/sh
/opt/SUNWspro/bin/f77 -c -PIC pythia*.f
/opt/SUNWspro/bin/f77 -c -PIC tpythia6_called_from_cc.F
echo 'void MAIN__() {}' >main.c
gcc -c main.c
gcc -c pythia6_common_address.c
/opt/SUNWspro/bin/f77 -G -o libPythia6.so  \
 main.o tpythia6_called_from_cc.o pythia*.o  -lsunmath -lfsu

