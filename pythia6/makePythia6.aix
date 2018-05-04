# /bin/sh
xlf -O2 -qextname  -qcharlen=32767 -c pythia*.f
xlf -O2 -qextname  -qcharlen=32767 -c tpythia6_called_from_cc.F
echo 'void MAIN__() {}' >main.c
gcc -c main.c
gcc -c pythia6_common_address.c
/usr/vacpp/bin/makeC++SharedLib -o libPythia6.a -p 0 \
 main.o tpythia6_called_from_cc.o pythia*.o  -lxlf -lxlf90
