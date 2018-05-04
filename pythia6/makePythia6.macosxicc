# /bin/sh
ifort -c pythia*.f
ifort -fno-second-underscore -c tpythia6_called_from_cc.F
echo 'void MAIN__() {}' > main.c
icc -c main.c
icc -c pythia6_common_address.c
icc -dynamiclib -undefined dynamic_lookup -single_module  \
-o libPythia6.dylib  main.o tpythia6_called_from_cc.o pythia*.o \
-L/opt/intel/fc/9.1.018/lib -lifcore -lirc -limf -lifport
ln -sf  libPythia6.dylib libPythia6.so
