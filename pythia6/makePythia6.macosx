# /bin/sh
g77 -c pythia*.f
g77 -fno-second-underscore -c tpythia6_called_from_cc.F
echo 'void MAIN__() {}' > main.c
gcc -c main.c
gcc -c pythia6_common_address.c
g77 -v -bundle -flat_namespace -undefined dynamic_lookup -o libPythia6.so \
 main.o tpythia6_called_from_cc.o pythia*.o
g77 -dynamiclib -flat_namespace -single_module -undefined dynamic_lookup \
 -install_name $HOME/pythia6/libPythia6.dylib -o libPythia6.dylib \
 main.o tpythia6_called_from_cc.o pythia*.o
rm -f libPythia6.a
ar rv libPythia6.a main.o tpythia6_called_from_cc.o pythia*.o
ranlib libPythia6.a
