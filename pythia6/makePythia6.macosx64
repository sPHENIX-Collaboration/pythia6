# /bin/sh

# Pythia6 compiled with gfortran on MacOS X >= 10.4 

# Look for a fink installation
FINKDIR=`which fink 2>&1 | sed -ne "s/\/bin\/fink//p"`
macosx_minor=`sw_vers | sed -n 's/ProductVersion://p' | cut -d . -f 2`

gfortran -m64 -c pythia*.f
gfortran -m64 -fno-second-underscore -c tpythia6_called_from_cc.F
echo 'void MAIN__() {}' > main.c
gcc -m64 -c main.c
gcc -m64 -c pythia6_common_address.c
rm -f libPythia6.dylib libPythia6.so
gcc -m64 -dynamiclib -single_module -undefined dynamic_lookup \
   -install_name $HOME/pythia6/libPythia6.dylib -o libPythia6.dylib \
   main.o tpythia6_called_from_cc.o pythia*.o \
   `gfortran -m64 -print-file-name=libgfortran.dylib` \
   `gfortran -m64 -print-file-name=libgfortran.a`
ln -s libPythia6.dylib libPythia6.so
