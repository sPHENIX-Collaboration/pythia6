g77 -O -c -foptimize-sibling-calls -finline-functions  *.f 
g77 -O -c  -fno-second-underscore -foptimize-sibling-calls -finline-functions  *.F 
echo 'void MAIN__() {}' >main.c
gcc -O -c -foptimize-sibling-calls -finline-functions  *.c 
$ROOTSYS/build/unix/wingcc_ld.sh -shared \
-Wl,--export-all-symbols -Wl,-soname=libPythia6.dll -O \
-foptimize-sibling-calls -finline-functions \
-Wl,--enable-auto-import -Wl,--enable-runtime-pseudo-reloc \
main.o  pythia*.o tpythia6_called_from_cc.o \
-o libPythia6.dll  -lg2c -lm  
ln -s libPythia6.dll libPythia6.dll.a


