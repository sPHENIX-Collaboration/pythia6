# pythiaeRHIC


--- Overview ---

This is the eRHIC implementation of PYTHIA, with radiative
corrections by Radgen. Output is generated as both an ASCII
file and a ROOT file (using the eic-smear EventPythia structure).
The ROOT file is equivalent to running the EIC BuildTree routine
on the ASCII file.

--- Dependencies ---

1) cernlib: http://cernlib.web.cern.ch/cernlib/
2) LHAPDF: http://lhapdf.hepforge.org/
3) ROOT: http://root.cern.ch/drupal/
4) eic-smear: https://wiki.bnl.gov/eic/index.php/Eic-smear

--- Building ---

pythiaeRHIC is built using configure and make.
Autotools must be available on your system.

1) In the top-level source directory (containing configure.ac) run
      autoreconf -f -i
2) In a build directory (this does not need to be the top-level
   source directory) run
      /path/to/configure [arguments]
   To see a full list of arguments run configure --help.
   Note particularly the arguments to specify the locations of
   external packages, all of which are required by pythiaeRHIC:
      --with-cern-libdir
      --with-lhapdf-incdir
      --with-lhapdf-libdir
      --with-eic-incdir
      --with-eic-libdir
3) To build pythiaeRHIC run
      make

A note on Fortran compiler:

A modern Fortran compiler (like gfortran) is recommended.
If using an older compiler additional flags may need to be
passed to configure. For example, if using g77 you may need
to pass -fno-second-underscore for successful linking:

configure [--with...] FFLAGS=-fno-second-underscore

--- Running ---

To run the program, use

pythiaeRHIC < input

where 'input' is a configuration file. See the example input
bundled with the distribution for an example.

A ROOT file is always generated. The name can be specified via
the option --out=<file name>. ASCII output is also enabled
by default, but can be suppressed with the option --noascii.
Note that the ASCII file name is set via the input file.

For example:

pythiaeRHIC --out=pythia.root --noascii < input

will generate a ROOT file named "pythia.root" and no ASCII file,
while

pythiaeRHIC --out=pythia.root < input

will generate the ROOT file, plus an ASCII file, named via input.
