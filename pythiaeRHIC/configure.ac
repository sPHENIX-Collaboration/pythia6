# configure.ac
# Thomas Burton
# 15th October 2012
# configure template for pythiaeRHIC

# Run initialisation.
AC_INIT(PYTHIA 6 for eRHIC, 1.0, tpb@bnl.gov, pythiaeRHIC, https://wiki.bnl.gov/eic/index.php/PYTHIA)

# Safety check against incorrect --srcdir.
# Check for existence of a file in the source directory.
AC_CONFIG_SRCDIR(erhic/pythia6.f95)

# Don't require boilerplate files.
AM_INIT_AUTOMAKE(foreign)

# Check for all the programmes we need.
AC_PROG_CPP
AC_PROG_CXX
AC_PROG_F77(gfortran)
AC_PROG_FC(gfortran)
dnl AC_PROG_LIBTOOL
LT_INIT([disable-static])

# Determine the linker flags (e.g., -L and -l) for
# Fortran intrinsic and runtime libraries.
AC_F77_LIBRARY_LDFLAGS
AC_FC_LIBRARY_LDFLAGS

# Set compiler/linker flags.
CXXFLAGS="${CXXFLAGS}"
# Get the ROOT linker flags now via root-config.
LDFLAGS="`root-config --libs` ${LDFLAGS}"
LIBS="${LIBS}"
# Get the ROOT include flags now via root-config.
CPPFLAGS="-I`root-config --incdir` ${CPPFLAGS}"
FFLAGS="${FFLAGS}"
# The Fortran 95 code is the only Fortran entirely under our
# control in this project (most of the F77 is from PYTHIA
# or RADGEN) so we are stricter with the warnings for it than
# for Fortran 77.
# Impose fixed-format for compatibility with F77 include files,
# but ignore line length limits.
FCFLAGS="-ffixed-form -ffixed-line-length-none ${FCFLAGS}"
FLIBS="${FLIBS}"

# Allow the user to set the cernlib path.
AC_SUBST(CERNLIB)
AC_ARG_WITH(
  cern-libdir,
  AS_HELP_STRING(
    [--with-cern-libdir=DIR],
    [location of cernlib libraries (kernlib, packlib etc)]
  ),
  [CERNLIB="-L$withval -lmathlib -lkernlib -lpacklib_noshift"],
  [CERNLIB="-lmathlib -lkernlib -lpacklib_noshift"]
)

# Allow the user to set the LHAPDF paths (lib and include).
# First the library path.
AC_SUBST(LHAPDF_LIB)
AC_ARG_WITH(
  lhapdf-libdir,
  AS_HELP_STRING(
    [--with-lhapdf-libdir=DIR],
    [location of libLHAPDF]
  ),
  [LHAPDF_LIB="-L$withval -lLHAPDF"],
  [LHAPDF_LIB="-lLHAPDF"]
)
# Now the include path.
AC_SUBST(LHAPDF_INC)
AC_ARG_WITH(
  lhapdf-incdir,
  AS_HELP_STRING(
    [--with-lhapdf-incdir=DIR],
    [location of LHAPDF header directory]
  ),
  [LHAPDF_INCDIR="-I$withval -leicsmear"],
  [LHAPDF_INCDIR="-leicsmear"]
)

# Allow the user to set the eic-smear paths (lib and include).
# First the library path.
AC_SUBST(EIC_LIB)
AC_ARG_WITH(
   eic-libdir,
   AS_HELP_STRING(
      [--with-eic-libdir=DIR],
      [location of libeicsmear]
   ),
   [EIC_LIB="-L$withval -leicsmear"],
   [EIC_LIB="-leicsmear"]
)
# Now the include path.
AC_SUBST(EIC_INC)
AC_ARG_WITH(
   eic-incdir,
   AS_HELP_STRING(
      [--with-eic-incdir=DIR],
      [location of eicsmear header directory]
   ),
   [EIC_INC="-I$withval"],
   [EIC_INC=""]
)

# Make AC_OUTPUT create Makefile (from Makefile.in).
AC_CONFIG_FILES([Makefile])

AC_OUTPUT
