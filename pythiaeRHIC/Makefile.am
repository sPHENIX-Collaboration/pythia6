# Makefile.am
# Thomas Burton
# 15th October 2012

# Generate PYTHIA 6 and Radgen libraries.
lib_LTLIBRARIES = \
libpythia6eRHIC.la \
libradgen.la

# Generate pythiaeRHIC executable.
bin_PROGRAMS = pythiaeRHIC

# PYTHIA 6 source files (all Fortran 77).
libpythia6eRHIC_la_SOURCES = \
pythia/fhhiggscorr.f \
pythia/fhsetflags.f \
pythia/fhsetpara.f \
pythia/py1ent.f \
pythia/py2ent.f \
pythia/py2frm.f \
pythia/py3ent.f \
pythia/py4ent.f \
pythia/py4frm.f \
pythia/py4jet.f \
pythia/py4jts.f \
pythia/py4jtw.f \
pythia/py6frm.f \
pythia/pyadsh.f \
pythia/pyalem.f \
pythia/pyalps.f \
pythia/pyangl.f \
pythia/pyapps.f \
pythia/pybesq.f \
pythia/pybksb.f \
pythia/pyboei.f \
pythia/pybook.f \
pythia/pycba2.f \
pythia/pycbal.f \
pythia/pycdiv.f \
pythia/pycell.f \
pythia/pychge.f \
pythia/pycjdc.f \
pythia/pyckbd.f \
pythia/pyclus.f \
pythia/pycmq2.f \
pythia/pycmqr.f \
pythia/pycomp.f \
pythia/pycrth.f \
pythia/pycsrt.f \
pythia/pyct5l.f \
pythia/pyct5m.f \
pythia/pycteq.f \
pythia/pycttr.f \
pythia/pydata.f \
pythia/pydcyk.f \
pythia/pydecy.f \
pythia/pydiff.f \
pythia/pydisg.f \
pythia/pydocu.f \
pythia/pydump.f \
pythia/pyedit.f \
pythia/pyeevt.f \
pythia/pyeicg.f \
pythia/pyeig4.f \
pythia/pyerrm.f \
pythia/pyevnt.f \
pythia/pyevnw.f \
pythia/pyevol.f \
pythia/pyevwt.f \
pythia/pyexec.f \
pythia/pyfact.f \
pythia/pyfcmp.f \
pythia/pyfeyn.f \
pythia/pyfill.f \
pythia/pyfint.f \
pythia/pyfisb.f \
pythia/pyfowo.f \
pythia/pyfram.f \
pythia/pyfscr.f \
pythia/pygaga.f \
pythia/pygamm.f \
pythia/pygano.f \
pythia/pygau2.f \
pythia/pygaus.f \
pythia/pygbeh.f \
pythia/pygdir.f \
pythia/pygfxx.f \
pythia/pyggam.f \
pythia/pygive.f \
pythia/pyglui.f \
pythia/pygrvd.f \
pythia/pygrvl.f \
pythia/pygrvm.f \
pythia/pygrvs.f \
pythia/pygrvv.f \
pythia/pygrvw.f \
pythia/pygvmd.f \
pythia/pyh2xx.f \
pythia/pyhepc.f \
pythia/pyhext.f \
pythia/pyhfth.f \
pythia/pyhggm.f \
pythia/pyhist.f \
pythia/pyi3au.f \
pythia/pyinbm.f \
pythia/pyindf.f \
pythia/pyinit.f \
pythia/pyinki.f \
pythia/pyinom.f \
pythia/pyinpr.f \
pythia/pyinre.f \
pythia/pyjmas.f \
pythia/pyjoin.f \
pythia/pyjurf.f \
pythia/pyk.f \
pythia/pykcut.f \
pythia/pykfdi.f \
pythia/pykfin.f \
pythia/pyklim.f \
pythia/pykmap.f \
pythia/pylamf.f \
pythia/pyldcm.f \
pythia/pylhef.f \
pythia/pylist.f \
pythia/pylogo.f \
pythia/pymael.f \
pythia/pymass.f \
pythia/pymaxi.f \
pythia/pymemx.f \
pythia/pymewt.f \
pythia/pymign.f \
pythia/pymihg.f \
pythia/pymihk.f \
pythia/pymirm.f \
pythia/pymrun.f \
pythia/pymsin.f \
pythia/pymult.f \
pythia/pyname.f \
pythia/pynjdc.f \
pythia/pynmes.f \
pythia/pynull.f \
pythia/pyofsh.f \
pythia/pyonia.f \
pythia/pyonof.f \
pythia/pyoper.f \
pythia/pyp.f \
pythia/pypcmp.f \
pythia/pypdel.f \
pythia/pypdfl.f \
pythia/pypdfu.f \
pythia/pypdga.f \
pythia/pypdpi.f \
pythia/pypdpo.f \
pythia/pypdpr.f \
pythia/pypile.f \
pythia/pyplot.f \
pythia/pypole.f \
pythia/pyprep.f \
pythia/pyptdi.f \
pythia/pyptfs.f \
pythia/pyptis.f \
pythia/pyptmi.f \
pythia/pyqqbh.f \
pythia/pyr.f \
pythia/pyradk.f \
pythia/pyrand.f \
pythia/pyreco.f \
pythia/pyremn.f \
pythia/pyresd.f \
pythia/pyrget.f \
pythia/pyrghm.f \
pythia/pyrnm3.f \
pythia/pyrnmq.f \
pythia/pyrobo.f \
pythia/pyrset.f \
pythia/pyrvch.f \
pythia/pyrvg1.f \
pythia/pyrvg2.f \
pythia/pyrvg3.f \
pythia/pyrvg4.f \
pythia/pyrvgl.f \
pythia/pyrvgw.f \
pythia/pyrvi1.f \
pythia/pyrvi2.f \
pythia/pyrvi3.f \
pythia/pyrvne.f \
pythia/pyrvr.f \
pythia/pyrvs.f \
pythia/pyrvsb.f \
pythia/pyrvsf.f \
pythia/pysave.f \
pythia/pyscat.f \
pythia/pysfdc.f \
pythia/pysgex.f \
pythia/pysghf.f \
pythia/pysghg.f \
pythia/pysgqc.f \
pythia/pysgsu.f \
pythia/pysgtc.f \
pythia/pysgwz.f \
pythia/pyshow.f \
pythia/pysigh.f \
pythia/pysimp.f \
pythia/pyslha.f \
pythia/pyspen.f \
pythia/pysphe.f \
pythia/pyspli.f \
pythia/pysspa.f \
pythia/pystat.f \
pythia/pystbh.f \
pythia/pystop.f \
pythia/pystrf.f \
pythia/pysubh.f \
pythia/pysugi.f \
pythia/pytabu.f \
pythia/pytaud.f \
pythia/pytbbc.f \
pythia/pytbbn.f \
pythia/pytbdy.f \
pythia/pytbhb.f \
pythia/pytbhg.f \
pythia/pytbhq.f \
pythia/pytbhs.f \
pythia/pytecm.f \
pythia/pytest.f \
pythia/pythag.f \
pythia/pythrg.f \
pythia/pythru.f \
pythia/pytime.f \
pythia/pytune.f \
pythia/pyupda.f \
pythia/pyupev.f \
pythia/pyupin.f \
pythia/pyupre.f \
pythia/pyveto.f \
pythia/pywaux.f \
pythia/pywidt.f \
pythia/pywidx.f \
pythia/pyx2xg.f \
pythia/pyx2xh.f \
pythia/pyx3jt.f \
pythia/pyx4jt.f \
pythia/pyxdif.f \
pythia/pyxjet.f \
pythia/pyxkfl.f \
pythia/pyxtee.f \
pythia/pyxtot.f \
pythia/pyxxga.f \
pythia/pyxxz6.f \
pythia/pyzdis.f \
pythia/ssmssm.f \
pythia/sugra.f \
pythia/upevnt.f \
pythia/upinit.f \
pythia/upveto.f \
pythia/visaje.f

libpythia6eRHIC_la_FFLAGS = -fPIC -I$(top_srcdir)/include

# RADGEN source files (all Fortran 77).
libradgen_la_SOURCES = \
radgen/radgen.f \
radgen/radgen_event.f \
radgen/radgen_init.f

libradgen_la_FFLAGS = -fPIC -I$(top_srcdir)/include

# pythiaeRHIC source files.
# These are a mixture of Fortran 77, Fortran 95 and C++.
pythiaeRHIC_SOURCES = \
erhic/factory.cpp \
erhic/factory.h \
erhic/gmc_random.f \
erhic/pyth_xsec.f \
erhic/pythia6.f95 \
erhic/pythia_commons.cpp \
erhic/pythia_commons.h \
erhic/pythia_erhic.h \
erhic/pythia_main.cpp \
erhic/pythia_radgen_extras.f

# Set CPP flags for pythiaeRHIC.
# C++ headers from this project are in the erhic directory.
# It also needs to know where the LHAPDF and eic-smear headers are.
# ROOT headers were located using root-config (see condfigure.ac).
pythiaeRHIC_CPPFLAGS = \
-I$(top_srcdir)/erhic \
$(EIC_INC) \
$(LHAPDF_INC)

pythiaeRHIC_FFLAGS = -I$(top_srcdir)/include

pythiaeRHIC_FCFLAGS = -I$(top_srcdir)/include

# pythiaeRHIC linker flags.
# Links against libpythia6 and libradgen from this project,
# libeicsmear and a few CERNLIB libraries.
# Linker flags for ROOT were specified in configure.ac
# using root-config.
pythiaeRHIC_LDFLAGS = -L.libs -lpythia6eRHIC -lradgen $(CERNLIB) $(EIC_LIB) $(LHAPDF_LIB)

clean-local:
	rm -f $(BUILT_SOURCES)
