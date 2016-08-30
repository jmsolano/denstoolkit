#
#
#Default path for installation of DensToolKit:
#
INSTALL_PATH = /usr/local/bin
#
#
#Change the following line for Windows OS, choose the
#appropriate name of povray
#USERPOVCMD = pvengine.exe
#USERPOVCMD = pvengine64.exe
#
USERPOVCMD = povray

#
#Set the number of processors that will be used by DTK.
#
SETDTKNPROC=1

#
#*************************************************************
# Please, do not modify below this line
#*************************************************************
#

LCL_SRCDIR   = $(shell pwd)
TOP          = $(shell cd ..; pwd)
LCL_OBJDIR   = $(LCL_SRCDIR)/objs
LCL_SHLDIR   = $(LCL_SRCDIR)/scripts
LCL_DOCDIR   = $(TOP)/doc
LCL_BINDIR   = $(TOP)/bin
LCL_WFDIR    = $(TOP)/wavefiles
LCL_OUTDIR   = $(TOP)/outputs
LCL_DVDCDIR  = $(LCL_SRCDIR)/devdoc

DTKCOMMONDIR    = $(LCL_SRCDIR)/common
DTKPOINTDIR     = $(LCL_SRCDIR)/dtkpoint
DTKLINEDIR      = $(LCL_SRCDIR)/dtkline
DTKPLANEDIR     = $(LCL_SRCDIR)/dtkplane
DTKCUBEDIR      = $(LCL_SRCDIR)/dtkcube
DTKFINDCPDIR    = $(LCL_SRCDIR)/dtkfindcp
DTKMOMDDIR      = $(LCL_SRCDIR)/dtkmomd
DTKQDMOLDIR     = $(LCL_SRCDIR)/dtkqdmol
DTKDEMAT1DIR    = $(LCL_SRCDIR)/dtkdemat1
DTKBPDENSDIR    = $(LCL_SRCDIR)/dtkbpdens
DTKCHARCICPDIR  = $(LCL_SRCDIR)/dtkcharcicp

DTKALLDIRS      = $(DTKCOMMONDIR) $(DTKPOINTDIR) $(DTKLINEDIR) $(DTKPLANEDIR) \
                  $(DTKCUBEDIR) $(DTKFINDCPDIR) $(DTKMOMDDIR) $(DTKQDMOLDIR) \
                  $(DTKDEMAT1DIR) $(DTKBPDENSDIR) $(DTKCHARCICPDIR)

DTKCOMMONFILES    = $(DTKCOMMONDIR)/bondnetwork.cpp $(DTKCOMMONDIR)/bondnetwork.h \
                  $(DTKCOMMONDIR)/critptnetwork.cpp $(DTKCOMMONDIR)/critptnetwork.h \
                  $(DTKCOMMONDIR)/eig2-4.cpp $(DTKCOMMONDIR)/eig2-4.h \
                  $(DTKCOMMONDIR)/fldtypesdef.h \
                  $(DTKCOMMONDIR)/iofuncts-wfn.cpp $(DTKCOMMONDIR)/iofuncts-wfn.h \
                  $(DTKCOMMONDIR)/iofuncts-wfx.cpp $(DTKCOMMONDIR)/iofuncts-wfx.h \
                  $(DTKCOMMONDIR)/solcubetools.cpp $(DTKCOMMONDIR)/solcubetools.h \
                  $(DTKCOMMONDIR)/solfileutils.cpp $(DTKCOMMONDIR)/solfileutils.h \
                  $(DTKCOMMONDIR)/solmath.cpp $(DTKCOMMONDIR)/solmath.h \
                  $(DTKCOMMONDIR)/solmemhand.cpp $(DTKCOMMONDIR)/solmemhand.h \
                  $(DTKCOMMONDIR)/solpovtools.cpp $(DTKCOMMONDIR)/solpovtools.h \
                  $(DTKCOMMONDIR)/solscrutils.cpp $(DTKCOMMONDIR)/solscrutils.h \
                  $(DTKCOMMONDIR)/solstringtools.cpp  $(DTKCOMMONDIR)/solstringtools.h \
                  $(DTKCOMMONDIR)/gausswavefunction.cpp $(DTKCOMMONDIR)/gausswavefunction.h \
                  $(DTKCOMMONDIR)/wfgrid1d.cpp $(DTKCOMMONDIR)/wfgrid1d.h \
                  $(DTKCOMMONDIR)/wfgrid2d.cpp $(DTKCOMMONDIR)/wfgrid2d.h \
                  $(DTKCOMMONDIR)/wfgrid3d.cpp $(DTKCOMMONDIR)/wfgrid3d.h \
                  $(DTKCOMMONDIR)/iofuncts-cpx.cpp $(DTKCOMMONDIR)/iofuncts-cpx.h
DTKPOINTFILES     = $(DTKPOINTDIR)/dtkpoint.cpp $(DTKPOINTDIR)/optflags.h \
                  $(DTKPOINTDIR)/optflags.cpp $(DTKPOINTDIR)/crtflnms.h \
                  $(DTKPOINTDIR)/crtflnms.cpp $(DTKPOINTDIR)/soldefines.h
DTKLINEFILES      = $(DTKLINEDIR)/dtkline.cpp $(DTKLINEDIR)/optflags.h \
                  $(DTKLINEDIR)/optflags.cpp $(DTKLINEDIR)/crtflnms.h \
                  $(DTKLINEDIR)/crtflnms.cpp $(DTKLINEDIR)/soldefines.h
DTKPLANEFILES     = $(DTKPLANEDIR)/dtkplane.cpp $(DTKPLANEDIR)/optflags.h \
                  $(DTKPLANEDIR)/optflags.cpp $(DTKPLANEDIR)/crtflnms.h \
                  $(DTKPLANEDIR)/crtflnms.cpp $(DTKPLANEDIR)/soldefines.h
DTKCUBEFILES      = $(DTKCUBEDIR)/dtkcube.cpp $(DTKCUBEDIR)/optflags.h \
                  $(DTKCUBEDIR)/optflags.cpp $(DTKCUBEDIR)/crtflnms.h \
                  $(DTKCUBEDIR)/crtflnms.cpp $(DTKCUBEDIR)/soldefines.h
DTKFINDCPFILES    = $(DTKFINDCPDIR)/dtkfindcp.cpp $(DTKFINDCPDIR)/optflags.h \
                  $(DTKFINDCPDIR)/optflags.cpp $(DTKFINDCPDIR)/crtflnms.h \
                  $(DTKFINDCPDIR)/crtflnms.cpp $(DTKFINDCPDIR)/soldefines.h \
                  $(DTKFINDCPDIR)/custfmtmathfuncts.cpp $(DTKFINDCPDIR)/custfmtmathfuncts.h
DTKMOMDFILES      = $(DTKMOMDDIR)/dtkmomd.cpp $(DTKMOMDDIR)/optflags.h \
                  $(DTKMOMDDIR)/optflags.cpp $(DTKMOMDDIR)/crtflnms.h \
                  $(DTKMOMDDIR)/crtflnms.cpp $(DTKMOMDDIR)/soldefines.h
DTKQDMOLFILES     = $(DTKQDMOLDIR)/dtkqdmol.cpp $(DTKQDMOLDIR)/optflags.h \
                  $(DTKQDMOLDIR)/optflags.cpp $(DTKQDMOLDIR)/crtflnms.h \
                  $(DTKQDMOLDIR)/crtflnms.cpp $(DTKQDMOLDIR)/soldefines.h
DTKDEMAT1FILES    = $(DTKDEMAT1DIR)/dtkdemat1.cpp $(DTKDEMAT1DIR)/optflags.h \
                  $(DTKDEMAT1DIR)/optflags.cpp $(DTKDEMAT1DIR)/crtflnms.h \
                  $(DTKDEMAT1DIR)/crtflnms.cpp $(DTKDEMAT1DIR)/soldefines.h
DTKBPDENSFILES    = $(DTKBPDENSDIR)/dtkbpdens.cpp $(DTKBPDENSDIR)/optflags.h \
                  $(DTKBPDENSDIR)/optflags.cpp $(DTKBPDENSDIR)/crtflnms.h \
                  $(DTKBPDENSDIR)/crtflnms.cpp $(DTKBPDENSDIR)/soldefines.h
DTKCHARCICPFILES  = $(DTKCHARCICPDIR)/dtkcharcicp.cpp $(DTKCHARCICPDIR)/optflags.h \
                  $(DTKCHARCICPDIR)/optflags.cpp $(DTKCHARCICPDIR)/crtflnms.h \
                  $(DTKCHARCICPDIR)/crtflnms.cpp $(DTKCHARCICPDIR)/soldefines.h

DTKALLFILES     = $(DTKCOMMONFILES) $(DTKPOINTFILES) $(DTKLINEFILES) $(DTKPLANEFILES) \
                  $(DTKCUBEFILES) $(DTKFINDCPFILES) $(DTKMOMDFILES) $(DTKQDMOLFILES) \
                  $(DTKDEMAT1FILES) $(DTKBPDENSFILES) $(DTKCHARCICPFILES)

ifeq ($(SETDTKNPROC),1)
      MKFLAGS   = DTKUSEOPENMP=0
else
      MKFLAGS   = DTKUSEOPENMP=$(SETDTKNPROC)
endif

# Target rules
all: $(DTKALLDIRS) $(DTKALLFILES) bindir outdir objdir build

build: $(LCL_BINDIR)/dtkpoint $(LCL_BINDIR)/dtkline $(LCL_BINDIR)/dtkplane \
       $(LCL_BINDIR)/dtkcube $(LCL_BINDIR)/dtkfindcp $(LCL_BINDIR)/dtkmomd \
       $(LCL_BINDIR)/dtkqdmol $(LCL_BINDIR)/dtkdemat1 $(LCL_BINDIR)/dtkbpdens \
		 $(LCL_BINDIR)/dtkcharcicp

$(LCL_BINDIR)/dtkpoint: $(DTKPOINTFILES) $(DTKCOMMONFILES)
	@cd $(DTKPOINTDIR); make $(MKFLAGS); mv dtkpoint $(LCL_BINDIR)/; cd $(LCL_SRCDIR)

$(LCL_BINDIR)/dtkline: $(DTKLINEFILES) $(DTKCOMMONFILES)
	@cd $(DTKLINEDIR); make $(MKFLAGS); mv dtkline $(LCL_BINDIR)/; cd $(LCL_SRCDIR)

$(LCL_BINDIR)/dtkplane: $(DTKPLANE) $(DTKCOMMONFILES)
	@cd $(DTKPLANEDIR); make $(MKFLAGS); mv dtkplane $(LCL_BINDIR)/; cd $(LCL_SRCDIR)

$(LCL_BINDIR)/dtkcube: $(DTKCUBEFILES) $(DTKCOMMONFILES)
	@cd $(DTKCUBEDIR); make $(MKFLAGS); mv dtkcube $(LCL_BINDIR)/; cd $(LCL_SRCDIR)

$(LCL_BINDIR)/dtkfindcp: $(DTKFINDCPFILES) $(DTKCOMMONFILES)
	@cd $(DTKFINDCPDIR); make $(MKFLAGS) SETCMDPOVRAY=\"$(USERPOVCMD)\";\
        mv dtkfindcp $(LCL_BINDIR)/; cd $(LCL_SRCDIR)

$(LCL_BINDIR)/dtkmomd: $(DTKMOMDFILES) $(DTKCOMMONFILES)
	@cd $(DTKMOMDDIR); make $(MKFLAGS); mv dtkmomd $(LCL_BINDIR)/; cd $(LCL_SRCDIR)

$(LCL_BINDIR)/dtkqdmol: $(DTKQDMOLFILES) $(DTKCOMMONFILES)
	@cd $(DTKQDMOLDIR); make $(MKFLAGS); mv dtkqdmol $(LCL_BINDIR)/; cd $(LCL_SRCDIR)

$(LCL_BINDIR)/dtkdemat1: $(DTKDEMAT1FILES) $(DTKCOMMONFILES)
	@cd $(DTKDEMAT1DIR); make $(MKFLAGS); mv dtkdemat1 $(LCL_BINDIR)/; cd $(LCL_SRCDIR)

$(LCL_BINDIR)/dtkbpdens: $(DTKBPDENSFILES) $(DTKCOMMONFILES)
	@cd $(DTKBPDENSDIR); make $(MKFLAGS); mv dtkbpdens $(LCL_BINDIR)/; cd $(LCL_SRCDIR)

$(LCL_BINDIR)/dtkcharcicp: $(DTKCHARCICPFILES) $(DTKCOMMONFILES)
	@cd $(DTKCHARCICPDIR); make $(MKFLAGS); mv dtkcharcicp $(LCL_BINDIR)/; cd $(LCL_SRCDIR)

clean: 
	@if test -d $(DTKPOINTDIR); then (cd $(DTKPOINTDIR); make clean; cd $(LCL_SRCDIR)) fi 
	@if test -d $(DTKLINEDIR); then (cd $(DTKLINEDIR); make clean; cd $(LCL_SRCDIR)) fi 
	@if test -d $(DTKPLANEDIR); then (cd $(DTKPLANEDIR); make clean; cd $(LCL_SRCDIR)) fi
	@if test -d $(DTKCUBEDIR); then (cd $(DTKCUBEDIR); make clean; cd $(LCL_SRCDIR)) fi 
	@if test -d $(DTKFINDCPDIR); then (cd $(DTKFINDCPDIR); make clean; cd $(LCL_SRCDIR)) fi 
	@if test -d $(DTKMOMDDIR); then (cd $(DTKMOMDDIR); make clean; cd $(LCL_SRCDIR)) fi 
	@if test -d $(DTKQDMOLDIR); then (cd $(DTKQDMOLDIR); make clean; cd $(LCL_SRCDIR)) fi 
	@if test -d $(DTKDEMAT1DIR); then (cd $(DTKDEMAT1DIR); make clean; cd $(LCL_SRCDIR)) fi 
	@if test -d $(DTKBPDENSDIR); then (cd $(DTKBPDENSDIR); make clean; cd $(LCL_SRCDIR)) fi 
	@if test -d $(DTKCHARCICPDIR); then (cd $(DTKCHARCICPDIR); make clean; cd $(LCL_SRCDIR)) fi 
	rm -f $(LCL_OBJDIR)/*.o
	rm -f $(LCL_BINDIR)/dtk*

install: build
	cp $(LCL_BINDIR)/dtk* $(INSTALL_PATH)
	cp $(LCL_SHLDIR)/dtk* $(INSTALL_PATH)
	@echo "\033[32m\nInstalled programs:\n"
	@if test -e $(INSTALL_PATH)/dtkpoint; then (dtkpoint -V; dtkline -V; \
		dtkplane -V; dtkcube -V; dtkfindcp -V; dtkmomd -V; dtkdemat1 -V; \
		dtkqdmol -V; dtkbpdens -V; dtkcharcicp -V; \
          cd $(INSTALL_PATH); ls -1 dtk*2*) fi
	@echo "\033[0m"

distclean: clean cleandeveldocs
	rm -f $(LCL_BINDIR)/*
	rm -f $(LCL_OUTDIR)/*
	rm -f $(INSTALL_PATH)/dtkpoint
	rm -f $(INSTALL_PATH)/dtkline
	rm -f $(INSTALL_PATH)/dtkplane
	rm -f $(INSTALL_PATH)/dtkcube
	rm -f $(INSTALL_PATH)/dtkfindcp
	rm -f $(INSTALL_PATH)/dtkmomd
	rm -f $(INSTALL_PATH)/dtkqdmol
	rm -f $(INSTALL_PATH)/dtkdemat1
	rm -f $(INSTALL_PATH)/dtkbpdens
	rm -f $(INSTALL_PATH)/dtkcharcicp
	rm -f $(INSTALL_PATH)/dtkeps2pdf
	rm -f $(INSTALL_PATH)/dtkpov2png

runtest: build
	cd $(LCL_OUTDIR); cp $(LCL_WFDIR)/benzene.wfn .; \
          $(LCL_BINDIR)/dtkpoint benzene.wfn -a 1; rm -f benzene.wfn
	cd $(LCL_OUTDIR); cp $(LCL_WFDIR)/ethanol.wfn .; \
          $(LCL_BINDIR)/dtkline ethanol.wfn -a 2 5 -p K -P -z; rm -f ethanol.wfn
	cd $(LCL_OUTDIR); cp $(LCL_WFDIR)/cubano_sto3g.wfx .; \
          $(LCL_BINDIR)/dtkplane cubano_sto3g.wfx -a 1 2 3 -p g -P -c -z -n 200; \
          rm -f cubano_sto3g.wfx
	cd $(LCL_OUTDIR); cp $(LCL_WFDIR)/ch4.wfx .; \
          $(LCL_BINDIR)/dtkcube ch4.wfx -s 100 -p d -z; rm ch4.wfx
	cd $(LCL_OUTDIR); cp $(LCL_WFDIR)/cyclopropane.wfx .; \
          $(LCL_BINDIR)/dtkfindcp cyclopropane.wfx -P -g -T -k -m -a; \
          rm -f cyclopropane.wfx
	cd $(LCL_OUTDIR); cp $(LCL_WFDIR)/benzene.wfn .; \
          $(LCL_BINDIR)/dtkmomd benzene.wfn -2 xy -P; rm benzene.wfn
	cd $(LCL_OUTDIR); cp $(LCL_WFDIR)/cubano_sto3g.wfx .; \
          $(LCL_BINDIR)/dtkqdmol cubano_sto3g.wfx; \
          rm -f cubano_sto3g.wfx
	@echo "Note: The actual gnp file (dtkqdmol) is $(LCL_OUTDIR)/cubano_sto3g*.gnp"
	cd $(LCL_OUTDIR); cp $(LCL_WFDIR)/cubano_sto3g.wfx .; \
          $(LCL_BINDIR)/dtkdemat1 cubano_sto3g.wfx -P -l -c -a 5 6 -s 0.02e0; \
          rm -f cubano_sto3g.wfx
	cd $(LCL_OUTDIR); cp $(LCL_WFDIR)/cyclopropane.wfx .; \
          $(LCL_BINDIR)/dtkbpdens cyclopropane.wfx -P -p M -n 160 -s 0.02e0 -l; \
          $(LCL_BINDIR)/dtkbpdens cyclopropane.wfx -P -p M -n 145 -l -L; \
          rm -f cyclopropane.wfx
	cd $(LCL_OUTDIR); cp $(LCL_WFDIR)/benzene.wfn.; \
          $(LCL_BINDIR)/dtkcharcicp cubano_sto3g.wfx; \
          rm -f benzene.wfn

develdocs: devdochtmldir
	@cd $(LCL_DVDCDIR); doxygen Doxyfile
	@if test -L $(LCL_DVDCDIR)/DeveloperDocumentation.html; \
        then : ; \
        else (ln -s $(LCL_DVDCDIR)/html/index.html $(LCL_DVDCDIR)/DeveloperDocumentation.html) \
        fi

latexdeveldocs: develdocs
	@cd $(LCL_DVDCDIR)/latex; make; mv $(LCL_DVDCDIR)/latex/refman.pdf $(LCL_DVDCDIR)

cleandeveldocs:
	@rm -rf $(LCL_DVDCDIR)/html/*
	@rm -rf $(LCL_DVDCDIR)/latex/*
	@rm -f $(LCL_DVDCDIR)/DeveloperDocumentation.html

bindir:
	@mkdir -p $(LCL_BINDIR)

outdir:
	@mkdir -p $(LCL_OUTDIR)

objdir:
	@mkdir -p $(LCL_OBJDIR)

devdochtmldir:
	@mkdir -p $(LCL_DVDCDIR)/html


