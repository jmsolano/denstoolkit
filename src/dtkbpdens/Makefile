TOP      = $(shell pwd)

# OS Name (Linux or Darwin)
OSUPPER = $(shell uname -s 2>/dev/null | tr [:lower:] [:upper:])
OSLOWER = $(shell uname -s 2>/dev/null | tr [:upper:] [:lower:])

# Flags to detect 32-bit or 64-bit OS platform
OS_SIZE = $(shell uname -m | sed -e "s/i.86/32/" -e "s/x86_64/64/")
OS_ARCH = $(shell uname -m | sed -e "s/i386/i686/")

# Flags to detect either a Linux system (linux) or Mac OSX (darwin)
DARWIN = $(strip $(findstring DARWIN, $(OSUPPER)))

GCC             ?= g++

# Extra user flags
OPTIM_CCFLAGS   ?= -O2 -funroll-loops -falign-loops=16
WARNING_CCFLAGS ?=  -W -Wall -pedantic -Wno-long-long

# OS-specific build flags
ifneq ($(DARWIN),) 
      CCFLAGS   := -arch $(OS_ARCH) 
else
  ifeq ($(OS_SIZE),32)
      CCFLAGS   := -m32
  else
      CCFLAGS   := -m64
  endif
endif

# Debug build flags
ifeq ($(DEBUGVERSION),1)
      CCFLAGS   += -DDEBUG=1
else
      CCFLAGS   += -DDEBUG=0
endif

# Use OpenMP
ifneq ($(DTKUSEOPENMP),0)
      CCFLAGS   += -DPARALLELISEDTK=$(DTKUSEOPENMP) -fopenmp
else
      CCFLAGS   += -DPARALLELISEDTK=0
endif

OBJDIR        := ../objs
INCSOLDEFS    := -include soldefines.h
PARTDEPS      := crtflnms.h crtflnms.cpp optflags.h optflags.cpp ../common/fldtypesdef.h

# Target rules
all: build

build: dtkbpdens

$(OBJDIR)/solmemhand.o: ../common/solmemhand.cpp soldefines.h
	$(GCC) $(CCFLAGS) $(WARNING_CCFLAGS) $(INCSOLDEFS) -o $@ -c $<

$(OBJDIR)/solscrutils.o: ../common/solscrutils.cpp soldefines.h
	$(GCC) $(CCFLAGS) $(WARNING_CCFLAGS) $(INCSOLDEFS) -o $@ -c $<

$(OBJDIR)/solstringtools.o: ../common/solstringtools.cpp soldefines.h
	$(GCC) $(CCFLAGS) $(WARNING_CCFLAGS) $(INCSOLDEFS) -o $@ -c $<

$(OBJDIR)/solfileutils.o: ../common/solfileutils.cpp soldefines.h
	$(GCC) $(CCFLAGS) $(WARNING_CCFLAGS) $(INCSOLDEFS) -o $@ -c $<

$(OBJDIR)/iofuncts-wfx.o: ../common/iofuncts-wfx.cpp ../common/solmemhand.cpp soldefines.h
	$(GCC) $(CCFLAGS) $(WARNING_CCFLAGS) $(INCSOLDEFS) -o $@ -c $<

$(OBJDIR)/iofuncts-wfn.o: ../common/iofuncts-wfn.cpp ../common/solmemhand.cpp soldefines.h
	$(GCC) $(CCFLAGS) $(WARNING_CCFLAGS) $(INCSOLDEFS) -o $@ -c $<

$(OBJDIR)/iofuncts-cpx.o: ../common/iofuncts-cpx.cpp ../common/solmemhand.cpp \
                          soldefines.h ../common/fldtypesdef.h ../common/solscrutils.h \
                          ../common/solstringtools.h
	$(GCC) $(CCFLAGS) $(WARNING_CCFLAGS) $(INCSOLDEFS) -o $@ -c $<

$(OBJDIR)/solmath.o: ../common/solmath.cpp soldefines.h
	$(GCC) $(CCFLAGS) $(WARNING_CCFLAGS) $(INCSOLDEFS) $(OPTIM_CCFLAGS) -o $@ -c $<

$(OBJDIR)/atomradiicust.o: ../common/atomradiicust.cpp ../common/atomradiicust.h soldefines.h
	$(GCC) $(CCFLAGS) $(WARNING_CCFLAGS) $(INCSOLDEFS) $(OPTIM_CCFLAGS) -o $@ -c $<

$(OBJDIR)/eig2-4.o: ../common/eig2-4.cpp soldefines.h
	$(GCC) $(CCFLAGS) $(WARNING_CCFLAGS) $(INCSOLDEFS) $(OPTIM_CCFLAGS) -o $@ -c $<

$(OBJDIR)/wavefunctionclass.o: ../common/wavefunctionclass.cpp soldefines.h
	$(GCC) $(CCFLAGS) $(WARNING_CCFLAGS) $(INCSOLDEFS) $(OPTIM_CCFLAGS) -o $@ -c $<

$(OBJDIR)/solpovtools.o: ../common/solpovtools.cpp soldefines.h
	$(GCC) $(CCFLAGS) $(WARNING_CCFLAGS) $(INCSOLDEFS) -o $@ -c $<

$(OBJDIR)/bondnetwork.o: ../common/bondnetwork.cpp soldefines.h #$(OBJDIR)/solpovtools.o
	$(GCC) $(CCFLAGS) $(WARNING_CCFLAGS) $(INCSOLDEFS) $(OPTIM_CCFLAGS) -o $@ -c $<

$(OBJDIR)/critptnetwork.o: ../common/critptnetwork.cpp soldefines.h 
	$(GCC) $(CCFLAGS) $(WARNING_CCFLAGS) $(INCSOLDEFS) $(OPTIM_CCFLAGS) -o $@ -c $<

$(OBJDIR)/wfgrid2d.o: ../common/wfgrid2d.cpp soldefines.h
	$(GCC) $(CCFLAGS) $(WARNING_CCFLAGS) $(INCSOLDEFS) $(OPTIM_CCFLAGS) -o $@ -c $<

$(OBJDIR)/atomcolschjmol.o: ../common/atomcolschjmol.cpp ../common/atomcolschjmol.h soldefines.h
	$(GCC) $(CCFLAGS) $(WARNING_CCFLAGS) $(INCSOLDEFS) $(OPTIM_CCFLAGS) -o $@ -c $<

optflags.o: optflags.cpp soldefines.h
	$(GCC) $(CCFLAGS) $(WARNING_CCFLAGS) $(INCSOLDEFS) -o $@ -c $<

crtflnms.o: crtflnms.cpp soldefines.h ../common/fldtypesdef.h
	$(GCC) $(CCFLAGS) $(WARNING_CCFLAGS) $(INCSOLDEFS) -o $@ -c $<

dtkbpdens.o: dtkbpdens.cpp soldefines.h $(PARTDEPS)
	$(GCC) $(CCFLAGS) $(WARNING_CCFLAGS) $(INCSOLDEFS) -o $@ -c $<

dtkbpdens: dtkbpdens.o optflags.o crtflnms.o $(OBJDIR)/solfileutils.o \
  $(OBJDIR)/iofuncts-wfx.o $(OBJDIR)/iofuncts-wfn.o $(OBJDIR)/iofuncts-cpx.o \
  $(OBJDIR)/wavefunctionclass.o $(OBJDIR)/solscrutils.o $(OBJDIR)/solmemhand.o \
  $(OBJDIR)/bondnetwork.o $(OBJDIR)/wfgrid2d.o $(OBJDIR)/solpovtools.o \
  $(OBJDIR)/solmath.o $(OBJDIR)/eig2-4.o $(OBJDIR)/critptnetwork.o \
  $(OBJDIR)/solstringtools.o $(OBJDIR)/atomradiicust.o $(OBJDIR)/atomcolschjmol.o
	$(GCC) $(CCFLAGS) $(OPTIM_CCFLAGS) $(INCSOLDEFS) -o $@ $+
clean:
	rm -f *.o dtkbpdens

distclean: clean
	rm -f $(OBJDIR)/*.o

runben: build
	./dtkbpdens  benzene.wfn -P -l -p M

rtestost: build
	./dtkbpdens testosterona_sto3g.wfx -a 33 45 -P -l


