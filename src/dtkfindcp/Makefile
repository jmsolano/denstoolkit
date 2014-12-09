TOP      = $(shell pwd)

# OS Name (Linux or Darwin)
OSUPPER = $(shell uname -s 2>/dev/null | tr [:lower:] [:upper:])
OSLOWER = $(shell uname -s 2>/dev/null | tr [:upper:] [:lower:])

# Flags to detect 32-bit or 64-bit OS platform
OS_SIZE = $(shell uname -m | sed -e "s/i.86/32/" -e "s/x86_64/64/")
OS_ARCH = $(shell uname -m | sed -e "s/i386/i686/")

# Flags to detect either a Linux system (linux) or Mac OSX (darwin)
DARWIN = $(strip $(findstring DARWIN, $(OSUPPER)))

#GCC             ?= g++-mp-4.8 -std=c++11
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

ifeq ($(SETCMDPOVRAY),)
      ADDPOVCMD = ""
else
      ADDPOVCMD = -D"CMD_POVRAY=\"${SETCMDPOVRAY}\""
endif

OBJDIR        := ../objs
INCSOLDEFS    := -include soldefines.h
PARTDEPS      := crtflnms.h crtflnms.cpp optflags.h optflags.cpp

# Target rules
all: build

build: dtkfindcp

$(OBJDIR)/solmemhand.o: ../common/solmemhand.cpp soldefines.h
	$(GCC) $(CCFLAGS) $(WARNING_CCFLAGS) $(INCSOLDEFS) -o $@ -c $<

$(OBJDIR)/solscrutils.o: ../common/solscrutils.cpp soldefines.h
	$(GCC) $(CCFLAGS) $(WARNING_CCFLAGS) $(INCSOLDEFS) -o $@ -c $<

$(OBJDIR)/solfileutils.o: ../common/solfileutils.cpp soldefines.h
	$(GCC) $(CCFLAGS) $(WARNING_CCFLAGS) $(INCSOLDEFS) -o $@ -c $<

$(OBJDIR)/iofuncts-wfx.o: ../common/iofuncts-wfx.cpp ../common/solmemhand.cpp soldefines.h
	$(GCC) $(CCFLAGS) $(WARNING_CCFLAGS) $(INCSOLDEFS) -o $@ -c $<

$(OBJDIR)/iofuncts-wfn.o: ../common/iofuncts-wfn.cpp ../common/solmemhand.cpp soldefines.h
	$(GCC) $(CCFLAGS) $(WARNING_CCFLAGS) $(INCSOLDEFS) -o $@ -c $<

$(OBJDIR)/solmath.o: ../common/solmath.cpp soldefines.h
	$(GCC) $(CCFLAGS) $(WARNING_CCFLAGS) $(INCSOLDEFS) $(OPTIM_CCFLAGS) -o $@ -c $<

$(OBJDIR)/eig2-4.o: ../common/eig2-4.cpp soldefines.h
	$(GCC) $(CCFLAGS) $(WARNING_CCFLAGS) $(INCSOLDEFS) $(OPTIM_CCFLAGS) -o $@ -c $<

$(OBJDIR)/solstringtools.o: ../common/solstringtools.cpp soldefines.h
	$(GCC) $(CCFLAGS) $(WARNING_CCFLAGS) $(INCSOLDEFS) -o $@ -c $<

$(OBJDIR)/solpovtools.o: ../common/solpovtools.cpp soldefines.h
	$(GCC) $(CCFLAGS) $(WARNING_CCFLAGS) $(INCSOLDEFS) -o $@ -c $<

$(OBJDIR)/wavefunctionclass.o: ../common/wavefunctionclass.cpp soldefines.h
	$(GCC) $(CCFLAGS) $(WARNING_CCFLAGS) $(INCSOLDEFS) $(OPTIM_CCFLAGS) -o $@ -c $<

$(OBJDIR)/bondnetwork.o: ../common/bondnetwork.cpp soldefines.h #$(OBJDIR)/solpovtools.o
	$(GCC) $(CCFLAGS) $(WARNING_CCFLAGS) $(INCSOLDEFS) $(OPTIM_CCFLAGS) -o $@ -c $<

$(OBJDIR)/solcubetools.o: ../common/solcubetools.cpp soldefines.h
	$(GCC) $(CCFLAGS) $(WARNING_CCFLAGS) $(INCSOLDEFS) -o $@ -c $<

$(OBJDIR)/critptnetwork.o: ../common/critptnetwork.cpp soldefines.h
	$(GCC) $(CCFLAGS) $(WARNING_CCFLAGS) $(INCSOLDEFS) $(OPTIM_CCFLAGS) -o $@ -c $<

$(OBJDIR)/iofuncts-cpx.o: ../common/iofuncts-cpx.cpp soldefines.h
	$(GCC) $(CCFLAGS) $(WARNING_CCFLAGS) $(INCSOLDEFS) $(OPTIM_CCFLAGS) -o $@ -c $<

$(OBJDIR)/atomradiicust.o: ../common/atomradiicust.cpp ../common/atomradiicust.h soldefines.h
	$(GCC) $(CCFLAGS) $(WARNING_CCFLAGS) $(INCSOLDEFS) $(OPTIM_CCFLAGS) -o $@ -c $<

$(OBJDIR)/atomcolschjmol.o: ../common/atomcolschjmol.cpp ../common/atomcolschjmol.h soldefines.h
	$(GCC) $(CCFLAGS) $(WARNING_CCFLAGS) $(INCSOLDEFS) $(OPTIM_CCFLAGS) -o $@ -c $<

optflags.o: optflags.cpp soldefines.h
	$(GCC) $(CCFLAGS) $(WARNING_CCFLAGS) $(INCSOLDEFS) -o $@ -c $<

crtflnms.o: crtflnms.cpp soldefines.h
	$(GCC) $(CCFLAGS) $(WARNING_CCFLAGS) $(INCSOLDEFS) -o $@ -c $<

custfmtmathfuncts.o: custfmtmathfuncts.cpp soldefines.h
	$(GCC) $(CCFLAGS) $(WARNING_CCFLAGS) $(INCSOLDEFS) -o $@ -c $<

dtkfindcp.o: dtkfindcp.cpp soldefines.h $(PARTDEPS)
	$(GCC) $(CCFLAGS) $(WARNING_CCFLAGS) $(ADDPOVCMD) $(INCSOLDEFS) -o $@ -c $<

dtkfindcp: dtkfindcp.o optflags.o crtflnms.o custfmtmathfuncts.o $(OBJDIR)/solfileutils.o \
  $(OBJDIR)/iofuncts-wfx.o $(OBJDIR)/iofuncts-wfn.o $(OBJDIR)/eig2-4.o \
  $(OBJDIR)/wavefunctionclass.o $(OBJDIR)/solscrutils.o $(OBJDIR)/solmemhand.o \
  $(OBJDIR)/bondnetwork.o $(OBJDIR)/solpovtools.o $(OBJDIR)/critptnetwork.o \
  $(OBJDIR)/iofuncts-cpx.o $(OBJDIR)/solstringtools.o $(OBJDIR)/atomradiicust.o \
  $(OBJDIR)/atomcolschjmol.o $(OBJDIR)/solmath.o
	$(GCC) $(CCFLAGS) $(OPTIM_CCFLAGS) $(INCSOLDEFS) -o $@ $+
clean:
	rm -f *.o dtkfindcp

distclean: clean
	rm -f $(OBJDIR)/*.o

runben: build
	./dtkfindcp benzene.wfn -P -g -T -k -m

runete: build
	./dtkfindcp ethane.wfn -P -g -T -k

runc60: build
	./dtkfindcp c60-321g.wfx -P -g -T -k

runh2ocomp: build
	./dtkfindcp 4h2o-anion.wfn -P -g -k -a

runli2: build
	./dtkfindcp li2.ccsd.augccpvtz.wfn -P -g -k -a -c 100

