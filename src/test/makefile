
SETDTKNPROC := 1

TOP      = $(shell pwd)

# OS Name (Linux or Darwin)
OSUPPER = $(shell uname -s 2>/dev/null | tr [:lower:] [:upper:])
OSLOWER = $(shell uname -s 2>/dev/null | tr [:upper:] [:lower:])

# Flags to detect 32-bit or 64-bit OS platform
OS_SIZE = $(shell uname -m | sed -e "s/i.86/32/" -e "s/x86_64/64/")
OS_ARCH = $(shell uname -m | sed -e "s/i386/i686/")

# Flags to detect either a Linux system (linux) or Mac OSX (darwin)
DARWIN = $(strip $(findstring DARWIN, $(OSUPPER)))

# OS-specific build flags
ifneq ($(DARWIN),)
    CXXFLAGS   := -arch $(OS_ARCH)
    DARWINVERSION=$(shell sw_vers -productVersion | tr '.' ' ' | awk '{print $$2}')
	 GTTN_14 = $(shell [ $(DARWINVERSION) -gt 14 ] && echo true)
    ifeq ($(GTTN_14),true)
	     CXXFLAGS += -isysroot $(shell xcrun --show-sdk-path)
    endif
else
  ifeq ($(OS_SIZE),32)
    CXXFLAGS   := -m32
  else
    ifeq ($(OS_ARCH),aarch64)
      CXXFLAGS   := -march=armv8-a
    else
      CXXFLAGS   := -m64
    endif
  endif
endif

# Use OpenMP
DTKUSEOPENMP := 0
ifneq ($(SETDTKNPROC),1)
      DTKUSEOPENMP=$(SETDTKNPROC)
      CXXFLAGS   += -DPARALLELISEDTK=$(DTKUSEOPENMP) -fopenmp
else
      CXXFLAGS   += -DPARALLELISEDTK=0
endif

MYMAKEFLAGS = SETDTKNPROC=$(SETDTKNPROC)

# COMPILERS OPTIONS
ifneq ($(DARWIN),)
  CXX          := g++-12
else
  CXX          := g++
endif

CXXFLAGS       += -std=c++11 -fPIC -Wall -pedantic
OPTIMFLAGS      = -O2 -funroll-loops -falign-loops=8

ARCHIVE      := ar
ARCHFLAG     := -rc

# Debug build flags
ifeq ($(DEBUGVERSION),1)
  CXXFLAGS     += -DDEBUG=1 -g -W -Wno-long-long
else
  CXXFLAGS     += -DDEBUG=0 -Wall $(OPTIMFLAGS)
endif

INCDEFS        :=  -include ../common/globaldefs.h -include localdefs.h

# FILES
CPPOBJS=$(patsubst %.cpp,%.o,$(wildcard *.cpp))
CCOBJS=$(patsubst %.cc,%.o,$(wildcard *.cc))
OBJS= $(CPPOBJS) $(CCOBJS)
HEADERS=$(wildcard *.h)
LIBOBJS=$(shell echo $(OBJS))
COMMONSTATICLIBNAME=$(shell cd ../..; pwd | sed -e 's;\(.*\)/\(.*\);\2;')
COMMONHEADERS=$(shell ls ../common/*.h)
COMMONSOURCES=$(shell ls ../common/*.cpp)
STATICLIB=$(TOP)/../common/lib$(COMMONSTATICLIBNAME).a
STATICOBJS=$(shell for i in $$(ls ../common/*cpp);do echo $${i%cpp}o;done)
TESTEXECS=$(patsubst %.cc,%.x,$(wildcard *.cc))

# -L: FOLDER LIBRARY
	LFLAGS+=-L$(TOP)/../common/
	LFLAGS+=#-L.

# -l: LIBRARY
	lLIBS+=-l$(COMMONSTATICLIBNAME)
ifneq ($(SETDTKNPROC),1)
      lLIBS   += -fopenmp
endif

# -I: INCLUDES
	IFLAGS+=-I. $(INCDEFS) -I$(TOP)/../common/

TARGET=target

$(TARGET): localdefs.h $(OBJS) $(TESTEXECS) $(STATICLIB)
	$(info LINKING: $@)

$(CPPOBJS): %.o: %.cpp %.h
	@echo "\033[32m   Compiling $@\\033[m"
	$(CXX) $(CXXFLAGS) -c $< -o $@  $(IFLAGS) 

$(CCOBJS): %.o: %.cc
	@echo "\033[32m   Compiling $@\\033[m"
	$(CXX) $(CXXFLAGS) -c $<  -o $@  $(IFLAGS) 

$(TESTEXECS): %.x: %.o $(STATICLIB)
	$(info Linking: $@)
	$(CXX) $(CXXFLAGS) $(LFLAGS) $< $(lLIBS) -o $@

$(STATICLIB): $(COMMONHEADERS) $(COMMONSOURCES)
	@cd $(TOP)/../common/;$(MAKE)

.PHONY: clean
clean:
	$(info CLEANING ALL)
	@$(RM) -f $(TARGET) $(OBJS) $(STATICLIB) $(TESTEXECS) 2>/dev/null || true

.PHONY: help
help: ## Shows this help
	@echo '=================================================================================='
	@echo 'usage: make [target] ...'
	@echo '=================================================================================='
	@echo ''
	@echo 'targets:'
	@egrep -e '^(.+)\:\ .*##\ (.+)' ${MAKEFILE_LIST} | sed 's/:.*##/#/' | column -t -c 2 -s '#'
	@echo ''
	@echo '=================================================================================='

localdefs.h: 
	@echo "//Here you can define the macros that are used WITHIN"  > localdefs.h
	@echo "//this directory. Notice that the macros defined here" >> localdefs.h
	@echo "//affect ALL and every single *cc and *cpp source." >> localdefs.h
	@echo "#ifndef _LOCALDEFS_H_" >> localdefs.h
	@echo "#define _LOCALDEFS_H_" >> localdefs.h
	@echo "#include \"globaldefs.h\"" >> localdefs.h
	@echo "" >> localdefs.h
	@echo "#endif /* _LOCALDEFS_H_ */" >> localdefs.h

fullclean: clean
	$(info TOTAL CLEANING)
	@cd ../common/;$(MAKE) clean

list:
	@echo $(OBJS)
	@echo $(HEADERS)

listobjs:
	@echo $(LIBOBJS)

