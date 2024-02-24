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

SETDTKNPROC := 1

# Use OpenMP
DTKUSEOPENMP := 0
ifneq ($(SETDTKNPROC),1)
      DTKUSEOPENMP=$(SETDTKNPROC)
      CXXFLAGS   += -DPARALLELISEDTK=$(DTKUSEOPENMP) -fopenmp
      LFLAGS     += -fopenmp
else
      CXXFLAGS   += -DPARALLELISEDTK=0
endif

MYMAKEFLAGS := 

# COMPILERS OPTIONS
ifneq ($(DARWIN),)
  CXX          := g++-12
else
  CXX          := g++
endif

CXXFLAGS       += -std=c++11 -fPIC -Wall -pedantic
OPTIMFLAGS      = -O2 -funroll-loops -falign-loops=8 -ffast-math

ARCHIVE      := ar
ARCHFLAG     := -rc

# Debug build flags
DEBUGVERSION := 0
ifeq ($(DEBUGVERSION),1)
  CXXFLAGS     += -DDEBUG=1 -ggdb -W -Wno-long-long
  MYMAKEFLAGS += DEBUGVERSION=1
else
   ifeq ($(PROFILEVERSION),1)
      CXXFLAGS     += -DDEBUG=1 -g -pg -Wall -pedantic $(OPTIMFLAGS)
      MYMAKEFLAGS += PROFILEVERSION=1
   else
      CXXFLAGS     += -DDEBUG=0 -Wall -pedantic $(OPTIMFLAGS)
   endif
endif

INCDEFS        :=  -include ../common/globaldefs.h -include localdefs.h

# FILES
CPPOBJS=$(patsubst %.cpp,%.o,$(wildcard *.cpp))
CCOBJS=$(patsubst %.cc,%.o,$(wildcard *.cc))
OBJS= $(CPPOBJS) $(CCOBJS)
HEADERS=$(wildcard *.h)
SOURCES=$(wildcard *.cpp *.cc)
LIBOBJS=$(shell echo $(OBJS))
COMMONSTATICLIBNAME=$(shell cd ../..; pwd | sed -e 's;\(.*\)/\(.*\);\2;')
COMMONHEADERS=$(shell ls ../common/*.h)
COMMONSOURCES=$(shell ls ../common/*.cpp)
STATICLIB=$(TOP)/../common/lib$(COMMONSTATICLIBNAME).a
STATICOBJS=$(shell for i in $$(ls ../common/*cpp);do echo $${i%cpp}o;done)
TESTEXECS=$(patsubst %.cc,%.x,$(wildcard *.cc))

# -L: FOLDER LIBRARY
	LFLAGS+=#-L$(TOP)/../common/
	LFLAGS+=#-L.

# -l: LIBRARY
	lLIBS+=-l$(COMMONSTATICLIBNAME)

# -I: INCLUDES
	IFLAGS+=-I. $(INCDEFS) -I$(TOP)/../common/

TARGET=target

$(TARGET): $(OBJS) $(TESTEXECS) $(STATICLIB)

$(CPPOBJS): %.o: %.cpp %.h
	@echo "\033[32m   Compiling $@\\033[m"
	$(CXX) $(CXXFLAGS) -c $< -o $@  $(IFLAGS) 

$(CCOBJS): %.o: %.cc $(COMMONHEADERS) $(HEADERS)
	@echo "\033[32m   Compiling $@\\033[m"
	$(CXX) $(CXXFLAGS) -c $<  -o $@  $(IFLAGS) 

$(TESTEXECS): %.x: %.o $(STATICLIB) $(CCOBJS) $(CPPOBJS) $(SOURCES) $(HEADERS) $(COMMONHEADERS) $(COMMONSOURCES)
	@echo "\033[32m   Linking $@\\033[m"
	$(CXX) $(LFLAGS) $(OPTIMFLAGS) $< $(CPPOBJS) $(STATICLIB) -o $@

$(STATICLIB): $(COMMONHEADERS) $(COMMONSOURCES)
	@echo "\033[32mBuilding lib$(COMMONSTATICLIBNAME).a\\033[m"
	@cd $(TOP)/../common/;make $(MYMAKEFLAGS)

.PHONY: clean
clean:
	$(info CLEANING ALL)
	@$(RM) -f $(TARGET) $(OBJS) $(TESTEXECS) 2>/dev/null || true

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

fullclean: clean
	@cd ../common/;make clean

list:
	@echo $(OBJS)
	@echo $(HEADERS)

listobjs:
	@echo $(LIBOBJS)

liststatic:
	@echo $(STATICOBJS)

