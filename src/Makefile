
SETDTKNPROC=4

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
else
  ifeq ($(OS_SIZE),32)
    CXXFLAGS   := -m32
  else
    CXXFLAGS   := -m64
  endif
endif

INSTALL_PATH := /usr/local/bin

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
  CXX          := g++-mp-4.9
else
  CXX          := g++
endif

# Debug build flags
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

BINDIR = $(TOP)/bin
#EXECUTABLE_DIRECTORIES
DIRSWITHEXES=$(shell ls */*.cc | grep -v test | sed -e 's;\(.*\)/\(.*\).cc;\1;g' | uniq)



# FILES
CCSOURCES=$(shell ls */*.cc | grep -v test)
HEADERS=$(shell ls */*.h | grep -v test | grep -v common)
CPPSOURCES=$(shell ls */*.cpp | grep -v test | grep -v common)
#CCSOURCES=$(shell MyList="$(DIRSWITHEXES)"; for i in $$MyList;do ls $$i/*.cc ; done)
EXECS=$(shell MyList="$(CCSOURCES)"; for i in $$MyList;do echo "$$i" | sed -e 's;\(.*\)/\(.*\).cc;bin/\2;g';done)
COMMONSTATICLIBNAME=$(shell cd ..; pwd | sed -e 's;\(.*\)/\(.*\);\2.a;')
LIBHEADERS=$(shell ls common/*.h)
LIBSOURCES=$(shell ls common/*.cpp)
ALLSOURCES= $(CCSOURCES) $(CPPSOURCES)
ALLHEADERS= $(HEADERS) $(LIBHEADERS)

.PHONY: all
all: $(BINDIR) compileexecs ## Compile every executable '*cc'

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

.PHONY: info
info: ## Shows information about headers, sources, etc.
	@echo '=================================================================================='
	@echo 'Info: '
	@echo '=================================================================================='
	@echo 'Current headers: $(HEADERS)'
	@echo 'Current sources (cc): $(CCSOURCES)'
	@echo 'Current sources (cpp): $(CPPSOURCES)'
	@echo 'Current executables: $(EXECS)'
	@echo ''
	@echo 'Current headers (all): $(ALLHEADERS)'
	@echo 'Current sources (all): $(ALLSOURCES)'
	@echo ''
	@echo 'Current libname: $(COMMONSTATICLIBNAME)'
	@echo 'Current libh: $(LIBHEADERS)'
	@echo 'Current libc: $(LIBSOURCES)'
	@echo '=================================================================================='

compileexecs: $(ALLHEADERS) $(ALLSOURCES)
	@echo "\033[32mCompiling all execs...\\033[m"
	@MyList="$(DIRSWITHEXES)";for i in $$MyList;do \
		echo "\033[32mCompiling in dir $$i\033[m"; cd $$i; make $(MYMAKEFLAGS); for j in $$(ls *.x);do cp $$j ../bin/$${j%.x};done; cd $(TOP) ;done

$(BINDIR):
	@mkdir -p bin

.PHONY: clean
clean: ## Cleans all binaries. Does not remove the installed executables, nor the static library (common/*)
	@echo "\033[33mCleaning bin dir...\\033[m"
	@$(RM) -r bin 2>/dev/null || true

.PHONY: fullclean
fullclean: clean ## Removes executables and static library (common/*)
	@echo "\033[33mCleaning subdirectories...\\033[m"
	@MyList="$(DIRSWITHEXES)";for i in $$MyList;do \
		cd $$i; make clean; cd $(TOP) ;done
	@echo "\033[33mCleaning library (common)...\033[m"
	@cd $(TOP)/common/;make clean
	@echo "\033[33mCleaning devdoc...\033[m"
	@cd $(TOP)/devdoc/;rm -rf html

.PHONY: install
install: all ## Installs the binaries (after compiling). Usually requires sudo mode.
	@cp bin/* $(INSTALL_PATH)/

.PHONY: doc
doc: ## Compiles the documentation (using doxygen)
	@cd devdoc; doxygen; cd $(TOP)

