# COMPILERS OPTIONS

# Aditional optimization flags
LOCAL_OPTIMFLAGS = 

# Replace flag
NEW_FLAG = 
OLD_FLAG = 

# Aditional Debug and Profile build flags
ifeq ($(DEBUGVERSION),1)
   LOCAL_CXXFLAGS     +=
else
   ifeq ($(PROFILEVERSION),1)
      LOCAL_CXXFLAGS     += $(LOCAL_OPTIMFLAGS)
   else
      LOCAL_CXXFLAGS     += $(LOCAL_OPTIMFLAGS)
   endif
endif

INCDEFS        := -include globaldefs.h
# Additional linking and include flags
# -L: FOLDER LIBRARY
LFLAGS=#-L/usr/lib
LFLAGS=#-L.

# -l: LIBRARY
lLIBS=#-lgfnl_util
lLIBS= #-lapr-1 -laprutil-1 -llog4cxx

# -I: INCLUDES
IFLAGS= $(INCDEFS)

#################################################################################
#################################################################################

# FILES
CPPOBJS=$(patsubst %.cpp,%.o,$(wildcard *.cpp))
OBJS= $(CPPOBJS)
HEADERS=$(wildcard *.h)
SOURCES=$(wildcard *.cpp)
LIBOBJS=$(shell echo $(OBJS))
COMMONSTATICLIBNAME=$(shell cd ../..; pwd | sed -e 's;\(.*\)/\(.*\);\2;')
STATICLIB="lib$(COMMONSTATICLIBNAME).a"

ARCHIVE      := ar
ARCHFLAG     := -rc

ifndef MAKE_SRC
    CXX = $(shell $(MAKE) -s -C .. loadVarCXX )
    CXXFLAGS = $(shell $(MAKE) -s -C .. loadVarCXXFLAGS )
endif

CXXFLAGS += $(LOCAL_CXXFLAGS)

# Replace the old flag with the new
TMPFLAGS := $(filter-out $(OLD_FLAG),$(CXXFLAGS))
CXXFLAGS := $(TMPFLAGS) $(NEW_FLAG)

TARGET=target

$(TARGET): $(OBJS) $(HEADERS) $(SOURCES) 
	@echo "\033[32mArchiving $(STATICLIB)\\033[m"
	@$(ARCHIVE) $(ARFLAGS) $(STATICLIB) $(LIBOBJS)

$(CPPOBJS): %.o: %.cpp
	@echo "\033[32m   Compiling $@\\033[m"
	$(CXX) $(CXXFLAGS) $(IFLAGS) -c $<  -o $@  

.PHONY: clean
clean:
	$(info CLEANING ALL)
	@$(RM) -f $(TARGET) $(OBJS) $(STATICLIB) 2>/dev/null || true

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

list:
	@echo $(OBJS)
	@echo $(HEADERS)

listobjs:
	@echo $(LIBOBJS)

