# Makefile is a part of the PYTHIA event generator.
# Copyright (C) 2017 Torbjorn Sjostrand.
# PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
# Please respect the MCnet Guidelines, see GUIDELINES for details.
# Author: Philip Ilten, September 2014.
#
# This is is the Makefile used to build PYTHIA examples on POSIX systems.
# Example usage is:
#     make main01
# For help using the make command please consult the local system documentation,
# i.e. "man make" or "make --help".

################################################################################
# VARIABLES: Definition of the relevant variables from the configuration script.
################################################################################

# Set the shell.
SHELL=/usr/bin/env bash

# Include the configuration.
-include Makefile.inc

# Handle GZIP support.
ifeq ($(GZIP_USE),true)
  CXX_COMMON+= -DGZIPSUPPORT -I$(GZIP_INCLUDE)
  CXX_COMMON+= -L$(GZIP_LIB) -Wl,-rpath,$(GZIP_LIB) -lz
endif

# Check distribution (use local version first, then installed version).
ifneq ("$(wildcard ../lib/libpythia8.*)","")
  PREFIX_LIB=../lib
  PREFIX_INCLUDE=../include
endif
CXX_COMMON:=-I$(PREFIX_INCLUDE) $(CXX_COMMON)
CXX_COMMON+= -L$(PREFIX_LIB) -Wl,-rpath,$(PREFIX_LIB) -lpythia8 -ldl 

################################################################################
# RULES: Definition of the rules used to build the PYTHIA examples.
################################################################################

# Rules without physical targets (secondary expansion for specific rules).
.SECONDEXPANSION:
.PHONY: all clean

# All targets 
all: VjetsPythia8

# The Makefile configuration.
Makefile.inc:
	$(error Error: PYTHIA must be configured, please run "./configure"\
                in the top PYTHIA directory)

# PYTHIA libraries.
$(PREFIX_LIB)/libpythia8.a :
	$(error Error: PYTHIA must be built, please run "make"\
                in the top PYTHIA directory)



VjetsPythia8: VjetsPythia8.o ANA_utils.o TruthPart.o TruthJets.o $(PREFIX_LIB)/libpythia8.a
ifeq ($(FASTJET3_USE)$(HEPMC2_USE)$(ROOT_USE),truetruetrue)
	$(CXX) VjetsPythia8.o ANA_utils.o TruthPart.o TruthJets.o -o $@ -w -I$(ROOT_INCLUDE) -I$(FASTJET3_INCLUDE) -I$(HEPMC2_INCLUDE) $(CXX_COMMON)\
	 `$(ROOTBIN)root-config --cflags`\
	 -Wl,-rpath,$(ROOT_LIB) `$(ROOT_BIN)root-config --glibs`\
	 -L$(HEPMC2_LIB) -Wl,-rpath,$(HEPMC2_LIB) -lHepMC\
	 -L$(FASTJET3_LIB) -Wl,-rpath,$(FASTJET3_LIB) -lfastjet
else
	@echo "Error: $@ requires ROOT, FASTJET3 and HEPMC2"
endif


VjetsPythia8.o: VjetsPythia8.cc VjetsPythia8.h ANA_utils.h TruthPart.h TruthJets.h $(PREFIX_LIB)/libpythia8.a 
ifeq ($(FASTJET3_USE)$(HEPMC2_USE)$(ROOT_USE),truetruetrue)
	${CXX} -c ${CXXFLAGS} VjetsPythia8.cc -w -I$(ROOT_INCLUDE) -I$(FASTJET3_INCLUDE) -I$(HEPMC2_INCLUDE) $(CXX_COMMON)\
	 `$(ROOTBIN)root-config --cflags`\
	 -Wl,-rpath,$(ROOT_LIB) `$(ROOT_BIN)root-config --glibs`\
	 -L$(HEPMC2_LIB) -Wl,-rpath,$(HEPMC2_LIB) -lHepMC\
	 -L$(FASTJET3_LIB) -Wl,-rpath,$(FASTJET3_LIB) -lfastjet	
else
	@echo "Error: $@ requires ROOT, FASTJET3 and HEPMC2"
endif


ANA_utils.o: ANA_utils.cc ANA_utils.h
ifeq ($(FASTJET3_USE)$(HEPMC2_USE)$(ROOT_USE),truetruetrue)
	${CXX} -c ${CXXFLAGS} ANA_utils.cc -w -I$(ROOT_INCLUDE) -I$(FASTJET3_INCLUDE) -I$(HEPMC2_INCLUDE) $(CXX_COMMON)\
	 `$(ROOTBIN)root-config --cflags`\
	 -Wl,-rpath,$(ROOT_LIB) `$(ROOT_BIN)root-config --glibs`\
	 -L$(HEPMC2_LIB) -Wl,-rpath,$(HEPMC2_LIB) -lHepMC\
	 -L$(FASTJET3_LIB) -Wl,-rpath,$(FASTJET3_LIB) -lfastjet	
else
	@echo "Error: $@ requires ROOT, FASTJET3 and HEPMC2"
endif

#   Compile the physics objets classes
#   ----------------------------------

TruthPart.o:  TruthPart.cc TruthPart.h 
	${CXX} -c ${CXXFLAGS} TruthPart.cc $(CXX_COMMON)

TruthJets.o:  TruthJets.cc TruthJets.h 
	${CXX} -c ${CXXFLAGS} TruthJets.cc $(CXX_COMMON)

#———————————————————————

# Internally used tests, without external dependencies.
test% : test%.cc $(PREFIX_LIB)/libpythia8.a
	$(CXX) $< -o $@ $(CXX_COMMON)

# Clean.
clean:
	@rm -f main[0-9][0-9]; rm -f out[0-9][0-9];\
	rm -f main[0-9][0-9][0-9]; rm -f out[0-9][0-9][0-9];\
	rm -f mymain[0-9][0-9]; rm -f myout[0-9][0-9];\
	rm -f test[0-9][0-9][0-9]; rm -f *.dat\
	rm -f weakbosons.lhe; rm -f Pythia8.promc; rm -f hist.root;\
	rm -f *~; rm -f \#*; rm -f core*; rm -f *Dct.*; rm -f *.so;
