# PYTHIA configuration file.
# Generated on Mon Oct 16 13:41:24 EDT 2017 with the user supplied options:
# --with-lhapdf6=/usr/local/opt/lhapdf

# Install directory prefixes.
PREFIX_BIN=/Users/HugoBeauchemin/Documents/beauche/Work/TheoryWork/QCDSystematics/Pythia/pythia8230/bin
PREFIX_INCLUDE=/Users/HugoBeauchemin/Documents/beauche/Work/TheoryWork/QCDSystematics/Pythia/pythia8230/include
PREFIX_LIB=/Users/HugoBeauchemin/Documents/beauche/Work/TheoryWork/QCDSystematics/Pythia/pythia8230/lib
PREFIX_SHARE=/Users/HugoBeauchemin/Documents/beauche/Work/TheoryWork/QCDSystematics/Pythia/pythia8230/share/Pythia8

# Compilation flags (see ./configure --help for further documentation).
ENABLE_SHARED=false
CXX=g++
CXX_COMMON=-O2  -pedantic -W -Wall -Wshadow -fPIC
CXX_SHARED=-dynamiclib
CXX_SONAME=-Wl,-dylib_install_name,@rpath/
LIB_SUFFIX=.dylib

# EVTGEN configuration.
EVTGEN_USE=false
EVTGEN_BIN=
EVTGEN_INCLUDE=./
EVTGEN_LIB=./

# FASTJET3 configuration.
FASTJET3_USE=true
FASTJET3_BIN=/usr/local/opt/fastjet/bin/
FASTJET3_INCLUDE=/usr/local/opt/fastjet/include
FASTJET3_LIB=/usr/local/opt/fastjet/lib

# HEPMC2 configuration.
HEPMC2_USE=true
HEPMC2_BIN=
HEPMC2_INCLUDE=/usr/local/opt/hepmc/include
HEPMC2_LIB=/usr/local/opt/hepmc/lib

# HEPMC3 configuration.
HEPMC3_USE=false
HEPMC3_BIN=
HEPMC3_INCLUDE=./
HEPMC3_LIB=./

# LHAPDF5 configuration.
LHAPDF5_USE=false
LHAPDF5_BIN=
LHAPDF5_INCLUDE=./
LHAPDF5_LIB=./

# LHAPDF6 configuration.
LHAPDF6_USE=true
LHAPDF6_BIN=/usr/local/opt/lhapdf/bin/
LHAPDF6_INCLUDE=/usr/local/opt/lhapdf/include
LHAPDF6_LIB=/usr/local/opt/lhapdf/lib

# POWHEG configuration.
POWHEG_USE=false
POWHEG_BIN=
POWHEG_INCLUDE=./
POWHEG_LIB=./

# PROMC configuration.
PROMC_USE=false
PROMC_BIN=
PROMC_INCLUDE=./
PROMC_LIB=./

# ROOT configuration.
ROOT_USE=true
ROOT_BIN=/usr/local/opt/root/bin/
ROOT_INCLUDE=/usr/local/opt/root/include
ROOT_LIB=/usr/local/opt/root/lib

# GZIP configuration.
GZIP_USE=false
GZIP_BIN=
GZIP_INCLUDE=./
GZIP_LIB=./

# BOOST configuration.
BOOST_USE=true
BOOST_BIN=
BOOST_INCLUDE=/usr/local/opt/boost/include
BOOST_LIB=/usr/local/opt/boost/lib

# PYTHON configuration.
PYTHON_USE=false
PYTHON_BIN=
PYTHON_INCLUDE=./
PYTHON_LIB=./
