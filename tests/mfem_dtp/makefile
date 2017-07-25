# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory. LLNL-CODE-734707. All Rights
# reserved. See files LICENSE and NOTICE for details.
#
# This file is part of CEED, a collection of benchmarks, miniapps, software
# libraries and APIs for efficient high-order finite element and spectral
# element discretizations for exascale applications. For more information and
# source code availability see http://github.com/ceed.
#
# The CEED research is supported by the Exascale Computing Project (17-SC-20-SC)
# a collaborative effort of two U.S. Department of Energy organizations (Office
# of Science and the National Nuclear Security Administration) responsible for
# the planning and preparation of a capable exascale ecosystem, including
# software, applications, hardware, advanced system engineering and early
# testbed platforms, in support of the nation's exascale computing imperative.

# Use the MFEM build directory
MFEM_DIR = ../../mfem
SRC =
BLD =
CONFIG_MK = $(MFEM_DIR)/config/config.mk

MFEM_LIB_FILE = mfem_is_not_built
-include $(CONFIG_MK)

# EXTRA_CXXFLAGS - extra flags that will be added to MFEM_CXXFLAGS. Some options
# for these flags are set below, based on MFEM_CXX. However, they can be
# overwritten at the make command line.
EXTRA_CXXFLAGS =

# Add any EXTRA_CXXFLAGS to MFEM_CXXFLAGS, which is part of MFEM_FLAGS.
MFEM_CXXFLAGS += $(EXTRA_CXXFLAGS)

DTP = dtp_baseline
ifeq ($(MFEM_USE_OCCA),YES)
DTP += dtp_occa
endif

.SUFFIXES:
.SUFFIXES: .o .cpp .mk
.PHONY: all clean clean-build clean-exec

# Remove built-in rule
%: %.cpp

# Replace the default implicit rule for *.cpp files
%: $(SRC)%.cpp $(MFEM_LIB_FILE) $(CONFIG_MK)
	$(MFEM_CXX) $(MFEM_FLAGS) $< -o $@ $(MFEM_LIBS)

all: $(DTP)

$(BLD)dtp_%: $(SRC)dtp_%.cpp $(MFEM_LIB_FILE) $(CONFIG_MK)
	$(MFEM_CXX) $(MFEM_FLAGS) $< -o $@ $(MFEM_LIBS)

# Generate an error message if the MFEM library is not built and exit
$(MFEM_LIB_FILE):
	$(error The MFEM library is not built)

clean: clean-build clean-exec

clean-build:
	rm -f *.o *~ $(SEQ_EXAMPLES) $(PAR_EXAMPLES)
	rm -rf *.dSYM *.TVD.*breakpoints

#---[ OCCA ]----------------------------
ifeq ($(MFEM_USE_OCCA),YES)
MFEM_REAL_DIR      := $(realpath $(MFEM_SOURCE_DIR))
OCCA_CACHE_DIR     ?= ${HOME}/.occa
OCCA_LIB_CACHE_DIR := $(OCCA_CACHE_DIR)/libraries

OKL_KERNELS        := $(realpath $(shell find $(MFEM_REAL_DIR) -type f -name '*.okl'))
OKL_CACHED_KERNELS := $(subst kernels/,,$(subst $(MFEM_REAL_DIR)/,$(OCCA_LIB_CACHE_DIR)/mfem/,$(OKL_KERNELS)))

# Cache kernels in the OCCA cache directory
.PHONY: cache-kernels
cache-kernels: $(OKL_CACHED_KERNELS)

.PHONY: clear-cache
clear-kernels: clear-mfem-kernels $(OKL_CACHED_KERNELS)

.PHONY: clear-cache
clear-mfem-kernels:
	@echo y | occa clear -l mfem

$(OCCA_LIB_CACHE_DIR)/mfem/fem/%.okl: $(MFEM_REAL_DIR)/fem/kernels/%.okl
	@echo "Caching: $(subst $(MFEM_REAL_DIR)/,,$<)"
	@occa cache mfem/$(subst $(OCCA_LIB_CACHE_DIR)/mfem/,,$(dir $@)) $<

$(OCCA_LIB_CACHE_DIR)/mfem/linalg/%.okl: $(MFEM_REAL_DIR)/linalg/kernels/%.okl
	@echo "Caching: $(subst $(MFEM_REAL_DIR)/,,$<)"
	@occa cache mfem/$(subst $(OCCA_LIB_CACHE_DIR)/mfem/,,$(dir $@)) $<
endif
#=======================================
