# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory. LLNL-CODE-XXXXXX. All Rights
# reserved. See file LICENSE for details.
#
# This file is part of CEED, a collection of benchmarks, miniapps, software
# libraries and APIs for efficient high-order finite element and spectral
# element discretizations for exascale applications. For more information and
# source code availability see http://github.com/ceed.
#
# The CEED research is supported by the Exascale Computing Project
# (17-SC-20-SC), a collaborative effort of two U.S. Department of Energy
# organizations (Office of Science and the National Nuclear Security
# Administration) responsible for the planning and preparation of a capable
# exascale ecosystem, including software, applications, hardware, advanced
# system engineering and early testbed platforms, in support of the nation's
# exascale computing imperative.

# Default build parameters
PROBLEM = 1
sol_p =
ir_order =
USE_MPI_WTIME =

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

ifeq ($(MFEM_USE_MPI),NO)
   $(error A parallel MFEM build is required.)
endif

.SUFFIXES:
.SUFFIXES: .o .cpp .mk
.PHONY: all clean clean-build clean-exec

# Remove built-in rule
%: %.cpp

# Replace the default implicit rule for *.cpp files
$(BLD)%: $(SRC)%.cpp $(MFEM_LIB_FILE) $(CONFIG_MK)
	$(MFEM_CXX) $(if $(PROBLEM),-DPROBLEM=$(PROBLEM)) $(if $(SOL_P),-DSOL_P=$(SOL_P)) $(if $(IR_ORDER),-DIR_ORDER=$(IR_ORDER)) $(MFEM_FLAGS) $< -o $@ $(MFEM_LIBS)

# Create rules for building executables
exe_list :=
define make_bp_rule
exe := $(BLD)bp$(PROBLEM)$(if $(1),_solp$(1))$(if $(2),_irorder$(2))
exe_list += $(exe)
$(exe): $(SRC)bp.cpp $(MFEM_LIB_FILE) $(CONFIG_MK)
	$(MFEM_CXX) \
	$(if $(PROBLEM),-DPROBLEM=$(PROBLEM)) \
	$(if $(1),-DSOL_P=$(1)) \
	$(if $(2),-DIR_ORDER=$(2)) \
	$(if $(USE_MPI_WTIME),-DUSE_MPI_WTIME) \
	$(MFEM_FLAGS) $$< -o $$@ $(MFEM_LIBS)
endef
comma = ,
args_list := $(join $(addsuffix /,$(sol_p)),$(ir_order))
$(foreach args,$(args_list),\
	$(eval $$(eval $$(call make_bp_rule,$(subst /,$(comma),$(args))))))

all: $(exe_list)

# Generate an error message if the MFEM library is not built and exit
$(MFEM_LIB_FILE):
	$(error The MFEM library is not built)

# Clean targets
clean: clean-build clean-exec
clean-build:
	rm -f *.o *~ $(BLD)*
	rm -rf *.dSYM *.TVD.*breakpoints
clean-exec:
	@rm -f mesh.* sol.*
