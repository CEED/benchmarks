# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory. LLNL-CODE-734707. All Rights
# reserved. See files LICENSE and NOTICE for details.
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

# Use the MFEM build directory
MFEM_DIR = ../../mfem
SRC =
BLD =
CONFIG_MK = $(or $(wildcard\
 $(MFEM_DIR)/config/config.mk),$(MFEM_DIR)/share/mfem/config.mk)

MFEM_LIB_FILE = mfem_is_not_built
-include $(CONFIG_MK)

EX1 = $(BLD)ex1p_amgx_precon
EX2 = $(BLD)ex1p_amgx_solver
EX3 = $(BLD)ex1p_amgx_agg_solver
EX4 = $(BLD)ex1p_amgx_agg_precon
EX5 = $(BLD)ex1p_amgx_precon_l1_jacobi
EX6 = $(BLD)ex1p_amgx_solver_l1_jacobi
EX7 = $(BLD)ex1p_amgx_agg_solver_l1_jacobi
EX8 = $(BLD)ex1p_amgx_agg_precon_l1_jacobi
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
	$(MFEM_CXX) $(MFEM_FLAGS) $< -o $@ $(MFEM_LIBS)

all: $(EX1) $(EX2) $(EX3) $(EX4) $(EX5) $(EX6) $(EX7) $(EX8)

# Generate an error message if the MFEM library is not built and exit
$(MFEM_LIB_FILE):
	$(error The MFEM library is not built)

clean:
	rm -f *.o *~ $(EX1)
	rm -rf *.dSYM *.TVD.*breakpoints
