# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory. LLNL-CODE-734707. All Rights
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

BPS = bp1_v1
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

all: $(BPS)

# Special set of rules for building bp1_v1
bp1_v1_DEF := $(if $(problem),-DPROBLEM=$(problem),)
bp1_v1_DEF += $(if $(geom),-DGEOM=$(geom),)
bp1_v1_DEF += $(if $(mesh_p),-DMESH_P=$(mesh_p),)
bp1_v1_DEF += $(if $(vdim),-DVDIM=$(vdim),)
bp1_v1_DEF += $(if $(vec_layout),-DVEC_LAYOUT=$(vec_layout),)
bp1_v1_DEF += $(if $(use_mpi_wtime),-DUSE_MPI_WTIME,)
bp1_v1_DEF := $(strip $(bp1_v1_DEF))
define make_bp1_v1_rule
$(BLD)bp1_v1$(3): $(SRC)bp1_v1.cpp $(MFEM_LIB_FILE) $(CONFIG_MK)
	cp -fp $(SRC)bp1_v1.cpp $(BLD)bp1_v1$(3).cpp
	$(MFEM_CXX) $(bp1_v1_DEF) $(if $(1),-DSOL_P=$(1),) \
	$(if $(2),-DIR_ORDER=$(2),) $(MFEM_FLAGS) \
	$(BLD)bp1_v1$(3).cpp -o $(BLD)bp1_v1$(3) $(MFEM_LIBS)
endef
ir_ord := $(wordlist 1,$(words $(sol_p)),$(ir_order))
exe_sfx := $(wordlist 1,$(words $(sol_p)),$(exe_suffix))
comma = ,
colon = :
list := $(join $(addsuffix /,$(sol_p)),$(ir_ord))
list := $(join $(addsuffix /,$(list)),$(exe_sfx))
$(foreach i,$(if $(list),$(list),/),\
   $(eval $$(eval $$(call make_bp1_v1_rule,$(subst /,$(comma),$(i))))))
$(if $(exe_sfx),$(eval .PHONY$(colon) bp1_v1))
$(if $(exe_sfx),$(eval bp1_v1$(colon) $(addprefix $(BLD)bp1_v1,$(exe_sfx))))


# Generate an error message if the MFEM library is not built and exit
$(MFEM_LIB_FILE):
	$(error The MFEM library is not built)

clean: clean-build clean-exec

clean-build:
	rm -f *.o *~ $(BPS)
	rm -rf *.dSYM *.TVD.*breakpoints

clean-exec:
	@rm -f refined.mesh mesh.* sol.*
