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

BPS = mass
ifeq ($(MFEM_USE_MPI),NO)
   $(error A parallel MFEM build is required.)
endif

.SUFFIXES:
.SUFFIXES: .o .cpp .mk
.PHONY: all clean clean-build clean-exec

# Remove built-in rules
%: %.cpp
%.o: %.cpp

# Rule for building mass-lib.o
mass-lib-src = mass-quad.c mass-hex.c diffusion-quad.c
$(BLD)mass-lib.o: $(SRC)mass-lib.cpp $(addprefix $(SRC),$(mass-lib-src))
	$(MFEM_CXX) -c $(MFEM_FLAGS) $< -o $@

# Replace the default implicit rule for *.cpp files
$(BLD)%: $(SRC)%.cpp $(MFEM_LIB_FILE) $(CONFIG_MK)
	$(MFEM_CXX) $(MFEM_FLAGS) $< -o $@ $(MFEM_LIBS)

all: $(BPS)

# Special set of rules for building mass
mass_DEF := $(if $(problem),-DPROBLEM=$(problem),)
mass_DEF += $(if $(geom),-DGEOM=$(geom),)
mass_DEF += $(if $(mesh_p),-DMESH_P=$(mesh_p),)
mass_DEF += $(if $(use_mpi_wtime),-DUSE_MPI_WTIME,)
mass_DEF := $(strip $(mass_DEF) -I$(SRC).)
define make_rule
$(BLD)$(1)$(4): $(SRC)$(1).cpp $(BLD)$(1)-lib.o $(MFEM_LIB_FILE) $(CONFIG_MK)
	cp -fp $(SRC)$(1).cpp $(BLD)$(1)$(4).cpp
	$(MFEM_CXX) $($(1)_DEF) $(if $(2),-DSOL_P=$(2),) \
	$(if $(3),-DIR_ORDER=$(3),) $(MFEM_FLAGS) \
	$(BLD)$(1)$(4).cpp $(BLD)$(1)-lib.o -o $(BLD)$(1)$(4) $(MFEM_LIBS)
endef
ir_ord := $(wordlist 1,$(words $(sol_p)),$(ir_order))
exe_sfx := $(wordlist 1,$(words $(sol_p)),$(exe_suffix))
comma = ,
colon = :
list := $(join $(addsuffix /,$(sol_p)),$(ir_ord))
list := $(join $(addsuffix /,$(list)),$(exe_sfx))
$(foreach i,$(if $(list),$(list),/),\
   $(eval $$(eval $$(call make_rule,mass,$(subst /,$(comma),$(i))))))
$(if $(exe_sfx),$(eval .PHONY$(colon) mass))
$(if $(exe_sfx),$(eval mass$(colon) $(addprefix $(BLD)mass,$(exe_sfx))))


# Generate an error message if the MFEM library is not built and exit
$(MFEM_LIB_FILE):
	$(error The MFEM library is not built)

clean: clean-build clean-exec

clean-build:
	rm -f *.o *~ $(BPS)
	rm -rf *.dSYM *.TVD.*breakpoints

clean-exec:
	@rm -f refined.mesh mesh.* sol.*
