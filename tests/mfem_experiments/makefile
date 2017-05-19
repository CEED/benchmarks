# This file is part of CEED. For more details, see exascaleproject.org.

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

BPS = mass2d
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

# Special set of rules for building mass2d
mass2d_DEF := $(if $(problem),-DPROBLEM=$(problem),)
mass2d_DEF += $(if $(geom),-DGEOM=$(geom),)
mass2d_DEF += $(if $(mesh_p),-DMESH_P=$(mesh_p),)
mass2d_DEF += $(if $(use_mpi_wtime),-DUSE_MPI_WTIME,)
mass2d_DEF := $(strip $(mass2d_DEF))
define make_mass2d_rule
$(BLD)mass2d$(3): $(SRC)mass2d.cpp $(MFEM_LIB_FILE) $(CONFIG_MK)
	cp -fp $(SRC)mass2d.cpp $(BLD)mass2d$(3).cpp
	$(MFEM_CXX) $(mass2d_DEF) $(if $(1),-DSOL_P=$(1),) \
	$(if $(2),-DIR_ORDER=$(2),) $(MFEM_FLAGS) \
	$(BLD)mass2d$(3).cpp -o $(BLD)mass2d$(3) $(MFEM_LIBS)
endef
ir_ord := $(wordlist 1,$(words $(sol_p)),$(ir_order))
exe_sfx := $(wordlist 1,$(words $(sol_p)),$(exe_suffix))
comma = ,
colon = :
list := $(join $(addsuffix /,$(sol_p)),$(ir_ord))
list := $(join $(addsuffix /,$(list)),$(exe_sfx))
$(foreach i,$(if $(list),$(list),/),\
   $(eval $$(eval $$(call make_mass2d_rule,$(subst /,$(comma),$(i))))))
$(if $(exe_sfx),$(eval .PHONY$(colon) mass2d))
$(if $(exe_sfx),$(eval mass2d$(colon) $(addprefix $(BLD)mass2d,$(exe_sfx))))


# Generate an error message if the MFEM library is not built and exit
$(MFEM_LIB_FILE):
	$(error The MFEM library is not built)

clean: clean-build clean-exec

clean-build:
	rm -f *.o *~ $(BPS)
	rm -rf *.dSYM *.TVD.*breakpoints

clean-exec:
	@rm -f refined.mesh mesh.* sol.*
