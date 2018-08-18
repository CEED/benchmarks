# Copyright (c) 2010, Lawrence Livermore National Security, LLC. Produced at the
# Lawrence Livermore National Laboratory. LLNL-CODE-443211. All Rights reserved.
# See file COPYRIGHT for details.
#
# This file is part of the MFEM library. For more information and source code
# availability see http://mfem.org.
#
# MFEM is free software; you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License (as published by the Free
# Software Foundation) version 2.1 dated February 1999.

ifndef MFEM_DIR
   $(error MFEM_DIR is not defined)
endif
override SRC := $(dir $(lastword $(MAKEFILE_LIST)))
# BLD is where the executables will be placed; ensure it ends with /
override BLD := $(patsubst %//,%/,$(or $(BLD),$(SRC))/)


# The MFEM configuration file
CONFIG_MK = $(MFEM_DIR)/config/config.mk

MFEM_LIB_FILE = mfem_is_not_built
-include $(CONFIG_MK)

TESTS = bp_main

.SUFFIXES:
.SUFFIXES: .o .cpp .mk
.PHONY: all clean clean-build clean-exec

# Remove built-in rules
%.o: %.cpp
%: %.cpp
%: %.o

# Replace the default implicit rule for *.cpp files
$(BLD)%: $(SRC)%.cpp $(MFEM_LIB_FILE) $(CONFIG_MK)
	$(MFEM_CXX) $(MFEM_FLAGS) $< -o $@ $(MFEM_LIBS)

all: $(addprefix $(BLD),$(TESTS))

# Generate an error message if the MFEM library is not built and exit
$(MFEM_LIB_FILE):
	$(error The MFEM library is not built)

clean: clean-build clean-exec

clean-build:
	$(RM) -r $(addprefix $(BLD),$(TESTS) *.dSYM)

clean-exec:
	$(RM) $(BLD)mesh.* $(BLD)sol.*
