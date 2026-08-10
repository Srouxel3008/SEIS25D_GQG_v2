# Seis2D_GQG_FWI portable Intel oneAPI build 
#
# Run from the repository root:
#   make -C src CONFIG=release
#   make -C src CONFIG=debug
#   make -C src clean 
#
# The default build uses LU/PARDISO and must not require MUMPS. Optional MUMPS
# support is provided for HPC installations.
# ----------------------------- User options -----------------------------

FC            := mpiifx
CONFIG        ?= release
GF_STORAGE    ?= sp
USE_MUMPS     ?= 0
MUMPS_DIR     ?=
MUMPS_FFLAGS  ?=
MUMPS_LIBS    ?=

# Command-line assignments override these values, for example:
#   make -C src CONFIG=debug FC=/opt/intel/oneapi/mpi/latest/bin/mpiifx
#   make -C src USE_MUMPS=1 MUMPS_DIR=/opt/mumps \
#       MUMPS_LIBS='-L/opt/mumps/lib -lzmumps -lmumps_common -lpord ...'

# ---------------------------- Project layout ----------------------------

ROOT_DIR      := $(abspath $(CURDIR)/..)
MODULE_SRC    := $(CURDIR)/modules
BUILD_DIR     ?= $(ROOT_DIR)/build
OBJECT_DIR    := $(BUILD_DIR)/objects
MODULE_DIR    := $(BUILD_DIR)/modules
BIN_DIR       ?= $(ROOT_DIR)/bin
TARGET        := $(BIN_DIR)/Seis2D_GQG_FWI

# ------------------------- Configuration checks -------------------------

VALID_CONFIGS := release debug
VALID_STORAGE := sp dp
VALID_SWITCH  := 0 1

ifeq ($(filter $(CONFIG),$(VALID_CONFIGS)),)
  $(error CONFIG must be one of: $(VALID_CONFIGS))
endif

ifeq ($(filter $(GF_STORAGE),$(VALID_STORAGE)),)
  $(error GF_STORAGE must be one of: $(VALID_STORAGE))
endif

ifeq ($(filter $(USE_MUMPS),$(VALID_SWITCH)),)
  $(error USE_MUMPS must be 0 or 1)
endif

# ------------------------- Compiler/linker flags ------------------------

CPPFLAGS      := -DSEIS_BUILD='"$(shell echo $(CONFIG) | tr a-z A-Z)"'
COMMON_FLAGS  := -fpp -qopenmp -qmkl -align array64byte -heap-arrays \
                 -fstack-protector -traceback
MODULE_FLAGS  := -I$(MODULE_DIR) -module $(MODULE_DIR)
LINK_FLAGS    := -qopenmp -qmkl=parallel

ifeq ($(CONFIG),release)
  CONFIG_FLAGS := -O3 -g
endif

ifeq ($(CONFIG),debug)
  CONFIG_FLAGS := -O0 -g -warn all -check all -fpe0 -ftrapuv
endif

ifeq ($(GF_STORAGE),dp)
  CPPFLAGS += -DGF_STORAGE_DP
endif

ifeq ($(USE_MUMPS),1)
  CPPFLAGS += -DHAVE_MUMPS=1

  ifneq ($(strip $(MUMPS_DIR)),)
    MUMPS_FFLAGS += -I$(MUMPS_DIR)/include
    MUMPS_LIBS   += -L$(MUMPS_DIR)/lib
  endif

  ifeq ($(strip $(MUMPS_FFLAGS)),)
    $(error USE_MUMPS=1 requires MUMPS_DIR or MUMPS_FFLAGS)
  endif

  ifeq ($(strip $(MUMPS_LIBS)),)
    $(error USE_MUMPS=1 requires site-specific MUMPS_LIBS)
  endif
else
  CPPFLAGS += -DHAVE_MUMPS=0
endif

FFLAGS  := $(CONFIG_FLAGS) $(COMMON_FLAGS) $(CPPFLAGS) $(MODULE_FLAGS) \
           $(MUMPS_FFLAGS)
LDLIBS  := $(LINK_FLAGS) $(MUMPS_LIBS)

# ----------------------------- Source files -----------------------------

MODULE_NAMES := \
  err_mpi_mod \
  hardware_mod \
  constant_mod \
  shared_mod \
  gridtype_mod \
  input_mod \
  output_mod \
  ky_sampling_mod \
  GLL_mod \
  grid_mod \
  boundaries_mod \
  scalers_mod \
  misfit_mod \
  partial_derivatives_mod \
  Frechet_mod \
  gradient_mod \
  stiffness_assembly_mod \
  solvers_mod \
  F_modeling \
  lbfgs_mod \
  lbfgs_lsearch

MODULE_OBJECTS := $(addprefix $(OBJECT_DIR)/,$(addsuffix .o,$(MODULE_NAMES)))
MAIN_OBJECT    := $(OBJECT_DIR)/Seis2D_GQG_FWI.o
OBJECTS        := $(MODULE_OBJECTS) $(MAIN_OBJECT)

# --------------------------- Primary targets ----------------------------

.PHONY: all release debug clean distclean check-env show-config

all: check-env $(TARGET)

release:
	$(MAKE) CONFIG=release all

debug:
	$(MAKE) CONFIG=debug all

$(TARGET): $(OBJECTS) | $(BIN_DIR)
	$(FC) $(OBJECTS) $(LDLIBS) -o $@

$(OBJECT_DIR)/%.o: $(MODULE_SRC)/%.f90 | $(OBJECT_DIR) $(MODULE_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(MAIN_OBJECT): $(CURDIR)/Seis2D_GQG_FWI.f90 | $(OBJECT_DIR) $(MODULE_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJECT_DIR) $(MODULE_DIR) $(BIN_DIR):
	mkdir -p $@

# -------------------- Fortran module dependencies ----------------------
# Explicit prerequisites make `make -j` safe and document the active graph.

ERR       := $(OBJECT_DIR)/err_mpi_mod.o
HARDWARE  := $(OBJECT_DIR)/hardware_mod.o
CONSTANT  := $(OBJECT_DIR)/constant_mod.o
SHARED    := $(OBJECT_DIR)/shared_mod.o
GRIDTYPE  := $(OBJECT_DIR)/gridtype_mod.o
INPUT     := $(OBJECT_DIR)/input_mod.o
OUTPUT    := $(OBJECT_DIR)/output_mod.o
KY        := $(OBJECT_DIR)/ky_sampling_mod.o
GLL       := $(OBJECT_DIR)/GLL_mod.o
GRID      := $(OBJECT_DIR)/grid_mod.o
BOUNDARY  := $(OBJECT_DIR)/boundaries_mod.o
SCALERS   := $(OBJECT_DIR)/scalers_mod.o
MISFIT    := $(OBJECT_DIR)/misfit_mod.o
PARTIAL   := $(OBJECT_DIR)/partial_derivatives_mod.o
FRECHET   := $(OBJECT_DIR)/Frechet_mod.o
GRADIENT  := $(OBJECT_DIR)/gradient_mod.o
STIFFNESS := $(OBJECT_DIR)/stiffness_assembly_mod.o
SOLVERS   := $(OBJECT_DIR)/solvers_mod.o
FMODEL    := $(OBJECT_DIR)/F_modeling.o
LBFGS     := $(OBJECT_DIR)/lbfgs_mod.o
LSEARCH   := $(OBJECT_DIR)/lbfgs_lsearch.o

$(SHARED): $(CONSTANT)
$(GRIDTYPE): $(SHARED)
$(INPUT): $(ERR) $(CONSTANT)
$(OUTPUT): $(SHARED)
$(KY): $(CONSTANT)
$(GRID): $(HARDWARE) $(SHARED) $(OUTPUT) $(KY) $(GLL) $(GRIDTYPE)
$(BOUNDARY): $(CONSTANT) $(ERR)
$(SCALERS): $(CONSTANT) $(GRID) $(GRIDTYPE)
$(MISFIT): $(CONSTANT) $(GRID) $(GRIDTYPE) $(OUTPUT) $(SHARED)
$(FRECHET): $(HARDWARE) $(SHARED) $(OUTPUT) $(PARTIAL) $(GRIDTYPE) $(CONSTANT)
$(GRADIENT): $(FRECHET) $(MISFIT) $(SHARED) $(GRID) $(OUTPUT) $(SCALERS) $(GRIDTYPE)
$(STIFFNESS): $(HARDWARE) $(SHARED) $(BOUNDARY) $(OUTPUT) $(GRIDTYPE) $(ERR) $(CONSTANT)
$(SOLVERS): $(BOUNDARY) $(CONSTANT) $(ERR) $(GRIDTYPE) $(HARDWARE) $(OUTPUT) $(SHARED) $(STIFFNESS)
$(FMODEL): $(HARDWARE) $(SHARED) $(BOUNDARY) $(SOLVERS) $(OUTPUT) $(GRIDTYPE) $(ERR) $(CONSTANT) $(STIFFNESS) $(SCALERS)
$(LBFGS): $(CONSTANT) $(FRECHET) $(GRADIENT) $(GRID) $(GRIDTYPE) $(SHARED)
$(LSEARCH): $(FMODEL) $(GRADIENT) $(GRIDTYPE) $(HARDWARE) $(MISFIT) $(OUTPUT) $(SHARED) $(CONSTANT)
$(MAIN_OBJECT): $(ERR) $(HARDWARE) $(CONSTANT) $(SHARED) $(GRIDTYPE) $(INPUT) \
                $(OUTPUT) $(KY) $(GLL) $(GRID) $(BOUNDARY) $(SCALERS) $(MISFIT) \
                $(PARTIAL) $(FRECHET) $(GRADIENT) $(STIFFNESS) $(SOLVERS) \
                $(FMODEL) $(LBFGS) $(LSEARCH)

# -------------------------- Utility targets -----------------------------

check-env:
	@command -v $(FC) >/dev/null 2>&1 || { \
	  echo "ERROR: $(FC) was not found. Initialize Intel oneAPI and Intel MPI." >&2; \
	  exit 2; \
	}

show-config:
	@echo "FC=$(FC)"
	@echo "CONFIG=$(CONFIG)"
	@echo "GF_STORAGE=$(GF_STORAGE)"
	@echo "USE_MUMPS=$(USE_MUMPS)"
	@echo "BUILD_DIR=$(BUILD_DIR)"
	@echo "MODULE_DIR=$(MODULE_DIR)"
	@echo "TARGET=$(TARGET)"

clean:
	rm -rf $(OBJECT_DIR) $(MODULE_DIR)
	rm -f $(TARGET)

distclean:
	rm -rf $(BUILD_DIR)
	rm -rf $(BIN_DIR)
