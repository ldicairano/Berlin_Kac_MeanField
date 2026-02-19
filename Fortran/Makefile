FC      = gfortran
FFLAGS  = -O3 -march=native -Wall -Wextra -fcheck=all -frecursive
LDFLAGS =

MODDIR  = build/mod
OBJDIR  = build/obj
BINDIR  = build/bin

FFLAGS += -J$(MODDIR) -I$(MODDIR)

EXE = $(BINDIR)/BKMF_RATTLE

SRCS = \
  mod_kinds.f90 \
  mod_types.f90 \
  mod_input.f90 \
  mod_random.f90 \
  mod_mf_spherical_model.f90 \
  mod_rattle_spherical.f90 \
  mod_init_spherical.f90 \
  mod_io_restart.f90 \
  mod_io_observables.f90 \
  main_berlin_kac_rattle.f90

OBJS = $(patsubst %.f90,$(OBJDIR)/%.o,$(SRCS))

all: dirs $(EXE)

dirs:
	@mkdir -p $(MODDIR) $(OBJDIR) $(BINDIR)

$(EXE): $(OBJS)
	$(FC) $(FFLAGS) $^ -o $@ $(LDFLAGS)

$(OBJDIR)/%.o: %.f90 | dirs
	$(FC) $(FFLAGS) -c $< -o $@

clean:
	rm -rf build *.o *.mod

# Example run (override with CLI at runtime)
run: all
	$(EXE) -i input.inp 102 1 1

.PHONY: all clean run dirs
