#-*-makefile-*-

############### Specify the compiler you want
# COMPILER options --- gnu, pgi
# e.g. make COMPILER=gnu

NETCDFLIB  =  $(NETCDF)/lib
NETCDFINC  =  $(NETCDF)/include


GNUFLAGS  = -O -ffree-form -fopenmp -J.
PGIFLAGS = -fast -Mfreeform -mp -DDEBUG
PGIFLAGS = -fast -w -V -Mfreeform -mp

LIBRARY = -L$(NETCDFLIB) -L$(NETCDF_FORTRAN)/lib -lnetcdff -lnetcdf

INCLUDE = -I. -I$(NETCDFINC) -I$(NETCDF_FORTRAN)/include


RM_LIST = *.o *.exe *.mod *.a core

# Misc. commands
RM  = /bin/rm -f
AR  = /usr/bin/ar r
CPP = /usr/bin/cpp
LN  = /bin/ln -sf
CD  = cd

##############  You shouldn't need to change anything below here
#+++++++++++++++++++++++++++++++++++++++++++++++++++++

ifeq ($(COMPILER),pgi)
  FC=pgf90
  FFLAGS=$(PGIFLAGS)
  LDFLAGS = $(LIBRARY) -mp
endif

ifeq ($(COMPILER),gnu)
  FC=gfortran
  FFLAGS=$(GNUFLAGS)
  LDFLAGS=$(LIBRARY) -fopenmp
  CPPFLAGS=-DNDEBUG -DpgiFortran
endif

OBJS = wacomm_main.o wacomm_core.o comblk.o
MODS = comblk.mod

wacomm_$(COMPILER).exe:	$(OBJS)
		$(FC) -o $@ $(OBJS) $(LDFLAGS)

wacomm_main.o:	$(MODS)
wacomm_core.o:	$(MODS)

%.o:            %.F90
		$(FC) -c $(FFLAGS) $(INCLUDE) $<

%.mod:          %.F90
		$(FC) -c $(FFLAGS) $<

clean:
		$(RM) $(RM_LIST)

