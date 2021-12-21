FC=gfortran

##$ Extends the implicit support of the Makefile to .f90 files
.SUFFIXES: .f90

GLOB_INC:=$(shell pkg-config --cflags scifor)
GLOB_LIB:=$(shell pkg-config --libs scifor)


OBJS     = COMMON_VARS.o AUX_FUNX.o


ifeq ($(FC),ifort)
FFLAG=-O2 -ftz
DFLAG=-p -O0 -g -fpe0 -warn -warn errors -debug extended -traceback -check all,noarg_temp_created
endif
ifeq ($(FC),gfortran)
FFLAG = -O2 -ffree-line-length-none
DFLAG = -O2 -p -g -fimplicit-none -Wsurprising  -Waliasing -fwhole-file -fcheck=all -pedantic -fbacktrace -ffree-line-length-none
endif


all: FLAG:=${FFLAG}
all: $(OBJS)
	${FC} ${FLAG} ${OBJS} test_build_sector_ImpBath.f90 -o test/test_build_sector_ImpBath ${GLOB_INC} ${GLOB_LIB}

debug: FLAG:=${DFLAG}
debug: ${OBJS}
	${FC} ${FLAG} ${OBJS} test_build_sector_ImpBath.f90 -o test/test_build_sector_ImpBath ${GLOB_INC} ${GLOB_LIB}

clean: 
	@echo "Cleaning:"
	@rm -fv *.mod *.o *~ test/test_build_sector_ImpBath
	@echo ""

.f90.o:	
	$(FC) $(FLAG) -c $< ${GLOB_INC}
