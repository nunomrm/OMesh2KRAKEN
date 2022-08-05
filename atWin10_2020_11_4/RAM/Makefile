
# just can just type
# gfortran ram1.5.f
# mv a.out ram.exe

# looks like the rest here was for a Windows system

FC=C:/Fortran/G77/bin/g77.exe

FFLAGS=-Os

EXECUTABLES = ram.exe

RAM_OBJECTS = ram.o

all: $(EXECUTABLES)
	@echo "RAM built"
	@echo "*************"
	@echo " "

install: all
	for f in $(EXECUTABLES) ; do \
		echo "----- Installing $$f"; cp -p $$f ../bin; \
	done

clean:
	-rm -f *.o *.mod *.exe

# Declare implicit rules (for make utilities with missing or broken build rules)

.SUFFIXES: .a .o .mod .f .f90

.f.o:
	$(FC) -c $(FFLAGS) $< 

.f90.o:
	$(FC) -c $(FFLAGS) $< 

.f90.mod:
	$(FC) -c $(FFLAGS) $< 

# Declare explicit dependencies (and explicit build rules as needed)

ram.exe: $(RAM_OBJECTS)
	$(FC) -o $@ $(FFLAGS) $(RAM_OBJECTS)

