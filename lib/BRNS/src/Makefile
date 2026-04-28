# BRNS Shared Library Makefile
# Builds brns.so for use with OpenFOAM solvers

F77 = gfortran
FFLAGS = -O4 -x f77-cpp-input -I. -fPIC
LDFLAGS = -L. -O4 -shared -s

# Object files needed for brns.so
OBJS = basic.o biogeo.o boundaries.o drivervalues.o gaussj.o invokebrns.o \
       jacobian.o LUBKSB.o LUDCMP.o MPROVE.o newtonsub.o residual.o \
       switches.o parameters.o limits.o rates.o

# Default target
all: brns.so

# Build shared library
brns.so: $(OBJS)
	$(F77) $(LDFLAGS) $(OBJS) -o brns.so

# Pattern rules
%.o: %.f
	$(F77) $(FFLAGS) -c $< -o $@

%.o: %.F
	$(F77) $(FFLAGS) -c $< -o $@

# Special rules for files in BrnsDll subdirectory
parameters.o: BrnsDll/parameters.f
	$(F77) $(FFLAGS) -c BrnsDll/parameters.f -o parameters.o

invokebrns.o: BrnsDll/invokebrns.f
	$(F77) $(FFLAGS) -c BrnsDll/invokebrns.f -o invokebrns.o

# Clean build artifacts
clean:
	rm -f *.o brns.so

.PHONY: all clean
