# This is the makefile for DELTAC on DarwinIntel computers.

FFLAGS = -O2 

F90 = gfortran $(FFLAGS)

INTELLIBS= \
	-L/usr/lib \
	-llapack \
	-L/usr/lib \
	-lblas \
	-lpthread \
	-lm

.f.o:
	$(F90) -c $*.f

#Note: I don't have DarwinIntel, which is required for lsode
LIBDIR = \
	-L../lib/DarwinIntel \

LIBS = \
	$(INTELLIBS) \
	$(LIBDIR) #\
	#-llsode

OBJS = \
	io.o \
	local.o \
	gamma.o \
	jacobi.o \
        lsode.o \
	deltar.o \
	deltac.o \
	driver.o

deltac: $(OBJS)
	$(F90) -o deltac $(OBJS) $(LIBS)

# dependencies

local.o: io.o
jacobi.o: local.o
gamma.o: local.o
deltar.o: gamma.o lsode.o
deltac.o: deltar.o jacobi.o lsode.o 
driver.o: deltar.o deltac.o

clean:
	rm -f *.o *.mod *.out *.bin *~ deltac
