S       = .

CC      = mpi
CFLAGS  = -std=c99

FC      = mpif77
FCFLAGS = -fdefault-real-8 -fdefault-double-8 -x f77-cpp-input \
          -g -O0 -mcmodel=medium #-fbacktrace #-Wall -Og

NOBJS   = comm_mpi.o x2p.o

##############################################################################

x2p: $(NOBJS)
	$(F77) $(FCFLAGS) -o x2p $(NOBJS)
x2p.o: x2p.f MGRID
	$(F77) $(FCFLAGS) -c x2p.f
comm_mpi.o: comm_mpi.f
	$(F77) $(FCFLAGS) -c comm_mpi.f

##############################################################################

.PHONY: clean
clean:
	rm -f $(NOBJS) x2p
