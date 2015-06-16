# Variables implicitly used by GNU Make to autocreate targets.
# Conditionally use different parameters if running on cray compiler.
# Flags: DEBUG_OUT, NUMERICAL_OUT
ifeq ($(PE_ENV),CRAY)
        FC = ftn
        FFLAGS = -O3 -s real64
        LDFLAGS = -O3
else
	ifeq ($(USER),oychang)
		FC = mpif77.mpich2
	else
		FC = mpif77
	endif

        FFLAGS = -O3 -mcmodel=medium -fdefault-real-8 -fdefault-double-8 \
		-DDEBUG_OUT
        LDFLAGS = -O3 -mcmodel=medium
endif

OBJS = comm_mpi.o x2p.o

##############################################################################

# Explicitly spell this one out otherwise uses C linker
x2p: $(OBJS)
	$(FC) $(LDFLAGS) -o $@ $^
x2p.o: x2p.F
comm_mpi.o: comm_mpi.F MGRID

##############################################################################

.PHONY: clean deploy run
run: x2p
	mpiexec ./x2p
clean:
	$(RM) $(OBJS) x2p x.x *.mod
deploy: x2p
	cp x2p $(HOME)/scratch/
	cp in.dat $(HOME)/scratch/
