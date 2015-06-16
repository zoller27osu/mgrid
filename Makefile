# Variables implicitly used by GNU Make to autocreate targets.
# Conditionally use different parameters if running on cray compiler.
ifeq ($(PE_ENV),CRAY)
        FC = ftn
        FFLAGS = -O3 -s real64
        LDFLAGS = -O3
endif
ifeq ($(USER),oychang)
        FC = mpif77.mpich2
        FFLAGS = -O3 -mcmodel=medium -fdefault-real-8 -fdefault-double-8 \
        	-x f77-cpp-input
        LDFLAGS = -O3 -mcmodel=medium
endif
ifeq ($(USER),ochang3)
        FC = mpif77
        FFLAGS = -O3 -mcmodel=medium -fdefault-real-8 -fdefault-double-8 \
        	-x f77-cpp-input
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
