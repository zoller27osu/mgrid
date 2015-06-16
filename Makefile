# Variables implicitly used by GNU Make to autocreate targets.
# Conditionally use different parameters if running on cray compiler.

#HOST=$(shell hostname)

ifeq ($(PE_ENV),CRAY)
        FC = ftn
	FLAGS = -O3 -s real64 -M 124,1058 -hnocaf -hnopgas_runtime -hmpi1 \
		#-hvector3 #-hscalar3
		#-hnegmsgs #-fbacktrace #-hdevelop -eD #-Wall -Og #-eI
        FFLAGS = $(FLAGS) -e chmnF -dX -r d -J bin -Q bin #-hkeepfiles #-S
        LDFLAGS = $(FLAGS) #-dynamic -hpic
else # normal MPI, as on workstations
	ifeq ($(USER),oychang)
		FC = mpif77.mpich2
	else
		FC = mpif77
	endif
        FFLAGS = -O3 -mcmodel=medium -fdefault-real-8 -fdefault-double-8 \
		-DDEBUG_OUT #-fbacktrace #-Wall -Og
        LDFLAGS = -O3 -mcmodel=medium
endif

OBJS = comm_mpi.o x2p.o

##############################################################################

# Explicitly spell this one out otherwise uses C linker
x2p: $(OBJS)
	$(FC) $(LDFLAGS) -o $@ $^
comm_mpi.o: comm_mpi.F MGRID
x2p.o: x2p.F comm_mpi.o

##############################################################################

.PHONY: clean deploy run
run: x2p
	mpiexec -n 8 ./x2p
clean:
	rm -f $(OBJS) x2p x.x *.mod bin/*
deploy: x2p
	cp x2p $(HOME)/scratch/
