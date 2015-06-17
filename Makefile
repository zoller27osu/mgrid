# Variables implicitly used by GNU Make to autocreate targets.
# Conditionally use different parameters if running on cray compiler.

HOST=$(shell hostname)

ifeq ($(PE_ENV),CRAY)
        FC = ftn
	FLAGS = -O3 -s real64 -M 124,1058 -hnocaf -hnopgas_runtime -hmpi1 \
		-hpic #first half of mcmodel=medium equivalent
		#-hvector3 -hscalar3 \
		#-hnegmsgs \
		#-fbacktrace \
		#-hdevelop -eD \
		#-Wall -Og #-eI
        FFLAGS = $(FLAGS) -e chmnF -dX -r d -J bin -Q bin #-hkeepfiles #-S
	#-dX: 10,000-variable-module initialize-before-main thing
        LDFLAGS = -dynamic #second half of mcmodel=medium equivalent
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

all: x2p ping_pong

# Explicitly spell this one out otherwise uses C linker
x2p: $(OBJS)
	$(FC) $(LDFLAGS) -o $@ $^
comm_mpi.o: comm_mpi.F MGRID
x2p.o: x2p.F comm_mpi.o

ping_pong: ping_pong.o
	$(FC) $(LDFLAGS) -o $@ $^
ping_pong.o: ping_pong.F

##############################################################################

.PHONY: clean deploy run
run: x2p
ifeq ($(PE_ENV),CRAY)
    ifeq (,$(findstring nid,$(HOST)))
	qsub -I -l gres=ccm -l nodes=4:ppn=16:xk -l walltime=01:00:00
    else
	#cd $(HOME)/scratch
	aprun -n 64 $(HOME)/scratch/./x2p
    endif
else
	mpiexec -n 8 ./x2p
endif

clean:
	rm -rf $(OBJS) x2p x.x *.mod bin/*
deploy: x2p
ifeq ($(PE_ENV),CRAY)
	cp x2p $(HOME)/scratch/
	cp in.dat $(HOME)/scratch/
endif
