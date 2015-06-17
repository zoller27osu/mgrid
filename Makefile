# Variables implicitly used by GNU Make to autocreate targets.
# Conditionally use different parameters if running on cray compiler.
ifeq ($(PE_ENV),CRAY)
        FC = ftn
	FLAGS = -O3 -s real64 -M 124,1058 -hnocaf -hnopgas_runtime -hmpi1 \
		-hpic 
		#-hvector3 -hscalar3 \
		#-hnegmsgs \
		#-fbacktrace \
		#-hdevelop -eD \
		#-Wall -Og #-eI
        FFLAGS = $(FLAGS) -e chmnF -dX -r d -DDEBUG_OUT#-hkeepfiles #-S
	#-dX: 10,000-variable-module initialize-before-main thing
        LDFLAGS = $(FLAGS) -dynamic
else
	FC = mpif77
        FFLAGS = -O3 -mcmodel=medium -fdefault-real-8 -fdefault-double-8
        LDFLAGS = -O3 -mcmodel=medium
endif

HOST=$(shell hostname)
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
	$(RM) $(OBJS) x2p x.x *.mod *.cg *.opt
deploy: x2p
	cp x2p $(HOME)/scratch/
	cp in.dat $(HOME)/scratch/
