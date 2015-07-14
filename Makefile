# Variables implicitly used by GNU Make to autocreate targets.
# Conditionally use different parameters if running on cray compiler.
HOST=$(shell hostname)

@FC = none
FFLAGS = -D alt_timing -D use_tags
ifeq ($(PE_ENV),CRAY)
    #@echo "Using CRAY."
    FC = ftn
    #-M 1058
    FLAGS = -M 124 -O3 -s default64 -hnocaf -hnopgas_runtime -hmpi1 \
        -hpic #first half of mcmodel=medium equivalent
        #-hvector3 -hscalar3 \
        #-hnegmsgs \
        #-fbacktrace \
        #-hdevelop -eD \
        #-Wall -Og #-eI
    FFLAGS += -e chmnF -dX -r d -J bin -Q bin
        #-hkeepfiles #-S
        #-dX: 10,000-variable-module initialize-before-main thing
    LDFLAGS = -dynamic #second half of mcmodel=medium equivalent
else ifeq ($(PE_ENV),PGI)
    #@echo "Using PGI."
    FC = ftn
    # TODO: check the -h flags as well as -fbacktrace
    #-M 1058,124
    #-default64: instructs compiler wrappers to include 64-bit MPI library
    FLAGS = -O3 -r8 -mcmodel=medium -Mdalign \
        -Mllalign -Munroll -Kieee -fastsse -Mipa=fast \
        #-hnocaf -hnopgas_runtime -hmpi1
        #-hvector3 -hscalar3 \
        #-hnegmsgs \
        #-fbacktrace \
        #-hdevelop -eD \
        #-Wall -Og #-eI
    FFLAGS += #-e chmnF -dX -r d -J bin -Q bin
        #-hkeepfiles #-S
        #-dX: 10,000-variable-module initialize-before-main thing
    LDFLAGS =
else # normal MPI, as on workstations
    #@echo "Using Workstation compiler."
    ifeq ($(USER),oychang)
        FC = mpif77.mpich2
    else
        FC = mpif77
    endif
    FLAGS = -O3 -mcmodel=medium -fdefault-real-8 -fdefault-double-8 #-DDEBUG_OUT
        #-fbacktrace #-Wall -Og
    FFLAGS += #
    #FCOMP = $(shell $(FC) -show)
    ifeq ($(word 1, $(shell $(FC) -show)),gfortran)
        FFLAGS += -D gfortran
    endif
#endif
    LDFLAGS =
endif

##############################################################################

X2P_OBJS = comm_mpi.o x2p.o
PP_OBJS = ping_pong.o
OBJS = $(X2P_OBJS) $(PP_OBJS)

EXECS = x2p intraer ping_pong

all: $(EXECS)

# Explicitly spell this one out otherwise uses C linker
x2p: $(X2P_OBJS)
	$(FC) $(LDFLAGS) -o $@ $^
comm_mpi.o: MGRID comm_mpi.F
x2p.o: comm_mpi.o x2p.F

#x2pTwisty: comm_mpi.o x2pTwisty.o
#	$(FC) $(LDFLAGS) -o $@ $^
#x2pTwisty.o: x2pTwisty.F comm_mpi.o

intraer: intraer.F
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $^

ping_pong: ping_pong.F
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $^

##############################################################################

.PHONY: clean deploy debug runx2p cray pgi

cray:
	-module swap PrgEnv-pgi PrgEnv-cray
	make clean
	@echo ""
	@echo "Swap complete. Run make as normal from now on."

pgi:
	-module swap PrgEnv-cray PrgEnv-pgi
	make clean
	@echo ""
	@echo "Swap complete. Run make as normal from now on."

clean:
	rm -rf $(OBJS) $(EXECS) x.x *.mod *.out bin/*

deploy: x2p ping_pong
ifneq ($(PE_ENV),)
	cp in.dat $(HOME)/scratch/
	cp x2p $(HOME)/scratch/
	cp x2p.pbs $(HOME)/scratch/
	cp x2p_graph.sh $(HOME)/scratch/
	cp ping_pong $(HOME)/scratch/
	cp ping_pong.pbs $(HOME)/scratch/
	cp bw_ping_pong.gp $(HOME)/scratch/
	cp pgo.sh ~/nek5_svn/examples/timing/.
endif

debug: deploy
ifneq ($(PE_ENV),)
    ifeq (,$(findstring nid,$(HOST))) # if not in a node
	#interactive qsub (aka CCM) is not for computing!!!
	# it is, however, appropriate for small debug jobs.
    	qsub -I -l gres=ccm -l nodes=4:ppn=16:xk -l walltime=01:00:00
    else
	#aprun -n 64 $(HOME)/scratch/./x2p # 4*16 = 64
    endif
else
	mpiexec -n 8 ./x2p # you're on your own
endif

runx2p: deploy
ifneq ($(PE_ENV),)
	# Assume we are computing, thus use qsub (proper even in CCM).
	./watch_x2p.sh $(shell qsub $(HOME)/scratch/x2p.pbs) # qsub should output the PBS_JOBID
    ifneq (,$(findstring nid,$(HOST))) # if in a node (i.e. CCM)
	# do special graphing stuff
    endif
else
	mpiexec -n 8 ./x2p
endif

runPing: deploy
ifneq ($(PE_ENV),)
	# Assume we are computing.
	qsub $(HOME)/scratch/ping_pong.pbs # proper for computing, even in CCM!
	#gnuplot bw_ping_pong.gp
else
	mpiexec -n 8 ./ping_pong > wks_ping_pong.out
	# printf '%s\n\n' "$(tail -n +6 wks_ping_pong.out)" > wks_ping.out
	# sed '1,5d' wks_ping_pong.out | sort -k1,1n -k3n > wks_ping.out
	#grep partner wks_ping_pong.out | sort -k1,1n -k3n > wks_ping.out
	tail -n +6 wks_ping_pong.out | sort -k1,1n -k3n > wks_ping.out
	gnuplot wks_ping_pong.gp
endif
