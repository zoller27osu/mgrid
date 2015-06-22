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

X2P_OBJS = comm_mpi.o x2p.o
PP_OBJS = ping_pong.o
OBJS = $(X2P_OBJS) $(PP_OBJS)

EXECS = x2p intraer ping_pong

##############################################################################

all: $(EXECS)

# Explicitly spell this one out otherwise uses C linker
x2p: $(X2P_OBJS)
	$(FC) $(LDFLAGS) -o $@ $^
comm_mpi.o: comm_mpi.F MGRID
x2p.o: x2p.F comm_mpi.o

intraer: intraer.F
	$(FC) $(LDFLAGS) -o $@ $^

ping_pong: ping_pong.F
	$(FC) $(LDFLAGS) -o $@ $^
ping_pong.o: ping_pong.F

##############################################################################

.PHONY: clean deploy debug runx2p

clean:
	rm -rf $(OBJS) $(EXECS) x.x *.mod *.out bin/*

deploy: x2p ping_pong
ifeq ($(PE_ENV),CRAY)
	cp in.dat $(HOME)/scratch/
	cp x2p $(HOME)/scratch/
	cp x2p.pbs $(HOME)/scratch/
	cp ping_pong $(HOME)/scratch/
	cp ping_pong.pbs $(HOME)/scratch/
endif

debug: deploy
ifeq ($(PE_ENV),CRAY)
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
ifeq ($(PE_ENV),CRAY)
	# Assume we are computing.	
	qsub $(HOME)/scratch/x2p.pbs # proper for computing, even in CCM!
    ifneq (,$(findstring nid,$(HOST))) # if in a node (i.e. CCM)
	# do special graphing stuff
    endif
else
	mpiexec -n 8 ./x2p
endif

runPing: deploy
ifeq ($(PE_ENV),CRAY)
	# Assume we are computing.
	qsub $(HOME)/scratch/ping_pong.pbs # proper for computing, even in CCM!
	#gnuplot bw_ping_pong.gp
else
	mpiexec -n 8 ./ping_pong > wks_ping_pong.out
	# grep partner pingout.$PBS_JOBID > ping.dat
	# printf '%s\n\n' "$(tail -n +6 wks_ping_pong.out)" > wks_ping.out
	#tail -n +6 wks_ping_pong.out > wks_ping.out
	sed '1,5d' wks_ping_pong.out | sort -k1n,1 -k3n > wks_ping.out
	gnuplot wks_ping_pong.gp
endif
