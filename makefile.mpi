S     = .
L     = 
X     = 
CC    = mpicc
F77   = mpif77
P     = -fdefault-real-8 -fdefault-double-8 -x f77-cpp-input

FLAGS = -g -O3 -mcmodel=medium #-fbacktrace #-Wall -Og

CFLAGS   = -std=c99
FFLAGS   = #-c

NOBJS = comm_mpi.o x2p.o 
NOBJSC = comm_mpi.ox2pc.o

all: f

f:	$(NOBJS)
	$(F77) $(FFLAGS) -o x2p $(NOBJS) $(FLAGS)

c:	$(NOBJSC)
	$(CC) $(CFLAGS) -o x2pc $(NOBJSC) $(FLAGS)

#x2p:	$(NOBJS)
#	$(F77) -o x2p $(NOBJS) $(FLAGS)

#x2pc:	$(NOBJS)
#	$(CC) -o x2pc $(NOBJSC) $(FLAGS)

clean:
	-'rm' $(NOBJS) $(NOBJSC) x2p x2pc

x2p.o   	: x2p.F 	; $(F77) -c $P $(FLAGS)  x2p.F 
comm_mpi.o   	: comm_mpi.F 	; $(F77) -c $P $(FLAGS)  comm_mpi.F 
x2pc.o   	: x2p.c 	; $(CC) -c $P $(FLAGS)  x2p.c 
#comm_mpi.o   	: comm_mpi.c 	; $(CC) -c $P $(FLAGS)  comm_mpi.c 
