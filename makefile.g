S     = .
L     = 
X     = 
CC    = pgcc
F77   = pgf77
P     = -r8

FLAGS = -mcmodel=medium -O0 -g

NOBJS = x2p.o comm_mpi.o mpi_dummy.o

all: x2p


x2p:	$(NOBJS)
	$(F77) -o x2p $(NOBJS) $(FLAGS)

clean:
	'rm' *.o
	'rm' x2p


x2p.o   	: x2p.f 	; $(F77) -c $P $(FLAGS)  x2p.f 
comm_mpi.o   	: comm_mpi.f 	; $(F77) -c $P $(FLAGS)  comm_mpi.f 
mpi_dummy.o   	: mpi_dummy.f 	; $(F77) -c $P $(FLAGS)  mpi_dummy.f 
