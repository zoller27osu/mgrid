Compiling (Instrumented) Tau Code
=================================

If things don't compile, make sure the `/sw/xe/tau/<version_num>` is correct with regards to the version in `module show tau`.
In addition, the `module swap` is to do the same sort of version correctness.

1. module load tau
2. export TAU_MAKEFILE=/sw/xe/tau/2.21.4/cnl4.1_cray8.1.1/craycnl/lib/Makefile.tau-cray-mpi
3. module swap darshan darshan/2.3.0
4. Switch `ftn` to `tau_f90.sh -optCompInst` in Makefile
5. `make x2p`


File Summaries
==============

```bash
comm_mpi.F # Wrappers for MPI functions
in.dat # Input data for run(s)
Makefile
MGRID # submodule-ish grouping of multigrid variables
x2p.F # multigrid algorithm
ping_pong.F # Performs an intra-node, then inter-node ping-pong test
```

x2p:
Input
=====

in.dat format
-------------

col# | var     | type            | description
-----|---------|-----------------|--------------------------------------------
col1 | nlev    | 4-byte int      | number of multigrid levels
col2 | nsmooth | 4-byte int      | number of smoothings (unused?)
col3 | igs     | 4-byte int      | {0,1} where 0 is jacobian, 1 is gaussian
col4 | sigma   | word size float | ...


Variables
=========

var | description
----|-------------------------------------------------------------------------
ue  | exact solution
u   | multigrid method solution


Output
======

if line ends in "max error[A-Z]"
--------------------------------

col#  | var  | description
------|------|----------------------------------------------------------------
col1  | ii   | number of iterations/passes (over by-one)
col2  | dmax | some maximum in prolongation stages; heuristically? always 0
col3  | emx  | "error", the max difference u[i] - ue[i] between all ranks
col4  | time | delta = mpi_wtime()_end - mpi_wtime()_start
col5  | ldim | number of dimensions; always equal to 3 (`MGRID` parameter)
col6  | nlev | number of multigrid levels
col7  | nx   | global x-dimension of current multigrid level
col8  | ny   | global y-dimension of current multigrid level
col9  | nz   | global z-dimension of current multigrid level
col10 | n3   | number of cells, i.e. nx*ny[*nz]
col11 | mp   | local x-dim; np = mp*mq*mr where np is a multiple of 8 (2*2*2)
col12 | mq   | local y-dim of current multigrid level
col13 | mr   | local z-dim of current multigrid level
col14 | np   | number of processors; result of MPI_COMM_SIZE(_, np, _)
col15 | c1   | trial label; 'A' => ..., 'B' => ..., 'X' => ...

if line ends in "flops[A-Z]"
----------------------------

col#  | var    | description
------|--------|--------------------------------------------------------------
col1  | ii     | number of iterations/passes
col2  | emx    | "error", the max difference u[i] - ue[i] between all ranks
col3  | time   | delta = mpi_wtime()_end - mpi_wtime()_start
col4  | fflops | FP ops per second; (cflops + pflops + rflops + rrlops)/time
col5  | cflops | sum of FLOP in smoothing steps
col6  | pflops | sum of FLOP in 3D prolongation steps
col7  | rflops | sum of FLOP in 3D restrict steps
col8  | rrlops | sum of FLOP in `compute_r`
col9  | n3     | number of cells, i.e. nx*ny[*nz]
col10 | np     | number of processors; result of MPI_COMM_SIZE(_, np, _)
col11 | c1     | trial label; 'A' => ..., 'B' => ..., 'X' => ...
