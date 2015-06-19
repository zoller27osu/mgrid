c*
c* Translated from C by Michael Zoller, 2015.
c* Original copyright information follows:
c***********************************************************************
c* Copyright (C) 2002-2012 the Network-Based Computing Laboratory
c* (NBCL), The Ohio State University.
c*
c* Contact: Dr. D. K. Panda (panda@cse.ohio-state.edu)
c*
c* For detailed copyright and licensing information, please refer to the
c* copyright file COPYRIGHT in the top level OMB directory.
c*/

!#define mem_realignment

#ifdef mem_realignment
#  define MESSAGE_ALIGNMENT (64)
#else
#  define MESSAGE_ALIGNMENT (0)
#endif
!#define MAX_MSG_SIZE 100000
#define MAX_MSG_SIZE (2**20)
#define MYBUFSIZE (MAX_MSG_SIZE + MESSAGE_ALIGNMENT)
#define SKIP_LARGE  10
#define LOOP_LARGE  100
#define SKIP_SMALL  1000
#define LOOP_SMALL  10000
#define LARGE_MSG_SIZE  (8192)

#ifdef mem_realignment
#  define MAX_ALIGNMENT (16384)
#else
#  define MAX_ALIGNMENT (0)
#endif
!#define MAX_MSG_SIZE (1<<7)
#define MAX_STEPS    (22+1)
#define MAXBUFSIZE (MAX_MSG_SIZE + MAX_ALIGNMENT)
#define MULTI_SKIP_LARGE  10
#define MULTI_LOOP_LARGE  1000
#define MULTI_SKIP_SMALL  100
#define MULTI_LOOP_SMALL  10000

!#ifndef FIELD_WIDTH
!#  define FIELD_WIDTH 20
!#endif

!#ifndef FLOAT_PRECISION
!#  define FLOAT_PRECISION 2
!#endif

      program ping_pong
      IMPLICIT NONE
      include 'mpif.h'

      integer nid, np, wdsize, ierr
      real :: eps, oeps

      call MPI_Init(ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, np, ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, nid, ierr)

      if (MOD(np,2).ne.0) then
        if (nid.eq.0) then
          !fprintf(stderr, "This test requires exactly two processes\n")
          write(6,*) "This test requires an even number of processes!"
        endif

        !call MPI_Finalize(ierr)

        stop
      endif

      wdsize = 4
      eps = 1.0e-12
      oeps = eps + 1
      if (oeps.ne.(1.0)) then
        wdsize = 8
      endif

      if (nid.eq.0) write(6,*) nid,": wordsize in bytes =", wdsize

      if (nid.eq.0) write(6,*) nid,": MYBUFSIZE, MAXBUFSIZE =", 
     &      MYBUFSIZE, MAXBUFSIZE, __LINE__
      
      call single_lat(nid, np)
      if (nid.eq.0) write(6,*) " "
      call multi_lat(nid, np)

      call MPI_Finalize(ierr)
      return !EXIT_SUCCESS
      end program

c-----------------------------------------------------------------------
      subroutine single_lat(nid, np)
      IMPLICIT NONE
      include 'mpif.h'
      integer :: nid, np
      
      integer :: i, size, skip, loop, partner, ierr
#ifdef mem_realignment
      integer :: align_size
#endif
      integer, dimension(MPI_STATUS_SIZE) :: reqstat
      logical :: send_first
      !char *s_buf, *r_buf
      character*1, dimension(MYBUFSIZE) :: s_buf, r_buf
      double precision :: t_start, t_end, latency, total_lat, avg_lat

      t_start = 0.0
      t_end = 0.0
      
#ifdef mem_realignment
      align_size = MESSAGE_ALIGNMENT
c**************Allocating Memory*********************
      ! TODO: translate to Fortran if possible, else remove
      s_buf =
        (char *) (((unsigned long) s_buf + (align_size - 1)) /
     &    align_size * align_size)

      r_buf =
        (char *) (((unsigned long) r_buf + (align_size - 1)) /
     &    align_size * align_size)
c**************Memory Allocation Done*********************
#endif

      if (nid.eq.0) then
        !write(6,*) HEADER ! where does HEADER come from?
        !fprintf(stdout, "%-*s%*s\n", 10, "# Size", FIELD_WIDTH, "Latency (us)")
        write(6,*) "Size (bytes)        Latency (us)"
        !call fflush(stdout)
      endif
      
      loop = LOOP_SMALL
      skip = SKIP_SMALL
      partner = nid + 1 - 2*MOD(nid,2) ! nid+1 for even, nid-1 for odd
      send_first = (MOD(nid,2).eq.0)
      
      size = 0
      do !size = 0,MAX_MSG_SIZE !size = (size ? size * 2 : 1)) {
        !* touch the data */
        do i = 1,size
          s_buf(i) = 'a'
          r_buf(i) = 'b'
        enddo

        if (size.gt.LARGE_MSG_SIZE) then
          loop = LOOP_LARGE
          skip = SKIP_LARGE
        endif
        !if (nid.eq.0) write(6,*) nid, ": size =", size, __LINE__
        do i = 0,loop + skip - 1
          if (i.eq.skip) then
            !write(6,*) nid, ": ", __LINE__
            call MPI_Barrier(MPI_COMM_WORLD, ierr)
            t_start = MPI_Wtime()
          endif
          !if (i.eq.0) write(6,*) nid, ": i =", i, __LINE__
          if (send_first) then
            call MPI_Send(s_buf, size, MPI_BYTE, partner, 1,
     &                    MPI_COMM_WORLD, ierr)
            call MPI_Recv(r_buf, size, MPI_BYTE, partner, 1,
     &                    MPI_COMM_WORLD, reqstat, ierr)
          else
            call MPI_Recv(r_buf, size, MPI_BYTE, partner, 1,
     &                    MPI_COMM_WORLD, reqstat, ierr)
            call MPI_Send(s_buf, size, MPI_BYTE, partner, 1,
     &                    MPI_COMM_WORLD, ierr)
          endif
        enddo
        t_end = MPI_Wtime()
        
        latency = (t_end - t_start) * 1e6 / (2.0 * loop)
        
        call MPI_Reduce(latency, total_lat, 1, MPI_DOUBLE, MPI_SUM, 0, 
     &                  MPI_COMM_WORLD, ierr)

        avg_lat = total_lat/np !(pairs * 2.0)

        if (nid.eq.0) then
          !fprintf(stdout,"%-*d%*.*f\n",10,size,FIELD_WIDTH,FLOAT_PRECISION,latency)
          write(6,*) size, avg_lat
          !call fflush(stdout)
        endif
        
        if (size.eq.0) then
          size = 1
        else
          size = size * 2
        endif
        
        if (size.gt.MAX_MSG_SIZE) return
      enddo
      
      !call MPI_Barrier(MPI_COMM_WORLD, ierr)
      
      return
      end subroutine
c-----------------------------------------------------------------------
      subroutine multi_lat(nid, np)
      IMPLICIT NONE
      include 'mpif.h'
      integer :: nid, np

      integer :: i, size, skip, loop, pairs, partner, ierr
#ifdef mem_realignment
      integer :: align_size
#endif
      integer, dimension(MPI_STATUS_SIZE) :: reqstat
      logical :: send_first
      character*1, dimension(MYBUFSIZE) :: s_buf, r_buf
      double precision :: t_start, t_end, latency, total_lat, avg_lat
      t_start = 0.0
      t_end = 0.0
      latency = 0.0
      total_lat = 0.0
      avg_lat = 0.0

#ifdef mem_realignment
      ! TODO: translate to Fortran if possible, else remove
      align_size = getpagesize()
      s_buf =
        (char *) (((unsigned long) s_buf + (align_size - 1)) /
                  align_size * align_size)
      r_buf =
        (char *) (((unsigned long) r_buf + (align_size - 1)) /
                  align_size * align_size)
#endif

      !memset(s_buf, 0, MAX_MSG_SIZE)
      !memset(r_buf, 0, MAX_MSG_SIZE)
      do i = 1,MAX_MSG_SIZE
        s_buf(i) = '\0'
        r_buf(i) = '\0'
      enddo

      pairs = np/2
      if (0.eq.0) then
        partner = MOD(nid + pairs, np)
      else
        if (nid.lt.pairs) then
          partner = nid + pairs
        else
          partner = nid - pairs
        endif
      endif
      write(6,*) nid, ": partner =", partner, __LINE__
      send_first = (nid.lt.pairs)
      loop = MULTI_LOOP_SMALL
      skip = MULTI_SKIP_SMALL

      if (nid.eq.0) then
        !fprintf(stdout, HEADER)
        !fprintf(stdout, "%-*s%*s\n", 10, "# Size", FIELD_WIDTH, "Latency (us)")
        write(6,*) "Size (bytes)        Latency (us)"
  !*hjw: */
        !fprintf(stdout, "# [ pairs: %d ] \n", pairs)
        !write(6,*) "# [ pairs:", pairs, "] "
        !fflush(stdout)
      endif
      
      size = 0
      do !size = 0,MAX_MSG_SIZE !size = (size ? size * 2 : 1)) {
        !call MPI_Barrier(MPI_COMM_WORLD, ierr)

        if (size.gt.LARGE_MSG_SIZE) then
          loop = MULTI_LOOP_LARGE
          skip = MULTI_SKIP_LARGE
        endif

        do i = 0,loop + skip - 1
          if (i.eq.skip) then
            call MPI_Barrier(MPI_COMM_WORLD, ierr)
            t_start = MPI_Wtime()
          endif
     
          if (send_first) then
            call MPI_Send(s_buf, size, MPI_BYTE, partner, 1,
     &                    MPI_COMM_WORLD, ierr)
            call MPI_Recv(r_buf, size, MPI_BYTE, partner, 1,
     &                    MPI_COMM_WORLD, reqstat, ierr)
          else
            call MPI_Recv(r_buf, size, MPI_BYTE, partner, 1,
     &                    MPI_COMM_WORLD, reqstat, ierr)
            call MPI_Send(s_buf, size, MPI_BYTE, partner, 1,
     &                    MPI_COMM_WORLD, ierr)
          endif
        enddo
        t_end = MPI_Wtime()

        latency = (t_end - t_start) * 1.0e6 / (2.0 * loop)

        call MPI_Reduce(latency, total_lat, 1, MPI_DOUBLE, MPI_SUM, 0, 
     &                  MPI_COMM_WORLD, ierr)

        avg_lat = total_lat/np !(pairs * 2.0)

        if (0.eq.nid) then
          !fprintf(stdout, "%-*d%*.*f\n", 10, size, FIELD_WIDTH, FLOAT_PRECISION, avg_lat)
          write(6,*) size, avg_lat
          !fflush(stdout)
        endif
        
        if (size.eq.0) then
          size = 1
        else
          size = size * 2
        endif
        
        if (size.gt.MAX_MSG_SIZE) exit
      enddo
      
      !call MPI_Barrier(MPI_COMM_WORLD, ierr)

      return !EXIT_SUCCESS
      end subroutine

c* vi: set sw=4 sts=4 tw=80: */
