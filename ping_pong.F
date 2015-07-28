!#define mem_realignment
!#define statuses

!probably doesn't matter, but has to match
#define TAG 1

#define MASTER_RANK 0

!how many different processes we are aiming to ping
#define TARGET_NP 2000.0

!#define ALL_SIZES

! 1 hr
!#define BASE 1.41421356237
! 2 hr
!#define BASE 1.189207115
! 4 hr
!#define BASE 1.09050773267
! 8 hr
!#define BASE 1.04427378243
! 16 hr
!#define BASE 1.02189714865
! "32 hr"
#define BASE 1.01088928605

!#define MAX_MSG_SIZE 100000
#define MAX_MSG_SIZE (2**20)
!#define MAX_MSG_SIZE (800000)
#define LARGE_MSG_SIZE  (8192)

#ifdef mem_realignment
#  define MAX_ALIGNMENT (16384)
#else
#  define MAX_ALIGNMENT (0)
#endif
!#define MAX_MSG_SIZE (1<<7)
#define MAX_STEPS    (22+1)
#define MAXBUFSIZE (MAX_MSG_SIZE + MAX_ALIGNMENT)
#define SKIP_LARGE  10
!#define LOOP_LARGE  100
#define LOOP_LARGE  1000
!#define SKIP_SMALL  1000
#define SKIP_SMALL  100
#define LOOP_SMALL  10000

!#ifndef FIELD_WIDTH
!#  define FIELD_WIDTH 20
!#endif

!#ifndef FLOAT_PRECISION
!#  define FLOAT_PRECISION 2
!#endif

      program ping_pong
      IMPLICIT NONE
      include 'mpif.h'

      integer nid, np, wdsize, ierr, step
      real :: eps, one_eps

      call MPI_Init(ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, np, ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, nid, ierr)

c      if (MOD(np,2).ne.0) then
c        if (nid.eq.0) then
c          write(6,*) "This test requires an even number of processes!"
c        endif

        !call MPI_Finalize(ierr)

c        stop
c      endif

c      wdsize = 4
c      eps = 1.0e-12
c      one_eps = eps + 1
c      if (one_eps.ne.(1.0)) then
c        wdsize = 8
c      endif
      
#ifdef SIZEOF
      wdsize = SIZEOF(eps)
#else
      !if (nid.eq.MASTER_RANK) write(6,*) nid, "SIZEOF not detected!"
      wdsize = 4
      if ((KIND(eps).eq.2).OR.(KIND(eps).eq.8)) then
        wdsize = 8
      endif
#endif

      if (nid.eq.MASTER_RANK) then
        write(6,*) nid,": wordsize in bytes =", wdsize
        write(6,*) nid,": MAXBUFSIZE =", MAXBUFSIZE, __LINE__
      endif
      
      step = NINT(np/TARGET_NP)
      if (nid.eq.MASTER_RANK) then
        write(6,*) nid,": calculated step, np=", step, np, __LINE__
      endif
      step = MAX(1, step)
c      if (nid.eq.MASTER_RANK) then
c        write(6,*) nid,": step, np=", step, np, __LINE__
c      endif
      
      call single_sender(nid, np, MASTER_RANK, step, real(wdsize))
      !if (nid.eq.0) write(6,*) " "
      !call multi_lat(nid, np)

      call MPI_Finalize(ierr)
      end program

c-----------------------------------------------------------------------
      subroutine single_sender(nid, np, sender, step, wdsize)
      IMPLICIT NONE
      include 'mpif.h'
      integer :: nid, np, step, sender
      real :: wdsize
      
      integer :: i, size, skip, loop, partner, ierr
#ifdef mem_realignment
      integer :: align_size
#endif
#ifdef statuses
      integer, dimension(MPI_STATUS_SIZE) :: reqstat
#endif
      real :: rsize
      character*1, dimension(MAXBUFSIZE) :: s_buf, r_buf
      double precision :: t_start, t_end, latency
      character(len=8) :: i_char

      t_start = 0.0
      t_end = 0.0
      
      ! Only nodes that sender will talk to need continue
      if (MOD(nid-sender, step).ne.0) then
        write(6,*) nid, ": returning! sender, step =", sender, step
        return
      endif
      
#ifdef mem_realignment
      align_size = MESSAGE_ALIGNMENT
c**************Allocating Memory*********************
      s_buf =
        (char *) (((unsigned long) s_buf + (align_size - 1)) /
     &    align_size * align_size)

      r_buf =
        (char *) (((unsigned long) r_buf + (align_size - 1)) /
     &    align_size * align_size)
c**************Memory Allocation Done*********************
#endif

      if (nid.eq.sender) then
        !write(6,*) HEADER ! where does HEADER come from?
        !fprintf(stdout, "%-*s%*s\n", 10, "# Size", FIELD_WIDTH, "Latency (us)")
        write(6,*) "Node with id", nid, "is sender."
        write(6,*) "Size (bytes)        Latency (us)"
        write(i_char, '(i8)') np
        write(6,*) "0 P=",adjustl(i_char)," time(us) P=",adjustl(i_char)
        !call fflush(stdout)
      endif
      
      loop = LOOP_SMALL
      skip = SKIP_SMALL
#ifdef ALL_SIZES
      do size = 0,MAX_MSG_SIZE !size = (size ? size * 2 : 1)) {
#else
      rsize = 0
      size = 0
      do
#endif
        
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
        
        if (nid.eq.sender) then
        
          do partner = MOD(nid,step),np-1,step
            if (partner.eq.nid) cycle
            !write(6,*) nid, ": partner =", partner, __LINE__
            do i = 0,loop + skip - 1
              if (i.eq.skip) then
                !write(6,*) nid, ": ", __LINE__
                !call MPI_Barrier(MPI_COMM_WORLD, ierr)
                t_start = MPI_Wtime()
              endif
              !if (i.eq.0) write(6,*) nid, ": i =", i, __LINE__
              call MPI_Send(s_buf, size, MPI_BYTE, partner, TAG,
     &                      MPI_COMM_WORLD, ierr)
#ifdef statuses
              call MPI_Recv(r_buf, size, MPI_BYTE, partner, TAG,
     &                      MPI_COMM_WORLD, reqstat, ierr)
#else
              call MPI_Recv(r_buf, size, MPI_BYTE, partner, TAG,
     &                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
#endif
            enddo
            t_end = MPI_Wtime()
            
            latency = (t_end - t_start) * 1e6 / (2.0 * loop)
            
            !fprintf(stdout,"%-*d%*.*f\n",10,size,FIELD_WIDTH,FLOAT_PRECISION,latency)
            write(6,*) partner, "partner", size/wdsize, latency
            !call fflush(stdout)
          enddo
        
        else
          ! Not the sender; wait until we hear from them and then respond, (loop + skip) times
          do i = 0,loop + skip - 1
            !write(6,*) nid, ": i =", i, __LINE__
#ifdef statuses
            call MPI_Recv(r_buf, size, MPI_BYTE, sender, TAG,
     &                    MPI_COMM_WORLD, reqstat, ierr)
#else
            call MPI_Recv(r_buf, size, MPI_BYTE, sender, TAG,
     &                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
#endif
            call MPI_Send(s_buf, size, MPI_BYTE, sender, TAG,
     &                    MPI_COMM_WORLD, ierr)
          enddo
        
        endif

#ifndef ALL_SIZES
        if (rsize.eq.0) then
          rsize = 1
        else
          do
            rsize = rsize * BASE
            if (int(rsize).ne.size) exit
            !if (nint(rsize).ne.size) exit
          enddo
        endif
        
        size = int(rsize)
        size = nint(rsize)
        
        if (size.gt.MAX_MSG_SIZE) return
#endif
      enddo
      
      !call MPI_Barrier(MPI_COMM_WORLD, ierr)
      
      return
      end subroutine

c* vi: set sw=4 sts=4 tw=80: */
