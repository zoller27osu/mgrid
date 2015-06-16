c-----------------------------------------------------------------------
      subroutine init_proc
      include 'mpif.h'
      include 'MGRID'

      call MPI_INIT(ierr)

      call MPI_COMM_SIZE(MPI_COMM_WORLD, np, ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, nid, ierr)

c     set mpi real type-
      wdsize=4
      eps=1.0e-12
      oneeps = 1.0+eps
      if (oneeps.ne.1.0) then
         wdsize=8
      endif
      mgreal = mpi_real
      if (wdsize.eq.8) mgreal = mpi_double_precision

      return
      end
c-----------------------------------------------------------------------
      subroutine gop( x, w, op, n)
c
c     Global vector commutative operation
c
      include 'mpif.h'
      common /cmgmpi/ nid,np,mgreal,wdsize
      integer wdsize

      real x(n), w(n)
      character*3 op

      if (op.eq.'+  ') then
         call mpi_allreduce (x,w,n,mgreal,mpi_sum ,mpi_comm_world,ierr)
      elseif (op.EQ.'M  ') then
         call mpi_allreduce (x,w,n,mgreal,mpi_max ,mpi_comm_world,ierr)
      elseif (op.EQ.'m  ') then
         call mpi_allreduce (x,w,n,mgreal,mpi_min ,mpi_comm_world,ierr)
      elseif (op.EQ.'*  ') then
         call mpi_allreduce (x,w,n,mgreal,mpi_prod,mpi_comm_world,ierr)
      else
         write(6,*) nid,' OP ',op,' not supported.  ABORT in GOP.'
         call exitt
      endif

      call copy(x,w,n)

      return
      end
c-----------------------------------------------------------------------
      subroutine csend(msgtag,buf,len,jnid,jpid)
      include 'mpif.h'
      common /cmgmpi/ nid,np,mgreal,wdsize
      integer wdsize
      real*4 buf(1)

      call mpi_send (buf,len,mpi_byte,jnid,msgtag,mpi_comm_world,ierr)

      return
      end
c-----------------------------------------------------------------------
      subroutine crecv(msgtag,buf,lenm)
      include 'mpif.h'
      common /cmgmpi/ nid,np,mgreal,wdsize
      integer wdsize
      integer status(mpi_status_size)

      real*4 buf(1)
      len = lenm
      jnid = mpi_any_source

      call mpi_recv (buf,len,mpi_byte
     $              ,jnid,msgtag,mpi_comm_world,status,ierr)

      if (len.gt.lenm) then
          write(6,*) nid,'long message in mpi_crecv:',len,lenm
          call exitt
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine gsync()

      include 'mpif.h'
      common /cmgmpi/ nid,np,mgreal,wdsize
      integer wdsize

      call mpi_barrier(mpi_comm_world,ierr)

      return
      end
c-----------------------------------------------------------------------
      subroutine msgwait(imsg)
      include 'mpif.h'
      common /cmgmpi/ nid,np,mgreal,wdsize
      integer wdsize

      integer status(mpi_status_size)

      if (imsg.eq.mpi_request_null) return

c     write(6,*) nid,' msgwait:',imsg
      call mpi_wait (imsg,status,ierr)

      return
      end
c-----------------------------------------------------------------------
      function isend(msgtag,x,len,jnid,jpid)
c
c     Note: len in bytes
c
      integer x(1)

      include 'mpif.h'
      common /cmgmpi/ nid,np,mgreal,wdsize
      integer wdsize

      call mpi_isend (x,len,mpi_byte,jnid,msgtag
     $       ,mpi_comm_world,imsg,ierr)
      isend = imsg
c     write(6,*) nid,' isend:',imsg,msgtag,len,jnid,(x(k),k=1,len/4)

      return
      end
c-----------------------------------------------------------------------
      function irecv(msgtag,x,len)
c
c     Note: len in bytes
c
      integer x(1)

      include 'mpif.h'
      common /cmgmpi/ nid,np,mgreal,wdsize
      integer wdsize

      call mpi_irecv (x,len,mpi_byte,mpi_any_source,msgtag
     $       ,mpi_comm_world,imsg,ierr)
      irecv = imsg
c     write(6,*) nid,' irecv:',imsg,msgtag,len


      return
      end
c-----------------------------------------------------------------------
      function irecv0(msgtag,x,len) ! return if msgtag < 0
      include 'mpif.h'
      common /cmgmpi/ nid,np,mgreal,wdsize
      integer wdsize

c     Note: len in bytes

      integer x(1)


      if (msgtag.lt.0) then ! here we assume x is real, based on current usage
         n = len/wdsize
         call rzero(x,n)
         irecv0 = mpi_request_null
      else
         call mpi_irecv (x,len,mpi_byte,mpi_any_source,msgtag
     $       ,mpi_comm_world,imsg,ierr)
         irecv0 = imsg
      endif

c     write(6,*) nid,' irecv:',imsg,msgtag,len

      return
      end
c-----------------------------------------------------------------------
      subroutine csend0(msgtag,buf,len,jnid,jpid)

!     return if msgtag < 0 or jnid < 0 or jnid > np

      include 'mpif.h'
      common /cmgmpi/ nid,np,mgreal,wdsize
      integer wdsize
      real*4 buf(1)

      if (msgtag.lt.0)  return
      if (jnid  .lt.0)  return
      if (jnid  .ge.np) return

      call mpi_send (buf,len,mpi_byte,jnid,msgtag,mpi_comm_world,ierr)

      return
      end
c-----------------------------------------------------------------------
      subroutine exitt

      call gsync
      call mpi_finalize (ierr)
      call exit(0)

      return
      end
c-----------------------------------------------------------------------
      subroutine bcast(buf,len)
      include 'mpif.h'
      common /cmgmpi/ nid,np,mgreal,wdsize
      integer wdsize

      real*4 buf(1)

      call mpi_bcast (buf,len,mpi_byte,0,mpi_comm_world,ierr)

      return
      end
c-----------------------------------------------------------------------
      function dclock()
      include 'mpif.h'

      dclock=mpi_wtime()

      return
      end
c-----------------------------------------------------------------------
