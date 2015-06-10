c-----------------------------------------------------------------------
      subroutine init_proc
      include 'mpif.h'
      include 'MGRID'

      double precision :: wtick

      call MPI_INIT(ierr)

      call MPI_COMM_SIZE(MPI_COMM_WORLD, np, ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, nid, ierr)

c     set mpi real type-
      wdsize=4
      eps=1.0e-12
      one_eps = 1.0+eps
      if (one_eps.ne.1.0) then
        wdsize=8
        if (nid.eq.0) write(6,*) "Wdsize = 8!"
      endif
      mgreal = mpi_real
      if (wdsize.eq.8) mgreal = mpi_double_precision

      if (nid.eq.0) then
        wtick = MPI_Wtick()
        write(6, *) "MPI_Wtick: ", wtick
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine make_3d_types(mx1,my1,mz1,iOldType)!,LRtype,UDtype)
      include 'MGRID'

      FBtype = 0
      call MPI_TYPE_VECTOR(my1*mz1, 1, mx1, iOldType, LRtype, ierr)
      call MPI_TYPE_COMMIT(LRtype, ierr)
      call MPI_TYPE_VECTOR(mz1, mx1, mx1*my1, iOldType, UDtype, ierr)
      call MPI_TYPE_COMMIT(UDtype, ierr)
      if (0.eq.1) then ! make false to safely disable FBtype
        call MPI_TYPE_CONTIGUOUS(mx1*my1, iOldType, FBtype, ierr)
        call MPI_TYPE_COMMIT(FBtype, ierr)
      endif
      !if (nid.eq.0) write(6,*) LRtype, UDtype, FBtype

      return
      end
c-----------------------------------------------------------------------
      subroutine erase_3d_types()!LRtype,UDtype,FBtype)
      include 'MGRID'

      !integer :: LRtype, UDtype

      call MPI_TYPE_FREE(LRtype, ierr)
      call MPI_TYPE_FREE(UDtype, ierr)
      if (FBtype.ne.0) call MPI_TYPE_FREE(FBtype, ierr)
      LRtype = 0
      UDtype = 0
      FBtype = 0

      return
      end
c-----------------------------------------------------------------------
      subroutine gop( x, w, op, n)
c
c     Global vector commutative operation
c
      include 'mpif.h'
      common /cmgmpi/ nid,np,mgreal,wdsize !Why? These are already in MGRID...
      integer wdsize
 
      real x(n), w(n)
      character*3 op
 
      if (op.eq.'+  ') then
         call mpi_allreduce (x,w,n,mgreal,mpi_sum ,MPI_COMM_WORLD,ierr)
      elseif (op.EQ.'M  ') then
         call mpi_allreduce (x,w,n,mgreal,mpi_max ,MPI_COMM_WORLD,ierr)
      elseif (op.EQ.'m  ') then
         call mpi_allreduce (x,w,n,mgreal,mpi_min ,MPI_COMM_WORLD,ierr)
      elseif (op.EQ.'*  ') then
         call mpi_allreduce (x,w,n,mgreal,mpi_prod,MPI_COMM_WORLD,ierr)
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

      call mpi_send (buf,len,mpi_byte,jnid,msgtag,MPI_COMM_WORLD,ierr)

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
     $              ,jnid,msgtag,MPI_COMM_WORLD,status,ierr)

      if (len.gt.lenm) then 
          write(6,*) nid,'long message in mpi_crecv:',len,lenm
          call exitt
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine gsync() !Why does this method even exist?

      include 'mpif.h'
      common /cmgmpi/ nid,np,mgreal,wdsize !Why?
      integer wdsize

      !write (6,*) nid, ": entering gsync"
      call mpi_barrier(MPI_COMM_WORLD,ierr)
      !write (6,*) nid, ": leaving gsync"

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
      subroutine csend0(msgtag,buf,len,jnid,jpid)
      include 'mpif.h'
      common /cmgmpi/ nid,np,mgreal,wdsize
      integer wdsize
      real*4 buf(1)

!     return if msgtag < 0 or jnid < 0 or jnid > np
      if (msgtag.lt.0)  return
      if (jnid  .lt.0)  return
      if (jnid  .ge.np) return

      call mpi_send (buf,len,mpi_byte,jnid,msgtag,MPI_COMM_WORLD,ierr)

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

      call mpi_bcast (buf,len,mpi_byte,0,MPI_COMM_WORLD,ierr)

      return
      end
c-----------------------------------------------------------------------
      function isend(msgtag,x,len,jnid,jpid)
c
c     Note: len in bytes
c
      !integer x(1)
   
      include 'mpif.h'
      common /cmgmpi/ nid,np,mgreal,wdsize
      integer wdsize
   
      call mpi_isend (x,len,mpi_byte,jnid,msgtag
     $       ,MPI_COMM_WORLD,imsg,ierr)
      isend = imsg
c     write(6,*) nid,' isend:',imsg,msgtag,len,jnid,(x(k),k=1,len/4)
   
      return
      end function
c-----------------------------------------------------------------------
      function irecv(msgtag,x,len)
c
c     Note: len in bytes
c
      !integer x(1)
   
      include 'mpif.h'
      common /cmgmpi/ nid,np,mgreal,wdsize
      integer wdsize
   
      call mpi_irecv (x,len,mpi_byte,mpi_any_source,msgtag
     $                ,MPI_COMM_WORLD,imsg,ierr)
      irecv = imsg
c     write(6,*) nid,' irecv:',imsg,msgtag,len
   
      return
      end function
c-----------------------------------------------------------------------
      function irecv0(msgtag,x,len) ! return if msgtag < 0
      include 'mpif.h'
      common /cmgmpi/ nid,np,mgreal,wdsize
      integer wdsize

c     Note: len in bytes

      !integer x(1)

      if (msgtag.lt.0) then ! here we assume x is real, based on current usage
         n = len/wdsize
         call rzero(x,n)
         irecv0 = mpi_request_null
      else
         call mpi_irecv (x,len,mpi_byte,mpi_any_source,msgtag
     $       ,MPI_COMM_WORLD,imsg,ierr)
         irecv0 = imsg
      endif

c     write(6,*) nid,' irecv:',imsg,msgtag,len

      return
      end function

c-----------------------------------------------------------------------
      subroutine csend1(msgtag,x,itype,num,jnid) 
      include 'mpif.h'
      common /cmgmpi/ nid,np,mgreal,wdsize
      integer wdsize
      !real*4 buf(1)

      if ((msgtag.lt.0).or.(jnid.lt.0).or.(jnid.ge.np)) return

      call mpi_send(x,num,itype,jnid,msgtag,MPI_COMM_WORLD,ierr)

      return
      end
c-----------------------------------------------------------------------
      function isend1(msgtag,x,itype,num,jnid)
      include 'mpif.h'
      common /cmgmpi/ nid,np,mgreal,wdsize
      integer wdsize
      !real*4 buf(1)

      if ((msgtag.lt.0).or.(jnid.lt.0).or.(jnid.ge.np)) return

      call mpi_isend(x,num,itype,jnid,msgtag,MPI_COMM_WORLD,isend1,ierr)

      return
      end function
c-----------------------------------------------------------------------
      function irecv1(msgtag,x,itype,num) ! return if msgtag < 0
      include 'mpif.h'
      common /cmgmpi/ nid,np,mgreal,wdsize
      integer wdsize

c     Note: len in bytes

      if (msgtag.lt.0) then
         irecv1 = mpi_request_null
      else
         !write(6,*) nid, ": posting recv with tag = ", msgtag
         call mpi_irecv (x,num,itype,mpi_any_source,msgtag
     $                    ,MPI_COMM_WORLD,irecv1,ierr)
         !irecv1 = imsg
      endif

c     write(6,*) nid,' irecv:',imsg,msgtag,len

      return
      end function
c-----------------------------------------------------------------------
      subroutine msgwaitall(imsgs,num)
      include 'mpif.h'
      common /cmgmpi/ nid,np,mgreal,wdsize
      integer wdsize

      integer imsgs(num)
      integer statuses(MPI_STATUS_SIZE,num)
      character str(MPI_MAX_ERROR_STRING)

      !if (imsg.eq.mpi_request_null) return

c     write(6,*) nid,' msgwait:',imsg
      call MPI_Comm_set_errhandler(MPI_COMM_WORLD,MPI_ERRORS_RETURN,ier)
      call MPI_Waitall(num,imsgs,statuses,ierr)
      if (ierr.ne.MPI_SUCCESS) then
        call MPI_Error_string(ierr,str,reslen,ierr)
        write(6,*) nid, ": ", str
        do i = 1,num
          if (statuses(MPI_ERROR,i).ne.MPI_SUCCESS) then
            write(6,*) nid, ": status(", i, ") error = "
            call MPI_Error_string(statuses(MPI_ERROR,i),str,reslen,ierr)
            write(6,*) nid, ": ", str
          endif
        enddo

        stop

      endif
      call MPI_Comm_set_errhandler(MPI_COMM_WORLD
     &                             ,MPI_ERRORS_ARE_FATAL,ierr)

      return
      end
c-----------------------------------------------------------------------
      MODULE Comm_Funcs
      CONTAINS
        function dclock()                               ! returns SECONDS
        double precision :: dclock
        include 'mpif.h'

        !if (nid.eq.0) write(6,*) "MPI_WTIME(): ", MPI_WTIME()
        dclock=mpi_wtime()
        !if (nid.eq.0) write(6,*) "dclock: ", dclock
        !call cpu_time(dclock)

        !call system_clock(itime, iresolution)
        !write(6,*) itime, iresolution
        !dclock=double precision(itime)/double precision(iresolution)

        return
        end function
c-----------------------------------------------------------------------
        function arreq(a,b,n)
        logical :: arreq
        real a(1),b(1)
        arreq = .TRUE.
        do idx=1,n
          if (a(idx).ne.b(idx)) then
            if (nid.eq.6) then
              write(6,*) "different at",idx,"(",a(idx),"vs",b(idx),")"
            endif
            arreq = .FALSE.
            return
          endif
        enddo
        return
        end function arreq

        function mateq(a,b,mnb)!mx,my,mz,mnb)
        include 'MGRID'
        logical :: mateq
        real, dimension(0:mx,0:my,0:mz) :: a,b
        integer atbound = 0
        mateq = .TRUE.
        do k=0,mz
          if ((k.eq.0).OR.(k.eq.mz)) atbound = atbound + 1
          if (k.eq.1) atbound = atbound - 1
          if (atbound.gt.mnb) cycle
        do j=0,my
          if ((j.eq.0).OR.(j.eq.my)) atbound = atbound + 1
          if (j.eq.1) atbound = atbound - 1
          if (atbound.gt.mnb) cycle
        do i=0,mx
          if ((i.eq.0).OR.(i.eq.mx)) atbound = atbound + 1
          if (i.eq.1) atbound = atbound - 1
          if (atbound.gt.mnb) cycle
          write(6,1) nid,mnb,i,mx,j,my,k,mz
    1 format(i3,": mnb=",i2," atbound=",i2," i,j,k=",i4,"/",i4,", ",i4,"/",i4,", ",i4,"/",i4)
          if (nid.eq.6) write(6,*) nid, ":", mnb, i, j, k
          if (a(i,j,k).ne.b(i,j,k)) then
            if (nid.eq.6) then
              write(6,*) "different at",i,j,k
     $                    ,"(",a(i,j,k),"vs",b(i,j,k),")"
            endif
            mateq = .FALSE.
            return
          endif
        enddo
        atbound = atbound - 1
        enddo
        atbound = atbound - 1
        enddo
        return
        end function mateq
      END MODULE Comm_Funcs
