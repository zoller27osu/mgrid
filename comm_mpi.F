!#define catch_mpi_errs
!#define waitall_statuses

!#define use_tags !if commented, sources will be used

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

#ifdef catch_mpi_errs
      if (nid.eq.0) call MPI_Comm_set_errhandler(MPI_COMM_WORLD
     &                                  ,MPI_ERRORS_RETURN,ier)
      call gsync
#endif

      LRtype = 0
      UDtype = 0
      FBtype = 0

      return
      end
c-----------------------------------------------------------------------
      subroutine make_3d_types(mx1,my1,mz1,iOldType,enable_FBtype)!,LRtype,UDtype)
      include 'MGRID'
      logical :: enable_FBtype

      !FBtype = 0
      call MPI_TYPE_VECTOR(my1*mz1, 1, mx1, iOldType, LRtype, ierr)
      call MPI_TYPE_COMMIT(LRtype, ierr)
      call MPI_TYPE_VECTOR(mz1, mx1, mx1*my1, iOldType, UDtype, ierr)
      call MPI_TYPE_COMMIT(UDtype, ierr)
      if (enable_FBtype) then
        call MPI_TYPE_CONTIGUOUS(mx1*my1, iOldType, FBtype, ierr)
        call MPI_TYPE_COMMIT(FBtype, ierr)
      endif
      !if (nid.eq.0) write(6,*) LRtype, UDtype, FBtype

      return
      end
c-----------------------------------------------------------------------
      subroutine erase_3d_types()!LRtype,UDtype,FBtype)
      include 'MGRID'

      if (UDtype.ne.0) then
        call MPI_TYPE_FREE(LRtype, ierr)
        call MPI_TYPE_FREE(UDtype, ierr)
      endif
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
     &              ,jnid,msgtag,MPI_COMM_WORLD,status,ierr)

      if (len.gt.lenm) then 
          write(6,*) nid,'long message in mpi_crecv:',len,lenm
          call exitt
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine gsync()
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

      integer status(MPI_STATUS_SIZE)

      if (imsg.eq.MPI_REQUEST_NULL) return

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

      !call erase_3d_types(LRtype, UDtype) ! erases the MPI datatypes used in xchange
      call erase_3d_types()

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

      integer buf(1)

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
     &       ,MPI_COMM_WORLD,imsg,ierr)
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
     &                ,MPI_COMM_WORLD,imsg,ierr)
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
     &       ,MPI_COMM_WORLD,imsg,ierr)
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

      if ((msgtag.lt.0).or.(jnid.eq.MPI_PROC_NULL)) return!(jnid.lt.0).or.(jnid.ge.np)) return

      call mpi_send(x,num,itype,jnid,msgtag,MPI_COMM_WORLD,ierr)

      return
      end
c-----------------------------------------------------------------------
      function isend1(msgtag,x,itype,num,jnid)
      include 'mpif.h'
      common /cmgmpi/ nid,np,mgreal,wdsize
      integer wdsize
      !real*4 buf(1)

      if ((msgtag.lt.0).or.(jnid.eq.MPI_PROC_NULL)) then !(jnid.lt.0).or.(jnid.ge.np)) then
        isend1 = MPI_REQUEST_NULL
      else
        call mpi_isend(x,num,itype,jnid,msgtag,
     &                  MPI_COMM_WORLD,isend1,ierr)
      endif

      return
      end function
c-----------------------------------------------------------------------
      function irecv1(msgtag,x,itype,num) ! return if msgtag < 0
      include 'mpif.h'
      common /cmgmpi/ nid,np,mgreal,wdsize
      integer wdsize
#ifdef catch_mpi_errs
      character str(MPI_MAX_ERROR_STRING)
#endif

c     Note: len in bytes

      if (msgtag.lt.0) then
        irecv1 = MPI_REQUEST_NULL
      else
        !write(6,*) nid, ": posting recv with tag = ", msgtag
        call mpi_irecv (x,num,itype,mpi_any_source,msgtag
     &                    ,MPI_COMM_WORLD,irecv1,ierr)
#ifdef catch_mpi_errs
        if (ierr.ne.MPI_SUCCESS) then
          call MPI_Error_string(ierr,str,lenres,ierr)
          write(6,*) nid, ": ", str
        endif
#endif
      endif

c     write(6,*) nid,' irecv:',imsg,msgtag,len

      return
      end function
c-----------------------------------------------------------------------
      function irecv2(msgtag,x,itype,num,jnid) ! return if msgtag < 0
      include 'mpif.h'
      common /cmgmpi/ nid,np,mgreal,wdsize
      integer wdsize
#ifdef catch_mpi_errs
      character str(MPI_MAX_ERROR_STRING)
#endif

c     Note: len in bytes

#ifdef catch_mpi_errs
      ierr = MPI_SUCCESS
#endif
#ifdef use_tags
      if (msgtag.lt.0) then
        irecv2 = MPI_REQUEST_NULL
      else
        !write(6,*) nid, ": posting recv with tag = ", msgtag
        call mpi_irecv (x,num,itype,MPI_ANY_SOURCE,msgtag
     &                    ,MPI_COMM_WORLD,irecv2,ierr)
      endif
#else
c     ======  mpi_irecv can handle MPI_PROC_NULL just fine!  =======
        call mpi_irecv (x,num,itype,jnid,MPI_ANY_TAG
     &                    ,MPI_COMM_WORLD,irecv2,ierr)
#endif
#ifdef catch_mpi_errs
      if (ierr.ne.MPI_SUCCESS) then
        call MPI_Error_string(ierr,str,lenres,ierr)
        write(6,*) nid, ": ", str
      endif
#endif

c     write(6,*) nid,' irecv:',imsg,msgtag,len

      return
      end function
c-----------------------------------------------------------------------
      subroutine msgwaitall(imsgs,num)
      include 'mpif.h'
      common /cmgmpi/ nid,np,mgreal,wdsize
      integer wdsize

      integer imsgs(num)
#ifdef catch_mpi_errs
      character str(MPI_MAX_ERROR_STRING)
#endif
#ifdef waitall_statuses
      integer statuses(MPI_STATUS_SIZE,num)
#endif

       !write(6,*) nid, ": "
#ifdef waitall_statuses
      call MPI_Waitall(num,imsgs,statuses,ierr)
#else
      call MPI_Waitall(num,imsgs,MPI_STATUSES_IGNORE,ierr)
#endif
#ifdef catch_mpi_errs
      if (ierr.ne.MPI_SUCCESS) then
        call MPI_Error_string(ierr,str,lenres,ierr)
        write(6,*) nid, ": ", str
#ifdef waitall_statuses
        do i = 1,num
          if (statuses(MPI_ERROR,i).ne.MPI_SUCCESS) then
            write(6,*) nid, ": status(", i, ") error = "
            call MPI_Error_string(statuses(MPI_ERROR,i),str,lenres,ierr)
            write(6,*) nid, ": ", str
          endif
        enddo
#endif
        stop

      endif
#endif

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
        include 'MGRID'
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
c-----------------------------------------------------------------------
        function mateq(a,b,mnb)!mx,my,mz,mnb)
        include 'MGRID'
        logical :: mateq, ib, jb, kb
        real, dimension(0:mx,0:my,0:mz) :: a,b
        mateq = .TRUE.

        do k=0,mz
          kb = (k.eq.0).OR.(k.eq.mz)
          if ((mnb.le.0).AND.kb) cycle
        do j=0,my
          jb = (j.eq.0).OR.(j.eq.my)
          if ((mnb.le.1).AND.jb.AND.((mnb.le.0).OR.kb)) cycle
        do i=0,mx
          ib = (i.eq.0).OR.(i.eq.mx)
          if ((mnb.le.2).AND.ib.AND.((jb.AND.kb).OR.(mnb.le.0).OR.
     &        ((mnb.le.1).AND.(jb.OR.kb)) ) ) cycle

          !if (nid.eq.6) then
          !  write(6,17) nid,mnb,i,mx,j,my,k,mz,nz
  !17     !  format(i3,": mnb=",i2," i,j,k= "
     &    !     ,i4,"/",i4,", ",i4,"/",i4,", ",i4,"/",i4,i4)
          !endif

          if ((a(i,j,k).ne.b(i,j,k)).OR.ISNAN(a(i,j,k))
     &                              .OR.ISNAN(b(i,j,k))) then
            if (nid.eq.6) then
              write(6,*) "different at",i,j,k
     &                    ,"(",a(i,j,k),"vs",b(i,j,k),")"
            endif
            mateq = .FALSE.
            return
          endif
        enddo
        enddo
        enddo
        !if (nid.eq.6) write(6,*) nid, ": leaving mateq"

        return
        end function mateq
      END MODULE Comm_Funcs
