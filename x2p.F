c-----------------------------------------------------------------------
c
c     wdsize   - working precision
c
c     nx,ny,nz - global dimensions, current MG level
c     mx,my,mz - local dimensions,  current MG level (=nx,ny,nz if np=1)
c     mm       - local point count, inc. ghosts,      current MG level
c
c     st     - global start index for each proc. inc. ghosts
c     istx   - global start index for X on each proc. inc. ghosts
c
c     nf/n   - number of (fine) points on each side, total
c     kf/kc  - offset of each level of the 'v-cycle' (fine/coarse)
c
c     nsmooth - number of times to smooth
c
c-----------------------------------------------------------------------

#define XCHANGE_JGX 1
#define JAC3D_JGX 0
!#define enable_FBtype

#define verify_xchange
#define verify_jacobi_alts
#define allow_math_differences

!#define conservative_mpi_dtype_recreate

#define NPASS_DEFAULT 100
#define NSMOOTH_DEFAULT 21
#define IGS_DEFAULT 2
#define SIGMA_DEFAULT 0.66666666666666667

      program mgrid_test
      include 'MGRID'

      parameter(ml=2*lpow)
      integer tflag(ml),dmflag(ml)
      common /cflag/ tflag,dmflag
      common /cinit/ icalld

      integer nlev_last
      save    nlev_last
      data    nlev_last /0/

      call init_proc   ! Initialize parallel mode

      if (nid.eq.0) write(6,*) nid, ": num MPI procs = ", np

      call get_proc_mpqr (mp,mq,mr,np,ldim)             ! processor decomp
      call get_proc_ijk  (ip,iq,ir,nid,mp,mq,mr,0,ldim) ! nid = (ip,iq,ir)

      do itest=1,2
        icalld = 0
        do i = 1,90 ! parses up to 90 test lines from in.dat
           call mgrid_tester(nlev_last)
        enddo
        close(9)
      enddo

      call exitt
      end

c-----------------------------------------------------------------------
      subroutine mgrid_tester(nlev_last)
      USE Comm_Funcs, ONLY: dclock
      include 'MGRID'

      double precision :: time0, time1

      common /xcdata/ u(0:lt),r(0:lt)
      common /xddata/ ue(0:lt),rf(0:lt)
      real sigma

      parameter(ml=2*lpow)
      integer tflag(ml),dmflag(ml)
      common /cflag/ tflag,dmflag
      common /cinit/ icalld

      !character ig_labels(5)

      call init_mg(npass,nsmooth,igs,sigma)       ! Get current mg parameters

      if (nlev.le.nlev_last) return
      nlev_last = nlev

      call reset_flops ! Reset flop counters

      lev = nlev
      nf  = 2**nlev    ! Global number of points on finest grid in each direction
      n   = nf
      h   = 1./n                 
      kf  = 0
      
      call izero     ( tflag,ml)
      call izero     (dmflag,ml)

      call dmesh     (mm,n) ! calculates mx, my, mz

      call rzero     (u ,mm)                    ! clear u
      call rzero     (r ,mm)                    ! clear r
      call rzero     (rf,mm)                    ! clear rf
      call rzero     (ue,mm)                    ! clear ue

      call initu     (ue,n,h,ldim)              ! Exact solution (Get exact solution in ue?)
      call xchange0  (ue,XCHANGE_JGX)           ! share ue

      call compute_r (r,ue,rf,mx,my,mz,h,ldim)  ! b = 0 - A ue
      call chsign    (r,mm)                     ! b = A ue
      call xchange0  (r,XCHANGE_JGX)            ! share r? (Why?)

      tol   = 1.e-10

      lev = nlev
      nf  = 2**lev
      n   = nf
      h   = 1./n
      kf  = 0

      call dmesh(mm,n) ! calculates mx, my, mz

      iz=1

      dmax = 0.

c     ======== Test Xchange Alts=========
#ifdef verify_xchange
      call test_xchange(u,1,2) ! last 2 parameters are min/max "jgx"
#endif

c     ======== Test Jacobi Alts =========
#ifdef verify_jacobi_alts
      call test_alts(u(kf),r(kf),h,sigma,nsmooth,rf,iz,2,4) ! last 2 parameters are min/max "igs"
#endif

c     ========   Time Xchange  =========
      call gsync
      time0 = dclock()
      !if (nid.eq.1) write(6,*) time0
      !if (nid.eq.1) write(6,"(F15.3)") time0
      do ipass=1,npass*200
         call xchange0(u,1)
      enddo
      call gsync
      time1 = dclock()
      if (nid.eq.1) write(6,*) time0, time1
      call errchk (emx,u,ue,ipass-1,time1-time0,dmax,'D') !xchange_3d

      call gsync
      time0 = dclock()
      !if (nid.eq.1) write(6,*) time0
      !if (nid.eq.1) write(6,"(F15.3)") time0
      do ipass=1,npass*200
         call xchange0(u,2)
      enddo
      call gsync
      time1 = dclock()
      !if (nid.eq.1) write(6,*) time0, time1
      call errchk (emx,u,ue,ipass-1,time1-time0,dmax,'Y') !xchange_3d_async

      call gsync
      time0 = dclock()
      !if (nid.eq.1) write(6,*) time0
      !if (nid.eq.1) write(6,"(F15.3)") time0
      do ipass=1,npass*200
         call xchange0(u,0)
      enddo
      call gsync
      time1 = dclock()
      !if (nid.eq.1) write(6,*) time0, time1
      call errchk (emx,u,ue,ipass-1,time1-time0,dmax,'O') !xchange_3d_old
      if (nid.eq.1) write(6,*) " "

c     ========    Time Alts    =========
      call gsync
      time0 = dclock()
      !if (nid.eq.1) write(6,*) time0
      !if (nid.eq.1) write(6,"(F15.3)") time0
      do ipass=1,npass
         call smooth_q (u(kf),r(kf),mx,my,mz,h,sigma,nsmooth,rf,iz)
      enddo
      call gsync
      time1 = dclock()
      if (nid.eq.1) write(6,*) time0, time1
      call errchk (emx,u,ue,ipass-1,time1-time0,dmax,'Q') !smooth_q

      call gsync
      time0 = dclock()
      !if (nid.eq.1) write(6,*) time0
      do ipass=1,npass
         call smooth_m  (u(kf),r(kf),mx,my,mz,h
     &                          ,sigma,nsmooth,rf,iz,2,ldim)
         call xchange0  (u(kf),XCHANGE_JGX)
      enddo
      call gsync
      time1 = dclock()
      !if (nid.eq.1) write(6,*) time0, time1
      call errchk (emx,u,ue,ipass-1,time1-time0,dmax,'A') !jac3alt

      call gsync
      time0 = dclock()
      !if (nid.eq.1) write(6,*) time0
      do ipass=1,npass
         call smooth_m  (u(kf),r(kf),mx,my,mz,h
     &                          ,sigma,nsmooth,rf,iz,3,ldim)
         call xchange0  (u(kf),XCHANGE_JGX)
      enddo
      call gsync
      time1 = dclock()
      !if (nid.eq.1) write(6,*) time0, time1
      call errchk (emx,u,ue,ipass-1,time1-time0,dmax,'C') !jac3cnd

      call gsync
      time0 = dclock()
      !if (nid.eq.1) write(6,*) time0
      do ipass=1,npass
         call smooth_m  (u(kf),r(kf),mx,my,mz,h
     &                          ,sigma,nsmooth,rf,iz,4,ldim)
         call xchange0  (u(kf),XCHANGE_JGX)
      enddo
      call gsync
      time1 = dclock()
      !if (nid.eq.1) write(6,*) time0, time1
      call errchk (emx,u,ue,ipass-1,time1-time0,dmax,'F') !jac3flg

      call gsync
      time0 = dclock()
      !if (nid.eq.1) write(6,*) time0
      do ipass=1,npass
         call smooth_m  (u(kf),r(kf),mx,my,mz,h
     &                          ,sigma,nsmooth,rf,iz,0,ldim)
         call xchange0  (u(kf),XCHANGE_JGX)
      enddo
      call gsync
      time1 = dclock()
      !if (nid.eq.1) write(6,*) time0, time1
      call errchk (emx,u,ue,ipass-1,time1-time0,dmax,'J') !jac3d

      call gsync
      time0 = dclock()
      do ipass=1,NPASS_DEFAULT/5 ! Weird NaN stuff in rf happens when this number gets too high!

        lev = nlev
        nf  = 2**lev
        n   = nf
        h   = 1./n
        kf  = 0
        call dmesh(mm,n)

        do ilev=1,nlev-1

            if (ilev.eq.1) then ! retain boundary info
               iz=1
            else               ! homogeneous boundary data for defect correction
               iz=0
               call rzero(u(kf),mm)
            endif

            call smooth_m  (u(kf),r(kf),mx,my,mz,h
     &                          ,sigma,nsmooth,rf,iz,igs,ldim)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !if (nid.eq.6) write(6, *) nid, ": xchanging u(kf)",__LINE__
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            call xchange0  (u(kf),XCHANGE_JGX)

            call compute_r (rf,u(kf),r(kf),mx,my,mz,h,ldim)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !if (nid.eq.6) then
            !  write(6, *) nid,": xchanging rf",__LINE__,"kf,lt=",kf,lt
            !  call print_3d(rf, 5, 2, 2)
            !endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            call xchange0  (rf,XCHANGE_JGX)

            lev     = lev-1
            kf      = kf + mm
            n       = n/2
            h       = 1./n

            mxf    = mx      ! mxf is fine element count
            myf    = my
            mzf    = mz

            ieven = 2-mod(istx,2)
            jeven = 2-mod(isty,2)
            keven = 2-mod(istz,2)
            
            call dmesh(mm,n)  ! find mesh decomposition & msg tags
            call restrict_i (r(kf),rf,mxf,myf,mzf
     &                               ,ieven,jeven,keven,ldim)
        enddo
    
        one = 1.
        iz  = 0
        call smooth_m(u(kf),r(kf),mx,my,mz,h
     &                              ,one,1,rf,iz,igs,ldim)
        call reset_mtags (lev+1)  ! Exchange on prior level processor set
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !if (nid.eq.6) write(6, *) nid, ": xchanging u(kf)", __LINE__
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call xchange0    (u(kf),XCHANGE_JGX)

        kc = kf
        nf = n*2
        h  = h/2

        do ilev=nlev-1,1,-1

           mxc   = mx        ! mxc is coarse element count
           myc   = my
           mzc   = mz

           lev = lev+1
           call dmesh(mm,nf) !Why?
           kf     = kc-mm

           ieven = 2-mod(istx,2)
           jeven = 2-mod(isty,2)
           keven = 2-mod(istz,2)
        
           ! prolongate and add

           call prolong_a (u(kf),dmax,u(kc),mxc,myc,mzc
     &                                ,ieven,jeven,keven,ldim) 
           call reset_mtags (lev+1)  ! Exchange on prior level processor set
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           !if (nid.eq.6) write(6, *) nid, ": xchanging u(kf)", __LINE__
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           call xchange0    (u(kf),XCHANGE_JGX)

           n  = n*2
           h  = 1./n
           kc = kf
           nf = n*2

        enddo

c       time1 = dclock()-time0
c       call errchk (emx,u,ue,ipass,time1,dmax,'I')

c       if (ipass.gt.10) then
c          dmax = glmax(dmax,1)  ! Check for finish
c          dmax = glmax(dm,1)
c          if (dmax.le.tol) goto 100
c        endif

      enddo
  100 continue

      time1 = dclock()    
      call errchk (emx,u,ue,ipass-1,time1-time0,dmax,'X')
c     call outmp(u ,ipass,lev,'SOLN',rf,r)
c     call outmp(ue,ipass,lev,'EXAC',rf,r)

      return
      end
c-----------------------------------------------------------------------
      subroutine errchk(emx,u,ue,ii,time,dmax,c1)
      include 'MGRID'
      real u(1),ue(1)
      integer*8 n3
      character*1 c1
      double precision :: time

      common /ccflops/ cflops,pflops,rflops,rrlops

      mm = (mx+1)*(my+1)*(mz+1)

      emx = 0
      do k=1,mm
         diff = u(k)-ue(k)
         emx = max(emx,abs(diff))
      enddo
      emx = glmax(emx,1)

      n3 = nx*ny
      if (ldim.eq.3) n3=n3*nz
      if (nid.eq.1) then
        write(6,*) 
     &"   ii   dmax        emx         time    ldim  nlev  nx    ny   ",
     &" nz          n3       mp   mq   mr       np          c1"
        write(6,1) 
     &    ii,dmax,emx,time,ldim,nlev,nx,ny,nz,n3,mp,mq,mr,np,c1

    1   format(i7,1p,3e12.4,i2,i3,3i6,i13,'=n ',3i5,i9,'  max error',a1)
      endif

      cflops = cflops*np
      pflops = pflops*np
      rflops = rflops*np
      rrlops = rrlops*np
      fflops = (cflops+pflops+rflops+rrlops)
      if (time.gt.0) fflops = fflops/(time*1.e6)

      if (nid.eq.1) then
        write(6,*)
     &"   ii   emx         time        fflops      cflops      pflops ",
     &"     rflops      rrlops               n3           np      c1"
        write(6,2) 
     &    ii,emx,time,fflops,cflops,pflops,rflops,rrlops,n3,np,c1

    2   format(i7,1p,7e12.4,i13,'=n ',i9,'  flops',a1)
      endif

      call reset_flops ! Reset flop counters

      return
      end
c-----------------------------------------------------------------------
      function glmin(a,n)
      real a(1)
      
      tmin = 99.0e20
      do i=1,n
         tmin=min(tmin,a(i))
      enddo
      call gop(tmin,work,'m  ',1)
      glmin = tmin
      return
      end
c-----------------------------------------------------------------------
      function glmax(a,n)
      real a(1)
      tmax=-99.0e20
      do i=1,n
         tmax=max(tmax,a(i))
      enddo
      call gop(tmax,work,'M  ',1)
      glmax=tmax
      return
      end
c-----------------------------------------------------------------------
      function mod1(i,n)

c     Yields MOD(I,N) with the exception that if I=K*N, result is N.

      mod1=0
      if (i.eq.0) return

      if (n.eq.0) then
         write(6,*)
     &   'WARNING:  Attempt to take MOD(I,0) in function mod1.'
         return
      endif
      ii = i+n-1
      mod1 = mod(ii,n)+1
      return
      end
c-----------------------------------------------------------------------
      integer function log2(k)
      rk=(k)
      rlog=log10(rk)
      rlog2=log10(2.0)
      rlog=rlog/rlog2+0.5
      log2=int(rlog)
      return
      end
c-----------------------------------------------------------------------
      subroutine rzero(a,n)
      real a(1)
      do i=1,n
         a(i)=0
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine izero(a,n)
      integer a(1)
      do i=1,n
         a(i)=0
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine sub3(a,b,c,n)
      real a(1),b(1),c(1)
      do i=1,n
         a(i)=b(i)-c(i)
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine copy(a,b,n)
      real a(1),b(1)
      do i=1,n
         a(i)=b(i)
         !write(6, *) nid, ": i = ", i
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine chsign(a,n)
      real a(1)
      do i=1,n
         a(i)=-a(i)
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine checkit(s10,n_used,n_max)

      character*10 s10

      include 'MGRID'

c     n_used_m = iglmax(n_used,1)
      if (n_used.le.n_max) return
      
      write(6,*) nid,n_used,n_max,' ERROR checkit: ',s10
      call exitt

      stop
      end
c-----------------------------------------------------------------------
      subroutine reset_flops ! Reset flop counters
      common /ccflops/ cflops,pflops,rflops,rrlops

      cflops = 0.
      pflops = 0.
      rflops = 0.
      rrlops = 0.

      return
      end
c-----------------------------------------------------------------------
      subroutine init_mg(npass,nsmooth,igs,sigma)! Get current mg parameters
      include 'MGRID'
      common /cinit/ icalld
     
      nlev_max = (lpow-1) + log2(np)/ldim

      if (nid.eq.0) then
        if (icalld.eq.0) open(unit=9,file='in.dat')
        nlev = 0
        npass = NPASS_DEFAULT
        nsmooth = NSMOOTH_DEFAULT
        igs = IGS_DEFAULT
        sigma = SIGMA_DEFAULT
        read (9,*,end=999)  nlev,npass,nsmooth,igs,sigma
        nlev = min(nlev_max,nlev)
        write(6,*) 'Input nlev, npass, nsmooth, igs, sigma:'
        write(6,*)  nlev,npass,nsmooth,igs,sigma
  999   if (nlev.eq.0) write(6,*) 'nlev = 0, stopping!'
      endif

      icalld = 1

      call bcast(nlev   ,     4) ! 4-byte integer
      call bcast(npass  ,     4)
      call bcast(nsmooth,     4)
      call bcast(igs    ,     4)
      call bcast(sigma  ,wdsize)
      if (nlev.eq.0) stop

      return
      end
c-----------------------------------------------------------------------
      subroutine get_proc_mpqr(mp,mq,mr,np,ldim)

c     Given np processors, find a decomposition such that
c     mp*mq*mr = np, with mp,mq,mr, and np integer
c     not just used for power of 2, but best decomp is for pow2

      if (ldim.eq.3) then 
        rp = np**(1./3.)
        mp = rp*(1.01)
        do ip=mp,1,-1
           mqr = np/ip
           nrem = np-mqr*ip
           if (nrem.eq.0) goto 10
        enddo
   10   mp = ip

        rq = mqr**(1./2.)
        mq = rq*(1.01)
        do iq=mq,1,-1
           mr = mqr/iq
           nrem = mqr-mr*iq
           if (nrem.eq.0) goto 20
        enddo
   20   mq = iq

        mr = np/(mp*mq)
      else
        rp = np**(1./2.)
        mp = rp*(1.01)
        do ip=mp,1,-1
           mqr = np/ip
           nrem = np-mqr*ip
           if (nrem.eq.0) goto 30
        enddo
   30   mp = ip
  
        mq = np/mp
        mr = 1
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine get_proc_ijk(ip,iq,ir,nid,mp,mq,mr,i0,ldim) 
c     i0 ==> zero-based index
c
c     Gets the (i,j,k) coordinate of the nid
c     In symmetric corner ordering
c

      mpq = mp*mq
      n1  = nid + (1-i0)    ! i0=0 --> zero-based index

      ir = 1 +  (n1-1)/mpq
      if (ldim.eq.2) ir = 1 

      iq = mod1 (n1,mpq)
      iq = 1 +  (iq-1)/mp
      ip = mod1 (n1,mp)

      if (i0.eq.0) ip = ip-1
      if (i0.eq.0) iq = iq-1
      if (i0.eq.0) ir = ir-1
 
      return
      end
c-----------------------------------------------------------------------
      subroutine mtag
      include 'MGRID'
      integer gid

      parameter(ml=2*lpow)
      integer tflag(ml),dmflag(ml)
      common /cflag/ tflag,dmflag

      if (tflag(lev).eq.0) then
c
c      Assigns destination nid/mtag for xchanges where nid = (ip,iq,ir)
c      But neighbor is not adjacent processor
c 
      
c      if (nid.eq.1) write(6,*) lev,nlev,ipass,nx,' level'

       call find_neighb (jpx,ip,mp,nx, 1) ! Right Neighbor
       call get_gid     (gid,jpx, iq, ir,mp,mq,mr)
       call get_mtag    (mtagsr,mtagrr,gid,nid)
       nidr = gid

       call find_neighb (jmx,ip,mp,nx,-1) ! Left Neighbor
       call get_gid     (gid,jmx, iq, ir,mp,mq,mr)
       call get_mtag    (mtagsl,mtagrl,gid,nid)
       nidl = gid

       call find_neighb (jpy,iq,mq,ny, 1) ! Top Neighbor
       call get_gid     (gid, ip,jpy, ir,mp,mq,mr)
       call get_mtag    (mtagsu,mtagru,gid,nid)
       nidu = gid

       call find_neighb (jmy,iq,mq,ny,-1) ! Bottom Neighbor
       call get_gid     (gid, ip,jmy, ir,mp,mq,mr)
       call get_mtag    (mtagsd,mtagrd,gid,nid)
       nidd = gid

       call find_neighb (jpz,ir,mr,nz, 1) ! Front Neighbor
       call get_gid     (gid, ip, iq,jpz,mp,mq,mr)
       call get_mtag    (mtagsf,mtagrf,gid,nid)
       nidf = gid

       call find_neighb (jmz,ir,mr,nz,-1) ! Back Neighbor
       call get_gid     (gid, ip, iq,jmz,mp,mq,mr)
       call get_mtag    (mtagsb,mtagrb,gid,nid)
       nidb = gid

       ilev = min(lev,nlev)
       if (lev.eq.0) ilev = nlev

       nnid(1,ilev)  = nidr
       nnid(2,ilev)  = nidl
       nnid(3,ilev)  = nidu
       nnid(4,ilev)  = nidd
       nnid(5,ilev)  = nidf
       nnid(6,ilev)  = nidb
       
       ntags(1,ilev) = mtagsr
       ntags(2,ilev) = mtagsl
       ntags(3,ilev) = mtagsu
       ntags(4,ilev) = mtagsd
       ntags(5,ilev) = mtagsf
       ntags(6,ilev) = mtagsb

       ntagr(1,ilev) = mtagrr
       ntagr(2,ilev) = mtagrl
       ntagr(3,ilev) = mtagru
       ntagr(4,ilev) = mtagrd
       ntagr(5,ilev) = mtagrf
       ntagr(6,ilev) = mtagrb

       tflag(lev) = 1

      else

        call reset_mtags(lev)

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine reset_mtags(jlev)
      include 'MGRID'

      ilev = min(jlev,nlev)
      if (ilev.eq.0) ilev = nlev

      nidr = nnid(1,ilev)
      nidl = nnid(2,ilev)
      nidu = nnid(3,ilev)
      nidd = nnid(4,ilev)
      nidf = nnid(5,ilev)
      nidb = nnid(6,ilev)
      
      mtagsr = ntags(1,ilev)
      mtagsl = ntags(2,ilev)
      mtagsu = ntags(3,ilev)
      mtagsd = ntags(4,ilev)
      mtagsf = ntags(5,ilev)
      mtagsb = ntags(6,ilev)

      mtagrr = ntagr(1,ilev)
      mtagrl = ntagr(2,ilev)
      mtagru = ntagr(3,ilev)
      mtagrd = ntagr(4,ilev)
      mtagrf = ntagr(5,ilev)
      mtagrb = ntagr(6,ilev)

c     write(6,1) nid,lev,ilev,jlev
c    &   ,nidu,nidd,(ntags(k,ilev),ntagr(k,ilev),k=3,4)
c  1  format(i3,'TAGS:',9i5)

      return
      end
c-----------------------------------------------------------------------
      subroutine find_neighb(jp,ip,mp,n,jdir)
c
c     jp is the ip coordinate for the neighbor in the jdir direction
c

      call get_mx(mx,jstx,jedx,n,ip,mp)

      jp = -1
      if (mx.le.1) return                     ! This node has no data

      if (jdir.gt.0) then

         if (ip.eq.mp-1) return   ! No Neighbor
         do jp=ip+1,mp-1
            call get_mx(mx,jstx,jedx,n,jp,mp)
            if (mx.gt.1) return
         enddo
         jp = -1

      else

         if (ip.eq.0) return   ! No Neighbor
         do jp=ip-1,0,-1
            call get_mx(mx,jstx,jedx,n,jp,mp)
            if (mx.gt.1) return
         enddo
         jp = -1

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine get_mtag(mtags,mtagr,gid,nid)
      integer gid

c     Assigns destination nid/mtag for xchanges where nid = (ip,iq,ir)
 
      mtags = nid       ! Send 
      mtagr = gid       ! Receive
      if (mtagr.lt.0) mtags = -1

      return
      end
c-----------------------------------------------------------------------
      subroutine get_gid(gid,ip,iq,ir,mp,mq,mr)
      include 'mpif.h'
      integer gid
c     ===========  Assign global numbering =======
      mpq = mp*mq
      
      gid = MPI_PROC_NULL !-1
      if (0.le.ip.and.ip.lt.mp.and.    ! Inside Box
     &    0.le.iq.and.iq.lt.mq.and.
     &    0.le.ir.and.ir.lt.mr      )  gid = ip + iq*mp + ir*mpq

      return
      end
c-----------------------------------------------------------------------
      subroutine dmesh(mm,n)
      include 'MGRID'
      parameter(ml=2*lpow)
      integer tflag(ml),dmflag(ml)
      common /cflag/ tflag,dmflag
      integer lmx(10,ml)
      save    lmx
c
c     This routine calculates mx,my,mz and their global start indx
c  
c     Partitions nx points across mp procs
c     mx is number of element per proc (0---mx)

#ifdef conservative_mpi_dtype_recreate
      imx = mx
      imy = my
      imz = mz
#endif      

      nx   = n ! Total number of active points in x-dir
      ny   = n !                                  y-dir
      nz   = n !                                  z-dir

      if (dmflag(lev).eq.0) then

        call get_mx(mx,istx,iedx,nx,ip,mp)
        call get_mx(my,isty,iedy,ny,iq,mq)
        call get_mx(mz,istz,iedz,nz,ir,mr)

        if (ldim.eq.2) then
           mz   = 0
           istz = 0
           nz   = 0
        endif

        mm = (mx+1)*(my+1)*(mz+1) ! Per-proc number of points, inc. ghosts

        dmflag(lev) = 1

        lmx(1,lev)  = mx
        lmx(2,lev)  = istx
        lmx(3,lev)  = iedx
        lmx(4,lev)  = my
        lmx(5,lev)  = isty
        lmx(6,lev)  = iedy
        lmx(7,lev)  = mz
        lmx(8,lev)  = istz
        lmx(9,lev)  = iedz
        lmx(10,lev) = mm

      else

        mx   = lmx(1,lev)
        istx = lmx(2,lev)
        iedx = lmx(3,lev)
        my   = lmx(4,lev)
        isty = lmx(5,lev)
        iedy = lmx(6,lev)
        mz   = lmx(7,lev)
        istz = lmx(8,lev)
        iedz = lmx(9,lev)
        mm   = lmx(10,lev)

      endif


      call mtag  ! Recompute message tags for case mx < mp, etc


#ifdef conservative_mpi_dtype_recreate
      if ((imx.ne.mx).OR.(imy.ne.my).OR.(imz.ne.mz)) then
        !if (nid.eq.1) write(6,*) nid, ": imx, imy, imz = ", imx, imy, imz
        !if (nid.eq.1) write(6,*) nid, ": mx, my, mz = ", mx, my, mz
#endif
        call erase_3d_types()
#ifdef enable_FBtype
        call make_3d_types(mx+1, my+1, mz+1, mgreal, .TRUE.)  ! makes MPI datatypes for xchange!
#else
        call make_3d_types(mx+1, my+1, mz+1, mgreal, .FALSE.) ! makes MPI datatypes for xchange!
#endif
#ifdef conservative_mpi_dtype_recreate
      endif
#endif

c     call out_mgrid

      return
      end
c-----------------------------------------------------------------------
      subroutine get_mx(mx,istx,iedx,nx,ip,mp)


      istx = nx*ip/mp
      iedx = nx*(ip+1)/mp + 1

      if (ip.eq.mp-1) iedx = iedx-1

      mx = iedx-istx

      return
      end
c-----------------------------------------------------------------------
      subroutine initu(u,n,h,ndim)  ! Fill grid
      real u(1) !Why?

      if (ndim.eq.2) call initu_2d(u,n,h)
      if (ndim.eq.3) call initu_3d(u,n,h)

      call xchange0  (u,XCHANGE_JGX)

      return
      end
c-----------------------------------------------------------------------
      subroutine initu_2d(u,n,h)  ! Fill grid
      include 'MGRID'
      real u(0:mx,0:my)

      dx = h
      dy = dx

      mm = (mx+1)*(my+1)
      call rzero(u,mm)

      one = 1.
      pi  = 4.*atan(one)
      pik = 5*pi
      pij = 15*pi

c     write(6,*) nid,istx,iedx,mx,' istx'
c     write(6,*) nid,isty,iedy,my,' isty'

      do j=1,my-1
      do i=1,mx-1
         x = (istx+i)*dx
         y = (isty+j)*dy
         u(i,j) = sin(pi*x)*sin(pi*y)
     &            + .25*sin(pik*x)*sin(pik*y)
     &            + .15*sin(pij*x)*sin(pij*y)
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine initu_3d(u,n,h)  ! Fill grid
      include 'MGRID'
      real u(0:mx,0:my,0:mz)

      dx = h
      dy = dx
      dz = dx

      mm = (mx+1)*(my+1)*(mz+1)
      call rzero(u,mm)

      one = 1.
      pi  = 4.*atan(one)
      pik = 5*pi
      pij = 3*pi

      do k=1,mz-1
      do j=1,my-1
      do i=1,mx-1
         x = (istx+i)*dx
         y = (isty+j)*dy
         z = (istz+k)*dz
         u(i,j,k) = sin(pi*x)*sin(pi*y)*sin(pi*z)
     &            + .25*sin(pik*x)*sin(pik*y)
     &            + .15*sin(pij*x)*sin(pij*z)
     &            + .15*sin(pij*z)*sin(pij*y)
         u(i,j,k) = u(i,j,k) + x + x*y + x*y*z
c        write(6,1) i,j,k,x,y,z,u(i,j,k)
c 1      format(3i5,4f12.5,' uinit')
      enddo
      enddo
      enddo
c     call exitt

      return
      end
c-----------------------------------------------------------------------
      subroutine compute_r (r,u,ro,mx,my,mz,h,ndim)
      !real r(1),u(1),ro(1) !Why?

      if (ndim.eq.2) call compute_r2d (r,u,ro,mx,my,h)
      if (ndim.eq.3) call compute_r3d (r,u,ro,mx,my,mz,h)

      return
      end
c-----------------------------------------------------------------------
      subroutine compute_r2d (r,u,ro,mx,my,h)
      !real r(0:mx,0:my),u(0:mx,0:my),ro(0:mx,0:my)
      real, dimension(0:mx,0:my) :: r, u, ro

      h2     = h*h
      h2i    = 1./h2

      do j=1,my-1
      do i=1,mx-1
         r(i,j)=h2i
     &       *(h2*ro(i,j)+u(i-1,j)+u(i,j-1)+u(i+1,j)+u(i,j+1)-4*u(i,j))
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine compute_r3d (r,u,ro,mx,my,mz,h)
      !real r(0:mx,0:my,0:mz)
      !real u(0:mx,0:my,0:mz),ro(0:mx,0:my,0:mz)
      real, dimension(0:mx,0:my,0:mz) :: r, u, ro

      common /ccflops/ cflops,pflops,rflops,rrlops

      h2     = h*h
      h2i    = 1./h2
      

      do k=1,mz-1
      do j=1,my-1
      do i=1,mx-1
         r(i,j,k)=h2i*(h2*ro(i,j,k)-6*u(i,j,k)
     &       + u(i-1,j,k) + u(i+1,j,k)
     &       + u(i,j-1,k) + u(i,j+1,k) 
     &       + u(i,j,k-1) + u(i,j,k+1))
      enddo
      enddo
      enddo
      rrlops = rrlops + (mx-1)*(my-1)*(mz-1)*10

      return
      end
c-----------------------------------------------------------------------
      subroutine restrict_i(rc,rf,mif,mjf,mkf,ie,je,ke,ndim)
      real rc(1),rf(1)

      if (ndim.eq.2) call restrict_2di(rc,rf,mif,mjf,ie,je)
      if (ndim.eq.3) call restrict_3di(rc,rf,mif,mjf,mkf,ie,je,ke)

      return
      end
c-----------------------------------------------------------------------
      subroutine restrict_2di(rc,rf,mif,mjf,ievn,jevn)
      include 'MGRID'
      real rc(0:mx,0:1),rf(0:mif,0:1)


      j=jevn           ! fine index         ! INTERIOR ONLY
      do jc=1,my-1     ! coarse index

         jm1=j-1
         jp1=j+1

         i=ievn        ! fine index
         do ic=1,mx-1  ! coarse

            im1=i-1
            ip1=i+1

            rc(ic,jc) = .25*rf(i,j)
     &        +(rf(im1,j)+rf(ip1,j)+rf(i,jm1)+rf(i,jp1))*.125
     &        +(rf(im1,jm1)+rf(ip1,jm1)+rf(im1,jp1)+rf(ip1,jp1))*.0625

            i=i+2
         enddo
         j=j+2
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine restrict_3di(rc,rf,mif,mjf,mkf,ievn,jevn,kevn)
      include 'MGRID'
      real rc(0:mx,0:my,0:1),rf(0:mif,0:mjf,0:1)

      common /ccflops/ cflops,pflops,rflops,rrlops

      k=kevn             ! fine index         ! INTERIOR ONLY
      do kc=1,mz-1       ! coarse index

        km1=k-1
        kp1=k+1

        j=jevn           ! fine index         ! INTERIOR ONLY
        do jc=1,my-1     ! coarse index

          jm1=j-1
          jp1=j+1

          i=ievn         ! fine index
          do ic=1,mx-1   ! coarse

            im1=i-1
            ip1=i+1
              
            rc(ic,jc,kc) = .125*rf(i,j,k)
     &      + .0625 * ( rf(im1,j,k)+rf(i,jm1,k)+rf(i,j,km1)
     &                + rf(ip1,j,k)+rf(i,jp1,k)+rf(i,j,kp1) ) ! +8 flops

     &      + .03125* ( rf(i,  jm1,km1)+rf(i,  jp1,km1)
     &                + rf(i,  jm1,kp1)+rf(i,  jp1,kp1)
     &                + rf(im1,j,  km1)+rf(ip1,j,  km1)
     &                + rf(im1,j,  kp1)+rf(ip1,j,  kp1)
     &                + rf(im1,jm1,k  )+rf(ip1,jm1,k  )
     &                + rf(im1,jp1,k  )+rf(ip1,jp1,k  ) )    ! +12 flops

     &      + .015625*( rf(im1,jm1,km1)+rf(ip1,jm1,km1)
     &                + rf(im1,jp1,km1)+rf(ip1,jp1,km1)
     &                + rf(im1,jm1,kp1)+rf(ip1,jm1,kp1)
     &                + rf(im1,jp1,kp1)+rf(ip1,jp1,kp1) )    ! +8 flops

            i=i+2
          enddo
          j=j+2
        enddo
        k=k+2
      enddo
      rflops = rflops + (mx-1)*(my-1)*(mz-1)*20

      return
      end
c-----------------------------------------------------------------------
      subroutine prolong_a(uf,dmax,uc,mic,mjc,mkc,ievn,jevn,kevn,ndim) 

c     uf = uf + P*uc

      !real uf(1),uc(1) !Why?

      if (ndim.eq.2) then
         call prolong_2da(uf,dmax,uc,mic,mjc,ievn,jevn)
      else
         call prolong_3da(uf,dmax,uc,mic,mjc,mkc,ievn,jevn,kevn)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine prolong_2da(uf,dmax,uc,mic,mjc,ievn,jevn) ! uf = uf + P*uc
      include 'MGRID'
      real uf(0:mx,0:1),uc(0:mic,0:1)
 
      dmax = 0

      j1   = jevn-1
      j2   = jevn
      jc   = 1              ! coarse index         ! INTERIOR ONLY
      jcm1 = jc-1
      do jfine=1,my,2       ! fine index

         i1   = ievn-1
         i2   = ievn
         ic   = 1           ! coarse index
         icm1 = ic-1

         do ifine=1,mx,2  ! fine index

            dmax = max(dmax,abs(uc(ic,jc)))    ! really cheap local estimate
            uf(i2,j2) = uf(i2,j2) + uc(ic,jc)
            uf(i2,j1) = uf(i2,j1) + .5*(uc(ic,jc)+uc(ic,jcm1))
            uf(i1,j2) = uf(i1,j2) + .5*(uc(ic,jc)+uc(icm1,jc))
            uf(i1,j1) = uf(i1,j1) + .25*
     &          (uc(icm1,jcm1)+uc(icm1,jc)+uc(ic,jcm1)+uc(ic,jc))
   
            i1   = i1  +2
            i2   = i2  +2
            ic   = ic  +1
            icm1 = icm1+1

         enddo

         j1   = j1  +2
         j2   = j2  +2
         jc   = jc  +1
         jcm1 = jcm1+1

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine prolong_3da(uf,dmax,uc,mic,mjc,mkc,ievn,jevn,kevn) ! uf = uf + P*uc
      include 'MGRID'
      real uf(0:mx,0:my,0:1),uc(0:mic,0:mjc,0:1)
      common /ccflops/ cflops,pflops,rflops,rrlops

      dcrn_mx = 0  ! Cheap estimate of max variation (local to this proc)
      dctr_mx = 0


      k1   = kevn-1
      k2   = kevn
      kc   = 1                 ! coarse index         ! INTERIOR ONLY
      kcm1 = kc-1
      do kfine=1,mz,2        ! fine index

        j1   = jevn-1
        j2   = jevn
        jc   = 1          
        jcm1 = jc-1
        do jfine=1,my,2      ! fine index

          i1   = ievn-1
          i2   = ievn
          ic   = 1       
          icm1 = ic-1
          do ifine=1,mx,2    ! fine index

            ! Corner  !           
            uf(i2,j2,k2)=uf(i2,j2,k2) + uc(ic,jc,kc)                  
            dcrn_mx = max(dcrn_mx,abs(uc(ic,jc,kc)))

            ! Edges   !
            uf(i1,j2,k2)=uf(i1,j2,k2) +.5*(uc(icm1,jc,kc)+uc(ic,jc,kc))
            uf(i2,j1,k2)=uf(i2,j1,k2) +.5*(uc(ic,jcm1,kc)+uc(ic,jc,kc))
            uf(i2,j2,k1)=uf(i2,j2,k1) +.5*(uc(ic,jc,kcm1)+uc(ic,jc,kc))

            ! Faces   !
            uf(i2,j1,k1)=uf(i2,j1,k1) + 0.25* ( uc(ic,jc,kc) +
     &                  uc(ic,jcm1,kc)+uc(ic,jc,kcm1)+uc(ic,jcm1,kcm1)) 
            uf(i1,j2,k1) = uf(i1,j2,k1) + 0.25* ( uc(ic,jc,kc) +
     &                  uc(icm1,jc,kc)+uc(ic,jc,kcm1)+uc(icm1,jc,kcm1))
            uf(i1,j1,k2) = uf(i1,j1,k2) + 0.25* ( uc(ic,jc,kc) +
     &                  uc(icm1,jc,kc)+uc(ic,jcm1,kc)+uc(icm1,jcm1,kc))

            ! Center !
            dctr = .125* ( uc(icm1,jcm1,kcm1) + uc(ic  ,jcm1,kcm1) +
     &                     uc(icm1,jc  ,kcm1) + uc(ic  ,jc  ,kcm1) +
     &                     uc(icm1,jcm1,kc  ) + uc(ic  ,jcm1,kc  ) +
     &                     uc(icm1,jc  ,kc  ) + uc(ic  ,jc  ,kc  ) )  

            uf(i1,j1,k1) = uf(i1,j1,k1) + dctr

            dctr_mx = max(dctr_mx,abs(dctr))

            i1   = i1  +2
            i2   = i2  +2
            ic   = ic  +1
            icm1 = icm1+1


          enddo
          j1   = j1  +2
          j2   = j2  +2
          jc   = jc  +1
          jcm1 = jcm1+1
        enddo
        k1   = k1  +2
        k2   = k2  +2
        kc   = kc  +1
        kcm1 = kcm1+1
      enddo
      dmax = max(dctr_mx,dcrn_mx)

      pflops = pflops + 34*(mx*my*mz)/8

      return
      end
c-----------------------------------------------------------------------
      subroutine smooth_m (u,r,mx,my,mz,h,sigma,m,uo,iz,igs,ndim)! m sweep
      real u(1),r(1),uo(1)

      if (ndim.eq.2) then
        if (igs.eq.0) then
          call smooth_m_jac2d(u,r,mx,my,h,sigma,m,uo,iz)
        elseif (igs.eq.1) then
          call smooth_m_gs_2d(u,r,mx,my,h,sigma,m,iz)
        else!if (nid.eq.1) then
          write(6, *) "no smooth_m func ",
     &      "corresponding to igs = ", igs, "!"
          call exitt
        endif
      else
        if (igs.eq.0) then
          call smooth_m_jac3d(u,r,mx,my,mz,h,sigma,m,uo,iz)
        elseif (igs.eq.1) then
          call smooth_m_gs_3d(u,r,mx,my,mz,h,sigma,m,iz)
        elseif (igs.eq.2) then
          call smooth_m_jac3alt(u,r,mx,my,mz,h,sigma,m,uo,iz)
        elseif (igs.eq.4) then
          call smooth_m_jac3cnd(u,r,mx,my,mz,h,sigma,m,uo,iz)
        elseif (igs.eq.3) then
          call smooth_m_jac3flg(u,r,mx,my,mz,h,sigma,m,uo,iz)
        elseif (igs.eq.5) then
          call smooth_m_jac3idx(u,r,mx,my,mz,h,sigma,m,uo,iz)
       ! elseif (igs.eq.6) then
       !   call smooth_m_jac3swp(u,r,mx,my,mz,h,sigma,m,uo,iz)
        else!if (nid.eq.1) then
          write(6, *) "no smooth_m func ",
     &      "corresponding to igs = ", igs, "!"
          call exitt
        endif
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine smooth_m_jac2d (u,r,mx,my,h,sigma,m,o,izero) !  Jacobi

      !real u(0:mx,0:my),r(0:mx,0:my),o(0:mx,0:my)
      real, dimension(0:mx,0:my) :: u, r, o

      h2       = h*h
      sa       = sigma/4

      if (izero.eq.0) then
        do j=1,my-1
        do i=1,mx-1
          u(i,j)=sa*h2*r(i,j)
        enddo
        enddo
      endif

      mm1 = (my+1)*(mx+1)
      do i_sweep = 2,m+izero
        call copy     (o,u,mm1) ! Old copy, w/ boundary data
        call xchange0 (o,XCHANGE_JGX)

        do j=1,my-1
        do i=1,mx-1
           u(i,j)=o(i,j) + sa
     &         *(h2*r(i,j)+o(i-1,j)+o(i,j-1)+o(i+1,j)+o(i,j+1)-4*o(i,j))
        enddo
        enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
c     ======= Calculation is independent of u =======
      subroutine jac3d_helper (ix,iy,iz,u,o,r,h2,sa,jgx) ! Jacobi helper
      !include 'MGRID'
      real, dimension(0:ix,0:iy,0:iz) :: u, r, o
      
      call xchange0(o,jgx) ! tile prev u

      do k=1,iz-1
      do j=1,iy-1
      do i=1,ix-1
        u(i,j,k) = o(i,j,k) + sa*( h2*r(i,j,k)-6*o(i,j,k)
     &            + o(i-1,j,k) + o(i+1,j,k)
     &            + o(i,j-1,k) + o(i,j+1,k)
     &            + o(i,j,k-1) + o(i,j,k+1) ) ! calculate new u
      enddo
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
c     ======= In this original version, u never contains any boundary data! =======
      subroutine smooth_m_jac3d (u,r,ix,iy,iz,h,sigma,m,o,izero) !  Jacobi
      !include 'MGRID'
      common /ccflops/ cflops,pflops,rflops,rrlops

      !real u(0:ix,0:iy,0:iz),r(0:ix,0:iy,0:iz),o(0:ix,0:iy,0:iz)
      real, dimension(0:ix,0:iy,0:iz) :: u, r, o

	    !if (nid.eq.6) write(6,*) nid,": ix,iy,iz=", ix, iy, iz

      h2       = h*h
      sa       = sigma/6

      if (izero.eq.0) then
        do k=1,iz-1
        do j=1,iy-1
        do i=1,ix-1
          u(i,j,k)=sa*h2*r(i,j,k) ! calculates u as if o were initially all 0's
        enddo
        enddo
        enddo
        cflops = cflops + 2*(ix-1)*(iy-1)*(iz-1)
      endif

      !write(6,*) nid, " m, izero, ix, iy, iz = ", m, izero, ix, iy, iz

      !iters = 0
      mm1 = (ix+1)*(iy+1)*(iz+1)
      do i_sweep = 2,m+izero
        call copy(o,u,mm1) ! Old copy, w/ old (EMPTY!) boundary data (o gets prev "u")

        !if (nid.eq.6) then
          !write(6,*) nid, ": u = "
          !call print_3d(u, 11, 2, 2)
          !write(6,*) nid, ": o = "
          !call print_3d(o, 11, 2, 2)
        !endif

        call jac3d_helper(ix,iy,iz,u,o,r,h2,sa,JAC3D_JGX) ! u gets new "u"

        !if (nid.eq.6) then
          !write(6,*) nid, ": xchanged o = "
          !call print_3d(o, 11, 2, 2)
          !write(6,*) nid, ": new u = "
          !call print_3d(u, 11, 2, 2)
        !endif

        !iters = iters + 1
      enddo
      cflops = cflops + (11*(ix-1)*(iy-1)*(iz-1))*(m+izero-2+1)

      !write(6,*) "alt - iters: ", iters, " m+izero-1: ", m+izero-1

      return
      end
c-----------------------------------------------------------------------
      subroutine smooth_m_jac3alt (u,r,ix,iy,iz,h,sigma,m,o,izero) !  Jacobi
      real, dimension(0:ix,0:iy,0:iz) :: u, r, o

      common /ccflops/ cflops,pflops,rflops,rrlops

      h2       = h*h
      sa       = sigma/6

      if (izero.eq.0) then
        do k=1,iz-1
        do j=1,iy-1
        do i=1,ix-1
          u(i,j,k)=sa*h2*r(i,j,k)
        enddo
        enddo
        enddo
        cflops = cflops + 2*(ix-1)*(iy-1)*(iz-1)
      endif

      !iters = 0
      i_end = (m+izero-2+1)/2+1
      do i_sweep = 2,i_end
        call jac3d_helper(ix,iy,iz,o,u,r,h2,sa,XCHANGE_JGX) ! prev "u" = u, o will get new "u"
        call jac3d_helper(ix,iy,iz,u,o,r,h2,sa,XCHANGE_JGX) ! prev "u" = o, u will get new "u"
        !iters = iters+2
      enddo

      if (MOD(m+izero-2+1,2).ne.0) then
        call jac3d_helper(ix,iy,iz,o,u,r,h2,sa,XCHANGE_JGX) ! prev "u" = u, o will get new "u"
        mm1 = (ix+1)*(iy+1)*(iz+1)
        call copy (u,o,mm1)                     ! u needs newest "u", which is in o
        !iters = iters+1
      endif

      cflops = cflops + (11*(ix-1)*(iy-1)*(iz-1))*(m+izero-2+1)

      !write(6,*) "alt - iters: ", iters, " m+izero-1: ", m+izero-1

      return
      end
c-----------------------------------------------------------------------
      subroutine smooth_m_jac3flg (u,r,ix,iy,iz,h,sigma,m,o,izero) !  Jacobi
      real, dimension(0:ix,0:iy,0:iz) :: u, r, o
      LOGICAL :: f

      common /ccflops/ cflops,pflops,rflops,rrlops

      h2       = h*h
      sa       = sigma/6

      if (izero.eq.0) then
        do k=1,iz-1
        do j=1,iy-1
        do i=1,ix-1
          u(i,j,k)=sa*h2*r(i,j,k)
        enddo
        enddo
        enddo
        cflops = cflops + 2*(ix-1)*(iy-1)*(iz-1)
      endif

      f = .TRUE. ! means u was last updated
      do i_sweep = 2,m+izero
        if (f) then
          call jac3d_helper(ix,iy,iz,o,u,r,h2,sa,XCHANGE_JGX)
        else
          call jac3d_helper(ix,iy,iz,u,o,r,h2,sa,XCHANGE_JGX)
        end if
        f = .NOT.f
      enddo
      cflops = cflops + (11*(ix-1)*(iy-1)*(iz-1))*(m+izero-2+1)

      if (.NOT.f) then
        mm1 = (ix+1)*(iy+1)*(iz+1)
        call copy (u,o,mm1) ! o was last updated
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine smooth_m_jac3cnd (u,r,ix,iy,iz,h,sigma,m,o,izero) !  Jacobi
      real, dimension(0:ix,0:iy,0:iz) :: u, r, o

      common /ccflops/ cflops,pflops,rflops,rrlops

      h2       = h*h
      sa       = sigma/6

      if (izero.eq.0) then
        do k=1,iz-1
        do j=1,iy-1
        do i=1,ix-1
          u(i,j,k)=sa*h2*r(i,j,k)
        enddo
        enddo
        enddo
        cflops = cflops + 2*(ix-1)*(iy-1)*(iz-1)
      endif

      do i_sweep = 2,m+izero
        if (MOD(i_sweep,2).eq.0) then
          call jac3d_helper(ix,iy,iz,o,u,r,h2,sa,XCHANGE_JGX)
        else
          call jac3d_helper(ix,iy,iz,u,o,r,h2,sa,XCHANGE_JGX)
        end if
      enddo
      cflops = cflops + (11*(ix-1)*(iy-1)*(iz-1))*(m+izero-2+1)

      if (MOD(m+izero,2).eq.0) then
        mm1 = (ix+1)*(iy+1)*(iz+1)
        call copy (u,o,mm1) ! o was updated last
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine smooth_m_jac3idx (u,r,ix,iy,iz,h,sigma,m,o,izero) !  Jacobi
      common /ccflops/ cflops,pflops,rflops,rrlops

      ! Requires explicit, constant bounds... inviable
      !type realmat
      !  real, dimension(:,:,:) :: m
      !end type realmat
      
      !type realptr
      !  type(realmat), pointer :: p
      !end type realptr

      real, dimension(0:ix,0:iy,0:iz) :: u, o, r
      !real, dimension(0:ix,0:iy,0:iz), target :: u, o
      !real r(0:ix,0:iy,0:iz)

      !type(realptr), dimension(0:1) :: uo
      !uo(0)%p => u; uo(1)%p = o
      real, dimension(0:1,0:ix,0:iy,0:iz) :: uo
      uo(0,:,:,:) = u; uo(1,:,:,:) = o

      h2       = h*h
      sa       = sigma/6

      if (izero.eq.0) then
        do k=1,iz-1
        do j=1,iy-1
        do i=1,ix-1
          uo(0,i,j,k)=sa*h2*r(i,j,k)
        enddo
        enddo
        enddo
        cflops = cflops + 2*(ix-1)*(iy-1)*(iz-1)
      endif

      io = 0
      !mm1 = (ix+1)*(iy+1)*(iz+1)
      do i_sweep = 2,m+izero

        !call copy    (o,u,mm1) ! Old copy, w/ boundary data
        call xchange0(uo(io,:,:,:),XCHANGE_JGX)

        do k=1,iz-1
        do j=1,iy-1
        do i=1,ix-1
          uo(1-io,i,j,k) = uo(io,i,j,k) + sa*(
     &              h2*r(i,j,k)-6*uo(io,i,j,k)
     &              + uo(io,i-1,j,k) + uo(io,i+1,j,k)
     &              + uo(io,i,j-1,k) + uo(io,i,j+1,k)
     &              + uo(io,i,j,k-1) + uo(io,i,j,k+1) )
        enddo
        enddo
        enddo
        
        io = 1 - io
      enddo
      cflops = cflops + (11*(ix-1)*(iy-1)*(iz-1))*(m+izero-2+1)

      u = uo(io,:,:,:)
      o = uo(1-io,:,:,:)
      return
      end
c-----------------------------------------------------------------------
      subroutine smooth_m_gs_2d (u,r,mx,my,h,sigma,m,izero) !  GS
      !real u(0:mx,0:my),r(0:mx,0:my)
      real, dimension(0:mx,0:my) :: u, r

      call zero_bdry(u,mx,my,1,2)

      h2       = h*h
      sa       = sigma/4

      if (izero.eq.0) then
         do j=1,my-1
         do i=1,mx-1
            u(i,j)=sa*h2*r(i,j)
         enddo
         enddo
      endif

      do i_sweep = 2,m+izero
        do j=1,my-1
        do i=1,mx-1
           u   (i,j)=u(i,j) + sa
     &         *(h2*r(i,j)+u(i-1,j)+u(i,j-1)+u(i+1,j)+u(i,j+1)-4*u(i,j))
        enddo
        enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine smooth_m_gs_3d (u,r,mx,my,mz,h,sigma,m,izero) !  GS
      !real u(0:mx,0:my,0:mz),r(0:mx,0:my,0:mz)
      real, dimension(0:mx,0:my,0:mz) :: u, r

      h2       = h*h
      sa       = sigma/6

      s        = h2*sigma/6
      sa2      = 2*s - 6*s*s/h2
      sh2      = s/h2
      s22      = s*s/h2


      if (izero.eq.0) then
         do k=1,mz-1
         do j=1,my-1
         do i=1,mx-1
           u(i,j,k) =  sa2*r(i,j,k)
     &               + sh2*(u(i-1,j,k) + u(i,j-1,k) + u(i,j,k-1))
     &               + s22*(r(i+1,j,k) + r(i,j+1,k) + r(i,j,k+1))
         enddo
         enddo
         enddo
      endif

      do i_sweep = 3,m+2*izero
        do k=1,mz-1
        do j=1,my-1
        do i=1,mx-1
           u(i,j,k) = u(i,j,k) + sa*( h2*r(i,j,k)-6*u(i,j,k)
     &              + u(i-1,j,k) + u(i+1,j,k)
     &              + u(i,j-1,k) + u(i,j+1,k)
     &              + u(i,j,k-1) + u(i,j,k+1) )
        enddo
        enddo
        enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine outm(a,kx,n,m,name4)
      include 'MGRID'

      real a(0:mx,0:my,0:mz)
      character*4 name4
      integer p

      my10 = min(my,18)

      do p=np-1,0,-1
         call gsync()
         if (p.eq.nid) then
            open(unit=12,file='x.x') 
            write(12,*) 'open ',p
            close(12)
         endif
         call gsync()

         if (p.eq.nid) then
            write(6,2) p,name4,n,m,mx,my,mz,nx,ny,nz,kx
  2         format(i2,'oo ',a4,10i4)
            do k=mz,0,-1
               write(6,1) p,kx,k,(a(kx,j,k),j=0,my10)
  1            format(i2,'oo ',2i2,19f9.4)
            enddo
         endif

         call gsync()
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine zero_bdry(a,mx,my,mz,ndim)
      real a(1)

      if (ndim.eq.2) call zero_bdry2(a)
      if (ndim.eq.3) call zero_bdry3(a)

      return
      end
c-----------------------------------------------------------------------
      subroutine zero_bdry2(a)
      include 'MGRID'
      real a(0:mx,0:my)
c     zero out only the outer edge of mesh

      if (ip.eq.0) then
         do j=0,my
            a(0,j)  = 0
         enddo
      endif
      if (ip.eq.mp) then
         do j=0,my
            a(mx,j)  = 0
         enddo
      endif
      if (iq.eq.0) then
         do i=0,mx
            a(i,0)  = 0
         enddo
      endif
      if (iq.eq.mq) then
         do i=0,mx
            a(i,my)  = 0
         enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine zero_bdry3(a)
      include 'MGRID'
      real a(0:mx,0:my,0:mz)

      if (ip.eq.0) then
         do k=0,mz
         do j=0,my
            a(0,j,k)  = 0
         enddo
         enddo
      endif
      if (ip.eq.mp) then
         do k=0,mz
         do j=0,my
           a(mx,j,k)  = 0
         enddo
         enddo
      endif
      if (iq.eq.0) then
         do k = 0,mz
         do i = 0,mx
            a(i,0,k)  = 0
         enddo
         enddo
      endif
      if (iq.eq.mq) then
         do k = 0,mz
         do i = 0,mx
            a(i,my,k)  = 0
         enddo
         enddo
      endif
      if (ir.eq.0) then
         do j = 0,my
         do i = 0,mx
            a(i,j,0)  = 0
         enddo
         enddo
      endif
      if (ir.eq.mr) then
         do j = 0,my
         do i = 0,mx
            a(i,j,mz)  = 0
         enddo
         enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine print_3d(u, ix, iy, iz)
      include 'MGRID'
      real u(0:mx,0:my,0:mz)

      do k=0,iz
        write(6, *) "k = ", k
        do j=0,iy
          write(6, '(1000F16.4)') real(u(0:ix, j, k))
          !write(6,*) real(u(:mx-1, j, k))
        enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
#ifdef verify_jacobi_alts
      subroutine test_alts(ukf,rkf,h,sigma,nsmooth,rf,iz,minIg,maxIg)
      USE Comm_Funcs, ONLY : mateq
      include 'MGRID'
      real ukf(1), rkf(1), rf(1)
      real, dimension(0:mx,0:my,0:mz) :: rfc, v, w

      !if (nid.eq.1) write(6,*) nid, ": mx, my, mz = ", mx, my, mz

      mm1 = (mx+1)*(my+1)
      if (ldim.eq.3) then
        mm1 = mm1*(mz+1)
      endif

      !call smooth_m(ukf,rkf,mx,my,mz,h,sigma,nsmooth,rf,iz,minIg,ldim)

      !if (nid.eq.6) then
        !write(6,*) nid, ": original u:"
        !call print_3d(ukf, 11, 2, 2)
      !endif

      do nsmoth = nsmooth, nsmooth+1
        call copy(rfc,rf,mm1)
        call copy(v,ukf,mm1)
        !call gsync
        call smooth_m(v,rkf,mx,my,mz,h,sigma,nsmoth,rfc,iz,0,ldim)
        if (nid.eq.6) then
          !write(6,*) nid, ": v = "
          !call print_3d(v, 11, 2, 2)
        endif
        do ig = minIg, maxIg
          call copy(rfc,rf,mm1)
          call copy(w,ukf,mm1)
          !call gsync
          call smooth_m(w,rkf,mx,my,mz,h,sigma,nsmoth,rfc,iz,ig,ldim)

          if (nid.eq.6) then
            if (mateq(v,w,0)) then ! the 0 excludes the faces from the comparison
              !write(6,*) nid, ": igs = ", ig, " matches result of "
     &       !   , "igs = 0 for nsmooth = ", nsmoth, "!"
            else
              write(6,*) nid, ": igs = ", ig, " gave a different"
     &          , " result than igs = 0 for nsmooth = ", nsmoth
     &          , " (not including faces!):"
              write(6,*) nid, ": result of igs = 0:"
              call print_3d(v(1:mx-1,1:my-1,1:mz-1), 5, 2, 2)
              write(6,*) nid, ": result of igs = ", ig, ":"
              call print_3d(w(1:mx-1,1:my-1,1:mz-1), 5, 2, 2)

#ifndef allow_math_differences
              stop
#endif
            endif
          endif
        enddo
        !if (nid.eq.6) write(6,*) " "
      enddo
      end
#endif
c-----------------------------------------------------------------------
#ifdef verify_xchange
      subroutine test_xchange(u,minJg,maxJg)
      USE Comm_Funcs, ONLY : mateq

      include 'MGRID'
      real u(1)
      real, dimension(0:mx,0:my,0:mz) :: v, w

      !if (nid.eq.1) write(6,*) nid, ": mx, my, mz = ", mx, my, mz

      mm1 = (mx+1)*(my+1)
      if (ldim.eq.3) then
        mm1 = mm1*(mz+1)
      endif

      !call smooth_m(ukf,rkf,mx,my,mz,h,sigma,nsmooth,rf,iz,minIg,ldim)

      if (nid.eq.6) then
        !write(6,*) nid, ": original u:"
        !call print_3d(u, 11, 2, 2)
      endif

      call copy(v,u,mm1)
      !call gsync
      call xchange0(v,0)
      !if (nid.eq.6) then
        !write(6,*) nid, ": v = "
        !call print_3d(v, 11, 2, 2)
      !endif
      do jg = minJg, maxJg
        call copy(w,u,mm1)
        !call gsync
        call xchange0(w,jg)

        if ((nid.eq.6)) then!.or.(nid.eq.7)) then
          if (mateq(v,w,1)) then ! the 1 means we only want to compare elements that are at the boundary in at most 1 dimensions (e.g. cube faces, but no edges)
            !write(6,*) nid,": jgx = ",jg," matches result of jgx = 0!"
          else
            write(6,*) nid, ": jgx = ", jg, " gave a different"
     &        , " result than jgx = 0 (not including edges!):"
            write(6,*) nid, ": result of jgx = 0:"
            call print_3d(v, 5, 2, 2)
            write(6,*) nid, ": result of jgx =", jg, ":"
            call print_3d(w, 5, 2, 2)

            stop

          endif
        endif
      enddo
      end
#endif
c-----------------------------------------------------------------------
      subroutine xchange0(u, jgx) ! Zeros out boundary data
      include 'MGRID'
      parameter (lb=(lx+2)**(ldim-1))
      common /cbuff/ b1(lb),b2(lb),bo(lb)

      real u(1) ! essentially "real *u" from c, allows typed passing

      if (ldim.eq.2) then
        call xchange_2d(u,b1,b2,bo)
      else
        if (jgx.eq.0) then
          call xchange_3d_old(u,b1,b2,bo)
        elseif (jgx.eq.1) then
          call xchange_3d(u)
        elseif (jgx.eq.2) then
          call xchange_3d_async(u)
        elseif (nid.eq.1) then
          write(6, *) nid, ": no xchange0 func ",
     &      "corresponding to jgx = ", jgx, "!"
          call exitt
        endif
      endif
 
      return
      end
c-----------------------------------------------------------------------
      subroutine xchange_3d_async(u)
      include 'MGRID'
      real u(0:mx,0:my,0:mz)
      integer reqs(12)
      !do i = 1,12
      !  reqs(i) = MPI_Request_Null
      !enddo

      !write(6,*) nid, ":", nidl, nidr, nidd, nidu, nidb, nidf
      !call gsync

      !        ^ y
      !        |
      !        +---> x
      !      z/
c     ======   Exchange x-faces   =======
c     ======    Recv/Send Left    =======
      reqs(1) = irecv2(mtagrl,u(0,0,0)   ,LRtype,1,nidl)      ! mtags from MGRID
      reqs(4) = isend1(mtagsl,u(1,0,0)   ,LRtype,1,nidl)      ! nidl from MGRID, null if mtag < 0
c     ======    Recv/Send Right   =======
      reqs(3) = irecv2(mtagrr,u(mx  ,0,0),LRtype,1,nidr)      ! no-op leaving 0's if mtag < 0
      reqs(2) = isend1(mtagsr,u(mx-1,0,0),LRtype,1,nidr)      ! nidr from MGRID
      !write(6,*) nid, ": LR done"
      !call msgwaitall(reqs,4)

c     ======   Exchange y-faces   =======
c     ======    Recv/Send Down    =======
      reqs(5) = irecv2(mtagrd,u(0,0,0)   ,UDtype,1,nidd)      ! mtags from MGRID
      reqs(8) = isend1(mtagsd,u(0,1,0)   ,UDtype,1,nidd)      ! nidd from MGRID
c     ======     Recv/Send Up     =======
      reqs(7) = irecv2(mtagru,u(0,my  ,0),UDtype,1,nidu)
      reqs(6) = isend1(mtagsu,u(0,my-1,0),UDtype,1,nidu)      ! nidu from MGRID

      !call msgwaitall(reqs(5:8),4)

c     ======   Exchange z-faces, no buffers (alignment already matches)   =======
#ifdef enable_FBtype
c     ======    Recv/Send Back    =======
      reqs(9 ) = irecv2(mtagrb,u(0,0,0)   ,FBtype,1,nidb)   ! mtags from MGRID
      reqs(12) = isend1(mtagsb,u(0,0,1)   ,FBtype,1,nidb)   ! nidb from MGRID
c     ======    Recv/Send Front   =======
      reqs(11) = irecv2(mtagrf,u(0,0,mz  ),FBtype,1,nidf)
      reqs(10) = isend1(mtagsf,u(0,0,mz-1),FBtype,1,nidf)   ! nidf from MGRID
#else
      isz = (mx+1)*(my+1)
c     ======    Recv/Send Back    =======
      reqs(9 ) = irecv2(mtagrb,u(0,0,0)   ,mgreal,isz,nidb) ! mtags from MGRID
      reqs(12) = isend1(mtagsb,u(0,0,1)   ,mgreal,isz,nidb) ! nidb from MGRID
c     ======    Recv/Send Front   =======
      reqs(11) = irecv2(mtagrf,u(0,0,mz  ),mgreal,isz,nidf)
      reqs(10) = isend1(mtagsf,u(0,0,mz-1),mgreal,isz,nidf) ! nidf from MGRID
#endif

c     ======       Wait All       =======
      !call msgwaitall(reqs(9:12),4)
      call msgwaitall(reqs,12)
      !do i = 1,12
      !  if (reqs(i).ne.MPI_REQUEST_NULL) call msgwait1(reqs(i))
      !enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine xchange_3d(u)!,buff1,buff2,buffo)
      include 'MGRID'
!     Zeros out boundary data !
      !USE Comm_Funcs, ONLY: irecv0

      real u(0:mx,0:my,0:mz)!,buff1(0:1),buff2(0:1),buffo(0:1)

      !if (nid.eq.1) write(6,*) nid, ": ", LRtype, UDtype
      !do
      !enddo

c     ======   Exchange x-faces   =======
      !myz  = (my+1)*(mz+1)
      !myz1 = myz-1

c     ======       Recv Left      =======
      !len = myz*wdsize                         ! wdsize from MGRID
      !imsg_left  = irecv0(mtagrl,buff1,len)    ! mtags from MGRID
      imsg_left = irecv2(mtagrl,u(0 ,0,0),LRtype,1,nidl)

c     ======      Recv Right      =======
      !imsg_right = irecv0(mtagrr,buff2,len)    ! Set buff=0 if mtag < 0
      imsg_right = irecv2(mtagrr,u(mx,0,0),LRtype,1,nidr)
      !call gsync()                             ! make sure all are posted

c     ======       Send Left      =======
      !len = myz*wdsize                         ! wdsize from MGRID
c      do j = 0,myz1; buffo(j) = u(1   ,j,0);! load left face into outbound buff      enddo
      !buffo(0:myz1) = u(1,0:myz1,0)
      !call csend0(mtagsl,buffo,len,nidl,0) ! nidl from MGRID, null if mtag < 0
      call csend1(mtagsl,u(1  ,0,0),LRtype,1,nidl)

c     ======      Send Right      =======
c      do j = 0,myz1; buffo(j) = u(mx-1,j,0);! load right face into outbound buff      enddo
      !buffo(0:myz1) = u(mx-1,0:myz1,0)
      !call csend0(mtagsr,buffo,len,nidr,0) ! nidr from MGRID
      call csend1(mtagsr,u(mx-1,0,0),LRtype,1,nidr)

c     ======      Wait Right      =======
      call msgwait(imsg_right)
c      do j = 0,myz1; u(mx,j,0)= buff2(j);! unpack inbound right face      enddo
      !u(mx,0:myz1,0) = buff2(0:myz1)

c     ======       Wait Left      =======
      call msgwait(imsg_left)
c      do j = 0,myz1; u(0,j,0)= buff1(j);! unpack inbound left face      enddo
      !u(0,0:myz1,0) = buff1(0:myz1)

      !write(6,*) nid, ": LR done"



c     ======   Exchange y-faces   =======
      !mxz = (mx+1)*(mz+1)

      !len = mxz*wdsize                        ! wdsize from MGRID
c     ======       Recv Down      =======
      !imsg_down  = irecv0(mtagrd,buff1,len)   ! mtags from MGRID
      imsg_down  = irecv2(mtagrd,u(0,0 ,0),UDtype,1,nidd)
c     ======        Recv Up       =======
      !imsg_up    = irecv0(mtagru,buff2,len)
      imsg_up    = irecv2(mtagru,u(0,my,0),UDtype,1,nidu)
      !call gsync()                            ! make sure all are posted

c     ======       Send Down      =======
      !l=0                        !        ^ y
      !do k = 0,mz                !        |
      !do i = 0,mx                !        +---> x
      !   buffo(l) = u(i,1,k)     !      z/
      !   l=l+1                   ! load down face into outbound buff
      !enddo
      !enddo
      !len = l*wdsize
      !call csend0(mtagsd,buffo,len,nidd,0) ! nidd from MGRID
      call csend1(mtagsd,u(0,1,0),UDtype,1,nidd)

c     ======        Send Up       =======
      !l=0                        !        ^ y
      !do k = 0,mz                !        |
      !do i = 0,mx                !        +---> x
      !   buffo(l) = u(i,my-1,k)  !      z/
      !   l=l+1                   ! load up face into outbound buff
      !enddo
      !enddo
      !len = l*wdsize
      !call csend0(mtagsu,buffo,len,nidu,0) ! nidu from MGRID
      call csend1(mtagsu,u(0,my-1,0),UDtype,1,nidu)

c     ======        Wait Up       =======
      call msgwait(imsg_up)
      !l=0                        !        ^ y
      !do k = 0,mz                !        |
      !do i = 0,mx                !        +---> x
      !   u(i,my,k)= buff2(l)     !      z/
      !   l=l+1                   ! unpack onto up face 
      !enddo
      !enddo

c     ======       Wait Down      =======
      call msgwait(imsg_down)
      !l=0                        !        ^ y
      !do k = 0,mz                !        |
      !do i = 0,mx                !        +---> x
      !   u(i,0,k)= buff1(l)      !      z/ 
      !   l=l+1                   ! unpack onto down face 
      !enddo
      !enddo



c     ======   Exchange z-faces, no buffers (alignment already matches)   =======
      isz = (mx+1)*(my+1)!*wdsize                       ! wdsize from MGRID
      !len = size*wdsize

#ifdef enable_FBtype
        imsg_back  = irecv2(mtagrb,u(0,0,0 ),FBtype,1,nidb)   ! mtags from MGRID
        imsg_front = irecv2(mtagrf,u(0,0,mz),FBtype,1,nidf)

        call csend1(mtagsb,u(0,0,   1),FBtype,1,nidb)         ! nidb from MGRID
        call csend1(mtagsf,u(0,0,mz-1),FBtype,1,nidf)         ! nidf from MGRID
#else
        !imsg_back  = irecv0(mtagrb,u(0,0,0 ),len)             ! mtags from MGRID
        !imsg_front = irecv0(mtagrf,u(0,0,mz),len) 
        imsg_back  = irecv2(mtagrb,u(0,0,0 ),mgreal,isz,nidb) ! mtags from MGRID
        imsg_front = irecv2(mtagrf,u(0,0,mz),mgreal,isz,nidf)
        !call gsync()                                          ! make sure all are posted

        !call csend0(mtagsb,u(0,0,   1),len,nidb,0)            ! nidb from MGRID
        !call csend0(mtagsf,u(0,0,mz-1),len,nidf,0)            ! nidf from MGRID
        call csend1(mtagsb,u(0,0,   1),mgreal,isz,nidb)       ! nidb from MGRID
        call csend1(mtagsf,u(0,0,mz-1),mgreal,isz,nidf)       ! nidf from MGRID
#endif

      call msgwait(imsg_front)
      call msgwait(imsg_back)

      return
      end
c-----------------------------------------------------------------------
      subroutine xchange_3d_old(u,buff1,buff2,buffo) 
!     Zeros out boundary data !
      include 'MGRID'

      real u(0:mx,0:my,0:mz),buff1(0:1),buff2(0:1),buffo(0:1)

c     ======    Exchange x-faces  =======

      myz  = (my+1)*(mz+1)
      myz1 = myz-1

      len = myz*wdsize                         ! wdsize from MGRID
      imsg_left  = irecv0(mtagrl,buff1,len)    ! mtags from MGRID
      imsg_right = irecv0(mtagrr,buff2,len)    ! Set buff=0 if mtag < 0
      call gsync()                             ! make sure all are posted

      len = myz*wdsize                         ! wdsize from MGRID
      do j = 0,myz1
         buffo(j) = u(1   ,j,0)     ! load left face into outbound buff
      enddo
      call csend0(mtagsl,buffo,len,nidl,0) ! nidl from MGRID, null if mtag < 0

      do j = 0,myz1
         buffo(j) = u(mx-1,j,0)     ! load right face into outbound buff
      enddo
      call csend0(mtagsr,buffo,len,nidr,0) ! nidr from MGRID

      call msgwait(imsg_right)
      do j = 0,myz1
         u(mx,j,0)= buff2(j)        ! unpack inbound right face
      enddo

      call msgwait(imsg_left)
      do j = 0,myz1
         u(0,j,0)= buff1(j)         ! unpack inbound left face
      enddo


c     ======   Exchange y-faces  =======

      mxz = (mx+1)*(mz+1)

      len = mxz*wdsize                        ! wdsize from MGRID
      imsg_bot   = irecv0(mtagrd,buff1,len)   ! mtags from MGRID
      imsg_top   = irecv0(mtagru,buff2,len)
      call gsync()                            ! make sure all are posted


      l=0                        !        ^ y
      do k = 0,mz                !        |
      do i = 0,mx                !        +---> x
         buffo(l) = u(i,1,k)     !      z/
         l=l+1                   ! load bottom face into outbound buff
      enddo
      enddo
      len = l*wdsize
      call csend0(mtagsd,buffo,len,nidd,0) ! nidd from MGRID

      l=0                        !        ^ y
      do k = 0,mz                !        |
      do i = 0,mx                !        +---> x
         buffo(l) = u(i,my-1,k)  !      z/
         l=l+1                   ! load top face into outbound buff
      enddo
      enddo
      len = l*wdsize
      call csend0(mtagsu,buffo,len,nidu,0) ! nidu from MGRID


      call msgwait(imsg_top)
      l=0                        !        ^ y
      do k = 0,mz                !        |
      do i = 0,mx                !        +---> x
         u(i,my,k)= buff2(l)     !      z/
         l=l+1                   ! unpack onto top face 
      enddo
      enddo

      call msgwait(imsg_bot)
      l=0                        !        ^ y
      do k = 0,mz                !        |
      do i = 0,mx                !        +---> x
         u(i,0,k)= buff1(l)      !      z/ 
         l=l+1                   ! unpack onto bottom face 
      enddo
      enddo

c     ======    Exchange z-faces, no buffers  =======

      len = (mx+1)*(my+1)*wdsize                  ! wdsize from MGRID
      imsg_back  = irecv0(mtagrb,u(0,0,0 ),len)   ! mtags from MGRID
      imsg_front = irecv0(mtagrf,u(0,0,mz),len) 
      call gsync()                                ! make sure all are posted

      len = (mx+1)*(my+1)*wdsize                  ! wdsize from MGRID
      call csend0(mtagsb,u(0,0,   1),len,nidb,0)  ! nidb from MGRID
      call csend0(mtagsf,u(0,0,mz-1),len,nidf,0)  ! nidf from MGRID

      call msgwait(imsg_front)
      call msgwait(imsg_back)

      return
      end
c-----------------------------------------------------------------------
      subroutine xchange_2d(u,buff1,buff2,buffo)
      !USE Comm_Funcs, ONLY: irecv0

      include 'MGRID'

      real u(0:mx,0:my),buff1(0:1),buff2(0:1),buffo(0:1)


c     ======    Exchange x-faces  =======

      myz = (my+1)

      len = myz*wdsize                         ! wdsize from MGRID
      imsg_left  = irecv0(mtagrl,buff1,len)    ! mtags from MGRID
      imsg_right = irecv0(mtagrr,buff2,len)
      call gsync()                            ! make sure all are posted

      len = myz*wdsize
c      do j = 0,my
c         buffo(j) = u(1,j)          ! load left face into outbound buff
c      enddo
      buffo(0:my) = u(1,0:my)
      call csend0(mtagsl,buffo,len,nidl,0) ! nidl from MGRID

c      do j = 0,my
c         buffo(j) = u(mx-1,j)       ! load right face into outbound buff
c      enddo
      buffo(0:my) = u(mx-1,0:my)
      call csend0(mtagsr,buffo,len,nidr,0) ! nidr from MGRID

      call msgwait(imsg_right)
c      do j = 0,my
c         u(mx,j)= buff2(j)          ! unpack inbound right face
c      enddo
      u(mx,0:my) = buff2(0:my)

      call msgwait(imsg_left)
c      do j = 0,my
c         u(0,j)= buff1(j)           ! unpack inbound left face
c      enddo
      u(0,0:my) = buff1(0:my)
      

c     ======    Exchange y-faces, no buffers  =======

      len = (mx+1)*wdsize                     ! wdsize from MGRID
      imsg_bot = irecv0(mtagrd,u(0,0 ),len)    ! mtags from MGRID
      imsg_top = irecv0(mtagru,u(0,my),len)
      call gsync()                            ! make sure all are posted

      !len = (mx+1)*wdsize                     ! wdsize from MGRID
      call csend0(mtagsd,u(0,   1),len,nidd,0) ! nidd from MGRID
      call csend0(mtagsu,u(0,my-1),len,nidu,0) ! nidu from MGRID

      call msgwait(imsg_bot)
      call msgwait(imsg_top)

      return
      end
c-----------------------------------------------------------------------
      subroutine out_mgrid
      include 'MGRID'
      integer ic
      save    ic
      data    ic /0/

      do mid=0,np-1
         call gsync()
         if (mid.eq.nid) then
            write(6,*) nid
            write(6,1) ic,nid,ip,iq,ir,mp,mq,mr
            write(6,2) ic,nid,istx,isty,istz,iedx,iedy,iedz
            write(6,3) ic,nid,mx,my,mz,nx,ny,nz
            write(6,4) ic,nid,mtagsl,mtagsr,mtagsd,mtagsu,mtagsb,mtagsf
            write(6,5) ic,nid,mtagrl,mtagrr,mtagrd,mtagru,mtagrb,mtagrf
            write(6,6) ic,nid,nidl,nidr,nidd,nidu,nidb,nidf
            write(6,9) ic,nid,lev,istx,iedx,mx,ip,nx,np,(i,i=istx,iedx)
            write(6,9) ic,nid,lev,isty,iedy,my,iq,ny,np,(i,i=isty,iedy)
            write(6,9) ic,nid,lev,istz,iedz,mz,ir,nz,np,(i,i=istz,iedz)
            write(6,*) nid
         endif
         call gsync()
      enddo
    1 format(i4,i5,'qq ip,iq,...,mr:',8i5)
    2 format(i4,i5,'qq istx,isty,..:',8i5)
    3 format(i4,i5,'qq mx,my,...,nz:',8i5)
    4 format(i4,i5,'qq mtagsl,.....:',8i5)
    5 format(i4,i5,'qq mtagrl,.....:',8i5)
    6 format(i4,i5,'qq nidl,nidr,..:',8i5)
    9 format(i4,i5,'qqx',7i4,' :',12i4)    

      ic = ic+1

      return
      end
c-----------------------------------------------------------------------
      subroutine outmp(a,n,m,name4,u,w)
      include 'MGRID'

      real a(0:mx,0:my,0:mz),u(-1:nx,0:ny,0:mz),w(-1:nx,0:ny,0:mz)
      character*4 name4

c     return

      do kx=0,nx
         call outm(a,kx,n,m,name4)
      enddo

      mx10 = min(mx-1,18)

      nn = (nx+2)*(ny+1)*(nz+1)
      call rzero(u,nn)
      call rzero(w,nn)

      if (mx.gt.1.and.my.gt.1.and.mz.gt.1) then

        k=0
        do kk=istz+1,iedz-1

          k=k+1
          j=0
          do jj=isty+1,iedy-1

            j=j+1
            i=0

            do ii=istx+1,iedx-1
               i=i+1
               u(ii,jj,kk) = a(i,j,k)
            enddo
            u(-1,jj,kk) = nid

          enddo

        enddo

      endif

      call gop(u,w,'+  ',nn)

      if (nid.eq.1) then
         write(6,2) lev,name4,n,m,mp,mq,mr,nx,ny,nz
  2      format(i2,'pp ',a4,10i4)
         do k=nz-1,1,-1
          write(6,*)
          do j=ny-1,1,-1
            mid = u(-1,j,k)
            write(6,1) lev,k,j,mid,(u(i,j,k),i=1,mx10)
  1         format(i2,'pp ',3i2,19f9.4)
          enddo
         enddo
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine smooth_q (u,r,mx,my,mz,h,sigma,m,uo,iz)! m sweep
      real u(1),r(1),uo(1)

      call smooth_m_jac3d_q(u,r,mx,my,mz,h,sigma,m,uo,iz)

      return
      end
c-----------------------------------------------------------------------
      subroutine smooth_m_jac3d_q (u,r,mx,my,mz,h,sigma,m,o,izero) !  Jacobi

      real u(0:mx,0:my,0:mz),r(0:mx,0:my,0:mz),o(0:mx,0:my,0:mz)
      common /ccflops/ cflops,pflops,rflops,rrlops

      h2       = h*h
      sa       = sigma/6

      if (izero.eq.0) then
         do k=1,mz-1
         do j=1,my-1
         do i=1,mx-1
            u(i,j,k)=sa*h2*r(i,j,k)
         enddo
         enddo
         enddo
         cflops = cflops + 2*(mx-1)*(my-1)*(mz-1)
      endif

      mm1 = (mx+1)*(my+1)*(mz+1)
      do i_sweep = 2,m+izero

        call copy    (o,u,mm1) ! Old copy, w/ boundary data
c       call xchange0(o,XCHANGE_JGX)

        do k=1,mz-1
        do j=1,my-1
        do i=1,mx-1
           u(i,j,k) = o(i,j,k) + sa*( h2*r(i,j,k)-6*o(i,j,k)
     &              + o(i-1,j,k) + o(i+1,j,k)
     &              + o(i,j-1,k) + o(i,j+1,k)
     &              + o(i,j,k-1) + o(i,j,k+1) )
        enddo
        enddo
        enddo
      enddo
      cflops = cflops + (11*(mx-1)*(my-1)*(mz-1))*(m+izero-1)

      return
      end
c-----------------------------------------------------------------------
