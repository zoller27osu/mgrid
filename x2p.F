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
c-----------------------------------------------------------------------
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

      call get_proc_mpqr (mp,mq,mr,np,ldim)             ! processor decomp
      call get_proc_ijk  (ip,iq,ir,nid,mp,mq,mr,0,ldim) ! nid = (ip,iq,ir)

      do itest=1,2
        icalld = 0
        do i = 1,90
           call mgrid_tester(nlev_last)
        enddo
        close(9)
      enddo

      call exitt
      end
c-----------------------------------------------------------------------
      subroutine mgrid_tester(nlev_last)
      include 'MGRID'

      common /xcdata/ u(0:lt),r(0:lt)
      common /xddata/ ue(0:lt),rf(0:lt)
      real sigma

      parameter(ml=2*lpow)
      integer tflag(ml),dmflag(ml)
      common /cflag/ tflag,dmflag
      common /cinit/ icalld


      call init_mg(nsmooth,igs,sigma)       ! Get current mg parameters

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

      call dmesh     (mm,n)

      call rzero     (u ,mm)
      call rzero     (r ,mm)
      call rzero     (rf,mm)
      call rzero     (ue,mm)

      call initu     (ue,n,h,ldim)             ! Exact solution
      call xchange0  (ue)

      call compute_r (r,ue,rf,mx,my,mz,h,ldim) ! b = 0 - A ue
      call chsign    (r,mm)                    ! b = A ue
      call xchange0  (r)

      tol   = 1.e-10


      lev = nlev
      nf  = 2**lev
      n   = nf
      h   = 1./n
      kf  = 0
      call dmesh(mm,n)
      iz=1

      dmax = 0.
      call gsync
      time0 = dclock()
      do ipass=1,100
         call smooth_q  (u(kf),r(kf),mx,my,mz,h
     &                          ,sigma,nsmooth,rf,iz,igs,ldim)
      enddo
#ifdef DEBUG_OUT
      time1 = dclock()-time0
      call errchk (emx,u,ue,ipass,time1,dmax,'A')
#endif

      call gsync
#ifdef DEBUG_OUT
      time0 = dclock()
#endif
      do ipass=1,100
         call smooth_m  (u(kf),r(kf),mx,my,mz,h
     &                          ,sigma,nsmooth,rf,iz,igs,ldim)
         call xchange0  (u(kf))
      enddo
#ifdef DEBUG_OUT
      time1 = dclock()-time0
      call errchk (emx,u,ue,ipass,time1,dmax,'B')
#endif

      call gsync
#ifdef DEBUG_OUT
      time0 = dclock()
#endif
      do ipass=1,20

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
            call xchange0  (u(kf))

            call compute_r (rf,u(kf),r(kf),mx,my,mz,h,ldim)
            call xchange0  (rf)

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
     $                               ,ieven,jeven,keven,ldim)
        enddo

        one = 1.
        iz  = 0
        call smooth_m(u(kf),r(kf),mx,my,mz,h
     &                              ,one,1,rf,iz,igs,ldim)
        call reset_mtags (lev+1)  ! Exchange on prior level processor set
        call xchange0    (u(kf))

        kc = kf
        nf = n*2
        h  = h/2

        do ilev=nlev-1,1,-1

           mxc   = mx        ! mxc is coarse element count
           myc   = my
           mzc   = mz

           lev = lev+1
           call dmesh(mm,nf)
           kf     = kc-mm

           ieven = 2-mod(istx,2)
           jeven = 2-mod(isty,2)
           keven = 2-mod(istz,2)

           ! prolongate and add

           call prolong_a (u(kf),dmax,u(kc),mxc,myc,mzc
     $                                ,ieven,jeven,keven,ldim)
           call reset_mtags (lev+1)  ! Exchange on prior level processor set
           call xchange0    (u(kf))

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

      time1 = dclock()-time0
      call errchk (emx,u,ue,ipass,time1,dmax,'X')
#ifdef NUMERICAL_OUT
c     call out_mgrid()
c     call outmp(u ,ipass,lev,'SOLN',rf,r)
      call outmp(ue,ipass,lev,'EXAC',rf,r)
#endif

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
     $   'WARNING:  Attempt to take MOD(I,0) in function mod1.'
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
      subroutine errchk(emx,u,ue,ii,time,dmax,c1)
      include 'MGRID'
      real u(1),ue(1)
      integer*8 n3
      character*1 c1

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
      if (nid.eq.0) then
        write(6,1)
     $    ii,dmax,emx,time,ldim,nlev,nx,ny,nz,n3,mp,mq,mr,np,c1
    1   format(i7,1p3e12.4,i2,i3,3i6,i13,'=n ',3i5,i9,'  max error',a1)
      endif

      cflops = cflops*np
      pflops = pflops*np
      rflops = rflops*np
      rrlops = rrlops*np
      fflops = (cflops+pflops+rflops+rrlops)
      if (time.gt.0) fflops = fflops/(time*1.e6)

      if (nid.eq.0) write(6,2)
     $  ii,emx,time,fflops,cflops,pflops,rflops,rrlops,n3,np,c1
    2 format(i7,1p7e12.4,i13,'=n ',i9,'  flops',a1)

      call reset_flops ! Reset flop counters

      return
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
      subroutine init_mg(nsmooth,igs,sigma)! Get current mg parameters

      include 'MGRID'
      common /cinit/ icalld

      nlev_max = (lpow-1) + log2(np)/ldim

      if (nid.eq.0) then
         if (icalld.eq.0) open(unit=9,file='in.dat')
#ifdef DEBUG_OUT
         write(6,*) 'Input nlev, nsmooth, igs, sigma:'
#endif
         nlev=0
         read (9,*,end=999)  nlev,nsmooth,igs,sigma
         nlev = min(nlev_max,nlev)
  999    write(6,*)  nlev,nsmooth,igs,sigma
      endif

      icalld = 1

      call bcast(nlev   ,     4) ! 4-byte integer
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

c      if (nid.eq.0) write(6,*) lev,nlev,ipass,nx,' level'

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
c    $   ,nidu,nidd,(ntags(k,ilev),ntagr(k,ilev),k=3,4)
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
      integer gid
c     ===========  Assign global numbering =======
      mpq = mp*mq

      gid = -1
      if (0.le.ip.and.ip.lt.mp.and.    ! Inside Box
     $    0.le.iq.and.iq.lt.mq.and.
     $    0.le.ir.and.ir.lt.mr      )  gid = ip + iq*mp + ir*mpq


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
c     This routine calculates mx.my,mz and their global start indx
c
c     Partitions nx points across mp procs
c     mx is number of element per proc (0---mx)

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
      real u(1)

      if (ndim.eq.2) call initu_2d(u,n,h)
      if (ndim.eq.3) call initu_3d(u,n,h)

      call xchange0  (u)

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
     $            + .25*sin(pik*x)*sin(pik*y)
     $            + .15*sin(pij*x)*sin(pij*y)
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
     $            + .25*sin(pik*x)*sin(pik*y)
     $            + .15*sin(pij*x)*sin(pij*z)
     $            + .15*sin(pij*z)*sin(pij*y)
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
      real r(1),u(1),ro(1)

      if (ndim.eq.2) call compute_r2d (r,u,ro,mx,my,h)
      if (ndim.eq.3) call compute_r3d (r,u,ro,mx,my,mz,h)

      return
      end
c-----------------------------------------------------------------------
      subroutine compute_r2d (r,u,ro,mx,my,h)
      real r(0:mx,0:my),u(0:mx,0:my),ro(0:mx,0:my)

      h2     = h*h
      h2i    = 1./h2

      do j=1,my-1
      do i=1,mx-1
         r(i,j)=h2i
     $       *(h2*ro(i,j)+u(i-1,j)+u(i,j-1)+u(i+1,j)+u(i,j+1)-4*u(i,j))
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine compute_r3d (r,u,ro,mx,my,mz,h)
      real r(0:mx,0:my,0:mz)
      real u(0:mx,0:my,0:mz),ro(0:mx,0:my,0:mz)

      common /ccflops/ cflops,pflops,rflops,rrlops

      h2     = h*h
      h2i    = 1./h2


      do k=1,mz-1
      do j=1,my-1
      do i=1,mx-1
         r(i,j,k)=h2i*(h2*ro(i,j,k)-6*u(i,j,k)
     $       + u(i-1,j,k) + u(i+1,j,k)
     $       + u(i,j-1,k) + u(i,j+1,k)
     $       + u(i,j,k-1) + u(i,j,k+1))
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
     $        +(rf(im1,j)+rf(ip1,j)+rf(i,jm1)+rf(i,jp1))*.125
     $        +(rf(im1,jm1)+rf(ip1,jm1)+rf(im1,jp1)+rf(ip1,jp1))*.0625

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
     $      + .0625 * ( rf(im1,j,k)+rf(i,jm1,k)+rf(i,j,km1)
     $                + rf(ip1,j,k)+rf(i,jp1,k)+rf(i,j,kp1) ) ! +8 flops

     $      + .03125* ( rf(i,  jm1,km1)+rf(i,  jp1,km1)
     $                + rf(i,  jm1,kp1)+rf(i,  jp1,kp1)
     $                + rf(im1,j,  km1)+rf(ip1,j,  km1)
     $                + rf(im1,j,  kp1)+rf(ip1,j,  kp1)
     $                + rf(im1,jm1,k  )+rf(ip1,jm1,k  )
     $                + rf(im1,jp1,k  )+rf(ip1,jp1,k  ) )    ! +12 flops

     $      + .015625*( rf(im1,jm1,km1)+rf(ip1,jm1,km1)
     $                + rf(im1,jp1,km1)+rf(ip1,jp1,km1)
     $                + rf(im1,jm1,kp1)+rf(ip1,jm1,kp1)
     $                + rf(im1,jp1,kp1)+rf(ip1,jp1,kp1) )    ! +8 flops

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

      real uf(1),uc(1)

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
     $          (uc(icm1,jcm1)+uc(icm1,jc)+uc(ic,jcm1)+uc(ic,jc))

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
     $                  uc(ic,jcm1,kc)+uc(ic,jc,kcm1)+uc(ic,jcm1,kcm1))
            uf(i1,j2,k1) = uf(i1,j2,k1) + 0.25* ( uc(ic,jc,kc) +
     $                  uc(icm1,jc,kc)+uc(ic,jc,kcm1)+uc(icm1,jc,kcm1))
            uf(i1,j1,k2) = uf(i1,j1,k2) + 0.25* ( uc(ic,jc,kc) +
     $                  uc(icm1,jc,kc)+uc(ic,jcm1,kc)+uc(icm1,jcm1,kc))

            ! Center !
            dctr = .125* ( uc(icm1,jcm1,kcm1) + uc(ic  ,jcm1,kcm1) +
     $                     uc(icm1,jc  ,kcm1) + uc(ic  ,jc  ,kcm1) +
     $                     uc(icm1,jcm1,kc  ) + uc(ic  ,jcm1,kc  ) +
     $                     uc(icm1,jc  ,kc  ) + uc(ic  ,jc  ,kc  ) )

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
       if (igs.eq.0) call smooth_m_jac2d(u,r,mx,my,h,sigma,m,uo,iz)
       if (igs.eq.1) call smooth_m_gs_2d(u,r,mx,my,h,sigma,m,iz)
      else
       if (igs.eq.0) call smooth_m_jac3d(u,r,mx,my,mz,h,sigma,m,uo,iz)
       if (igs.eq.1) call smooth_m_gs_3d(u,r,mx,my,mz,h,sigma,m,iz)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine smooth_m_gs_2d (u,r,mx,my,h,sigma,m,izero) !  GS
      real u(0:mx,0:my),r(0:mx,0:my)

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
     $         *(h2*r(i,j)+u(i-1,j)+u(i,j-1)+u(i+1,j)+u(i,j+1)-4*u(i,j))
        enddo
        enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine smooth_m_gs_3d (u,r,mx,my,mz,h,sigma,m,izero) !  Jacobi
      real u(0:mx,0:my,0:mz),r(0:mx,0:my,0:mz)

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
     $               + sh2*(u(i-1,j,k) + u(i,j-1,k) + u(i,j,k-1))
     $               + s22*(r(i+1,j,k) + r(i,j+1,k) + r(i,j,k+1))
         enddo
         enddo
         enddo
      endif

      do i_sweep = 3,m+2*izero
        do k=1,mz-1
        do j=1,my-1
        do i=1,mx-1
           u(i,j,k) = u(i,j,k) + sa*( h2*r(i,j,k)-6*u(i,j,k)
     $              + u(i-1,j,k) + u(i+1,j,k)
     $              + u(i,j-1,k) + u(i,j+1,k)
     $              + u(i,j,k-1) + u(i,j,k+1) )
        enddo
        enddo
        enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine smooth_m_jac2d (u,r,mx,my,h,sigma,m,o,izero) !  Jacobi
      real u(0:mx,0:my),r(0:mx,0:my),o(0:mx,0:my)

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
        call xchange0 (o)

        do j=1,my-1
        do i=1,mx-1
           u(i,j)=o(i,j) + sa
     $         *(h2*r(i,j)+o(i-1,j)+o(i,j-1)+o(i+1,j)+o(i,j+1)-4*o(i,j))
        enddo
        enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine smooth_m_jac3d (u,r,mx,my,mz,h,sigma,m,o,izero) !  Jacobi

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
        call xchange0(o)

        do k=1,mz-1
        do j=1,my-1
        do i=1,mx-1
           u(i,j,k) = o(i,j,k) + sa*( h2*r(i,j,k)-6*o(i,j,k)
     $              + o(i-1,j,k) + o(i+1,j,k)
     $              + o(i,j-1,k) + o(i,j+1,k)
     $              + o(i,j,k-1) + o(i,j,k+1) )
        enddo
        enddo
        enddo
      enddo
      cflops = cflops + (11*(mx-1)*(my-1)*(mz-1))*(m+izero-1)

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
  1            format(i2,'oo ',2i2,19f19.16)
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
      subroutine xchange0(u) ! Zeros out boundary data

      include 'MGRID'
      parameter (lb=(lx+2)**(ldim-1))
      common /cbuff/ b1(lb),b2(lb),bo(lb)

      real u(1)

      if (ldim.eq.2) call xchange_2d(u,b1,b2,bo)
      if (ldim.eq.3) call xchange_3d(u,b1,b2,bo)

      return
      end
c-----------------------------------------------------------------------
      subroutine xchange_3d(u,buff1,buff2,buffo)
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
      call csend0(mtagsd,buffo,len,nidd,0) ! nidb from MGRID

      l=0                        !        ^ y
      do k = 0,mz                !        |
      do i = 0,mx                !        +---> x
         buffo(l) = u(i,my-1,k)  !      z/
         l=l+1                   ! load top face into outbound buff
      enddo
      enddo
      len = l*wdsize
      call csend0(mtagsu,buffo,len,nidu,0) ! nidt from MGRID


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

      len = (mx+1)*(my+1)*wdsize                ! wdsize from MGRID
      imsg_back  = irecv0(mtagrb,u(0,0,0 ),len) ! mtags from MGRID
      imsg_front = irecv0(mtagrf,u(0,0,mz),len)
      call gsync()                              ! make sure all are posted

      len = (mx+1)*(my+1)*wdsize                ! wdsize from MGRID
      call csend0(mtagsb,u(0,0,   1),len,nidb,0) ! nid2 from MGRID
      call csend0(mtagsf,u(0,0,mz-1),len,nidf,0) ! nid1 from MGRID

      call msgwait(imsg_front)
      call msgwait(imsg_back)

      return
      end
c-----------------------------------------------------------------------
      subroutine xchange_2d(u,buff1,buff2,buffo)

      include 'MGRID'

      real u(0:mx,0:my),buff1(0:1),buff2(0:1),buffo(0:1)


c     ======    Exchange x-faces  =======

      myz = (my+1)

      len = myz*wdsize                         ! wdsize from MGRID
      imsg_left  = irecv0(mtagrl,buff1,len)    ! mtags from MGRID
      imsg_right = irecv0(mtagrr,buff2,len)
      call gsync()                            ! make sure all are posted

      len = myz*wdsize
      do j = 0,my
         buffo(j) = u(1,j)          ! load left face into outbound buff
      enddo
      call csend0(mtagsl,buffo,len,nidl,0) ! nidl from MGRID

      do j = 0,my
         buffo(j) = u(mx-1,j)       ! load right face into outbound buff
      enddo
      call csend0(mtagsr,buffo,len,nidr,0) ! nidr from MGRID

      call msgwait(imsg_right)
      do j = 0,my
         u(mx,j)= buff2(j)          ! unpack inbound right face
      enddo

      call msgwait(imsg_left)
      do j = 0,my
         u(0,j)= buff1(j)           ! unpack inbound left face
      enddo


c     ======    Exchange y-faces, no buffers  =======

      len = (mx+1)*wdsize                     ! wdsize from MGRID

      imsg_bot = irecv0(mtagrd,u(0,0 ),len)    ! mtags from MGRID
      imsg_top = irecv0(mtagru,u(0,my),len)
      call gsync()                            ! make sure all are posted

      len = (mx+1)*wdsize                     ! wdsize from MGRID
      call csend0(mtagsd,u(0,   1),len,nidd,0) ! nidb from MGRID
      call csend0(mtagsu,u(0,my-1),len,nidu,0) ! nidt from MGRID

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

      if (nid.eq.0) then
c        write(6,2) lev,name4,n,m,mp,mq,mr,nx,ny,nz
c 2      format(i2,'pp ',a4,10i4)
         do k=nz-1,1,-1
          write(6,*)
          do j=ny-1,1,-1
            mid = u(-1,j,k)
c           write(6,1) lev,k,j,mid,(u(i,j,k),i=1,mx10)
c 1         format(i2,'pp ',3i2,19f9.4)
          enddo
         enddo
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine smooth_q (u,r,mx,my,mz,h,sigma,m,uo,iz,igs,ndim)! m sweep
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
c       call xchange0(o)

        do k=1,mz-1
        do j=1,my-1
        do i=1,mx-1
           u(i,j,k) = o(i,j,k) + sa*( h2*r(i,j,k)-6*o(i,j,k)
     $              + o(i-1,j,k) + o(i+1,j,k)
     $              + o(i,j-1,k) + o(i,j+1,k)
     $              + o(i,j,k-1) + o(i,j,k+1) )
        enddo
        enddo
        enddo
      enddo
      cflops = cflops + (11*(mx-1)*(my-1)*(mz-1))*(m+izero-1)

      return
      end
c-----------------------------------------------------------------------
