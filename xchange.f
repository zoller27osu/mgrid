c-----------------------------------------------------------------------
      subroutine xchange0(u,nx,ny,nz,ndim)
!     Zeros out boundary data

      real u(1)

      if (ndim.eq.2) call xchange_2d(u,nx,ny)
      if (ndim.eq.3) call xchange_3d(u,nx,ny,nz)
 
      return
      end
c-----------------------------------------------------------------------
      subroutine xchangeb(u,nx,ny,nz,ndim,w)
!     Preserves boundary data

      real u(0:nx,0:ny,0:nz),w(0:nx,0:ny,0:nz)

      n = (nx+1)*(ny+1)
      if (ndim.eq.3) n = n*(nz+1)
      call copy(w,u,n)

      if (ndim.eq.2) call xchange_2d(w,nx,ny)
      if (ndim.eq.3) call xchange_3d(w,nx,ny,nz)

      if (ndim.eq.3) then
         do k=1,nz-1
         do j=1,ny-1
         do i=1,nx-1
            u(i,j,k)=w(i,j,k)
         enddo
         enddo
         enddo
      else
         do j=1,ny-1
         do i=1,nx-1
            u(i,j,0)=w(i,j,0)
         enddo
         enddo
      endif
 
      return
      end
c-----------------------------------------------------------------------
      subroutine xchange_3d(u,nx,ny,nz) ! Zeros out boundary data !
      include 'MGRID'

      common /cbuff/ buff1(nxzm),buff2(nxzm),buffo(nxzm)

      real u(0:nx,0:ny,0:nz)

      call checkit('nx nxm xchange_3d$',nx,nxm)

c     ======    Exchange x-faces  =======

      nxyz = (ny+1)*(nz+1)
      call rzero(buff1,nxyz)        ! zero out buffers         
      call rzero(buff2,nxyz)

      len = nxzm*wdsize                        ! wdsize from MGRID
      imsg_left  = irecv0(mtagrl,buff1,len)    ! mtags from MGRID
      imsg_right = irecv0(mtagrr,buff2,len)    ! Set buff=0 if mtag < 0
      call gsync()                             ! make sure all are posted

      nyz1 = (ny+1)*(nz+1)-1
      l=0
      do j = 0,nyz1
         l=l+1
         buffo(l) = u(1   ,j,0)     ! load left face into outbound buff
      enddo
      len = l*wdsize
      call csend0(mtagsl,buffo,len,nidl,0) ! nidl from MGRID, null if mtag < 0

      l=0
      do j = 0,nyz1
         l=l+1
         buffo(l) = u(nx-1,j,0)     ! load right face into outbound buff
      enddo
      len = l*wdsize
      call csend0(mtagsr,buffo,len,nidr,0) ! nidr from MGRID

      call msgwait(imsg_right)
      l=0
      do j = 0,nyz1
         l=l+1
         u(nx,j,0)= buff2(l)        ! unpack inbound right face
      enddo

      call msgwait(imsg_left)
      l=0
      do j = 0,nyz1
         l=l+1
         u(0,j,0)= buff1(l)         ! unpack inbound left face
      enddo


c     ======   Exchange y-faces  =======

      nxyz = (nx+1)*(nz+1)
      call rzero(buff1,nxyz)     ! zero out buffers         
      call rzero(buff2,nxyz)

      len = nxzm*wdsize                       ! wdsize from MGRID
      imsg_bot   = irecv0(mtagrd,buff1,len)   ! mtags from MGRID
      imsg_top   = irecv0(mtagru,buff2,len)
      call gsync()                            ! make sure all are posted


      l=0                        !        ^ y
      do k = 0,nz                !        |
      do i = 0,nx                !        +---> x
         l=l+1                   !      z/
         buffo(l) = u(i,1,k)     ! load bottom face into outbound buff
      enddo
      enddo
      len = l*wdsize
      call csend0(mtagsd,buffo,len,nidd,0) ! nidb from MGRID

      l=0                        !        ^ y
      do k = 0,nz                !        |
      do i = 0,nx                !        +---> x
         l=l+1                   !      z/
         buffo(l) = u(i,ny-1,k)  ! load top face into outbound buff
      enddo
      enddo
      len = l*wdsize
      call csend0(mtagsu,buffo,len,nidu,0) ! nidt from MGRID


      call msgwait(imsg_top)
      l=0                        !        ^ y
      do k = 0,nz                !        |
      do i = 0,nx                !        +---> x
         l=l+1                   !      z/
         u(i,ny,k)= buff2(l)     ! unpack onto top face 
      enddo
      enddo

      call msgwait(imsg_bot)
      l=0                        !        ^ y
      do k = 0,nz                !        |
      do i = 0,nx                !        +---> x
         l=l+1                   !      z/
         u(i,0,k)= buff1(l)      ! unpack onto bottom face 
      enddo
      enddo

c     ======    Exchange z-faces, no buffers  =======

      len = (nx+1)*(ny+1)*wdsize                ! wdsize from MGRID
      imsg_back  = irecv0(mtagrb,u(0,0,0 ),len) ! mtags from MGRID
      imsg_front = irecv0(mtagrf,u(0,0,nz),len) 
      call gsync()                              ! make sure all are posted

      len = (nx+1)*(ny+1)*wdsize                ! wdsize from MGRID
      call csend0(mtagsb,u(0,0,   1),len,nidb,0) ! nid2 from MGRID
      call csend0(mtagsf,u(0,0,nz-1),len,nidf,0) ! nid1 from MGRID

      call msgwait(imsg_front)
      call msgwait(imsg_back)

      return
      end
c-----------------------------------------------------------------------
      subroutine xchange_2d(u,nx,ny) ! Zeros out boundary data !
      include 'MGRID'

      common /cbuff/ buff1(nxzm),buff2(nxzm),buffo(nxzm)

      real u(0:nx,0:ny)

      call checkit('nx nxm xchange_2d$',nx,nxm)


c     ======    Exchange x-faces  =======

      nxyz = (ny+1)
      call rzero(buff1,nxyz)     ! zero out buffers         
      call rzero(buff2,nxyz)

      len = nxzm*wdsize                       ! wdsize from MGRID
      imsg_left  = irecv0(mtagrl,buff1,len)    ! mtags from MGRID
      imsg_right = irecv0(mtagrr,buff2,len)
      call gsync()                            ! make sure all are posted

      l=0
      do j = 0,ny
         l=l+1
         buffo(l) = u(1,j)          ! load left face into outbound buff
      enddo
      len = l*wdsize
      call csend0(mtagsl,buffo,len,nidl,0) ! nidl from MGRID

      l=0
      do j = 0,ny
         l=l+1
         buffo(l) = u(nx-1,j)       ! load right face into outbound buff
      enddo
      len = l*wdsize
      call csend0(mtagsr,buffo,len,nidr,0) ! nidr from MGRID

      call msgwait(imsg_right)
      l=0
      do j = 0,ny
         l=l+1
         u(nx,j)= buff2(l)          ! unpack inbound right face
      enddo

      call msgwait(imsg_left)
      l=0
      do j = 0,ny
         l=l+1
         u(0,j)= buff1(l)           ! unpack inbound left face
      enddo


c     ======    Exchange y-faces, no buffers  =======

      len = (nx+1)*wdsize                     ! wdsize from MGRID
      
      imsg_bot = irecv0(mtagrd,u(0,0 ),len)    ! mtags from MGRID
      imsg_top = irecv0(mtagru,u(0,ny),len)
      call gsync()                            ! make sure all are posted

      len = (nx+1)*wdsize                     ! wdsize from MGRID
      call csend0(mtagsd,u(0,   1),len,nidd,0) ! nidb from MGRID
      call csend0(mtagsu,u(0,ny-1),len,nidu,0) ! nidt from MGRID

      call msgwait(imsg_bot)
      call msgwait(imsg_top)

      return
      end
c-----------------------------------------------------------------------
