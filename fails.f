      subroutine swap(ix, iy)
        it = ix
        ix = iy
        iy = it
      end

c-----------------------------------------------------------------------
      call gsync
      time0 = dclock()
      do ipass=1,100
         call smooth_m  (u(kf),r(kf),mx,my,mz,h
     &                          ,sigma,nsmooth,rf,iz,5,ldim)
         call xchange0  (u(kf))
      enddo
      time1 = dclock()
      call errchk (emx,u,ue,ipass,time1-time0,dmax,'S')

c-----------------------------------------------------------------------
      subroutine smooth_m_jac3dsw (u,r,mx,my,mz,h,sigma,m,o,izero) !  Jacobi
      common /ccflops/ cflops,pflops,rflops,rrlops

      !real, target :: u(0:mx,0:my,0:mz), o(0:mx,0:my,0:mz)
      real, dimension(0:mx,0:my,0:mz), target :: u, o
      real r(0:mx,0:my,0:mz)
      real, dimension(:,:,:), pointer :: pu, po, temp
      pu => u; po => o

      h2       = h*h
      sa       = sigma/6

      if (izero.eq.0) then
         do k=1,mz-1
         do j=1,my-1
         do i=1,mx-1
            pu(i,j,k)=sa*h2*r(i,j,k)
         enddo
         enddo
         enddo
         cflops = cflops + 2*(mx-1)*(my-1)*(mz-1)
      endif

      mm1 = (mx+1)*(my+1)*(mz+1)
      do i_sweep = 2,m+izero
        !call copy    (o,u,mm1) ! Old copy, w/ boundary data
        temp => pu; pu => po; po => temp
        call xchange0(o)

        do k=1,mz-1
        do j=1,my-1
        do i=1,mx-1
           pu(i,j,k) = po(i,j,k) + sa*( h2*r(i,j,k)-6*po(i,j,k)
     $              + po(i-1,j,k) + po(i+1,j,k)
     $              + po(i,j-1,k) + po(i,j+1,k)
     $              + po(i,j,k-1) + po(i,j,k+1) )
        enddo
        enddo
        enddo
      enddo
      cflops = cflops + (11*(mx-1)*(my-1)*(mz-1))*(m+izero-1)

      if (((m+izero-1)%2).ne.0) then !pu points to o
        call copy (u,o,mm1)
      endif

      return
      end

c-----------------------------------------------------------------------
      subroutine smooth_m_jac3cnd (u,r,mx,my,mz,h,sigma,m,o,izero) !  Jacobi
      real, dimension(0:mx,0:my,0:mz) :: u, r, o

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

      do i_sweep = 2,m+izero
        !call copy    (o,u,mm1) ! Old copy, w/ boundary data
        if (MOD(i_sweep,2).eq.0) then
          call jac3d_helper(mx,my,mz,o,u,r,h2,sa)
          !call xchange0(u)

          !do k=1,mz-1
          !do j=1,my-1
          !do i=1,mx-1
          !  o(i,j,k) = u(i,j,k) + sa*( h2*r(i,j,k)-6*u(i,j,k)
     $    !            + u(i-1,j,k) + u(i+1,j,k)
     $    !            + u(i,j-1,k) + u(i,j+1,k)
     $    !            + u(i,j,k-1) + u(i,j,k+1) )
          !enddo
          !enddo
          !enddo
        else
          call jac3d_helper(mx,my,mz,u,o,r,h2,sa)
          !call xchange0(o)

          !do k=1,mz-1
          !do j=1,my-1
          !do i=1,mx-1
          !  u(i,j,k) = o(i,j,k) + sa*( h2*r(i,j,k)-6*o(i,j,k)
     $    !            + o(i-1,j,k) + o(i+1,j,k)
     $    !            + o(i,j-1,k) + o(i,j+1,k)
     $    !            + o(i,j,k-1) + o(i,j,k+1) )
          !enddo
          !enddo
          !enddo
        end if
      enddo
      cflops = cflops + (11*(mx-1)*(my-1)*(mz-1))*(m+izero-2+1)

      if (MOD(m+izero,2).eq.0) then
        mm1 = (mx+1)*(my+1)*(mz+1)
        call copy (u,o,mm1) ! o was updated last
      endif

      return
      end
c-----------------------------------------------------------------------
      !if (mtagrl.ge.0) then
        !call mpi_irecv(u(0,0,0),myz,realcol,mpi_any_source,mtagrl
     $  !     ,mpi_comm_world,imsg,ierr)
        !imsg_left = imsg
      !else
        !imsg_left = MPI_Request_Null
      !endif

      !if (mtagrr.ge.0) then
      !  call mpi_irecv(u(mx,0,0),myz,realcol,mpi_any_source,mtagrr
     $!       ,mpi_comm_world,imsg,ierr)
      !  imsg_right = imsg
      !else
      !  imsg_right = MPI_Request_Null
      !endif

      !if ((mtagsl.ge.0).and.(jnid.ge.0).and.(jnid.lt.np)) then
      !call mpi_send(u(1,0,0),1,realcol,nidl,mtagsl,MPI_COMM_WORLD,ierr)
      !endif
c-----------------------------------------------------------------------
      subroutine make_type(size, ioldtype, inewtype, ierr)
      call MPI_TYPE_CONTIGUOUS(size, ioldtype, inewtype, ierr)
      call MPI_TYPE_COMMIT(inewtype, ierr)
      
      return
      end

              ! previously in arreq
              ioff = -1
              if (n.gt.((mx-1)*(my-1)*(mz-1))) ioff = 1
              k = MIN(mz+ioff, idx/((mx+ioff)*(my+ioff)) + 1)
              i = idx - (mx+ioff)*(my+ioff)*(k-1)
              j = MIN(my+ioff, i/(mx+ioff) + 1)
              i = i - (mx+ioff) * (j-1)
              write(6,*) "different at",i,j,k,"(",a(idx),"vs",b(idx),")"

! previously in xchange0 -----------------------------------------------
      if (ldim.eq.2) call xchange_2d(u,b1,b2,bo)
      if (ldim.eq.3) then
#ifdef verify_xchange
        mxyz = (mx+1)*(my+1)*(mz+1)
        call copy(v,u,mxyz)             ! must be done by all processes
        call xchange_3d_old(v,b1,b2,bo) ! to xchange results successfully
#endif
#ifdef xchange_async
	      call xchange_3d_async(u)
#else
        call xchange_3d(u)
        !call xchange_3d_old(u,b1,b2,bo)
#endif
#ifdef verify_xchange
        if ((nid.eq.6)) then!.or.(nid.eq.7)) then
          if (mateq(u,v,1)) then ! the 1 means we only want to compare elements that are at the boundary in at most 1 dimensions
            !write(6,*) nid, ": xchange with datatypes is working!"
          else
            write(6,*) nid, ": new xchange is not working..."
            write(6,*) nid, ": old is "
            call print_3d(v, 5, 2, 2)
            write(6,*) nid, ": new is "
            call print_3d(u, 5, 2, 2)

            call exit(11) ! Y U NO SEGFAULT?!?!?!

          endif
        endif
#endif
      endif

      !if (nid.eq.6) then
        !call print_3d(u, 5, 2, 2)
        !write(6, '(8F8.4)') "u: ", real(u(0:mx, 0, 0))
        !write(6, '(8F8.4)') "v: ", real(v(0:mx, 0, 0))
        !write(6,*) "u:", real(u(0:3, 0, 1))
        !write(6,*) "v:", real(v(0:3, 0, 1))
        !write(6,*) " "
      !endif
 
      return
      end

