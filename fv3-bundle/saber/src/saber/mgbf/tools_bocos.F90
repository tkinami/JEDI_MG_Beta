!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        module tools_bocos
!***********************************************************************
!                                                                      !
!  Provide communication between subdomains and supply halos on        !
!  filter grid                                                         !
!                                                                      !
!                                                     M. Rancic (2020) !
!***********************************************************************
use tools_kinds, only: kind_real,kind_int
use fckit_configuration_module, only: fckit_configuration
use fckit_mpi_module, only: fckit_mpi_comm,fckit_mpi_sum, fckit_mpi_status
use type_mpl, only: mpl_type
use type_mgbf_grid

implicit none

public:: multiply_add

public:: bocosHn
public:: bocosHnT

public:: upsend
public:: downsend

public:: boco05
public:: v02v


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine multiply_add                         &
!***********************************************************************
!                                                                      *
!         Adjoint test                                                 *
!                                                                      *
!***********************************************************************
(mpl,grid,a,b,im,jm,hx,hy,res_glob)
!-----------------------------------------------------------------------
!use mpi
!use mg_domain1, only: grid%Fleast,grid%Flwest,grid%Flsouth,grid%Flnorth

implicit none
!-----------------------------------------------------------------------
type(mpl_type),intent(inout) :: mpl
type(mgbf_grid_type),intent(inout) :: grid
integer(kind_int), intent(in):: im,jm,hx,hy
real(kind_real), dimension(-hx:im+hx,-hy:jm+hy),intent(in):: a,b
real(kind_real), intent(out) :: res_glob

real(kind_real) :: res
integer(kind_int) i,j,imax,jmax,ierr
!-----------------------------------------------------------------------


      res=0.
      res_glob=0.


      if(grid%Fleast(1)) then
        imax=im
      else
        imax=im-1
      endif     
      if(grid%Flnorth(1)) then
        jmax=jm
      else
        jmax=jm-1
      endif     

      do j=0,jmax
      do i=0,imax
        res=res+a(i,j)*b(i,j)
      end do
      end do

!-----------------------------------------------------------------------
      call mpl%f_comm%allreduce(res,res_glob,fckit_mpi_sum())
!      call MPI_REDUCE(res,res_glob,1,dtype,MPI_SUM,0,mpi_comm_comp,ierr)

!-----------------------------------------------------------------------
                        endsubroutine multiply_add                         

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine bocosHn                              &
!**********************************************************************!
!                                                                      !
! Side sending subroutine:                                             !
! Supplies (nbx,nby) lines of halos in (x,y) directions, including     !
! values at the edges of the subdomains and assuming mirror boundary   !
! conditions                                                           !
!                                                                      !
!**********************************************************************!
(mpl,grid,Warray,km,im,jm,Lm,nbx,nby,Fimax,Fjmax,mygen_min,mygen_max)
!-----------------------------------------------------------------------
!use mg_domain1, only: grid%Fitarg_n,grid%Fitarg_s,grid%Fitarg_w,grid%Fitarg_e
!&
!                     ,grid%Flwest,grid%Fleast,grid%Flsouth,grid%Flnorth
!                     &
!                     ,grid%Fitarg_ne,grid%Fitarg_se,grid%Fitarg_sw,grid%Fitarg_nw
!                     &
!                     ,grid%Flcorner_sw,grid%Flcorner_nw,grid%Flcorner_se,grid%Flcorner_ne

!use mpi

implicit none

!-----------------------------------------------------------------------
type(mpl_type),intent(inout) :: mpl
type(mgbf_grid_type),intent(inout) :: grid
integer(kind_int), intent(in):: km,im,jm,Lm,nbx,nby,mygen_min,mygen_max
real(kind_real),dimension(-nbx:im+nbx,-nby:jm+nby,1:Lm,km),intent(inout):: &
                                  Warray
integer(kind_int), dimension(grid%gm), intent(in):: Fimax,Fjmax
!-----------------------------------------------------------------------

real(kind_real), allocatable, dimension(:,:,:,:)::                         &
                                  sBuf_N,sBuf_E,sBuf_S,sBuf_W           &
                                 ,rBuf_N,rBuf_E,rBuf_S,rBuf_W
real(kind_real), dimension(nbx,nby,LM,km)::                                &
                                  sBuf_NE,sBuf_SE,sBuf_SW,sBuf_NW       &
                                 ,rBuf_NE,rBuf_SE,rBuf_SW,rBuf_NW

integer(kind_int) itarg_n,itarg_s,itarg_w,itarg_e                         &
               ,itarg_ne,itarg_se,itarg_sw,itarg_nw                     &
               ,imax,jmax
logical:: lcorner_sw,lcorner_nw,lcorner_se,lcorner_ne                   &
         ,lwest,least,lsouth,lnorth

integer(kind_int) sHandle(4),rHandle(4)!,ISTAT(MPI_STATUS_SIZE)
integer(kind_int) iaerr,ierr,iderr,l,i,j
integer(kind_int) isend,irecv,nebpe
integer(kind_int) ndatax,ndatay,nbxy
integer(kind_int) g_ind,g
logical l_sidesend
type(fckit_mpi_status) :: istat
!-----------------------------------------------------------------------
!
! Limit communications to selected number of generations
!

         l_sidesend=.false.

       if(mygen_min==1.and.mygen_max==1) then
         g_ind=1
         g = 1
         l_sidesend=.true.
       else &
       if(mygen_min <= grid%my_hgen .and. grid%my_hgen <= mygen_max) then
         g_ind=2
         g = grid%my_hgen
         l_sidesend=.true.
       endif

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

FILT_GRID:    if(l_sidesend) then

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!
! from mg_domain
!
          itarg_n = grid%Fitarg_n(g_ind)
          itarg_s = grid%Fitarg_s(g_ind)
          itarg_w = grid%Fitarg_w(g_ind)
          itarg_e = grid%Fitarg_e(g_ind)

          lwest   = grid%Flwest(g_ind)
          least   = grid%Fleast(g_ind)
          lsouth  = grid%Flsouth(g_ind)
          lnorth  = grid%Flnorth(g_ind)

          itarg_ne = grid%Fitarg_ne(g_ind)
          itarg_se = grid%Fitarg_se(g_ind)
          itarg_sw = grid%Fitarg_sw(g_ind)
          itarg_nw = grid%Fitarg_nw(g_ind)

          lcorner_sw = grid%Flcorner_sw(g_ind)
          lcorner_nw = grid%Flcorner_nw(g_ind)
          lcorner_se = grid%Flcorner_se(g_ind)
          lcorner_ne = grid%Flcorner_ne(g_ind)

          if(least) then
            imax = Fimax(g)
          else
            imax = im       !   << Note that is not necesseraly im from
          endif             !      mg_parameter.  Could be also imL >>>
          if(lnorth) then
            jmax = Fjmax(g)
          else
            jmax = jm
          endif


!-----------------------------------------------------------------------
      ndatay = km*(imax+1)*nby*Lm
      ndatax = km*(jmax+1)*nbx*Lm
      nbxy   = km*nbx*nby*Lm


!
!     SEND boundaries
!

! --- toward SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s

            allocate( sBuf_S(0:imax,nby,1:Lm,1:km), stat = iaerr )

              do L=1,Lm
                do j=1,nby
                  do i=0,imax
                    sBuf_S(i,j,L,:) = Warray(i,j,L,:)
                  enddo
                enddo
              enddo

              sHandle(3) = mpl%f_comm%isend(sBuf_S,nebpe, grid%mype)
!              call MPI_ISEND( sBuf_S, ndatay, dtype, nebpe, mype,  &
!                              mpi_comm_comp, sHandle(3), isend)
      end if

! --- toward NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n

            allocate( sBuf_N(0:imax,nby,1:Lm,1:km), stat = iaerr )

              do L=1,Lm
                do j=1,nby
                  do i=0,imax
                    sBuf_N(i,j,L,:)=Warray(i,jmax-nby-1+j,L,:)
                  enddo
                enddo
              enddo

              sHandle(1) = mpl%f_comm%isend(sBuf_N, nebpe, grid%mype)
!              call MPI_ISEND( sBuf_N, ndatay, dtype, nebpe, mype,        &
!                              mpi_comm_comp, sHandle(1), isend)

      end if

! --- toward WEST ---

      if( itarg_w >= 0) then
        nebpe = itarg_w

              allocate( sBuf_W(nbx,0:jmax,1:Lm,1:km), stat = iaerr )

              do L=1,Lm
                do j=0,jmax
                  do i=1,nbx
                    sBuf_W(i,j,L,:) = Warray(i,j,L,:)
                  enddo
                enddo
              enddo

              sHandle(4) = mpl%f_comm%isend(sBuf_W,nebpe, grid%mype)
!              call MPI_ISEND( sBuf_W, ndatax, dtype, nebpe, mype, &
!                              mpi_comm_comp, sHandle(4), isend)

      end if

! --- toward EAST ---

      if( itarg_e >= 0 ) then
        nebpe = itarg_e

              allocate( sBuf_E(nbx,0:jmax,1:Lm,1:km), stat = iaerr )

              do L=1,Lm
                do j=0,jmax
                  do i=1,nbx
                    sBuf_E(i,j,L,:) = Warray(imax-nbx-1+i,j,L,:)
                  enddo
                enddo
              enddo

              sHandle(2) =  mpl%f_comm%isend(sBuf_E, nebpe, grid%mype)
!              call MPI_ISEND( sBuf_E, ndatax, dtype, nebpe, mype, &
!                              mpi_comm_comp, sHandle(2), isend)

      end if


!
!     RECEIVE boundaries
!

! --- from NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n

          allocate( rBuf_N(0:imax,nby,1:Lm,1:km), stat = iaerr )
          rHandle(1) = mpl%f_comm%ireceive(rBuf_N, nebpe, nebpe)
!          call MPI_IRECV( rBuf_N, ndatay, dtype, nebpe, nebpe, &
!                      mpi_comm_comp, rHandle(1), irecv)

          call mpl%f_comm%wait(rHandle(1), istat)
!          call MPI_WAIT( rHandle(1), istat, ierr )

      end if

! --- from SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s

          allocate( rBuf_S(0:imax,nby,1:Lm,1:km), stat = iaerr )
          rHandle(3) = mpl%f_comm%ireceive(rBuf_S, nebpe, nebpe)
!          call MPI_IRECV( rBuf_S, ndatay, dtype, nebpe, nebpe,  &
!                       mpi_comm_comp, rHandle(3), irecv)
          call mpl%f_comm%wait(rHandle(3), istat)
!          call MPI_WAIT( rHandle(3), istat, ierr )

      end if


! --- from EAST ---

      if( itarg_e >= 0 ) then
        nebpe = itarg_e

          allocate( rBuf_E(nbx,0:jmax,1:Lm,1:km), stat = iaerr )
          rHandle(2) = mpl%f_comm%ireceive(rBuf_E, nebpe, nebpe)
!          call MPI_IRECV( rBuf_E, ndatax, dtype, nebpe, nebpe,  &
!                       mpi_comm_comp, rHandle(2), irecv)
          call mpl%f_comm%wait(rHandle(2), istat)
!          call MPI_WAIT( rHandle(2), istat, ierr )

      end if

! --- from WEST ---

      if( itarg_w >= 0 ) then
        nebpe = itarg_w

          allocate( rBuf_W(nbx,0:jmax,1:Lm,1:km), stat = iaerr )
          rHandle(4) = mpl%f_comm%ireceive(rBuf_W, nebpe, nebpe)
!          call MPI_IRECV( rBuf_W, ndatax, dtype, nebpe, nebpe,  &
!                       mpi_comm_comp, rHandle(4), irecv)
          call mpl%f_comm%wait(rHandle(4), istat)
!          call MPI_WAIT( rHandle(4), istat, ierr )

      end if

!
!                           DEALLOCATE sBufferes
!


      if( itarg_n >= 0 ) then
         call mpl%f_comm%wait(sHandle(1), istat)
!         call MPI_WAIT( sHandle(1), istat, ierr )
         deallocate( sBuf_N, stat = ierr )
      end if
      if( itarg_e >= 0 ) then
         call mpl%f_comm%wait(sHandle(2), istat)
!         call MPI_WAIT( sHandle(2), istat, ierr )
         deallocate( sBuf_E, stat = ierr )
      end if
      if( itarg_s >= 0 ) then
         call mpl%f_comm%wait(sHandle(3), istat)
!         call MPI_WAIT( sHandle(3), istat, ierr )
         deallocate( sBuf_S, stat = ierr )
      end if
      if( itarg_w >= 0 ) then
         call mpl%f_comm%wait(sHandle(4), istat)
!         call MPI_WAIT( sHandle(4), istat, ierr )
         deallocate( sBuf_W, stat = ierr )
      end if


!-----------------------------------------------------------------------
!
!                           SEND corners
!
! --- toward SOUTH-WEST ---

      if( itarg_sw >= 0 ) then
        nebpe = itarg_sw

          do L=1,Lm
          do j=1,nby
          do i=1,nbx
            sBuf_SW(i,j,L,:)= Warray(i,j,L,:)
          enddo
          enddo
          enddo

        sHandle(3) = mpl%f_comm%isend(sBuf_SW,nebpe, grid%mype)
!        call MPI_ISEND( sBuf_SW, nbxy, dtype, nebpe, mype,  &
!                        mpi_comm_comp, sHandle(3), isend)
      end if

! --- toward SOUTH-EAST ---

      if( itarg_se >= 0 ) then
        nebpe = itarg_se

          do L=1,Lm
          do j=1,nby
          do i=im-nbx,im-1
            sBuf_SE(i-(im-nbx)+1,j,L,:)=Warray(i,j,L,:)
          enddo
          enddo
          enddo

        sHandle(2) = mpl%f_comm%isend(sBuf_SE,nebpe, grid%mype) 
!        call MPI_ISEND( sBuf_SE, nbxy, dtype, nebpe, mype,  &
!                        mpi_comm_comp, sHandle(2), isend)
      end if

! --- toward NORTH-EAST ---

      if( itarg_ne >= 0 ) then
        nebpe = itarg_ne

          do L=1,Lm
          do j=jm-nby,jm-1
          do i=im-nbx,im-1
            sBuf_NE(i-(im-nbx)+1,j-(jm-nby)+1,L,:) = Warray(i,j,L,:)
          enddo
          enddo
          enddo

        sHandle(1) = mpl%f_comm%isend(sBuf_NE,nebpe, grid%mype)
!        call MPI_ISEND( sBuf_NE, nbxy, dtype, nebpe, mype,  &
!                        mpi_comm_comp, sHandle(1), isend)
      end if


! --- toward NORTH-WEST ---

      if( itarg_nw >= 0 ) then
        nebpe = itarg_nw

          do L=1,Lm
          do j=jm-nby,jm-1
          do i=1,nbx
            sBuf_NW(i,j-(jm-nby)+1,L,:) = Warray(i,j,L,:)
          enddo
          enddo
          enddo

        sHandle(4) = mpl%f_comm%isend(sBuf_NW,nebpe, grid%mype)
!       call MPI_ISEND( sBuf_NW, nbxy, dtype, nebpe, mype,  &
!                       mpi_comm_comp, sHandle(4), isend)
      end if

!
!                           RECEIVE corners
!
! --- from NORTH-EAST  ---

      if( itarg_ne >= 0 ) then
        nebpe = itarg_ne
        rHandle(1) = mpl%f_comm%ireceive(rBuf_NE, nebpe, nebpe)
        call mpl%f_comm%wait(rHandle(1), istat)
!        call MPI_IRECV( rBuf_NE, nbxy, dtype, nebpe, nebpe,  &
!                        mpi_comm_comp, rHandle(1), irecv)
!        call MPI_WAIT( rHandle(1), istat, ierr )
      end if

! --- from NORTH-WEST ---

      if( itarg_nw >= 0 ) then
        nebpe = itarg_nw
        rHandle(4) = mpl%f_comm%ireceive(rBuf_NW, nebpe, nebpe)
        call mpl%f_comm%wait(rHandle(4), istat)
!        call MPI_IRECV( rBuf_NW, nbxy, dtype, nebpe, nebpe, &
!                        mpi_comm_comp, rHandle(4), irecv)
!        call MPI_WAIT( rHandle(4), istat, ierr )
      end if

! --- from SOUTH-EAST ---

      if( itarg_se >= 0 ) then
        nebpe = itarg_se
        rHandle(2) = mpl%f_comm%ireceive(rBuf_SE, nebpe, nebpe)
        call mpl%f_comm%wait(rHandle(2), istat)
!        call MPI_IRECV( rBuf_SE, nbxy, dtype, nebpe, nebpe, &
!                        mpi_comm_comp, rHandle(2), irecv)
!        call MPI_WAIT( rHandle(2), istat, ierr )
      end if

! --- from SOUTH-WEST ---

      if(  itarg_sw >= 0 ) then
        nebpe = itarg_sw
        rHandle(3) = mpl%f_comm%ireceive(rBuf_SW, nebpe, nebpe)
        call mpl%f_comm%wait(rHandle(3), istat)
!        call MPI_IRECV( rBuf_SW, nbxy, dtype, nebpe, nebpe,  &
!                        mpi_comm_comp, rHandle(3), irecv)
!        call MPI_WAIT( rHandle(3), istat, ierr )
      end if

!
! Assign received values from NORTH, SOUTH, EAST and WEST
!

! --- from NORTH ---

   if( lnorth) then

     do L=1,Lm
     do j=1,nby
     do i=0,imax
       Warray(i,jmax+j,L,:)=Warray(i,jmax-j,L,:)
     enddo
     enddo
     enddo

   else

     do L=1,Lm
     do j=1,nby
     do i=0,imax
       Warray(i,jmax+j,L,:)=rBuf_N(i,j,L,:)
     enddo
     enddo
     enddo

   endif

! From south

   if(lsouth) then

     do L=1,Lm
     do j=1,nby
     do i=0,imax
       Warray(i,-nby-1+j,L,:)=Warray(i,nby+1-j,L,:)
     end do
     end do
     end do

   else

     do L=1,Lm
     do j=1,nby
     do i=0,imax
       Warray(i,-nby-1+j,L,:)=rBuf_S(i,j,L,:)
     enddo
     enddo
     enddo

   endif

! From west

   if(lwest) then

     do L=1,Lm
     do j=0,jmax
     do i=1,nbx
       Warray(-nbx-1+i,j,L,:)= Warray(nbx+1-i,j,L,:)
     end do
     end do
     end do

   else

     do L=1,Lm
     do j=0,jmax
     do i=1,nbx
       Warray(-nbx-1+i,j,L,:)= rBuf_W(i,j,L,:)
     enddo
     enddo
     enddo


   endif

! From east

   if(least) then

     do L=1,Lm
     do j=0,jmax
     do i=1,nbx
       Warray(imax+i,j,L,:)=Warray(imax-i,j,L,:)
     end do
     end do
     end do

   else

     do L=1,Lm
     do j=0,jmax
     do i=1,nbx
       Warray(imax+i,j,L,:)=rBuf_E(i,j,L,:)
     enddo
     enddo
     enddo

   endif

! From South-West

   if(lcorner_sw)then

     do L=1,Lm
     do j=-1,-nby,-1
     do i=-1,-nbx,-1
       Warray(i,j,L,:)=Warray(-i,-j,L,:)
     end do
     end do
     end do

   else &
   if(lwest) then

     do L=1,Lm
     do j=-nby,-1
     do i=1,nbx
        Warray(-nbx-1+i,j,L,:)= Warray(nbx+1-i,j,L,:)
     end do
     end do
     end do

   else &
   if(lsouth) then

     do L=1,Lm
     do j=1,nby
     do i=-nbx,-1
       Warray(i,-nby-1+j,L,:)=Warray(i,nby+1-j,L,:)
     end do
     end do
     end do

   else

     do L=1,Lm
     do j=1,nby
     do i=1,nbx
       Warray(-nbx-1+i,-nby-1+j,L,:)=rBuf_SW(i,j,L,:)
     end do
     end do
     end do

   end if


! From North-West

   if(lcorner_nw)then

     do L=1,Lm
     do j=1,nby
     do i=-1,-nbx,-1
       Warray(i,jmax+j,L,:)=Warray(-i,jmax-j,L,:)
     end do
     end do
     end do

   else &
   if(lwest) then

     do L=1,Lm
     do j=jmax+1,jmax+nby
     do i=1,nbx
        Warray(-nbx-1+i,j,L,:)= Warray(nbx+1-i,j,L,:)
      end do
      end do
      end do

   else &
   if(lnorth) then

     do L=1,Lm
     do j=1,nby
     do i=-nbx,-1
       Warray(i,jmax+j,L,:)=Warray(i,jmax-j,L,:)
     enddo
     enddo
     enddo

   else

     do L=1,Lm
     do j=1,nby
     do i=1,nbx
       Warray(-nbx-1+i,jmax+j,L,:)=rBuf_NW(i,j,L,:)
     end do
     end do
     end do

   end if


! From South-East

   if(lcorner_se)then

     do L=1,Lm
     do j=-1,-nby,-1
     do i=1,nbx
       Warray(imax+i,j,L,:)=Warray(imax-i,-j,L,:)
     end do
     end do
     end do

   else &
   if(least) then

     do L=1,Lm
     do j=-nby,-1
     do i=1,nbx
       Warray(imax+i,j,L,:)=Warray(imax-i,j,L,:)
     end do
     end do
     end do

   else &
   if(lsouth) then

     do L=1,Lm
     do j=1,nby
     do i=1,nbx
       Warray(imax+i,-nby-1+j,L,:)=Warray(imax+i,nby+1-j,L,:)
     end do
     end do
     end do

   else

     do L=1,Lm
     do j=1,nby
     do i=1,nbx
       Warray(imax+i,-nby-1+j,L,:)=rBuf_SE(i,j,L,:)
     end do
     end do
     end do

   end if

! From North-East


   if(lcorner_ne)then

     do L=1,Lm
     do j=1,nby
     do i=1,nbx
       Warray(imax+i,jmax+j,L,:)=Warray(imax-i,jmax-j,L,:)
     end do
     end do
     end do

   else &
   if(least) then

     do L=1,Lm
     do j=jmax+1,jmax+nby
     do i=1,nbx
       Warray(imax+i,j,L,:)=Warray(imax-i,j,L,:)
     end do
     end do
     end do

   else &
   if(lnorth) then

     do L=1,Lm
     do j=1,nby
     do i=imax+1,imax+nbx
       Warray(i,jmax+j,L,:)=Warray(i,jmax-j,L,:)
     enddo
     enddo
     enddo

   else

     do L=1,Lm
     do j=1,nby
     do i=1,nbx
       Warray(imax+i,jmax+j,L,:)=rBuf_NE(i,j,L,:)
     end do
     end do
     end do

   end if


!-----------------------------------------------------------------------
!
!                           DEALLOCATE rBufferes
!

        deallocate( rBuf_W, stat = iderr)
        deallocate( rBuf_E, stat = iderr)
        deallocate( rBuf_S, stat = iderr)
        deallocate( rBuf_N, stat = iderr)

!
!                           DEALLOCATE sBufferes
!

      if( itarg_n >= 0 ) then
        call mpl%f_comm%wait(sHandle(1), istat)
!        call MPI_WAIT( sHandle(1), istat, ierr )
        deallocate( sBuf_N, stat = ierr )
      end if
      if( itarg_e >= 0 ) then
        call mpl%f_comm%wait(sHandle(2), istat)
!        call MPI_WAIT( sHandle(2), istat, ierr )
        deallocate( sBuf_E, stat = ierr )
      end if
      if( itarg_s >= 0 ) then
        call mpl%f_comm%wait(sHandle(3), istat)
!        call MPI_WAIT( sHandle(3), istat, ierr )
        deallocate( sBuf_S, stat = ierr )
      end if
      if( itarg_w >= 0 ) then
        call mpl%f_comm%wait(sHandle(4), istat)
!        call MPI_WAIT( sHandle(4), istat, ierr )
        deallocate( sBuf_W, stat = ierr )
      end if


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

     endif FILT_GRID

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

!-----------------------------------------------------------------------
                        endsubroutine bocosHn

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine bocosHnT                             &
!***********************************************************************
!                                                                      *
!  Supply n-lines inside of domains, including edges, with halos from  *
!  the surrounding domains.  Assume mirror boundary conditions at the  *
!  boundaries of the domain                                            *
!                                                                      *
!***********************************************************************
(mpl,grid,W,km,im,jm,Lm,nbx,nby,Fimax,Fjmax,mygen_min,mygen_max)
!-----------------------------------------------------------------------
!use mg_domain1, only: grid%Flwest,grid%Fleast,grid%Flsouth,grid%Flnorth
!&
!                     ,grid%Fitarg_n,grid%Fitarg_s,grid%Fitarg_w,grid%Fitarg_e
!&
!                     ,grid%Fitarg_ne,grid%Fitarg_se,grid%Fitarg_sw,grid%Fitarg_nw
!&
!                     ,grid%Flcorner_sw,grid%Flcorner_nw,grid%Flcorner_se,grid%Flcorner_ne
!use mpi

implicit none

!-----------------------------------------------------------------------
type(mpl_type),intent(inout) :: mpl
type(mgbf_grid_type),intent(inout) :: grid
integer(kind_int), intent(in):: km,im,jm,Lm,nbx,nby,mygen_min,mygen_max
real(kind_real), dimension(-nbx:im+nbx,-nby:jm+nby,Lm,km),intent(inout):: W
integer(kind_int), dimension(grid%gm), intent(in):: Fimax,Fjmax
!-----------------------------------------------------------------------

real(kind_real), allocatable, dimension(:,:,:,:)::                         &
                                        sBuf_N,sBuf_E,sBuf_S,sBuf_W     &
                                       ,rBuf_N,rBuf_E,rBuf_S,rBuf_W
real(kind_real), dimension(0:nbx,0:nby,1:Lm,1:km)::                        &
                                        sBuf_NE,sBuf_SE,sBuf_SW,sBuf_NW &
                                       ,rBuf_NE,rBuf_SE,rBuf_SW,rBuf_NW

integer(kind_int) itarg_n,itarg_s,itarg_w,itarg_e                         &
               ,itarg_ne,itarg_se,itarg_sw,itarg_nw                     &
               ,imax,jmax
logical lwest,least,lsouth,lnorth                                       &
       ,lcorner_sw,lcorner_nw,lcorner_se,lcorner_ne

integer(kind_int) sHandle(4),rHandle(4)!,ISTAT(MPI_STATUS_SIZE)
integer(kind_int) iaerr,ierr,iderr,L,i,j
integer(kind_int) isend,irecv,nebpe,nbxy
integer(kind_int) ndatax,ndatay
logical l_sidesend
integer(kind_int) g_ind,g
type(fckit_mpi_status) :: istat

         l_sidesend=.false.

       if(mygen_min==1.and.mygen_max==1) then
         g_ind=1
         g = 1
         l_sidesend=.true.
       else &
       if(mygen_min <= grid%my_hgen .and. grid%my_hgen <= mygen_max) then
         g_ind=2
         g = grid%my_hgen
         l_sidesend=.true.
       endif


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

FILT_GRID:    if(l_sidesend) then

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

!
! from mg_domain
!
          itarg_n = grid%Fitarg_n(g_ind)
          itarg_s = grid%Fitarg_s(g_ind)
          itarg_w = grid%Fitarg_w(g_ind)
          itarg_e = grid%Fitarg_e(g_ind)

          lwest   = grid%Flwest(g_ind)
          least   = grid%Fleast(g_ind)
          lsouth  = grid%Flsouth(g_ind)
          lnorth  = grid%Flnorth(g_ind)

          itarg_ne = grid%Fitarg_ne(g_ind)
          itarg_se = grid%Fitarg_se(g_ind)
          itarg_sw = grid%Fitarg_sw(g_ind)
          itarg_nw = grid%Fitarg_nw(g_ind)

          lcorner_sw = grid%Flcorner_sw(g_ind)
          lcorner_nw = grid%Flcorner_nw(g_ind)
          lcorner_se = grid%Flcorner_se(g_ind)
          lcorner_ne = grid%Flcorner_ne(g_ind)

          if(least) then
            imax = Fimax(g)
          else
            imax = im       !   << Note that is not necesseraly im from
          endif             !      mg_parameter.  Could be also imL >>>
          if(lnorth) then
            jmax = Fjmax(g)
          else
            jmax = jm
          endif


!----------------------------------------------------------------------
      ndatax =km*(jmax+1)*(nbx+1) *Lm
      ndatay =km*(imax+1)*(nby+1) *Lm
      nbxy   =km*(nbx+1)*(nby+1)*Lm

! --- toward SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s

              allocate( sBuf_S(0:imax,0:nby,1:Lm,1:km), stat = iaerr )

              do L=Lm,1,-1
              do j=-nby,0
              do i=0,imax
                sBuf_S(i,j+nby,L,:) = W(i,j,L,:)
              enddo
              enddo
              enddo

              sHandle(3) = mpl%f_comm%isend(sBuf_S,nebpe, grid%mype)
!              call MPI_ISEND( sBuf_S, ndatay, dtype, nebpe, mype,  &
!                              mpi_comm_comp, sHandle(3), isend)
      end if

! --- toward NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n

             allocate( sBuf_N(0:imax,0:nby,1:Lm,1:km), stat = iaerr )

              do L=Lm,1,-1
              do j=0,nby
              do i=0,imax
                sBuf_N(i,j,L,:)=W(i,jmax+j,L,:)
              enddo
              enddo
              enddo

              sHandle(1) = mpl%f_comm%isend(sBuf_N,nebpe, grid%mype)
!             call MPI_ISEND( sBuf_N, ndatay, dtype, nebpe, mype,        &
!                             mpi_comm_comp, sHandle(1), isend)

      end if

! --- toward WEST ---

      if( itarg_w >= 0) then
        nebpe = itarg_w

              allocate( sBuf_W(0:nbx,0:jmax,1:Lm,1:km), stat = iaerr )

              do L=Lm,1,-1
              do j=0,jmax
              do i=-nbx,0
                sBuf_W(i+nbx,j,Lm,:) = W(i,j,Lm,:)
              enddo
              enddo
              enddo

              sHandle(4) = mpl%f_comm%isend(sBuf_W,nebpe, grid%mype)
!              call MPI_ISEND( sBuf_W, ndatax, dtype, nebpe, mype,       &
!                              mpi_comm_comp, sHandle(4), isend)

      end if

! --- toward EAST ---

      if( itarg_e >= 0 ) then
        nebpe = itarg_e

              allocate( sBuf_E(0:nbx,0:jmax,1:Lm,1:km), stat = iaerr )

              do L=Lm,1,-1
              do j=0,jmax
              do i=0,nbx
                sBuf_E(i,j,L,:) = W(imax+i,j,L,:)
              enddo
              enddo
              enddo

              sHandle(2) = mpl%f_comm%isend(sBuf_E,nebpe, grid%mype)
!              call MPI_ISEND( sBuf_E, ndatax, dtype, nebpe, mype,       &
!                              mpi_comm_comp, sHandle(2), isend)

      end if

!----------------------------------------------------------------------
!
!                           RECEIVE boundaries
!
! --- from NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n


          allocate( rBuf_N(0:imax,0:nby,1:Lm,1:km), stat = iaerr )
          rHandle(1) = mpl%f_comm%ireceive(rBuf_N,nebpe, nebpe)
          call mpl%f_comm%wait(rHandle(1), istat)
!          call MPI_IRECV( rBuf_N, ndatay, dtype, nebpe, nebpe,          &
!                      mpi_comm_comp, rHandle(1), irecv)
!          call MPI_WAIT( rHandle(1), istat, ierr )

      end if

! --- from SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s


          allocate( rBuf_S(0:imax,0:nby,1:Lm,1:km), stat = iaerr )
          rHandle(3) = mpl%f_comm%ireceive(rBuf_S,nebpe, nebpe)
          call mpl%f_comm%wait(rHandle(3), istat)
!          call MPI_IRECV( rBuf_S, ndatay, dtype, nebpe, nebpe,          &
!                       mpi_comm_comp, rHandle(3), irecv)
!          call MPI_WAIT( rHandle(3), istat, ierr )


      end if

! --- from EAST ---

      if(  itarg_e >= 0 ) then
        nebpe = itarg_e


          allocate( rBuf_E(0:nbx,0:jmax,1:Lm,1:km), stat = iaerr )
          rHandle(2) = mpl%f_comm%ireceive(rBuf_E,nebpe, nebpe)
          call mpl%f_comm%wait(rHandle(2), istat)
!          call MPI_IRECV( rBuf_E, ndatax, dtype, nebpe, nebpe,          &
!                       mpi_comm_comp, rHandle(2), irecv)
!          call MPI_WAIT( rHandle(2), istat, ierr )

      end if

! --- from WEST ---

      if(  itarg_w >= 0 ) then
        nebpe = itarg_w


          allocate( rBuf_W(0:nbx,0:jmax,1:Lm,1:km), stat = iaerr )
          rHandle(4) = mpl%f_comm%ireceive(rBuf_W,nebpe, nebpe)
          call mpl%f_comm%wait(rHandle(4), istat)
!          call MPI_IRECV( rBuf_W, ndatax, dtype, nebpe, nebpe,          &
!                       mpi_comm_comp, rHandle(4), irecv)
!          call MPI_WAIT( rHandle(4), istat, ierr )


      end if



      if( itarg_w >= 0 ) then
         call mpl%f_comm%wait(sHandle(4), istat)
!         call MPI_WAIT( sHandle(4), istat, ierr )
         deallocate( sBuf_W, stat = ierr )
      end if
      if( itarg_e >= 0 ) then
         call mpl%f_comm%wait(sHandle(2), istat)
!         call MPI_WAIT( sHandle(2), istat, ierr )
         deallocate( sBuf_E, stat = ierr )
      end if
      if( itarg_s >= 0 ) then
         call mpl%f_comm%wait(sHandle(3), istat)
!         call MPI_WAIT( sHandle(3), istat, ierr )
         deallocate( sBuf_S, stat = ierr )
      end if
      if( itarg_n >= 0 ) then
         call mpl%f_comm%wait(sHandle(1), istat)
!         call MPI_WAIT( sHandle(1), istat, ierr )
         deallocate( sBuf_N, stat = ierr )
      end if

!-----------------------------------------------------------------------
!
!                           SEND corners
!
! --- toward SOUTH-WEST ---

      if(  itarg_sw >= 0 ) then
        nebpe = itarg_sw

          do L=Lm,1,-1
          do j=-nby,0
          do i=-nbx,0
            sBuf_SW(nbx+i,nby+j,L,:)= W(i,j,L,:)
          enddo
          enddo
          enddo

        sHandle(3) = mpl%f_comm%isend(sBuf_SW,nebpe, grid%mype)
!        call MPI_ISEND( sBuf_SW, nbxy, dtype, nebpe, mype,              &
!                        mpi_comm_comp, sHandle(3), isend)
      end if

! --- toward SOUTH-EAST ---

      if( itarg_se >= 0 ) then
        nebpe = itarg_se

          do L=Lm,1,-1
          do j=-nby,0
          do i=imax,imax+nbx
            sBuf_SE(i-imax,nby+j,L,:)=W(i,j,L,:)
          enddo
          enddo
          enddo

        sHandle(2) = mpl%f_comm%isend(sBuf_SE,nebpe, grid%mype)
!        call MPI_ISEND( sBuf_SE, nbxy, dtype, nebpe, mype,              &
!                       mpi_comm_comp, sHandle(2), isend)
      end if

! --- toward NORTH-WEST ---

      if( itarg_nw >= 0 ) then
        nebpe = itarg_nw

          do L=Lm,1,-1
          do j=jmax,jmax+nby
          do i=-nbx,0
            sBuf_NW(nbx+i,j-jmax,L,:) = W(i,j,L,:)
          enddo
          enddo
          enddo

        sHandle(4) = mpl%f_comm%isend(sBuf_NW,nebpe, grid%mype)
!       call MPI_ISEND( sBuf_NW, nbxy, dtype, nebpe, mype,               &
!                       mpi_comm_comp, sHandle(4), isend)
      end if

! --- toward NORTH-EAST ---

      if( itarg_ne >= 0 ) then
        nebpe = itarg_ne

          do L=Lm,1,-1
          do j=jmax,jmax+nby
          do i=imax,imax+nbx
            sBuf_NE(i-imax,j-jmax,L,:) = W(i,j,L,:)
          enddo
          enddo
          enddo

        sHandle(1) = mpl%f_comm%isend(sBuf_NE,nebpe, grid%mype)
!        call MPI_ISEND( sBuf_NE, nbxy, dtype, nebpe, mype,              &
!                        mpi_comm_comp, sHandle(1), isend)
      end if


!
!                           RECEIVE corners
!
! --- from NORTH-EAST  ---

      if( itarg_ne >= 0 ) then
        nebpe = itarg_ne
        rHandle(1) = mpl%f_comm%ireceive(rBuf_NE, nebpe, nebpe)
        call mpl%f_comm%wait(rHandle(1), istat)
!        call MPI_IRECV( rBuf_NE, nbxy, dtype, nebpe, nebpe,             &
!                      mpi_comm_comp, rHandle(1), irecv)
!        call MPI_WAIT( rHandle(1), istat, ierr )
      end if

! --- from NORTH-WEST ---

      if( itarg_nw >= 0 ) then
        nebpe = itarg_nw
        rHandle(4) = mpl%f_comm%ireceive(rBuf_NW, nebpe, nebpe)
        call mpl%f_comm%wait(rHandle(4), istat)
!        call MPI_IRECV( rBuf_NW, nbxy, dtype, nebpe, nebpe,             &
!                       mpi_comm_comp, rHandle(4), irecv)
!        call MPI_WAIT( rHandle(4), istat, ierr )
      end if

! --- from SOUTH-EAST ---

      if( itarg_se >= 0 ) then
        nebpe = itarg_se
        rHandle(2) = mpl%f_comm%ireceive(rBuf_SE, nebpe, nebpe)
        call mpl%f_comm%wait(rHandle(2), istat)
!        call MPI_IRECV( rBuf_SE, nbxy, dtype, nebpe, nebpe,             &
!                       mpi_comm_comp, rHandle(2), irecv)
!        call MPI_WAIT( rHandle(2), istat, ierr )
      end if

! --- from SOUTH-WEST ---

      if( itarg_sw >= 0 ) then
        nebpe = itarg_sw
        rHandle(3) = mpl%f_comm%ireceive(rBuf_SW, nebpe, nebpe)
        call mpl%f_comm%wait(rHandle(3), istat)
!        call MPI_IRECV( rBuf_SW, nbxy, dtype, nebpe, nebpe,             &
!                       mpi_comm_comp, rHandle(3), irecv)
!                       mpi_comm_comp, rHandle(3), irecv)
!        call MPI_WAIT( rHandle(3), istat, ierr )
      end if


!
! Assign received values
!

! From west

   if(lwest) then
     do L=1,lm
     do j=0,jmax
     do i=0,nbx
       W(i,j,L,:)= W(i,j,L,:)+W(-i,j,L,:)
     end do
     end do
     end do
   else
     do L=1,lm
     do j=0,jmax
     do i=0,nbx
      W(i,j,L,:)= W(i,j,L,:)+rBuf_W(i,j,L,:)
     end do
     end do
     end do
   endif

! From east

   if(least) then
     do L=1,lm
     do j=0,jmax
     do i=0,nbx
       W(imax-nbx+i,j,L,:)= W(imax-nbx+i,j,L,:)+W(imax+nbx-i,j,L,:)
     end do
     end do
     end do
   else
     do L=1,lm
     do j=0,jmax
     do i=0,nbx
       W(imax-nbx+i,j,L,:)= W(imax-nbx+i,j,L,:)+rBuf_E(i,j,L,:)
     end do
     end do
     end do
   endif

! From south

   if(lsouth) then
     do L=1,lm
     do j=0,nby
     do i=0,imax
       W(i,j,L,:)= W(i,j,L,:)+W(i,-j,L,:)
     end do
     end do
     end do
   else
     do L=1,lm
     do j=0,nby
     do i=0,imax
       W(i,j,L,:)= W(i,j,L,:)+rBuf_S(i,j,L,:)
     end do
     end do
     end do
   endif

!  From north

   if(lnorth) then
     do L=1,lm
     do j=0,nby
     do i=0,imax
       W(i,jmax-nby+j,L,:)= W(i,jmax-nby+j,L,:)+W(i,jmax+nby-j,L,:)
     enddo
     enddo
     enddo
   else
     do L=1,lm
     do j=0,nby
     do i=0,imax
       W(i,jmax-nby+j,L,:)= W(i,jmax-nby+j,L,:)+rBuf_N(i,j,L,:)
     enddo
     enddo
     enddo
   endif

! From South-West


   if(lcorner_sw) then

     do L=1,lm
     do j=1,nby
     do i=0,nbx
       W(i,j,L,:)= W(i,j,L,:)+W(-i,-j,L,:)
     enddo
     enddo
     enddo

   else &

   if(lwest) then
     do L=1,lm
     do j=0,nby
     do i=0,nbx
       W(i,j,L,:)=W(i,j,L,:)+rBuf_SW(nbx-i,j,L,:)
     enddo
     enddo
     enddo

   else &
   if(lsouth) then
     do L=1,lm
     do j=1,nby
     do i=0,nbx
       W(i,j,L,:)= W(i,j,L,:)+rBuf_SW(i,nby-j,L,:)
     enddo
     enddo
     enddo

   else

     do L=1,lm
     do j=0,nby
     do i=0,nbx
       W(i,j,L,:)= W(i,j,L,:)+ rBuf_SW(i,j,L,:)
     end do
     end do
     end do

    end if

! From South-East

    if(lcorner_se) then

     do L=1,lm
     do j=1,nby
     do i=0,nbx
       W(imax-i,j,L,:)=W(imax-i,j,L,:)+W(imax+i,-j,L,:)
     enddo
     enddo
     enddo

    else &
    if(least) then

      do L=1,lm
      do j=0,nby
      do i=0,nbx
        W(imax-i,j,L,:)=W(imax-i,j,L,:)+rBuf_SE(i,j,L,:)
      enddo
      enddo
      enddo

    else &
    if(lsouth) then

      do L=1,lm
      do j=1,nby
      do i=0,nbx
        W(imax-nby+i,j,L,:)=W(imax-nby+i,j,L,:)+rBuf_SE(i,nby-j,L,:)
      enddo
      enddo
      enddo

    else

      do L=1,lm
      do j=0,nby
      do i=0,nbx
        W(imax-nbx+i,j,L,:)= W(imax-nbx+i,j,L,:)+rBuf_SE(i,j,L,:)
      end do
      end do
      end do

    end if

! From North-West

    if(lcorner_nw) then

      do L=1,lm
      do j=1,nby
      do i=0,nbx
        W(i,jmax-j,L,:)=W(i,jmax-j,L,:)+W(-i,jmax+j,L,:)
      end do
      end do
      end do

    else &
    if(lwest) then

      do L=1,lm
      do j=0,nby
      do i=0,nbx
        W(i,jmax-nby+j,L,:)=W(i,jmax-nby+j,L,:)+rBuf_NW(nbx-i,j,L,:)
      end do
      end do
      end do

    else &
    if(lnorth) then

      do L=1,lm
      do j=0,nby-1
      do i=0,nbx
        W(i,jmax-nby+j,L,:)=W(i,jmax-nby+j,L,:)+rBuf_NW(i,nby-j,L,:)
      end do
      end do
      end do

     else

      do L=1,lm
      do j=0,nby
      do i=0,nbx
        W(i,jmax-nby+j,L,:)= W(i,jmax-nby+j,L,:)+rBuf_NW(i,j,L,:)
      end do
      end do
      end do

    end if


! From North-East

    if(lcorner_ne) then

      do L=1,lm
      do j=1,nby
      do i=0,nbx
        W(imax-i,jmax-j,L,:)=W(imax-i,jmax-j,L,:)+W(imax+i,jmax+j,L,:)
      end do
      end do
      end do

    else &
    if(least) then

      do L=1,lm
      do j=0,nby
      do i=0,nbx
        W(imax-i,jmax-nby+j,L,:)=W(imax-i,jmax-nby+j,L,:)+rBuf_NE(i,j,L,:)
      end do
      end do
      end do

    else &
    if(lnorth) then

      do L=1,lm
      do j=0,nby-1
      do i=0,nbx
        W(imax-nbx+i,jmax-nby+j,L,:)=W(imax-nbx+i,jmax-nby+j,L,:)+rBuf_NE(i,nby-j,L,:)
      end do
      end do
      end do

     else

      do L=1,lm
      do j=0,nby
      do i=0,nbx
        W(imax-nbx+i,jmax-nby+j,L,:)=W(imax-nbx+i,jmax-nby+j,L,:)+rBuf_NE(i,j,L,:)
      end do
      end do
      end do

    end if

        deallocate( rBuf_W, stat = iderr)
        deallocate( rBuf_E, stat = iderr)
        deallocate( rBuf_S, stat = iderr)
        deallocate( rBuf_N, stat = iderr)

!-----------------------------------------------------------------------
!
!                           DEALLOCATE sBufferes

      if( itarg_w  >= 0 ) then
         call mpl%f_comm%wait(sHandle(4), istat)
!         call MPI_WAIT( sHandle(4), istat, ierr )
      end if
      if( itarg_e  >= 0 ) then
         call mpl%f_comm%wait(sHandle(2), istat)
!         call MPI_WAIT( sHandle(2), istat, ierr )
      end if
      if( itarg_s  >= 0 ) then
         call mpl%f_comm%wait(sHandle(3), istat)
!         call MPI_WAIT( sHandle(3), istat, ierr )
      end if
      if( itarg_n  >= 0 ) then
         call mpl%f_comm%wait(sHandle(1), istat)
!         call MPI_WAIT( sHandle(1), istat, ierr )
      end if


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

     endif FILT_GRID

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

!-----------------------------------------------------------------------
                        endsubroutine bocosHnT

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine upsend                               &
!***********************************************************************
!                                                                      *
!         Upsending data from the one resolution pes                   *
!         to the next low-resolution pes                               *
!                                                                      *
!***********************************************************************
(mpl,grid,Harray,Warray,km,Lm,mygen_dn,mygen_up)
!-----------------------------------------------------------------------
!use mg_parameter1, only: im,jm,imL,jmL,hx,hy
!use mg_domain1, only: Flsendup_sw,Flsendup_se,Flsendup_nw,Flsendup_ne   &
!                     ,Fitarg_up                                         &
!                     ,itargdn_sw,itargdn_se,itargdn_ne,itargdn_nw
!use mpi

implicit none

!-----------------------------------------------------------------------
type(mpl_type),intent(inout) :: mpl
type(mgbf_grid_type),intent(inout) :: grid
integer(kind_int), intent(in):: km,Lm
real(kind_real), dimension(0:grid%imL,0:grid%jmL,Lm,km),intent(in):: Harray
real(kind_real), dimension(-grid%hx:grid%im+grid%hx,-grid%hy:grid%jm+grid%hy,Lm,km),intent(out):: Warray
integer(kind_int),intent(in):: mygen_dn,mygen_up

!-----------------------------------------------------------------------
real(kind_real), allocatable, dimension(:,:,:,:)::                          &
                                         sBuf_SW,sBuf_SE,sBuf_NW,sBuf_NE &
                                        ,rBuf_SW,rBuf_SE,rBuf_NW,rBuf_NE

real(kind_real),dimension(0:grid%imL,0:grid%jmL,1:Lm,1:km):: dBuf_SW
real(kind_real),dimension(0:grid%imL,0:grid%jmL,1:Lm,1:km):: dBuf_SE
real(kind_real),dimension(0:grid%imL,0:grid%jmL,1:Lm,1:km):: dBuf_NW
real(kind_real),dimension(0:grid%imL,0:grid%jmL,1:Lm,1:km):: dBuf_NE

integer(kind_int) sHandle(4),rHandle(4)!,ISTAT(MPI_STATUS_SIZE)
integer(kind_int) iaerr,ierr,iderr,ndata,i,j,L
integer(kind_int) isend,irecv,nebpe

logical:: lsendup_sw,lsendup_se,lsendup_nw,lsendup_ne,flag_up
integer(kind_int):: itarg_up
integer:: g_ind
type(fckit_mpi_status) :: istat

!-----------------------------------------------------------------------
!
! Define generational flags
!

     if(mygen_dn==1) then
         g_ind=1
       lsendup_sw=grid%Flsendup_sw(g_ind)
       lsendup_se=grid%Flsendup_se(g_ind)
       lsendup_nw=grid%Flsendup_nw(g_ind)
       lsendup_ne=grid%Flsendup_ne(g_ind)
     else
         g_ind=2
       lsendup_sw=grid%Flsendup_sw(g_ind).and.(grid%my_hgen==mygen_dn)
       lsendup_se=grid%Flsendup_se(g_ind).and.(grid%my_hgen==mygen_dn)
       lsendup_nw=grid%Flsendup_nw(g_ind).and.(grid%my_hgen==mygen_dn)
       lsendup_ne=grid%Flsendup_ne(g_ind).and.(grid%my_hgen==mygen_dn)
     endif


       itarg_up=grid%Fitarg_up(g_ind)


!-----------------------------------------------------------------------

   if(grid%my_hgen==mygen_up) then
      Warray(:,:,:,:)=0.
   endif

     ndata =km*(grid%imL+1)*(grid%jmL+1)*Lm


      if(  lsendup_sw ) then

        nebpe = itarg_up

        if(nebpe == grid%mype) then

             do L=1,Lm
             do j=0,grid%jmL
             do i=0,grid%imL
                dBuf_SW(i,j,L,:) = Harray(i,j,L,:)
             enddo
             enddo
             enddo

        else

        allocate( sBuf_SW(0:grid%imL,0:grid%jmL,1:Lm,1:km), stat = iaerr )

             do L=1,Lm
             do j=0,grid%jmL
             do i=0,grid%imL
                sBuf_SW(i,j,L,:) = Harray(i,j,L,:)
             enddo
             enddo
             enddo

        sHandle(1) = mpl%f_comm%isend(sBuf_SW,nebpe, grid%mype)
        call mpl%f_comm%wait(sHandle(1), istat)
!        call MPI_ISEND( sBuf_SW, ndata, dtype, nebpe, mype,  &
!                       mpi_comm_comp, sHandle(1), isend)
!        call MPI_WAIT( sHandle(1), istat, ierr )

        deallocate( sBuf_SW, stat = ierr )

        endif

      end if

!
! --- Receive SW portion of data at higher generation
!

      if( grid%my_hgen==mygen_up .and. grid%itargdn_sw >= 0 ) then

        nebpe = grid%itargdn_sw

        if(nebpe == grid%mype) then

             do L=1,Lm
             do j=0,grid%jmL
             do i=0,grid%imL
               Warray(i,j,L,:)=dBuf_SW(i,j,L,:)
             enddo
             enddo
             enddo

        else

        allocate( rBuf_SW(0:grid%imL,0:grid%jmL,1:Lm,1:km), stat = iaerr )

        rHandle(1) = mpl%f_comm%ireceive(rBuf_SW,nebpe,nebpe)
        call mpl%f_comm%wait(rHandle(1), istat)
!        call MPI_IRECV( rBuf_SW, ndata, dtype, nebpe, nebpe, &
!                       mpi_comm_comp, rHandle(1), irecv)
!        call MPI_WAIT( rHandle(1), istat, ierr )

             do L=1,Lm
             do j=0,grid%jmL
             do i=0,grid%imL
               Warray(i,j,L,:)=Rbuf_SW(i,j,L,:)
             enddo
             enddo
             enddo

        endif

      endif

      call mpl%f_comm%barrier()
!      call barrierMPI

!
! --- Send data to SE portion of processors at higher generation
!

      if( lsendup_se ) then
        nebpe = itarg_up

        if(nebpe == grid%mype) then

             do L=1,Lm
             do j=0,grid%jmL
             do i=0,grid%imL
                dBuf_SE(i,j,L,:) = Harray(i,j,L,:)
             enddo
             enddo
             enddo

        else

        allocate( sBuf_SE(0:grid%imL,0:grid%jmL,1:Lm,1:km), stat = iaerr )

             do L=1,Lm
             do j=0,grid%jmL
             do i=0,grid%imL
               sBuf_SE(i,j,L,:) = Harray(i,j,L,:)
             enddo
             enddo
             enddo

        sHandle(2) = mpl%f_comm%isend(sBuf_SE,nebpe, grid%mype)
        call mpl%f_comm%wait(sHandle(2), istat)
!        call MPI_ISEND( sBuf_SE, ndata, dtype, nebpe, mype, &
!                       mpi_comm_comp, sHandle(2), isend)

!        call MPI_WAIT( sHandle(2), istat, ierr )

        deallocate( sBuf_SE, stat = ierr )

        endif

      end if

!
! --- Receive SE portion of data at higher generation


      if( grid%my_hgen==mygen_up .and. grid%itargdn_se >= 0 ) then
        nebpe = grid%itargdn_se

        if(nebpe == grid%mype) then

             do L=1,Lm
             do j=0,grid%jmL
             do i=0,grid%imL
               Warray(grid%imL+i,j,L,:)=dBuf_SE(i,j,L,:)
             enddo
             enddo
             enddo

        else

        allocate( rBuf_SE(0:grid%imL,0:grid%jmL,1:Lm,1:km), stat = iaerr )

        rHandle(2) = mpl%f_comm%ireceive(rBuf_SE,nebpe,nebpe)
        call mpl%f_comm%wait(rHandle(2), istat)
!        call MPI_IRECV( rBuf_SE, ndata, dtype, nebpe, nebpe,  &
!                       mpi_comm_comp, rHandle(2), irecv)
!        call MPI_WAIT( rHandle(2), istat, ierr )

             do L=1,Lm
             do j=0,grid%jmL
             do i=0,grid%imL
               Warray(grid%imL+i,j,L,:)=Rbuf_SE(i,j,L,:)
             enddo
             enddo
             enddo



        endif

      endif

      call mpl%f_comm%barrier()
!      call barrierMPI

!
! --- Send data to NW portion of processors at higher generation
!

      if( lsendup_nw ) then
        nebpe = itarg_up

        if(nebpe == grid%mype) then

             do L=1,Lm
             do j=0,grid%jmL
             do i=0,grid%imL
               dBuf_NW(i,j,L,:) = Harray(i,j,L,:)
             enddo
             enddo
             enddo

        else

        allocate( sBuf_NW(0:grid%imL,0:grid%jmL,1:Lm,1:km), stat = iaerr )

             do L=1,Lm
             do j=0,grid%jmL
             do i=0,grid%imL
               sBuf_NW(i,j,L,:) = Harray(i,j,L,:)
             enddo
             enddo
             enddo

        sHandle(3) = mpl%f_comm%isend(sBuf_NW,nebpe, grid%mype)
        call mpl%f_comm%wait(sHandle(3), istat)
!         call MPI_ISEND( sBuf_NW, ndata, dtype, nebpe, mype,  &
!                        mpi_comm_comp, sHandle(3), isend)

!         call MPI_WAIT( sHandle(3), istat, ierr )

         deallocate( sBuf_NW, stat = ierr )

      end if

    end if

!
! --- Receive NW portion of data at higher generation
!

      if( grid%my_hgen==mygen_up .and. grid%itargdn_nw >= 0 ) then
        nebpe = grid%itargdn_nw

        if(nebpe == grid%mype) then

             do L=1,Lm
             do j=0,grid%jmL
             do i=0,grid%imL
               Warray(i,grid%jmL+j,L,:)=dBuf_NW(i,j,L,:)
             enddo
             enddo
             enddo

        else

        allocate( rBuf_NW(0:grid%imL,0:grid%jmL,1:Lm,1:km), stat = iaerr )

        rHandle(3) = mpl%f_comm%ireceive(rBuf_NW,nebpe,nebpe)
        call mpl%f_comm%wait(rHandle(3), istat)
!        call MPI_IRECV( rBuf_NW, ndata, dtype, nebpe, nebpe,  &
!                       mpi_comm_comp, rHandle(3), irecv)

!        call MPI_WAIT( rHandle(3), istat, ierr )

             do L=1,Lm
             do j=0,grid%jmL
             do i=0,grid%imL
               Warray(i,grid%jmL+j,L,:)=rBuf_NW(i,j,L,:)
             enddo
             enddo
             enddo

        deallocate( rBuf_NW, stat = iderr)

        end if

      end if

      call mpl%f_comm%barrier()
!      call barrierMPI
!
! --- Send data to NE portion of processors at higher generation
!

      if( lsendup_ne ) then
        nebpe = itarg_up

        if(nebpe == grid%mype) then

             do L=1,Lm
             do j=0,grid%jmL
             do i=0,grid%imL
               dBuf_NE(i,j,L,:) = Harray(i,j,L,:)
             enddo
             enddo
             enddo

        else

        allocate( sBuf_NE(0:grid%imL,0:grid%jmL,1:Lm,1:km), stat = iaerr )

             do L=1,Lm
             do j=0,grid%jmL
             do i=0,grid%imL
               sBuf_NE(i,j,L,:) = Harray(i,j,L,:)
             enddo
             enddo
             enddo

        sHandle(4) = mpl%f_comm%isend(sBuf_NE,nebpe, grid%mype)
        call mpl%f_comm%wait(sHandle(4), istat)
!        call MPI_ISEND( sBuf_NE, ndata, dtype, nebpe, mype, &
!                       mpi_comm_comp, sHandle(4), isend)

!         call MPI_WAIT( sHandle(4), istat, ierr )

         deallocate( sBuf_NE, stat = ierr )

        endif

      end if

!
! --- Receive NE portion of data at higher generation
!

      if( grid%my_hgen==mygen_up .and. grid%itargdn_ne >= 0 ) then
        nebpe = grid%itargdn_ne

        if(nebpe == grid%mype) then

             do L=1,Lm
             do j=0,grid%jmL
             do i=0,grid%imL
               Warray(grid%imL+i,grid%jmL+j,L,:)=dBuf_NE(i,j,L,:)
             enddo
             enddo
             enddo

        else

        allocate( rBuf_NE(0:grid%imL,0:grid%jmL,1:Lm,1:km), stat = iaerr )

        rHandle(4) = mpl%f_comm%ireceive(rBuf_NE,nebpe,nebpe)
        call mpl%f_comm%wait(rHandle(4), istat)
!        call MPI_IRECV( rBuf_NE, ndata, dtype, nebpe, nebpe,  &
!                       mpi_comm_comp, rHandle(4), irecv)

!        call MPI_WAIT( rHandle(4), istat, ierr )

             do L=1,Lm
             do j=0,grid%jmL
             do i=0,grid%imL
               Warray(grid%imL+i,grid%jmL+j,L,:)=rBuf_NE(i,j,L,:)
             enddo
             enddo
             enddo

          deallocate( rBuf_NE, stat = iderr)

        endif
      endif

      call mpl%f_comm%barrier()
!      call barrierMPI

!-----------------------------------------------------------------------
                        endsubroutine upsend

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine downsend                             &
!***********************************************************************
!                                                                      *
!         Downsending data from low resolution pes    (mygen_up)       *
!         to the concurent high-resolution pes        (mygen_dn)       *
!         and add the existing and the recevied values                 *
!                       (MPI version)                                  *
!                                                                      *
!***********************************************************************
(mpl,grid,Warray,Harray,km,Lm,mygen_up,mygen_dn)
!-----------------------------------------------------------------------
!use mg_parameter1, only: im,jm,imL,jmL,hx,hy
!use mg_domain1, only: Flsendup_sw,Flsendup_se,Flsendup_nw,Flsendup_ne    &
!                     ,Fitarg_up                                          &
!                     ,itargdn_sw,itargdn_se,itargdn_ne,itargdn_nw
!use mpi

implicit none
!-----------------------------------------------------------------------
type(mpl_type),intent(inout) :: mpl
type(mgbf_grid_type),intent(inout) :: grid
integer(kind_int), intent(in)::km,Lm
real(kind_real), dimension(0:grid%im,0:grid%jm,Lm,1:km),intent(in):: Warray
real(kind_real), dimension(-1:grid%imL+1,-1:grid%jmL+1,Lm,1:km),intent(out):: Harray
integer, intent(in):: mygen_up,mygen_dn
!-----------------------------------------------------------------------
real(kind_real), allocatable, dimension(:,:,:,:)::                          &
                            sBuf_SW,sBuf_SE,sBuf_NW,sBuf_NE              &
                           ,rBuf_SW,rBuf_SE,rBuf_NW,rBuf_NE

real(kind_real),dimension(0:grid%imL,0:grid%jmL,1:Lm,1:km):: dBuf_SW
real(kind_real),dimension(0:grid%imL,0:grid%jmL,1:Lm,1:km):: dBuf_SE
real(kind_real),dimension(0:grid%imL,0:grid%jmL,1:Lm,1:km):: dBuf_NW
real(kind_real),dimension(0:grid%imL,0:grid%jmL,1:Lm,1:km):: dBuf_NE

integer(kind_int) sHandle(4),rHandle(4)!,ISTAT(MPI_STATUS_SIZE)
integer(kind_int) iaerr,ierr,iderr,ndata,i,j,L
integer(kind_int) isend,irecv,nebpe

logical:: lsendup_sw,lsendup_se,lsendup_nw,lsendup_ne
integer(kind_int):: itarg_up
integer(kind_int):: g_ind
type(fckit_mpi_status) :: istat
!-----------------------------------------------------------------------
!
! Define generational flags
!

     if(mygen_dn==1) then
         g_ind=1
       lsendup_sw=grid%Flsendup_sw(g_ind)
       lsendup_se=grid%Flsendup_se(g_ind)
       lsendup_nw=grid%Flsendup_nw(g_ind)
       lsendup_ne=grid%Flsendup_ne(g_ind)
     else
         g_ind=2
       lsendup_sw=grid%Flsendup_sw(g_ind).and.(grid%my_hgen==mygen_dn)
       lsendup_se=grid%Flsendup_se(g_ind).and.(grid%my_hgen==mygen_dn)
       lsendup_nw=grid%Flsendup_nw(g_ind).and.(grid%my_hgen==mygen_dn)
       lsendup_ne=grid%Flsendup_ne(g_ind).and.(grid%my_hgen==mygen_dn)
     endif

       itarg_up=grid%Fitarg_up(g_ind)

!
      ndata =km*(grid%imL+1)*(grid%jmL+1)*Lm

!
! --- Send data from SW portion of processors at the higher generation
!     to corresponding  PE's at lower generation


  if(grid%my_hgen==mygen_up .and. grid%itargdn_sw >= 0 ) then
        nebpe = grid%itargdn_sw

        if(nebpe == grid%mype) then

             do L=1,Lm
             do j=0,grid%jmL
             do i=0,grid%imL
                dBuf_SW(i,j,L,:) = Warray(i,j,L,:)
             enddo
             enddo
             enddo

        else

        allocate( sBuf_SW(0:grid%imL,0:grid%jmL,1:Lm,1:km), stat = iaerr )

             do L=1,Lm
             do j=0,grid%jmL
             do i=0,grid%imL
                sBuf_SW(i,j,L,:) = Warray(i,j,L,:)
             enddo
             enddo
             enddo

        sHandle(1) = mpl%f_comm%isend(sBuf_SW,nebpe, grid%mype)
        call mpl%f_comm%wait(sHandle(1), istat)
!        call MPI_ISEND( sBuf_SW, ndata, dtype, nebpe, grid%mype,  &
!                        mpi_comm_comp, sHandle(1), isend)
!        call MPI_WAIT( sHandle(1), istat, ierr )
        deallocate( sBuf_SW, stat = ierr )

        endif

  endif
!
! --- Receive SW portion of data at lower generation


      if( lsendup_sw ) then

        nebpe = itarg_up

        if(nebpe == grid%mype) then

             do L=1,Lm
             do j=0,grid%jmL
             do i=0,grid%imL
               Harray(i,j,L,:)=dBuf_SW(i,j,L,:)
             enddo
             enddo
             enddo

        else


        allocate( rBuf_SW(0:grid%imL,0:grid%jmL,Lm,1:km), stat = iaerr )

        rHandle(1) = mpl%f_comm%ireceive(rBuf_SW,nebpe, nebpe)
        call mpl%f_comm%wait(rHandle(1), istat)
!        call MPI_IRECV( rBuf_SW, ndata, dtype, nebpe, nebpe, &
!                        mpi_comm_comp, rHandle(1), irecv)
!        call MPI_WAIT( rHandle(1), istat, ierr )

             do L=1,Lm
             do j=0,grid%jmL
             do i=0,grid%imL
               Harray(i,j,L,:)=Rbuf_SW(i,j,L,:)
             enddo
             enddo
             enddo

        deallocate( rBuf_SW, stat = iderr)

        endif

      endif

      call mpl%f_comm%barrier()
!      call barrierMPI

!
! --- Send data from SE portion of processors at the higher generation
!     to corresponding  PE's at lower generation

  if(grid%my_hgen==mygen_up .and.  grid%itargdn_se >= 0 ) then
        nebpe = grid%itargdn_se


        if(nebpe == grid%mype) then

             do L=1,Lm
             do j=0,grid%jmL
             do i=0,grid%imL
                dBuf_SE(i,j,L,:) = Warray(grid%imL+i,j,L,:)
             enddo
             enddo
             enddo

        else

        allocate( sBuf_SE(0:grid%imL,0:grid%jmL,1:Lm,1:km), stat = iaerr )

             do L=1,Lm
             do j=0,grid%jmL
             do i=0,grid%imL
               sBuf_SE(i,j,L,:) = Warray(grid%imL+i,j,L,:)
             enddo
             enddo
             enddo

        sHandle(2) = mpl%f_comm%isend(sBuf_SE,nebpe, grid%mype)
        call mpl%f_comm%wait(sHandle(2), istat)
!        call MPI_ISEND( sBuf_SE, ndata, dtype, nebpe, grid%mype,  &
!                       mpi_comm_comp, sHandle(2), isend)
!        call MPI_WAIT( sHandle(2), istat, ierr )
        deallocate( sBuf_SE, stat = ierr )

        endif

  endif
!
! --- Receive SE portion of data at lower generation


      if( lsendup_se ) then
        nebpe = itarg_up

        if(nebpe == grid%mype) then

             do L=1,Lm
             do j=0,grid%jmL
             do i=0,grid%imL
               Harray(i,j,L,:)=dBuf_SE(i,j,L,:)
             enddo
             enddo
             enddo

        else

        allocate( rBuf_SE(0:grid%imL,0:grid%jmL,Lm,1:km), stat = iaerr )

        rHandle(2) = mpl%f_comm%ireceive(rBuf_SE,nebpe, nebpe)
        call mpl%f_comm%wait(rHandle(2), istat)
!        call MPI_IRECV( rBuf_SE, ndata, dtype, nebpe, nebpe, &
!                        mpi_comm_comp, rHandle(2), irecv)
!        call MPI_WAIT( rHandle(2), istat, ierr )

             do L=1,Lm
             do j=0,grid%jmL
             do i=0,grid%imL
               Harray(i,j,L,:)=Rbuf_SE(i,j,L,:)
             enddo
             enddo
             enddo

       deallocate( rBuf_SE, stat = iderr)

       endif

     end if

      call mpl%f_comm%barrier()
!     call barrierMPI

! --- Send data from NW portion of processors at the higher generation
!     to corresponding  PE's at lower generantion

  if(grid%my_hgen==mygen_up .and. grid%itargdn_nw >= 0 ) then
        nebpe = grid%itargdn_nw


        if(nebpe == grid%mype) then

             do L=1,Lm
             do j=0,grid%jmL
             do i=0,grid%imL
                dBuf_NW(i,j,L,:) = Warray(i,grid%jmL+j,L,:)
             enddo
             enddo
             enddo

        else


        allocate( sBuf_NW(0:grid%imL,0:grid%jmL,1:Lm,1:km), stat = iaerr )

             do L=1,Lm
             do j=0,grid%jmL
             do i=0,grid%imL
                sBuf_NW(i,j,L,:) = Warray(i,grid%jmL+j,L,:)
             enddo
             enddo
             enddo

        sHandle(3) = mpl%f_comm%isend(sBuf_NW,nebpe, grid%mype)
        call mpl%f_comm%wait(sHandle(3), istat)
!        call MPI_ISEND( sBuf_NW, ndata, dtype, nebpe, grid%mype,  &
!                        mpi_comm_comp, sHandle(3), isend)
!        call MPI_WAIT( sHandle(3), istat, ierr )
        deallocate( sBuf_NW, stat = ierr )

        endif

  endif
!
! --- Receive NW portion of data at lower generation


      if( lsendup_nw ) then

        nebpe = itarg_up

        if(nebpe == grid%mype) then

             do L=1,Lm
             do j=0,grid%jmL
             do i=0,grid%imL
               Harray(i,j,L,:)=dBuf_NW(i,j,L,:)
             enddo
             enddo
             enddo

        else

        allocate( rBuf_NW(0:grid%imL,0:grid%jmL,1:Lm,1:km), stat = iaerr )

        rHandle(3) = mpl%f_comm%ireceive(rBuf_NW,nebpe, nebpe)
        call mpl%f_comm%wait(rHandle(3), istat)
!        call MPI_IRECV( rBuf_NW, ndata, dtype, nebpe, nebpe, &
!                       mpi_comm_comp, rHandle(3), irecv)
!        call MPI_WAIT( rHandle(3), istat, ierr )

             do L=1,Lm
             do j=0,grid%jmL
             do i=0,grid%imL
               Harray(i,j,L,:)=Rbuf_NW(i,j,L,:)
             enddo
             enddo
             enddo

        deallocate( rBuf_NW, stat = iderr)

        endif

      end if

      call mpl%f_comm%barrier()
!      call barrierMPI

! --- Send data from NE portion of processors at the higher generation
!     to corresponding  PE's at lower generation

  if(grid%my_hgen==mygen_up .and. grid%itargdn_ne >= 0 ) then
        nebpe = grid%itargdn_ne
        if(nebpe == grid%mype) then

             do L=1,Lm
             do j=0,grid%jmL
             do i=0,grid%imL
                dBuf_NE(i,j,L,:) = Warray(grid%imL+i,grid%jmL+j,L,:)
             enddo
             enddo
             enddo

        else

        allocate( sBuf_NE(0:grid%imL,0:grid%jmL,1:Lm,1:km), stat = iaerr )

             do L=1,Lm
             do j=0,grid%jmL
             do i=0,grid%imL
                sBuf_NE(i,j,L,:) = Warray(grid%imL+i,grid%jmL+j,L,:)
             enddo
             enddo
             enddo

        sHandle(4) = mpl%f_comm%isend(sBuf_NE,nebpe, grid%mype)
        call mpl%f_comm%wait(sHandle(4), istat)
!        call MPI_ISEND( sBuf_NE, ndata, dtype, nebpe, grid%mype,  &
!                        mpi_comm_comp, sHandle(4), isend)
!        call MPI_WAIT( sHandle(4), istat, ierr )
        deallocate( sBuf_NE, stat = ierr )

        endif

  endif
!
! --- Receive NE portion of data at lower generation
!

      if( lsendup_ne ) then
        nebpe = itarg_up

        if(nebpe == grid%mype) then

             do L=1,Lm
             do j=0,grid%jmL
             do i=0,grid%imL
               Harray(i,j,L,:)=dBuf_NE(i,j,L,:)
             enddo
             enddo
             enddo

        else

        allocate( rBuf_NE(0:grid%imL,0:grid%jmL,1:Lm,1:km), stat = iaerr )

        rHandle(4) = mpl%f_comm%ireceive(rBuf_NE,nebpe, nebpe)
        call mpl%f_comm%wait(rHandle(4), istat)
!        call MPI_IRECV( rBuf_NE, ndata, dtype, nebpe, nebpe, &
!                        mpi_comm_comp, rHandle(4), irecv)
!        call MPI_WAIT( rHandle(4), istat, ierr )

             do L=1,Lm
             do j=0,grid%jmL
             do i=0,grid%imL
               Harray(i,j,L,:)=rBuf_NE(i,j,L,:)
             enddo
             enddo
             enddo

        deallocate( rBuf_NE, stat = iderr)

        endif

      end if

      call mpl%f_comm%barrier()
!      call barrierMPI


!-----------------------------------------------------------------------
                        endsubroutine downsend

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine boco05                               &
!***********************************************************************
!                                                                      *
! Half values at the edges to accomodate adjoint of beta filter        *
!                                                                      *
!***********************************************************************
(grid,V,H,km,im,jm,Lm,nbx,nby)
!-----------------------------------------------------------------------
!use mpi

implicit none

!-----------------------------------------------------------------------
type(mgbf_grid_type),intent(inout) :: grid
integer(kind_int), intent(in):: km,im,jm,Lm,nbx,nby
real(kind_real), dimension(-nbx:im+nbx,-nby:jm+nby,Lm,1:km),intent(inout):: V,H
!-----------------------------------------------------------------------
logical lgen
integer(kind_int) L
!-----------------------------------------------------------------------
!
! Limit comminications to selected number of generations
!
!
! Define new boundaries
!

           do L=1,LM
              V(0 ,0:jm,L,:)=V(0 ,0:jm,L,:)*0.5
              V(im,0:jm,L,:)=V(im,0:jm,L,:)*0.5
              V(0:im,0 ,L,:)=V(0:im,0 ,L,:)*0.5
              V(0:im,jm,L,:)=V(0:im,jm,L,:)*0.5
            end do

        if(grid%l_hgen) then

           do L=1,LM
              H(0 ,0:jm,L,:)=H(0 ,0:jm,L,:)*0.5
              H(im,0:jm,L,:)=H(im,0:jm,L,:)*0.5
              H(0:im,0 ,L,:)=H(0:im,0 ,L,:)*0.5
              H(0:im,jm,L,:)=H(0:im,jm,L,:)*0.5
            end do

        endif
!-----------------------------------------------------------------------
                        endsubroutine boco05

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine v02v                                 &
!**********************************************************************!
!                                                                      !
!  Provide common edges for analysis subdomains                        !
!                                                                      !
!  Conversion from 'original' analysis array (which does not share     !
!  neighbors) to 'work' analysis array (which share edges with         !
!  neighbors                                                           !
!                                                                      !
!          o o o o o o          x x x x x x                            !
!          o x x x x o          x x x x x x                            !
!          o x x x x o    ->    x x x x x x                            !
!          o x x x x o          x x x x x x                            !
!          o x x x x o          x x x x x x                            !
!          o o o o o o          x x x x x x                            !
!                                                                      !
!**********************************************************************!
(mpl,grid,WA,nmax,mmax,lmax)
!-----------------------------------------------------------------------
!use mg_domain1, only: itarg_wA,itarg_eA,itarg_sA,itarg_nA               &
!                     ,lwestA,leastA,lsouthA,lnorthA
!use mpi

implicit none

!-----------------------------------------------------------------------
type(mpl_type),intent(inout) :: mpl
type(mgbf_grid_type),intent(inout) :: grid
integer(kind_int), intent(in):: nmax,mmax,lmax
real(kind_real), dimension(0:nmax,0:mmax,lmax),intent(inout):: WA
!-----------------------------------------------------------------------

real(kind_real), allocatable, dimension(:,:):: sBuf_N,rBuf_S
real(kind_real), allocatable, dimension(:,:):: sBuf_E,rBuf_W
integer(kind_int) sHandle(2),rHandle(2)!,ISTAT(MPI_STATUS_SIZE)
integer(kind_int) iaerr,ierr,iderr,l,n,m
integer(kind_int) isend,irecv,nebpe
integer(kind_int) ndatax,ndatay
logical lgen
type(fckit_mpi_status) :: istat
!-----------------------------------------------------------------------

!
! Define boundary conditions
!


!
!                           SEND boundaries toward North
!

      ndatay = nmax*lmax



      if( grid%itarg_nA >= 0 ) then
        nebpe = grid%itarg_nA

           allocate( sBuf_N(1:nmax,1:lmax), stat = iaerr )

              do n=1,nmax
                sBuf_N(n,:) = WA(n,mmax,:)
              enddo 
              sHandle(1) = mpl%f_comm%isend(sBuf_N, nebpe, grid%mype)
!              call MPI_ISEND( sBuf_N, ndatay, dtype, nebpe, mype,  &
!                              mpi_comm_comp, sHandle(1), isend)
      end if

!
!                           RECEIVE boundaries from South
!

      if( grid%itarg_sA >= 0 ) then
        nebpe = grid%itarg_sA
          allocate( rBuf_S(1:nmax,1:lmax), stat = iaerr )
          rHandle(1) = mpl%f_comm%ireceive(rBuf_S, nebpe, nebpe)
          call mpl%f_comm%wait(rHandle(1), istat)
!          call MPI_IRECV( rBuf_S, ndatay, dtype, nebpe, nebpe, &
!                      mpi_comm_comp, rHandle(1), irecv)
!          call MPI_WAIT( rHandle(1), istat, ierr )

      end if

!
!                           ASSIGN boundaries received from South
!

      if(  grid%itarg_sA >= 0 ) then

          do n=1,nmax
            WA(n,0,:)=rBuf_S(n,:)
          enddo

      endif

                                  call mpl%f_comm%barrier()
!                                  call barrierMPI

!----------------------------------------------------------------------

!
!                           SEND boundaries toward Eeast
!

      ndatax = (mmax+1)*lmax

      if( grid%itarg_eA >= 0) then
        nebpe = grid%itarg_eA

              allocate( sBuf_E(0:mmax,1:lmax), stat = iaerr )

              do m=0,mmax
                sBuf_E(m,1:lmax) = WA(nmax,m,1:lmax)
              enddo
              sHandle(2) = mpl%f_comm%isend(sBuf_E, nebpe, grid%mype)
!              call MPI_ISEND( sBuf_E, ndatax, dtype, nebpe, mype, &
!                             mpi_comm_comp, sHandle(2), isend)

      end if

!
!                           RECEIVE boundaries from Weast
!

      if( grid%itarg_wA >= 0 ) then
        nebpe = grid%itarg_wA
          allocate( rBuf_W(0:mmax,1:lmax), stat = iaerr )
          rHandle(2) = mpl%f_comm%ireceive(rBuf_W, nebpe, nebpe)
          call mpl%f_comm%wait(rHandle(2), istat)
!          call MPI_IRECV( rBuf_W, ndatax, dtype, nebpe, nebpe,  &
!                       mpi_comm_comp, rHandle(2), irecv)
!          call MPI_WAIT( rHandle(2), istat, ierr )

      end if
!
!                           ASSIGN boundaries from east

      if( grid%itarg_wA >= 0 ) then

          do m=0,mmax
            WA(0,m,1:lmax)=rBuf_W(m,1:lmax)
          enddo
      else
          do m=0,mmax                            !  Do not need that
            WA(0,m,1:lmax)=0.                    !
          enddo                                  !

      end if


                                  call mpl%f_comm%barrier()
!                                  call barrierMPI

!-----------------------------------------------------------------------
!
!                           DEALLOCATE Bufferes
!
      if(grid%itarg_eA >= 0) deallocate( sBuf_E)
      if(grid%itarg_wA >= 0) deallocate( rBuf_W)

      if(grid%itarg_nA >= 0) deallocate( sBuf_N)
      if(grid%itarg_sA >= 0) deallocate( rBuf_S)


!-----------------------------------------------------------------------
                        endsubroutine v02v



!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        endmodule tools_bocos
