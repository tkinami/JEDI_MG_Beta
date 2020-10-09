!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        module tools_generations
!***********************************************************************
!                                                                      !
!  Contains procedures that include differrent generations             !
!                                                                      !
!                                                     M. Rancic (2020) !
!***********************************************************************
use tools_kinds, only: kind_real,kind_int
use tools_bocos, only: bocosHn,bocosHnT,upsend,downsend
use type_mgbf_grid
use type_mgbf_state
use type_mpl

interface weighting1
module procedure weighting1_2d,weighting1_3d
endinterface weighting1

interface forward_k
module procedure forward_k2d,forward_k3d
endinterface forward_k

interface adjoint_k
module procedure adjoint_k2d,adjoint_k3d
endinterface adjoint_k

interface upsending1
module procedure upsending1_2d,upsending1_3d
endinterface upsending1

interface downsending1
module procedure downsending1_2d,downsending1_3d
endinterface downsending1

interface differencing1
module procedure differencing1_2d,differencing1_3d
endinterface differencing1

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        contains

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine weighting1_2d                        &
!***********************************************************************
!                                                                      !
!  Apply weight coefficients                                           !
!                                                                      !
!***********************************************************************
(grid,state,V,H)
!-----------------------------------------------------------------------
IMPLICIT NONE

type(mgbf_grid_type),intent(in) :: grid
type(mgbf_state_type),intent(inout) :: state
real(kind_real),dimension(-grid%hx:grid%im+grid%hx,-grid%hy:grid%jm+grid%hy,1:grid%km2),intent(inout):: V,H
!-----------------------------------------------------------------------

          V(:,:,:)=V(:,:,:)*state%weig_2d(:,:,:,1)

       if(grid%l_hgen) then
          H(:,:,:)=H(:,:,:)*state%weig_2d(:,:,:,2)
       endif

!-----------------------------------------------------------------------
                        endsubroutine weighting1_2d

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine weighting1_3d                        &
!***********************************************************************
!                                                                      !
!  Apply weight coefficients                                           !                                  
!                                                                      !
!***********************************************************************
(grid,state,V,H)
!-----------------------------------------------------------------------
IMPLICIT NONE

type(mgbf_grid_type),intent(in) :: grid
type(mgbf_state_type),intent(inout) :: state
real(kind_real),dimension(-grid%hx:grid%im+grid%hx,-grid%hy:grid%jm+grid%hy,grid%lm,grid%km3),intent(inout):: V,H
!-----------------------------------------------------------------------

          V(:,:,:,:)=V(:,:,:,:)*state%weig_3d(:,:,:,:,1)

       if(grid%l_hgen) then
          H(:,:,:,:)=H(:,:,:,:)*state%weig_3d(:,:,:,:,2)
       endif

!-----------------------------------------------------------------------
                        endsubroutine weighting1_3d

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine forward_k2d                          &
!***********************************************************************
!                                                                      !
!   Mapping from the low to high resolution grid                       !
!   using linearly squared interpolations   (2d case)                  !
!                                                                      !
!***********************************************************************
(grid,L02D,H2D,g)
!-----------------------------------------------------------------------
IMPLICIT NONE

type(mgbf_grid_type),intent(in) :: grid
real(kind_real), dimension(-1:grid%imL+1,-1:grid%jmL+1,grid%km2), intent(in):: L02D
real(kind_real), dimension(0:grid%im,0:grid%jm,grid%km2), intent(out):: H2D
integer(kind_int):: g
real(kind_real), dimension(0:grid%im,-1:grid%jmL+1,grid%km2):: HXLY
integer(kind_int):: i,iL,j,jL
!-----------------------------------------------------------------------
 
 if(g==grid%my_hgen) then

    do jL=-1,grid%jmL+1
    do i=0,grid%im,2
      iL=i/2
      HXLY(i,jL,:)=L02D(iL,jL,:)
    enddo
    enddo

    do jL=-1,grid%jmL+1
    do i=1,grid%im-1,2
      iL=(i+1)/2
      HXLY(i,jL,:)=-0.0625*L02D(iL-2,jL,:) &
                   +0.5625*L02D(iL-1,jL,:) &
                   +0.5625*L02D(iL  ,jL,:) &
                   -0.0625*L02D(iL+1,jL,:)
     enddo
     enddo

     do j=0,grid%jm,2
       jL = j/2
     do i=0,grid%im
       H2D(i,j,:)=HXLY(i,jL,:)
     enddo
     enddo 

     do j=1,grid%jm-1,2
       jL=(j+1)/2
     do i=0,grid%im
      H2D(i,j,:)=-0.0625*HXLY(i,jL-2,:) &
                 +0.5625*HXLY(i,jL-1,:) &
                 +0.5625*HXLY(i,jL  ,:) &
                 -0.0625*HXLY(i,jL+1,:)
     enddo
     enddo

 endif
!-----------------------------------------------------------------------
                        endsubroutine forward_k2d    

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine forward_k3d                          &
!***********************************************************************
!                                                                      !
!   Mapping from the low to high resolution grid                       !
!   using linearly squared interpolations   (3d case)                  !
!                                                                      !
!***********************************************************************
(grid,L02D,H2D,g)
!-----------------------------------------------------------------------
IMPLICIT NONE
type(mgbf_grid_type),intent(in) :: grid
real(kind_real), dimension(-1:grid%imL+1,-1:grid%jmL+1,grid%lm,grid%km3), intent(in):: L02D
real(kind_real), dimension(0:grid%im,0:grid%jm,grid%lm,grid%km3), intent(out):: H2D
integer(kind_int):: g
real(kind_real), dimension(0:grid%im,-1:grid%jmL+1,grid%km3):: HXLY
integer(kind_int):: i,iL,j,jL,L
!-----------------------------------------------------------------------

 if(g==grid%my_hgen) then

  do L=1,grid%lm

    do jL=-1,grid%jmL+1
    do i=0,grid%im,2
      iL=i/2
      HXLY(i,jL,:)=L02D(iL,jL,L,:)
    enddo
    enddo

    do jL=-1,grid%jmL+1
    do i=1,grid%im-1,2
      iL=(i+1)/2
      HXLY(i,jL,:)=-0.0625*L02D(iL-2,jL,L,:) &
                   +0.5625*L02D(iL-1,jL,L,:) &
                   +0.5625*L02D(iL  ,jL,L,:) &
                   -0.0625*L02D(iL+1,jL,L,:)
     enddo
     enddo

     do j=0,grid%jm,2
       jL = j/2
     do i=0,grid%im
       H2D(i,j,L,:)=HXLY(i,jL,:)
     enddo
     enddo 

     do j=1,grid%jm-1,2
       jL=(j+1)/2
     do i=0,grid%im
      H2D(i,j,L,:)=-0.0625*HXLY(i,jL-2,:) &
                   +0.5625*HXLY(i,jL-1,:) &
                   +0.5625*HXLY(i,jL  ,:) &
                   -0.0625*HXLY(i,jL+1,:)
     enddo
     enddo

  enddo

 endif
!-----------------------------------------------------------------------
                        endsubroutine forward_k3d    

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine adjoint_k2d                          &
!***********************************************************************
!                                                                      !
!   Mapping from the high to low resolution grid                       !
!   using linearly squared interpolations   (2d case)                  !
!                                                                      !
!***********************************************************************
(grid,H02D,L2D,g)
!-----------------------------------------------------------------------
IMPLICIT NONE
type(mgbf_grid_type),intent(in) :: grid
real(kind_real), dimension(0:grid%im,0:grid%jm,grid%km2), intent(in):: H02D
real(kind_real), dimension(-1:grid%imL+1,-1:grid%jmL+1,grid%km2), intent(out):: L2D
integer(kind_int), intent(in):: g
real(kind_real), dimension(0:grid%im,-1:grid%jmL+1,grid%km2):: HXLY
integer(kind_int):: i,j,iL,jL
!-----------------------------------------------------------------------

  if(g==1 .or. g==grid%my_hgen) then

      L2D(:,:,:) = 0.
      HXLY(:,:,:)= 0.

    do j=grid%jm-1,1,-2
      jL = (j+1)/2
    do i=grid%im,0,-1
       iL=i/2
      HXLY(i,jL-2,:)=HXLY(i,jL-2,:)-0.0625*H02D(i,j,:)
      HXLY(i,jL-1,:)=HXLY(i,jL-1,:)+0.5625*H02D(i,j,:)
      HXLY(i,jL  ,:)=HXLY(i,jL  ,:)+0.5625*H02D(i,j,:)
      HXLY(i,jL+1,:)=HXLY(i,jL+1,:)-0.0625*H02D(i,j,:)
    enddo
    enddo
  
    do j=grid%jm,0,-2
      jL = j/2
    do i=grid%im,0,-1
      HXLY(i,jL,:)=HXLY(i,jL,:)+H02D(i,j,:)
    enddo
    enddo

    do jL=grid%jmL+1,-1,-1
    do i=grid%im-1,1,-2
      iL = (i+1)/2
      L2D(iL-2,jL,:)=L2D(iL-2,jL,:)-0.0625*HXLY(i,jL,:)
      L2D(iL-1,jL,:)=L2D(iL-1,jL,:)+0.5625*HXLY(i,jL,:)
      L2D(iL  ,jL,:)=L2D(iL  ,jL,:)+0.5625*HXLY(i,jL,:)
      L2D(iL+1,jL,:)=L2D(iL+1,jL,:)-0.0625*HXLY(i,jL,:)
    enddo
    enddo

    do jL=grid%jmL+1,-1,-1
    do i =grid%im,0,-2
      iL = i/2
      L2D(iL,jL,:)=L2D(iL,jL,:)+HXLY(i,jL,:)
    enddo  
    enddo  

  endif
!-----------------------------------------------------------------------
                        endsubroutine adjoint_k2d    

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine adjoint_k3d                          &
!***********************************************************************
!                                                                      !
!   Mapping from the high to low resolution grid                       !
!   using linearly squared interpolations   (3d case)                  !
!                                                                      !
!***********************************************************************
(grid,H02D,L2D,g)
!-----------------------------------------------------------------------
IMPLICIT NONE
type(mgbf_grid_type),intent(in) :: grid
real(kind_real), dimension(0:grid%im,0:grid%jm,grid%lm,grid%km3), intent(in):: H02D
real(kind_real), dimension(-1:grid%imL+1,-1:grid%jmL+1,grid%lm,grid%km3), intent(out):: L2D
integer(kind_int), intent(in):: g
real(kind_real), dimension(0:grid%im,-1:grid%jmL+1,grid%km3):: HXLY
integer(kind_int):: i,j,iL,jL,l
!-----------------------------------------------------------------------

  if(g==1 .or. g==grid%my_hgen) then

      L2D(:,:,:,:) = 0.

    Do l=1,grid%lm

    do j=grid%jm-1,1,-2
      jL = (j+1)/2
    do i=grid%im,0,-1
       iL=i/2
      HXLY(i,jL-2,:)=HXLY(i,jL-2,:)-0.0625*H02D(i,j,L,:)
      HXLY(i,jL-1,:)=HXLY(i,jL-1,:)+0.5625*H02D(i,j,L,:)
      HXLY(i,jL  ,:)=HXLY(i,jL  ,:)+0.5625*H02D(i,j,L,:)
      HXLY(i,jL+1,:)=HXLY(i,jL+1,:)-0.0625*H02D(i,j,L,:)
    enddo
    enddo
  
    do j=grid%jm,0,-2
      jL = j/2
    do i=grid%im,0,-1
      HXLY(i,jL,:)=HXLY(i,jL,:)+H02D(i,j,L,:)
    enddo
    enddo

    do jL=grid%jmL+1,-1,-1
    do i=grid%im-1,1,-2
      iL = (i+1)/2
      L2D(iL-2,jL,L,:)=L2D(iL-2,jL,L,:)-0.0625*HXLY(i,jL,:)
      L2D(iL-1,jL,L,:)=L2D(iL-1,jL,L,:)+0.5625*HXLY(i,jL,:)
      L2D(iL  ,jL,L,:)=L2D(iL  ,jL,L,:)+0.5625*HXLY(i,jL,:)
      L2D(iL+1,jL,L,:)=L2D(iL+1,jL,L,:)-0.0625*HXLY(i,jL,:)
    enddo
    enddo

    do jL=grid%jmL+1,-1,-1
    do i =grid%im,0,-2
      iL = i/2
      L2D(iL,jL,L,:)=L2D(iL,jL,L,:)+HXLY(i,jL,:)
    enddo  
    enddo  

    enddo

  endif
!-----------------------------------------------------------------------
                        endsubroutine adjoint_k3d    

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine upsending1_2d                        &
!***********************************************************************
!                                                                      !
!  Adjoint interpolate and upsend:                                     !
!       First from g1->g2 (V -> H)                                     !
!       Then  from g2->...->gn  (H -> H)                               !
!                                                                      !
!***********************************************************************
(mpl,grid,V,H)
!-----------------------------------------------------------------------
IMPLICIT NONE

type(mpl_type),intent(inout) :: mpl
type(mgbf_grid_type),intent(inout) :: grid
real(kind_real),dimension(-grid%hx:grid%im+grid%hx,-grid%hy:grid%jm+grid%hy,1:grid%km2),intent(in):: V
real(kind_real),dimension(-grid%hx:grid%im+grid%hx,-grid%hy:grid%jm+grid%hy,1:grid%km2),intent(out):: H

real(kind_real),dimension(0:grid%im,0:grid%jm,1:grid%km2):: V_PROX,H_PROX
real(kind_real),dimension(-1:grid%imL+1,-1:grid%jmL+1,1:grid%km2):: V_INT,H_INT
integer(kind_int):: g,L
!-----------------------------------------------------------------------
!
! From generation 1 to generation 2
!
        call adjoint_k(grid,V(0:grid%im,0:grid%jm,1:grid%km2),V_INT,1) 
        call bocosHnT(mpl,grid,V_INT,grid%km2,grid%imL,grid%jmL,1,1,1,grid%FimaxL,grid%FjmaxL,1,1)
        call upsend(mpl,grid,V_INT(0:grid%imL,0:grid%jmL,1:grid%km2),H,grid%km2,1,1,2)
!
! From generation 2 sequentially to higher generations
!
      do g=2,grid%gm-1

        call adjoint_k(grid,H(0:grid%im,0:grid%jm,1:grid%km2),H_INT,g) 
        call bocosHnT(mpl,grid,H_INT,grid%km2,grid%imL,grid%jmL,1,1,1,grid%FimaxL,grid%FjmaxL,g,g)
        call upsend(mpl,grid,H_INT(0:grid%imL,0:grid%jmL,1:grid%km2),H,grid%km2,1,g,g+1)

      end do 

!-----------------------------------------------------------------------
                        endsubroutine upsending1_2d

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine upsending1_3d                        &
!***********************************************************************
!                                                                      !
!  Adjoint interpolate and upsend:                                     !
!       First from g1->g2 (V -> H)                                     !
!       Then  from g2->...->gn  (H -> H)                               !
!                                                                      !
!***********************************************************************
(mpl,grid,V,H)
!-----------------------------------------------------------------------
IMPLICIT NONE
type(mpl_type),intent(inout):: mpl
type(mgbf_grid_type),intent(inout) :: grid
real(kind_real),dimension(-grid%hx:grid%im+grid%hx,-grid%hy:grid%jm+grid%hy,grid%lm,grid%km3),intent(in):: V
real(kind_real),dimension(-grid%hx:grid%im+grid%hx,-grid%hy:grid%jm+grid%hy,grid%lm,grid%km3),intent(out):: H

real(kind_real),dimension(-1:grid%imL+1,-1:grid%jmL+1,grid%lm,grid%km3):: V_INT,H_INT
integer(kind_int):: g,L
!-----------------------------------------------------------------------
!
! From generation 1 to generation 2
!

        call adjoint_k(grid,V(0:grid%im,0:grid%jm,1:grid%lm,1:grid%km3),V_INT,1) 
        call bocosHnT(mpl,grid,V_INT,grid%km3,grid%imL,grid%jmL,grid%lm,1,1,grid%FimaxL,grid%FjmaxL,1,1)
        call upsend(mpl,grid,V_INT(0:grid%imL,0:grid%jmL,1:grid%lm,1:grid%km3),H,grid%km3,grid%lm,1,2)

!
! From generation 2 sequentially to higher generations
!
      do g=2,grid%gm-1

        call adjoint_k(grid,H(0:grid%im,0:grid%jm,1:grid%lm,1:grid%km3),H_INT,g) 
        call bocosHnT(mpl,grid,H_INT,grid%km3,grid%imL,grid%jmL,grid%lm,1,1,grid%FimaxL,grid%FjmaxL,g,g)
        call upsend(mpl,grid,H_INT(0:grid%imL,0:grid%jmL,1:grid%lm,1:grid%km3),H,grid%km3,grid%lm,g,g+1)

      end do 


!-----------------------------------------------------------------------
                        endsubroutine upsending1_3d

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine downsending1_2d                      &
!***********************************************************************
!                                                                      !
!  Downsend, interpolate and add:                                      !
!      First from grid%gm->g3...->g2                                        !
!      Then  from g2->g1                                               !
!                                                                      !
!***********************************************************************
(mpl,grid,H,V)
!-----------------------------------------------------------------------
IMPLICIT NONE
type(mpl_type),intent(inout):: mpl
type(mgbf_grid_type),intent(inout) :: grid
real(kind_real),dimension(-grid%hx:grid%im+grid%hx,-grid%hy:grid%jm+grid%hy,1:grid%km2),intent(inout):: H
real(kind_real),dimension(-grid%hx:grid%im+grid%hx,-grid%hy:grid%jm+grid%hy,1:grid%km2),intent(out):: V
real(kind_real),dimension(-1:grid%imL+1,-1:grid%jmL+1,1:grid%km2):: H_INT,V_INT
integer(kind_int):: g,l,k
!-----------------------------------------------------------------------

      do g=grid%gm,3,-1
        call downsend(mpl,grid,H(0:grid%im,0:grid%jm,1:grid%km2),H_INT,grid%km2,1,g,g-1)
          if(g==grid%my_hgen) H(:,:,:)=0.
        call bocosHn(mpl,grid,H_INT,grid%km2,grid%imL,grid%jmL,1,1,1,grid%FimaxL,grid%FjmaxL,g-1,g-1)
        call forward_k(grid,H_INT,H(0:grid%im,0:grid%jm,1:grid%km2),g-1)
      enddo


        call downsend(mpl,grid,H(0:grid%im,0:grid%jm,1:grid%km2),V_INT,grid%km2,1,2,1)
          H(:,:,:)=0.
        call bocosHn(mpl,grid,V_INT,grid%km2,grid%imL,grid%jmL,1,1,1,grid%FimaxL,grid%FjmaxL,1,1)
        call forward_k(grid,V_INT,V(0:grid%im,0:grid%jm,1:grid%km2),1)

!-----------------------------------------------------------------------
                        endsubroutine downsending1_2d

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine downsending1_3d                      &
!***********************************************************************
!                                                                      !
!  Downsend, interpolate and add:                                      !
!      First from grid%gm->g3...->g2                                        !
!      Then  from g2->g1                                               !
!                                                                      !
!***********************************************************************
(mpl,grid,H,V)
!-----------------------------------------------------------------------
IMPLICIT NONE

type(mpl_type),intent(inout) :: mpl
type(mgbf_grid_type),intent(inout) :: grid
real(kind_real),dimension(-grid%hx:grid%im+grid%hx,-grid%hy:grid%jm+grid%hy,grid%lm,grid%km3),intent(inout):: H
real(kind_real),dimension(-grid%hx:grid%im+grid%hx,-grid%hy:grid%jm+grid%hy,grid%lm,grid%km3),intent(out):: V
real(kind_real),dimension(-1:grid%imL+1,-1:grid%jmL+1,grid%lm,grid%km3):: H_INT,V_INT
integer(kind_int):: g,l,k
!-----------------------------------------------------------------------

      do g=grid%gm,3,-1
        call downsend(mpl,grid,H(0:grid%im,0:grid%jm,1:grid%lm,1:grid%km3),H_INT,grid%km3,grid%lm,g,g-1)
          if(g==grid%my_hgen) H(:,:,:,:)=0.
        call bocosHn(mpl,grid,H_INT,grid%km3,grid%imL,grid%jmL,grid%lm,1,1,grid%FimaxL,grid%FjmaxL,g-1,g-1)
        call forward_k(grid,H_INT,H(0:grid%im,0:grid%jm,1:grid%lm,1:grid%km3),g-1)
      enddo

        call downsend(mpl,grid,H(0:grid%im,0:grid%jm,1:grid%lm,1:grid%km3),V_INT,grid%km3,1,2,1)
          H(:,:,:,:)=0.
        call bocosHn(mpl,grid,V_INT,grid%km3,grid%imL,grid%jmL,grid%lm,1,1,grid%FimaxL,grid%FjmaxL,1,1)
        call forward_k(grid,V_INT,V(0:grid%im,0:grid%jm,1:grid%lm,1:grid%km3),1)

!-----------------------------------------------------------------------
                        endsubroutine downsending1_3d

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine differencing1_2d                     &
!***********************************************************************
!                                                                      !
!  Apply 2D differential operator                                      !
!                                                                      !
!***********************************************************************
(grid,state,V,H)
!-----------------------------------------------------------------------
IMPLICIT NONE

type(mgbf_grid_type),intent(in) :: grid
type(mgbf_state_type),intent(in) :: state
real(kind_real),dimension(-grid%hx:grid%im+grid%hx,-grid%hy:grid%jm+grid%hy,grid%km2),intent(inout):: V
real(kind_real),dimension(-grid%hx:grid%im+grid%hx,-grid%hy:grid%jm+grid%hy,grid%km2),intent(inout):: H
real(kind_real),dimension(-1:grid%im,0:grid%jm, grid%km2):: DIFX
real(kind_real),dimension( 0:grid%im,-1:grid%jm,grid%km2):: DIFY
integer(kind_int):: i,j,k
!-----------------------------------------------------------------------


     do j=0,grid%jm
     do i=-1,grid%im
       DIFX(i,j,:)=V(i+1,j,:)-V(i,j,:)
     enddo
     enddo
     do j=-1,grid%jm
     do i=0,grid%im
       DIFY(i,j,:)=V(i,j+1,:)-V(i,j,:)
     enddo
     enddo
     do j=0,grid%jm
     do i=0,grid%im
       V(i,j,:)=state%a2_diff(i,j,:,1)*V(i,j,:)                               &
               -state%b2_diff(i,j,:,1)*(DIFX(i,j,:)-DIFX(i-1,j,:)             &
                                 +DIFY(i,j,:)-DIFY(i,j-1,:))
     enddo
     enddo

 if(grid%l_hgen) then

     do j=0,grid%jm
     do i=-1,grid%im
       DIFX(i,j,:)=H(i+1,j,:)-H(i,j,:)
     enddo
     enddo
     do j=-1,grid%jm
     do i=0,grid%im
       DIFY(i,j,:)=H(i,j+1,:)-H(i,j,:)
     enddo
     enddo
     do j=0,grid%jm
     do i=0,grid%im
       H(i,j,:)=state%a2_diff(i,j,:,2)*H(i,j,:)                               &
               -state%b2_diff(i,j,:,2)*(DIFX(i,j,:)-DIFX(i-1,j,:)             &
                                 +DIFY(i,j,:)-DIFY(i,j-1,:))
     enddo
     enddo

 endif


!-----------------------------------------------------------------------
                        endsubroutine differencing1_2d

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine differencing1_3d                     &
!***********************************************************************
!                                                                      !
!  Apply 3D differential operator                                      !
!                                                                      !
!***********************************************************************
(grid,state,V,H)
!-----------------------------------------------------------------------
IMPLICIT NONE

type(mgbf_grid_type),intent(in) :: grid
type(mgbf_state_type),intent(in) :: state
real(kind_real),dimension(-grid%hx:grid%im+grid%hx,-grid%hy:grid%jm+grid%hy,grid%lm,grid%km3),intent(inout):: V
real(kind_real),dimension(-grid%hx:grid%im+grid%hx,-grid%hy:grid%jm+grid%hy,grid%lm,grid%km3),intent(inout):: H
real(kind_real),dimension(-1:grid%im,0:grid%jm ,grid%lm  ,grid%km3):: DIFX
real(kind_real),dimension(0:grid%im ,-1:grid%jm,grid%lm  ,grid%km3):: DIFY
real(kind_real),dimension(0:grid%im ,0:grid%jm ,0:grid%lm,grid%km3):: DIFZ
integer(kind_int):: i,j,l,k
!-----------------------------------------------------------------------

   do l=1,grid%lm
     do j=0,grid%jm
     do i=-1,grid%im
       DIFX(i,j,l,:)=V(i+1,j,l,:)-V(i,j,l,:)
     enddo
     enddo
     do j=-1,grid%jm
     do i=0,grid%im
       DIFY(i,j,l,:)=V(i,j+1,l,:)-V(i,j,l,:)
     enddo
     enddo
   enddo

   do l=1,grid%lm-1
     do j=0,grid%jm
     do i=0,grid%im
       DIFZ(i,j,l,:)=V(i,j,l+1,:)-V(i,j,l,:)
     enddo
     enddo
   enddo
     do j=0,grid%jm
     do i=0,grid%im
       DIFZ(i,j,0 ,:)=-DIFZ(i,j,1   ,:)
       DIFZ(i,j,grid%lm,:)=-DIFZ(i,j,grid%lm-1,:)
     enddo
     enddo

   do l=1,grid%lm
     do j=0,grid%jm
     do i=0,grid%im
        V(i,j,l,:)=state%a3_diff(i,j,l,:,1)*V(i,j,l,:)                        &
                  -state%b3_diff(i,j,l,:,1)*(DIFX(i,j,l,:)-DIFX(i-1,j,l,:)    &
                                      +DIFY(i,j,l,:)-DIFY(i,j-1,l,:)    &
                                      +DIFZ(i,j,l,:)-DIFZ(i,j,l-1,:))
     enddo
     enddo
   enddo

if(grid%l_hgen) then

   do l=1,grid%lm
     do j=0,grid%jm
     do i=-1,grid%im
       DIFX(i,j,l,:)=H(i+1,j,l,:)-H(i,j,l,:)
     enddo
     enddo
     do j=-1,grid%jm
     do i=0,grid%im
       DIFY(i,j,l,:)=H(i,j+1,l,:)-H(i,j,l,:)
     enddo
     enddo
   enddo

   do l=1,grid%lm-1
     do j=0,grid%jm
     do i=0,grid%im
       DIFZ(i,j,l,:)=H(i,j,l+1,:)-H(i,j,l,:)
     enddo
     enddo
   enddo
     do j=0,grid%jm
     do i=0,grid%im
       DIFZ(i,j,0 ,:)=-DIFZ(i,j,1   ,:)
       DIFZ(i,j,grid%lm,:)=-DIFZ(i,j,grid%lm-1,:)
     enddo
     enddo

   do l=1,grid%lm
     do j=0,grid%jm
     do i=0,grid%im
        H(i,j,l,:)=state%a3_diff(i,j,l,:,2)*H(i,j,l,:)                       &
                  -state%b3_diff(i,j,l,:,2)*(DIFX(i,j,l,:)-DIFX(i-1,j,l,:)   &
                                      +DIFY(i,j,l,:)-DIFY(i,j-1,l,:)   &
                                      +DIFZ(i,j,l,:)-DIFZ(i,j,l-1,:))
     enddo
     enddo
   enddo

endif

!-----------------------------------------------------------------------
                        endsubroutine differencing1_3d


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        endmodule tools_generations
