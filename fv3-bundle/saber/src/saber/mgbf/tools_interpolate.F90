!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        module tools_interpolate
!***********************************************************************
!                                                                      !
!    general mapping between 2d arrays using linerly squared           !
!    interpolations                                                    !
!                                                                      !
!                                                     M. Rancic (2020) !
!***********************************************************************
use tools_kinds
!use mg_parameter1, only: grid%x0,grid%y0,grid%dxa,grid%dxf,grid%dya,grid%dyf                          &
!                        ,grid%nm,grid%mm,grid%km,grid%km2,grid%km3,grid%lm                            &
!                        ,grid%im,grid%jm,grid%ib,grid%jb                                    
!use mg_intstate1, only: state%iref,state%jref                                       &
!                      ,state%cx0,state%cx1,state%cx2,state%cx3                                  &
!                      ,state%cy0,state%cy1,state%cy2,state%cy3                                  &
!                      ,state%cf00,state%cf01,state%cf02,state%cf03                              &
!                      ,state%cf10,state%cf11,state%cf12,state%cf13                              &
!                      ,state%cf20,state%cf21,state%cf22,state%cf23                              &
!                      ,state%cf30,state%cf31,state%cf32,state%cf33
!use mg_mppstuff1, only: finishMPI,mype
use type_mgbf_grid, only: mgbf_grid_type
use type_mgbf_state, only: mgbf_state_type
implicit none

public lsqr_mg_coef


interface lsqr_forward_xk
 module procedure lsqr_forward_x_k2d
 module procedure lsqr_forward_x_k3d
endinterface

interface lsqr_adjoint_xk
  module procedure lsqr_adjoint_x_k2d
  module procedure lsqr_adjoint_x_k3d
endinterface

interface lsqr_forward_xyk
 module procedure lsqr_forward_xy_k2d
 module procedure lsqr_forward_xy_k3d
endinterface

interface lsqr_adjoint_xyk
  module procedure lsqr_adjoint_xy_k2d
  module procedure lsqr_adjoint_xy_k3d
endinterface

interface lsqr_forward
  module procedure lsqr_forward_k2d
  module procedure lsqr_forward_k3d
endinterface

interface lsqr_adjoint
  module procedure lsqr_adjoint_k2d
  module procedure lsqr_adjoint_k3d
endinterface

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine lsqr_mg_coef(grid,state)                            
!***********************************************************************
!                                                                      !
!   Prepare coeficients for mapping between:                           !
!        filter grid on analysis decomposition: W(-grid%ib:grid%im+grid%ib,-grid%jb:grid%jm+grid%jb) !
!        and analysis grid:                     V(0:grid%nm,0:grid%mm)           !  
!                                                                      !
!              (  grid%im <= grid%nm  and  grid%jm < grid%mm   )                           !
!                                                                      !
!***********************************************************************
implicit none

type(mgbf_grid_type),intent(in) :: grid
type(mgbf_state_type),intent(inout) :: state

real(kind_real), dimension(0:grid%nm):: xa
real(kind_real), dimension(-grid%ib:grid%im+grid%ib):: xf
real(kind_real), dimension(0:grid%mm):: ya
real(kind_real), dimension(-grid%jb:grid%jm+grid%jb):: yf
integer(kind_int):: i,j,n,m
real(kind_real) x1,x2,x3,x4,x
real(kind_real) x1x,x2x,x3x,x4x
real(kind_real) rx2x1,rx3x1,rx4x1,rx3x2,rx4x2,rx4x3
real(kind_real) y1,y2,y3,y4,y
real(kind_real) y1y,y2y,y3y,y4y
real(kind_real) ry2y1,ry3y1,ry4y1,ry3y2,ry4y2,ry4y3
real(kind_real) cfl1,cfl2,cfl3,cll
real(kind_real) cfr1,cfr2,cfr3,crr
!-----------------------------------------------------------------------
!
! Initialize
!
 
   do n=0,grid%nm
     xa(n)=grid%x0+n*grid%dxa
   enddo

   do i=-grid%ib,grid%im+grid%ib
     xf(i)=grid%x0+i*grid%dxf
   enddo

   do m=0,grid%mm
     ya(m)=grid%y0+m*grid%dya
   enddo

   do j=-grid%jb,grid%jm+grid%jb
     yf(j)=grid%y0+j*grid%dyf
   enddo

!
! Find state%iref and state%jref
!
   do n=0,grid%nm
     do i=-grid%ib,grid%im+grid%ib-1
       if(xf(i)<=xa(n) .and. xa(n)<xf(i+1)) then
         state%iref(n)=i-1
         exit
       endif
     enddo
   enddo

   do m=0,grid%mm
     do j=-grid%jb,grid%jm+grid%jb-1
       if(yf(j)<=ya(m) .and. ya(m)<yf(j+1)) then
         state%jref(m)=j-1
         exit
       endif
     enddo
   enddo


   do n=0,grid%nm
     i=state%iref(n)
     x1=xf(i)
     x2=xf(i+1)
     x3=xf(i+2)
     x4=xf(i+3)
     x = xa(n)
       x1x = x1-x   
       x2x = x2-x   
       x3x = x3-x   
       x4x = x4-x   
       rx2x1 = 1./(x2-x1)
       rx3x1 = 1./(x3-x1)
       rx4x1 = 1./(x4-x1)
       rx3x2 = 1./(x3-x2)
       rx4x2 = 1./(x4-x2)
       rx4x3 = 1./(x4-x3)
     CFL1 = x2x*x3x*rx2x1*rx3x1
     CFL2 =-x1x*x3x*rx2x1*rx3x2
     CFL3 = x1x*x2x*rx3x1*rx3x2
     CLL = x3x*rx3x2
     CFR1 = x3x*x4x*rx3x2*rx4x2
     CFR2 =-x2x*x4x*rx3x2*rx4x3
     CFR3 = x2x*x3x*rx4x2*rx4x3
     CRR =-x2x*rx3x2
       state%cx0(n)=CFL1*CLL
       state%cx1(n)=CFL2*CLL+CFR1*CRR
       state%cx2(n)=CFL3*CLL+CFR2*CRR
       state%cx3(n)=CFR3*CRR
   enddo

   do m=0,grid%mm
     j=state%jref(m)
     y1=yf(j)
     y2=yf(j+1)
     y3=yf(j+2)
     y4=yf(j+3)
     y = ya(m)
       y1y = y1-y   
       y2y = y2-y   
       y3y = y3-y   
       y4y = y4-y   
       ry2y1 = 1./(y2-y1)
       ry3y1 = 1./(y3-y1)
       ry4y1 = 1./(y4-y1)
       ry3y2 = 1./(y3-y2)
       ry4y2 = 1./(y4-y2)
       ry4y3 = 1./(y4-y3)
     CFL1 = y2y*y3y*ry2y1*ry3y1
     CFL2 =-y1y*y3y*ry2y1*ry3y2
     CFL3 = y1y*y2y*ry3y1*ry3y2
     CLL = y3y*ry3y2
     CFR1 = y3y*y4y*ry3y2*ry4y2
     CFR2 =-y2y*y4y*ry3y2*ry4y3
     CFR3 = y2y*y3y*ry4y2*ry4y3
     CRR =-y2y*ry3y2
       state%cy0(m)=CFL1*CLL
       state%cy1(m)=CFL2*CLL+CFR1*CRR
       state%cy2(m)=CFL3*CLL+CFR2*CRR
       state%cy3(m)=CFR3*CRR
   enddo

   do m=0,grid%mm
   do n=0,grid%nm
     state%cf00(n,m)=state%cx0(n)*state%cy0(m)
     state%cf01(n,m)=state%cx0(n)*state%cy1(m)
     state%cf02(n,m)=state%cx0(n)*state%cy2(m)
     state%cf03(n,m)=state%cx0(n)*state%cy3(m)
     state%cf10(n,m)=state%cx1(n)*state%cy0(m)
     state%cf11(n,m)=state%cx1(n)*state%cy1(m)
     state%cf12(n,m)=state%cx1(n)*state%cy2(m)
     state%cf13(n,m)=state%cx1(n)*state%cy3(m)
     state%cf20(n,m)=state%cx2(n)*state%cy0(m)
     state%cf21(n,m)=state%cx2(n)*state%cy1(m)
     state%cf22(n,m)=state%cx2(n)*state%cy2(m)
     state%cf23(n,m)=state%cx2(n)*state%cy3(m)
     state%cf30(n,m)=state%cx3(n)*state%cy0(m)
     state%cf31(n,m)=state%cx3(n)*state%cy1(m)
     state%cf32(n,m)=state%cx3(n)*state%cy2(m)
     state%cf33(n,m)=state%cx3(n)*state%cy3(m)
   enddo
   enddo
 
!-----------------------------------------------------------------------
                        endsubroutine lsqr_mg_coef


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine lsqr_forward_k2d                     &
!***********************************************************************
!                                                                      !
! Given a source array  V(-grid%ib:grid%im+grid%ib,-grid%jb:grid%jm+grid%jb,grid%km2) perform             !
! forward  interpolations to get target array W(0:grid%nm,0,grid%mm,grid%km2)         ! 
!                                                                      !
!***********************************************************************
(grid,state,V,W)
!-----------------------------------------------------------------------
implicit none
type(mgbf_grid_type),intent(in) :: grid
type(mgbf_state_type),intent(in) :: state
real(kind_real), dimension(-grid%ib:grid%im+grid%ib,-grid%ib:grid%im+grid%ib,grid%km2), intent(in):: V
real(kind_real), dimension(0:grid%nm,0:grid%mm,grid%km2),intent(out):: W  
integer(kind_int):: i,j,n,m,l
real(kind_real),dimension(grid%km2):: v00,v01,v02,v03                            &
                             ,v10,v11,v12,v13                            &
                             ,v20,v21,v22,v23                            &
                             ,v30,v31,v32,v33                                            
!-----------------------------------------------------------------------

   do m=0,grid%mm
     j = state%jref(m)
   do n=0,grid%nm
     i = state%iref(n)
     v00(:)=V(i  ,j  ,:)
     v10(:)=V(i+1,j  ,:)
     v20(:)=V(i+2,j  ,:)
     v30(:)=V(i+3,j  ,:)
     v01(:)=V(i  ,j+1,:)
     v11(:)=V(i+1,j+1,:)
     v21(:)=V(i+2,j+1,:)
     v31(:)=V(i+3,j+1,:)
     v02(:)=V(i  ,j+2,:)
     v12(:)=V(i+1,j+2,:)
     v22(:)=V(i+2,j+2,:)
     v32(:)=V(i+3,j+2,:)
     v03(:)=V(i  ,j+3,:)
     v13(:)=V(i+1,j+3,:)
     v23(:)=V(i+2,j+3,:)
     v33(:)=V(i+3,j+3,:)
     W(n,m,:) = state%cf00(n,m)*v00(:)+state%cf10(n,m)*v10(:)+state%cf20(n,m)*v20(:)+state%cf30(n,m)*v30(:)   &
               +state%cf01(n,m)*v01(:)+state%cf11(n,m)*v11(:)+state%cf21(n,m)*v21(:)+state%cf31(n,m)*v31(:)   &
               +state%cf02(n,m)*v02(:)+state%cf12(n,m)*v12(:)+state%cf22(n,m)*v22(:)+state%cf32(n,m)*v32(:)   &
               +state%cf03(n,m)*v03(:)+state%cf13(n,m)*v13(:)+state%cf23(n,m)*v23(:)+state%cf33(n,m)*v33(:) 
   enddo
   enddo


!-----------------------------------------------------------------------
                        endsubroutine lsqr_forward_k2d   

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine lsqr_forward_k3d                     &
!***********************************************************************
!                                                                      !
! Given a source array  V(-grid%ib:grid%im+grid%ib,-grid%jb:grid%jm+grid%jb,grid%lm,grid%km3) perform forward  !
! interpolations to get target array W(0:grid%nm,0,grid%mm,grid%lm,grid%km3)               ! 
!                                                                      !
!***********************************************************************
(grid,state,V,W)
!-----------------------------------------------------------------------
implicit none
type(mgbf_grid_type),intent(in) :: grid
type(mgbf_state_type),intent(in) :: state
real(kind_real), dimension(-grid%ib:grid%im+grid%ib,-grid%jb:grid%jm+grid%jb,grid%lm,grid%km3), intent(in):: V
real(kind_real), dimension(0:grid%nm,0:grid%mm,grid%lm,grid%km3),intent(out):: W  
integer(kind_int):: i,j,n,m,l
real(kind_real),dimension(grid%km3):: v00,v01,v02,v03                             &
                             ,v10,v11,v12,v13                             &
                             ,v20,v21,v22,v23                             &
                             ,v30,v31,v32,v33                                            
!-----------------------------------------------------------------------

   do l=1,grid%lm
   do m=0,grid%mm
     j = state%jref(m)
   do n=0,grid%nm
     i = state%iref(n)
     v00(:)=V(i  ,j  ,l,:)
     v10(:)=V(i+1,j  ,l,:)
     v20(:)=V(i+2,j  ,l,:)
     v30(:)=V(i+3,j  ,l,:)
     v01(:)=V(i  ,j+1,l,:)
     v11(:)=V(i+1,j+1,l,:)
     v21(:)=V(i+2,j+1,l,:)
     v31(:)=V(i+3,j+1,l,:)
     v02(:)=V(i  ,j+2,l,:)
     v12(:)=V(i+1,j+2,l,:)
     v22(:)=V(i+2,j+2,l,:)
     v32(:)=V(i+3,j+2,l,:)
     v03(:)=V(i  ,j+3,l,:)
     v13(:)=V(i+1,j+3,l,:)
     v23(:)=V(i+2,j+3,l,:)
     v33(:)=V(i+3,j+3,l,:)
     W(n,m,l,:) = state%cf00(n,m)*v00(:)+state%cf10(n,m)*v10(:)+state%cf20(n,m)*v20(:)+state%cf30(n,m)*v30(:)   &
                 +state%cf01(n,m)*v01(:)+state%cf11(n,m)*v11(:)+state%cf21(n,m)*v21(:)+state%cf31(n,m)*v31(:)   &
                 +state%cf02(n,m)*v02(:)+state%cf12(n,m)*v12(:)+state%cf22(n,m)*v22(:)+state%cf32(n,m)*v32(:)   &
                 +state%cf03(n,m)*v03(:)+state%cf13(n,m)*v13(:)+state%cf23(n,m)*v23(:)+state%cf33(n,m)*v33(:) 
   enddo
   enddo
   enddo

!-----------------------------------------------------------------------
                        endsubroutine lsqr_forward_k3d   

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine lsqr_adjoint_k2d                      &
!***********************************************************************
!                                                                      !
! Given a target array W(0:grid%nm,0:grid%mm,grid%km2) perform adjoint interpolations !
! in order to get source array V(-grid%ib:grid%im+grid%ib,-grid%jb:grid%jm+grid%jb,grid%km2)              ! 
!                                                                      !
!***********************************************************************
(grid,state,W,V)
!-----------------------------------------------------------------------
implicit none
type(mgbf_grid_type),intent(in) :: grid
type(mgbf_state_type),intent(in) :: state
real(kind_real), dimension(0:grid%nm,0:grid%mm,grid%km2),intent(in):: W  
real(kind_real), dimension(-grid%ib:grid%im+grid%ib,-grid%jb:grid%jm+grid%jb,grid%km2), intent(out):: V
integer(kind_int):: i,j,n,m
integer(kind_int):: ip1,ip2,ip3,jp1,jp2,jp3
!-----------------------------------------------------------------------
     
    V(:,:,:) = 0.

!$OMP PARALLEL DO 
   do m=0,grid%mm
       j = state%jref(m)
       jp1=j+1
       jp2=j+2
       jp3=j+3
   do n=0,grid%nm
       i = state%iref(n) 
       ip1=i+1
       ip2=i+2
       ip3=i+3
     V(i  ,j  ,:) = V(i  ,j  ,:)+W(n,m,:)*state%cf00(n,m)
     V(ip1,j  ,:) = V(ip1,j  ,:)+W(n,m,:)*state%cf10(n,m)
     V(ip2,j  ,:) = V(ip2,j  ,:)+W(n,m,:)*state%cf20(n,m)
     V(ip3,j  ,:) = V(ip3,j  ,:)+W(n,m,:)*state%cf30(n,m)
     V(i  ,jp1,:) = V(i  ,jp1,:)+W(n,m,:)*state%cf01(n,m) 
     V(ip1,jp1,:) = V(ip1,jp1,:)+W(n,m,:)*state%cf11(n,m)
     V(ip2,jp1,:) = V(ip2,jp1,:)+W(n,m,:)*state%cf21(n,m)
     V(ip3,jp1,:) = V(ip3,jp1,:)+W(n,m,:)*state%cf31(n,m)
     V(i  ,jp2,:) = V(i  ,jp2,:)+W(n,m,:)*state%cf02(n,m)
     V(ip1,jp2,:) = V(ip1,jp2,:)+W(n,m,:)*state%cf12(n,m)
     V(ip2,jp2,:) = V(ip2,jp2,:)+W(n,m,:)*state%cf22(n,m)
     V(ip3,jp2,:) = V(ip3,jp2,:)+W(n,m,:)*state%cf32(n,m)
     V(i  ,jp3,:) = V(i  ,jp3,:)+W(n,m,:)*state%cf03(n,m)
     V(ip1,jp3,:) = V(ip1,jp3,:)+W(n,m,:)*state%cf13(n,m)
     V(ip2,jp3,:) = V(ip2,jp3,:)+W(n,m,:)*state%cf23(n,m)
     V(ip3,jp3,:) = V(ip3,jp3,:)+W(n,m,:)*state%cf33(n,m)
   enddo
   enddo

!$OMP END PARALLEL DO 
!-----------------------------------------------------------------------
                        endsubroutine lsqr_adjoint_k2d  

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine lsqr_adjoint_k3d                      &
!***********************************************************************
!                                                                      !
! Given a target array W(0:grid%nm,0:grid%mm,grid%lm,grid%km3) perform adjoint             !
! interpolations in order to get source V(-grid%ib:grid%im+grid%ib,-grid%jb:grid%jm+grid%jb,grid%lm,grid%km3)  ! 
!                                                                      !
!***********************************************************************
(grid,state,W,V)
!-----------------------------------------------------------------------
implicit none
type(mgbf_grid_type),intent(in) :: grid
type(mgbf_state_type),intent(in) :: state
real(kind_real), dimension(0:grid%nm,0:grid%mm,grid%lm,grid%km3),intent(in):: W  
real(kind_real), dimension(-grid%ib:grid%im+grid%ib,-grid%jb:grid%jm+grid%jb,grid%lm,grid%km3), intent(out):: V
integer(kind_int):: i,j,n,m,l
integer(kind_int):: ip1,ip2,ip3,jp1,jp2,jp3
!-----------------------------------------------------------------------
     
    V(:,:,:,:) = 0.

   do l=1,grid%lm
   do m=0,grid%mm
       j = state%jref(m)
       jp1=j+1
       jp2=j+2
       jp3=j+3
   do n=0,grid%nm
       i = state%iref(n) 
       ip1=i+1
       ip2=i+2
       ip3=i+3
     V(i  ,j  ,l,:) = V(i  ,j  ,l,:)+W(n,m,l,:)*state%cf00(n,m)
     V(ip1,j  ,l,:) = V(ip1,j  ,l,:)+W(n,m,l,:)*state%cf10(n,m)
     V(ip2,j  ,l,:) = V(ip2,j  ,l,:)+W(n,m,l,:)*state%cf20(n,m)
     V(ip3,j  ,l,:) = V(ip3,j  ,l,:)+W(n,m,l,:)*state%cf30(n,m)
     V(i  ,jp1,l,:) = V(i  ,jp1,l,:)+W(n,m,l,:)*state%cf01(n,m) 
     V(ip1,jp1,l,:) = V(ip1,jp1,l,:)+W(n,m,l,:)*state%cf11(n,m)
     V(ip2,jp1,l,:) = V(ip2,jp1,l,:)+W(n,m,l,:)*state%cf21(n,m)
     V(ip3,jp1,l,:) = V(ip3,jp1,l,:)+W(n,m,l,:)*state%cf31(n,m)
     V(i  ,jp2,l,:) = V(i  ,jp2,l,:)+W(n,m,l,:)*state%cf02(n,m)
     V(ip1,jp2,l,:) = V(ip1,jp2,l,:)+W(n,m,l,:)*state%cf12(n,m)
     V(ip2,jp2,l,:) = V(ip2,jp2,l,:)+W(n,m,l,:)*state%cf22(n,m)
     V(ip3,jp2,l,:) = V(ip3,jp2,l,:)+W(n,m,l,:)*state%cf32(n,m)
     V(i  ,jp3,l,:) = V(i  ,jp3,l,:)+W(n,m,l,:)*state%cf03(n,m)
     V(ip1,jp3,l,:) = V(ip1,jp3,l,:)+W(n,m,l,:)*state%cf13(n,m)
     V(ip2,jp3,l,:) = V(ip2,jp3,l,:)+W(n,m,l,:)*state%cf23(n,m)
     V(ip3,jp3,l,:) = V(ip3,jp3,l,:)+W(n,m,l,:)*state%cf33(n,m)
   enddo
   enddo
   enddo

!-----------------------------------------------------------------------
                        endsubroutine lsqr_adjoint_k3d  

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine lsqr_forward_x_k2d                   &
!***********************************************************************
!                                                                      !
! Given a source array  V(-grid%ib:grid%im+grid%ib,0:grid%mm,grid%km2) perform forward          !
! interpolations to get target array W(0:grid%nm,0:grid%mm,grid%km2)                  ! 
! only in x direction                                                  !
!                                                                      !
!***********************************************************************
(grid,state,V,W)
!-----------------------------------------------------------------------
implicit none
type(mgbf_grid_type),intent(in) :: grid
type(mgbf_state_type),intent(in) :: state
real(kind_real), dimension(-grid%ib:grid%im+grid%ib,0:grid%mm,grid%km2), intent(in):: V
real(kind_real), dimension(0:grid%nm,0:grid%mm,grid%km2),intent(out):: W  
integer(kind_int):: i,j,n,m,l
real(kind_real),dimension(grid%km2):: v0,v1,v2,v3     
!-----------------------------------------------------------------------

   do m=0,grid%mm
       j = m
   do n=0,grid%nm
       i = state%iref(n)
     v0(:)=V(i  ,j  ,:)
     v1(:)=V(i+1,j  ,:)
     v2(:)=V(i+2,j  ,:)
     v3(:)=V(i+3,j  ,:)
     W(n,m,:) = state%cx0(n)*v0(:)+state%cx1(n)*v1(:)+state%cx2(n)*v2(:)+state%cx3(n)*v3(:)
   enddo
   enddo


!-----------------------------------------------------------------------
                        endsubroutine lsqr_forward_x_k2d   

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine lsqr_forward_x_k3d                   &
!***********************************************************************
!                                                                      !
! Given a source array  V(-grid%ib:grid%im+grid%ib,-grid%jb:grid%jm+grid%jb,grid%lm,grid%km3) perform forward  !
! interpolations to get target array W(0:grid%nm,0:grid%mm,grid%lm,grid%km3)               ! 
! only in x direction                                                  !
!                                                                      !
!***********************************************************************
(grid,state,V,W)
!-----------------------------------------------------------------------
implicit none
type(mgbf_grid_type),intent(in) :: grid
type(mgbf_state_type),intent(in) :: state
real(kind_real), dimension(-grid%ib:grid%im+grid%ib,0:grid%mm,grid%lm,grid%km3), intent(in):: V
real(kind_real), dimension(0:grid%nm,0:grid%mm,grid%lm,grid%km3),intent(out):: W  
integer(kind_int):: i,j,n,m,l
real(kind_real),dimension(grid%km3):: v0,v1,v2,v3 
!-----------------------------------------------------------------------

   do l=1,grid%lm
   do m=0,grid%mm
       j = m
   do n=0,grid%nm
       i = state%iref(n)
     v0(:)=V(i  ,j  ,l,:)
     v1(:)=V(i+1,j  ,l,:)
     v2(:)=V(i+2,j  ,l,:)
     v3(:)=V(i+3,j  ,l,:)
     W(n,m,l,:) = state%cx0(n)*v0(:)+state%cx1(n)*v1(:)+state%cx2(n)*v2(:)+state%cx3(n)*v3(:)  
   enddo
   enddo
   enddo

!-----------------------------------------------------------------------
                        endsubroutine lsqr_forward_x_k3d   

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine lsqr_adjoint_x_k2d                   &
!***********************************************************************
!                                                                      !
! Given a target array W(0:grid%nm,0:grid%mm,grid%km2) perform adjoint                !
! interpolations in order to get source V(-grid%ib:grid%im+grid%ib,0:grid%mm,grid%km2)          ! 
! only in x direction                                                  !
!                                                                      !
!***********************************************************************
(grid,state,W,V)
!-----------------------------------------------------------------------
implicit none
type(mgbf_grid_type),intent(in) :: grid
type(mgbf_state_type),intent(in) :: state
real(kind_real), dimension(0:grid%nm,0:grid%mm,grid%km2),intent(in):: W  
real(kind_real), dimension(-grid%ib:grid%im+grid%ib,0:grid%mm,grid%km2), intent(out):: V
integer(kind_int):: i,j,n,m
integer(kind_int):: ip1,ip2,ip3
!-----------------------------------------------------------------------
     
    V(:,:,:) = 0.

   do m=0,grid%mm
       j = m
   do n=0,grid%nm
       i = state%iref(n) 
       ip1=i+1
       ip2=i+2
       ip3=i+3
       j = m
     V(i  ,j  ,:) = V(i  ,j  ,:)+W(n,m,:)*state%cx0(n)
     V(ip1,j  ,:) = V(ip1,j  ,:)+W(n,m,:)*state%cx1(n)
     V(ip2,j  ,:) = V(ip2,j  ,:)+W(n,m,:)*state%cx2(n)
     V(ip3,j  ,:) = V(ip3,j  ,:)+W(n,m,:)*state%cx3(n)
   enddo
   enddo

!-----------------------------------------------------------------------
                        endsubroutine lsqr_adjoint_x_k2d  

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine lsqr_adjoint_x_k3d                   &
!***********************************************************************
!                                                                      !
! Given a target array W(0:grid%nm,0:grid%mm,grid%lm,grid%km3) perform adjoint             !
! interpolations in order to get source V(-grid%ib:grid%im+grid%ib,0:grid%mm,grid%lm,grid%km3)       ! 
! only in x direction                                                  !
!                                                                      !
!***********************************************************************
(grid,state,W,V)
!-----------------------------------------------------------------------
implicit none
type(mgbf_grid_type),intent(in) :: grid
type(mgbf_state_type),intent(in) :: state
real(kind_real), dimension(0:grid%nm,0:grid%mm,grid%lm,grid%km3),intent(in):: W  
real(kind_real), dimension(-grid%ib:grid%im+grid%ib,0:grid%mm,grid%lm,grid%km3), intent(out):: V
integer(kind_int):: i,j,n,m,l
integer(kind_int):: ip1,ip2,ip3
!-----------------------------------------------------------------------
     
    V(:,:,:,:) = 0.

   do l=1,grid%lm
   do m=0,grid%mm
       j = m
   do n=0,grid%nm
       i = state%iref(n) 
       ip1=i+1
       ip2=i+2
       ip3=i+3
       j = m
     V(i  ,j  ,l,:) = V(i  ,j  ,l,:)+W(n,m,l,:)*state%cx0(n)
     V(ip1,j  ,l,:) = V(ip1,j  ,l,:)+W(n,m,l,:)*state%cx1(n)
     V(ip2,j  ,l,:) = V(ip2,j  ,l,:)+W(n,m,l,:)*state%cx2(n)
     V(ip3,j  ,l,:) = V(ip3,j  ,l,:)+W(n,m,l,:)*state%cx3(n)
   enddo
   enddo
   enddo

!-----------------------------------------------------------------------
                        endsubroutine lsqr_adjoint_x_k3d  

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine lsqr_forward_xy_k2d                  &
!***********************************************************************
!                                                                      !
!  Given a source array  V(-grid%ib:grid%im+grid%ib,-grid%jb:grid%jm+grid%jb,grid%km2) perform forward    !
!  interpolations to get target array W(0:grid%nm,0:grid%mm,grid%km2)                 ! 
!  using two passes of 1d interpolator                                 !
!                                                                      !
!***********************************************************************
(grid,state,V,W)
!-----------------------------------------------------------------------
implicit none
type(mgbf_grid_type),intent(in) :: grid
type(mgbf_state_type),intent(in) :: state
real(kind_real), dimension(-grid%ib:grid%im+grid%ib,-grid%jb:grid%jm+grid%jb,grid%km2), intent(in):: V
real(kind_real), dimension(0:grid%nm,0:grid%mm,grid%km2),intent(out):: W  

real(kind_real), dimension(0:grid%nm,-grid%jb:grid%jm+grid%jb,grid%km2):: VX
integer(kind_int):: i,j,n,m,l
real(kind_real),dimension(grid%km2):: v0,v1,v2,v3     
!-----------------------------------------------------------------------

   do j=-grid%jb,grid%jm+grid%jb
   do n=0,grid%nm
       i = state%iref(n)
     v0(:)=V(i  ,j  ,:)
     v1(:)=V(i+1,j  ,:)
     v2(:)=V(i+2,j  ,:)
     v3(:)=V(i+3,j  ,:)
     VX(n,j,:) = state%cx0(n)*v0(:)+state%cx1(n)*v1(:)+state%cx2(n)*v2(:)+state%cx3(n)*v3(:)
   enddo
   enddo

   do m=0,grid%mm
     j = state%jref(m)
   do n=0,grid%nm
     v0(:)=VX(n,j  ,:) 
     v1(:)=VX(n,j+1,:) 
     v2(:)=VX(n,j+2,:) 
     v3(:)=VX(n,j+3,:) 
     W(n,m,:) =  state%cy0(m)*v0(:)+state%cy1(m)*v1(:)+state%cy2(m)*v2(:)+state%cy3(m)*v3(:)
   enddo
   enddo


!-----------------------------------------------------------------------
                        endsubroutine lsqr_forward_xy_k2d   

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine lsqr_forward_xy_k3d                  &
!***********************************************************************
!                                                                      !
! Given a source array  V(-grid%ib:grid%im+grid%ib,-grid%jb:grid%jm+grid%jb,grid%lm,grid%km3) perform forward  !
! interpolations to get target array W(0:grid%nm,0:grid%mm,grid%lm,grid%km3)               ! 
! using two passes of 1d interpolator                                  !
!                                                                      !
!***********************************************************************
(grid,state,V,W)
!-----------------------------------------------------------------------
implicit none
type(mgbf_grid_type),intent(in) :: grid
type(mgbf_state_type),intent(in) :: state
real(kind_real), dimension(-grid%ib:grid%im+grid%ib,-grid%jb:grid%jm+grid%jb,grid%lm,grid%km3), intent(in):: V
real(kind_real), dimension(0:grid%nm,0:grid%mm,grid%lm,grid%km3),intent(out):: W  

real(kind_real), dimension(0:grid%nm,-grid%jb:grid%jm+grid%jb,grid%km3):: VX
integer(kind_int):: i,j,n,m,l
real(kind_real),dimension(grid%km3):: v0,v1,v2,v3     
!-----------------------------------------------------------------------

 do l=1,grid%lm

   do j=-grid%jb,grid%jm+grid%jb
   do n=0,grid%nm
       i = state%iref(n)
     v0(:)=V(i  ,j  ,l,:)
     v1(:)=V(i+1,j  ,l,:)
     v2(:)=V(i+2,j  ,l,:)
     v3(:)=V(i+3,j  ,l,:)
     VX(n,j,:) = state%cx0(n)*v0(:)+state%cx1(n)*v1(:)+state%cx2(n)*v2(:)+state%cx3(n)*v3(:)
   enddo
   enddo

   do m=0,grid%mm
     j = state%jref(m)
   do n=0,grid%nm
     v0(:)=VX(n,j  ,:) 
     v1(:)=VX(n,j+1,:) 
     v2(:)=VX(n,j+2,:) 
     v3(:)=VX(n,j+3,:) 
     W(n,m,l,:) =  state%cy0(m)*v0(:)+state%cy1(m)*v1(:)+state%cy2(m)*v2(:)+state%cy3(m)*v3(:)
   enddo
   enddo

 enddo

!-----------------------------------------------------------------------
                        endsubroutine lsqr_forward_xy_k3d   

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine lsqr_adjoint_xy_k2d                  &
!***********************************************************************
!                                                                      !
! Given a target array W(0:grid%nm,0:grid%mm,grid%km2) perform adjoint                !
! interpolations to get source array V(-grid%ib:grid%im+grid%ib,-grid%jb:grid%jm+grid%jb,grid%km2)        ! 
! using two passes of 1d interpolator                                  !
!                                                                      !
!***********************************************************************
(grid,state,W,V)
!-----------------------------------------------------------------------
implicit none
type(mgbf_grid_type),intent(in) :: grid
type(mgbf_state_type),intent(in) :: state
real(kind_real), dimension(0:grid%nm,0:grid%mm,grid%km2),intent(in):: W  
real(kind_real), dimension(-grid%ib:grid%im+grid%ib,-grid%jb:grid%jm+grid%jb,grid%km2), intent(out):: V
real(kind_real), dimension(0:grid%nm,-grid%jb:grid%jm+grid%jb,grid%km2):: VX
integer(kind_int):: i,j,n,m,l
integer(kind_int):: ip1,ip2,ip3
integer(kind_int):: jp1,jp2,jp3
!-----------------------------------------------------------------------
   
   VX(:,:,:)=0.

   do m=0,grid%mm
     j = state%jref(m)
     jp1=j+1
     jp2=j+2
     jp3=j+3
   do n=0,grid%nm
     VX(n,j  ,:) = VX(n,j  ,:)+W(n,m,:)*state%cy0(m)
     VX(n,jp1,:) = VX(n,jp1,:)+W(n,m,:)*state%cy1(m)
     VX(n,jp2,:) = VX(n,jp2,:)+W(n,m,:)*state%cy2(m)
     VX(n,jp3,:) = VX(n,jp3,:)+W(n,m,:)*state%cy3(m)
   enddo
   enddo
 
   V(:,:,:) = 0.

   do j=-grid%jb,grid%jm+grid%jb
   do n=0,grid%nm
     i = state%iref(n)
     ip1=i+1
     ip2=i+2
     ip3=i+3

     V(i  ,j,:) = V(i  ,j,:)+VX(n,j,:)*state%cx0(n)
     V(ip1,j,:) = V(ip1,j,:)+VX(n,j,:)*state%cx1(n)
     V(ip2,j,:) = V(ip2,j,:)+VX(n,j,:)*state%cx2(n)
     V(ip3,j,:) = V(ip3,j,:)+VX(n,j,:)*state%cx3(n)
   enddo
   enddo


!-----------------------------------------------------------------------
                        endsubroutine lsqr_adjoint_xy_k2d   

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine lsqr_adjoint_xy_k3d                  &
!***********************************************************************
!                                                                      !
! Given a target array W(0:grid%nm,0:grid%mm,grid%lm,grid%km3) perform adjoint             !
! interpolations to get source array V(-grid%ib:grid%im+grid%ib,-grid%jb:grid%jm+grid%jb,grid%lm,grid%km3)     ! 
! using two passes of 1d interpolator                                  !
!                                                                      !
!***********************************************************************
(grid,state,W,V)
!-----------------------------------------------------------------------
implicit none
type(mgbf_grid_type),intent(in) :: grid
type(mgbf_state_type),intent(in) :: state
real(kind_real), dimension(0:grid%nm,0:grid%mm,grid%lm,grid%km3),intent(in):: W  
real(kind_real), dimension(-grid%ib:grid%im+grid%ib,-grid%jb:grid%jm+grid%jb,grid%lm,grid%km3), intent(out):: V
real(kind_real), dimension(0:grid%nm,-grid%jb:grid%jm+grid%jb,grid%km3):: VX
integer(kind_int):: i,j,n,m,l,k
integer(kind_int):: ip1,ip2,ip3
integer(kind_int):: jp1,jp2,jp3
!-----------------------------------------------------------------------

   V(:,:,:,:) = 0.

!TEST
     write(60,*)'Bih ovde'
!TEST
   
 do l=1,grid%lm

   VX(:,:,:)=0.

 do k=1,grid%km3
   do m=0,grid%mm
     j = state%jref(m)
     jp1=j+1
     jp2=j+2
     jp3=j+3
   do n=0,grid%nm
     VX(n,j  ,k) = VX(n,j  ,k)+W(n,m,l,k)*state%cy0(m)
     VX(n,jp1,k) = VX(n,jp1,k)+W(n,m,l,k)*state%cy1(m)
     VX(n,jp2,k) = VX(n,jp2,k)+W(n,m,l,k)*state%cy2(m)
     VX(n,jp3,k) = VX(n,jp3,k)+W(n,m,l,k)*state%cy3(m)
   enddo
   enddo
 end do
 

 do k=1,grid%km3
   do j=-grid%jb,grid%jm+grid%jb
   do n=0,grid%nm
     i = state%iref(n)
     ip1=i+1
     ip2=i+2
     ip3=i+3

     V(i  ,j,l,k) = V(i  ,j,l,k)+VX(n,j,k)*state%cx0(n)
     V(ip1,j,l,k) = V(ip1,j,l,k)+VX(n,j,k)*state%cx1(n)
     V(ip2,j,l,k) = V(ip2,j,l,k)+VX(n,j,k)*state%cx2(n)
     V(ip3,j,l,k) = V(ip3,j,l,k)+VX(n,j,k)*state%cx3(n)
   enddo
   enddo
 enddo

 enddo
!-----------------------------------------------------------------------
                        endsubroutine lsqr_adjoint_xy_k3d   


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        endmodule tools_interpolate
