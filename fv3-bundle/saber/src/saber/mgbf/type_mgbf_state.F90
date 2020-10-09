module type_mgbf_state

use tools_kinds, only: kind_int,kind_real
use type_mgbf_grid
use tools_betafilter

implicit none

! state derived type
type mgbf_state_type
real(kind_real), allocatable,dimension(:,:,:):: V
!
! Composite control variable on first generation o filter grid
!
real(kind_real), allocatable,dimension(:,:,:):: VK2D
real(kind_real), allocatable,dimension(:,:,:,:):: VK3D
!
! Composite control variable on high generations of filter grid
!
real(kind_real), allocatable,dimension(:,:,:):: HK2D
real(kind_real), allocatable,dimension(:,:,:,:):: HK3D
real(kind_real), allocatable,dimension(:,:,:,:):: a2_diff
real(kind_real), allocatable,dimension(:,:,:,:):: b2_diff
real(kind_real), allocatable,dimension(:,:,:,:,:):: a3_diff
real(kind_real), allocatable,dimension(:,:,:,:,:):: b3_diff

real(kind_real), allocatable,dimension(:,:,:,:):: weig_2d
real(kind_real), allocatable,dimension(:,:,:,:,:):: weig_3d

real(kind_real), allocatable,dimension(:,:):: p_eps
real(kind_real), allocatable,dimension(:,:):: p_del
real(kind_real), allocatable,dimension(:,:):: p_sig
real(kind_real), allocatable,dimension(:,:):: p_rho

real(kind_real), allocatable,dimension(:,:,:):: pasp1
real(kind_real), allocatable,dimension(:,:,:,:):: pasp2
real(kind_real), allocatable,dimension(:,:,:,:,:):: pasp3

real(kind_real), allocatable,dimension(:):: ss1
real(kind_real), allocatable,dimension(:,:):: ss2
real(kind_real), allocatable,dimension(:,:,:):: ss3

! 2D analysis variables
!
real(kind_real), allocatable,dimension(:,:):: PA30

real(kind_real), allocatable,dimension(:,:):: PA1
real(kind_real), allocatable,dimension(:,:):: PA2
real(kind_real), allocatable,dimension(:,:):: PA3
real(kind_real), allocatable,dimension(:,:):: PA4
!

! 3D analysis variables
!
real(kind_real), allocatable,dimension(:,:,:):: WA20

real(kind_real), allocatable,dimension(:,:,:):: WA1
real(kind_real), allocatable,dimension(:,:,:):: WA2
real(kind_real), allocatable,dimension(:,:,:):: WA3
real(kind_real), allocatable,dimension(:,:,:):: WA4
real(kind_real), allocatable,dimension(:,:,:):: WA5
real(kind_real), allocatable,dimension(:,:,:):: WA6

integer(kind_int),allocatable,dimension(:):: iref,jref

real(kind_real),allocatable,dimension(:):: cx0,cx1,cx2,cx3
real(kind_real),allocatable,dimension(:):: cy0,cy1,cy2,cy3

real(kind_real),allocatable,dimension(:,:):: cf00,cf01,cf02,cf03           &
                                         ,cf10,cf11,cf12,cf13           &
                                         ,cf20,cf21,cf22,cf23           &
                                         ,cf30,cf31,cf32,cf33
contains
  procedure :: setup => mgbf_state_setup
  procedure :: allocate => allocate_mg_intstate
  procedure :: def_weights => def_mg_weights
!  procedure :: 
end type mgbf_state_type

contains

subroutine mgbf_state_setup(state,grid)
implicit none

class(mgbf_state_type),intent(inout):: state
type(mgbf_grid_type),intent(in):: grid

call state%allocate(grid)
call state%def_weights(grid)

end subroutine mgbf_state_setup

subroutine allocate_mg_intstate(state,grid)
implicit none

class(mgbf_state_type),intent(inout):: state
type(mgbf_grid_type),intent(in):: grid

allocate(state%V(-grid%hx:grid%im+grid%hx,-grid%hy:grid%jm+grid%hy,grid%lm))                          ; state%V=0.0_kind_real
allocate(state%VK3D(-grid%hx:grid%im+grid%hx,-grid%hy:grid%jm+grid%hy,grid%lm,grid%km3))                   ; state%VK3D=0.0_kind_real
allocate(state%VK2D(-grid%hx:grid%im+grid%hx,-grid%hy:grid%jm+grid%hy,grid%km2))                      ; state%VK2D=0.0_kind_real
allocate(state%HK3D(-grid%hx:grid%im+grid%hx,-grid%hy:grid%jm+grid%hy,grid%lm,grid%km3))                   ; state%HK3D=0.0_kind_real
allocate(state%HK2D(-grid%hx:grid%im+grid%hx,-grid%hy:grid%jm+grid%hy,grid%km2))                      ; state%HK2D=0.0_kind_real

!
!FOR ADJOINT TEST
!
!allocate(state%A(-grid%hx:grid%im+grid%hx,-grid%hy:grid%jm+grid%hy))                              ; state%A=0.0_kind_real
!allocate(state%B(-grid%hx:grid%im+grid%hx,-grid%hy:grid%jm+grid%hy))                              ; state%B=0.0_kind_real
!allocate(state%A0(-grid%hx:grid%im+grid%hx,-grid%hy:grid%jm+grid%hy))                             ; state%A0=0.0_kind_real
!allocate(state%B0(-grid%hx:grid%im+grid%hx,-grid%hy:grid%jm+grid%hy))                             ; state%B0=0.0_kind_real



allocate(state%a2_diff(-grid%hx:grid%im+grid%hx,-grid%hy:grid%jm+grid%hy,1:grid%km2,1:2))                ; state%a2_diff=0.0_kind_real
allocate(state%b2_diff(-grid%hx:grid%im+grid%hx,-grid%hy:grid%jm+grid%hy,1:grid%km2,1:2))                ; state%b2_diff=0.0_kind_real

allocate(state%a3_diff(-grid%hx:grid%im+grid%hx,-grid%hy:grid%jm+grid%hy,1:grid%lm,1:grid%km3,1:2))           ; state%a3_diff=0.0_kind_real
allocate(state%b3_diff(-grid%hx:grid%im+grid%hx,-grid%hy:grid%jm+grid%hy,1:grid%lm,1:grid%km3,1:2))           ; state%b3_diff=0.0_kind_real

allocate(state%weig_2d(-grid%hx:grid%im+grid%hx,-grid%hy:grid%jm+grid%hy,1:grid%km2,1:2))                ; state%weig_2d=0.0_kind_real
allocate(state%weig_3d(-grid%hx:grid%im+grid%hx,-grid%hy:grid%jm+grid%hy,1:grid%lm,1:grid%km3,1:2))           ; state%weig_3d=0.0_kind_real

allocate(state%p_eps(-grid%hx:grid%im+grid%hx,-grid%hy:grid%jm+grid%hy))                            ; state%p_eps=0.0_kind_real
allocate(state%p_del(-grid%hx:grid%im+grid%hx,-grid%hy:grid%jm+grid%hy))                            ; state%p_del=0.0_kind_real
allocate(state%p_sig(-grid%hx:grid%im+grid%hx,-grid%hy:grid%jm+grid%hy))                            ; state%p_sig=0.0_kind_real
allocate(state%p_rho(-grid%hx:grid%im+grid%hx,-grid%hy:grid%jm+grid%hy))                            ; state%p_rho=0.0_kind_real
allocate(state%pasp1(1,1,1:grid%lm))                                       ; state%pasp1=0.0_kind_real
allocate(state%pasp2(2,2,0:grid%im,0:grid%jm))                                  ; state%pasp2=0.0_kind_real
allocate(state%pasp3(3,3,0:grid%im,0:grid%jm,1:grid%lm))                             ; state%pasp3=0.0_kind_real
allocate(state%ss1(1:grid%lm))                                             ; state%ss1=0.0_kind_real
allocate(state%ss2(0:grid%im,0:grid%jm))                                        ; state%ss2=0.0_kind_real
allocate(state%ss3(0:grid%im,0:grid%jm,1:grid%lm))                                   ; state%ss3=0.0_kind_real
!
!
! for testing
!

allocate(state%PA30(1:grid%nm,1:grid%mm))                            ; state%PA30=0.0_kind_real
allocate(state%WA20(1:grid%nm,1:grid%mm,grid%lm))                         ; state%WA20=0.0_kind_real

allocate(state%PA1(0:grid%nm,0:grid%mm))                             ; state%PA1=0.0_kind_real
allocate(state%PA2(0:grid%nm,0:grid%mm))                             ; state%PA2=0.0_kind_real
allocate(state%PA3(0:grid%nm,0:grid%mm))                             ; state%PA2=0.0_kind_real
allocate(state%PA4(0:grid%nm,0:grid%mm))                             ; state%PA2=0.0_kind_real

allocate(state%WA1(0:grid%nm,0:grid%mm,grid%lm))                          ; state%WA1=0.0_kind_real
allocate(state%WA2(0:grid%nm,0:grid%mm,grid%lm))                          ; state%WA2=0.0_kind_real
allocate(state%WA3(0:grid%nm,0:grid%mm,grid%lm))                          ; state%WA3=0.0_kind_real
allocate(state%WA4(0:grid%nm,0:grid%mm,grid%lm))                          ; state%WA4=0.0_kind_real
allocate(state%WA5(0:grid%nm,0:grid%mm,grid%lm))                          ; state%WA5=0.0_kind_real
allocate(state%WA6(0:grid%nm,0:grid%mm,grid%lm))                          ; state%WA6=0.0_kind_real

!
! for re-decomposition
!

allocate(state%iref(0:grid%nm))                                     ; state%iref=0
allocate(state%jref(0:grid%mm))                                     ; state%jref=0

allocate(state%cx0(0:grid%nm))                                      ; state%cx0=0.0_kind_real
allocate(state%cx1(0:grid%nm))                                      ; state%cx1=0.0_kind_real
allocate(state%cx2(0:grid%nm))                                      ; state%cx2=0.0_kind_real
allocate(state%cx3(0:grid%nm))                                      ; state%cx3=0.0_kind_real

allocate(state%cy0(0:grid%mm))                                      ; state%cy0=0.0_kind_real
allocate(state%cy1(0:grid%mm))                                      ; state%cy1=0.0_kind_real
allocate(state%cy2(0:grid%mm))                                      ; state%cy2=0.0_kind_real
allocate(state%cy3(0:grid%mm))                                      ; state%cy3=0.0_kind_real

allocate(state%cf00(0:grid%nm,0:grid%mm))                            ; state%cf00=0.0_kind_real
allocate(state%cf01(0:grid%nm,0:grid%mm))                            ; state%cf01=0.0_kind_real
allocate(state%cf02(0:grid%nm,0:grid%mm))                            ; state%cf02=0.0_kind_real
allocate(state%cf03(0:grid%nm,0:grid%mm))                            ; state%cf03=0.0_kind_real
allocate(state%cf10(0:grid%nm,0:grid%mm))                            ; state%cf10=0.0_kind_real
allocate(state%cf11(0:grid%nm,0:grid%mm))                            ; state%cf11=0.0_kind_real
allocate(state%cf12(0:grid%nm,0:grid%mm))                            ; state%cf12=0.0_kind_real
allocate(state%cf13(0:grid%nm,0:grid%mm))                            ; state%cf13=0.0_kind_real
allocate(state%cf20(0:grid%nm,0:grid%mm))                            ; state%cf20=0.0_kind_real
allocate(state%cf21(0:grid%nm,0:grid%mm))                            ; state%cf21=0.0_kind_real
allocate(state%cf22(0:grid%nm,0:grid%mm))                            ; state%cf22=0.0_kind_real
allocate(state%cf23(0:grid%nm,0:grid%mm))                            ; state%cf23=0.0_kind_real
allocate(state%cf30(0:grid%nm,0:grid%mm))                            ; state%cf30=0.0_kind_real
allocate(state%cf31(0:grid%nm,0:grid%mm))                            ; state%cf31=0.0_kind_real
allocate(state%cf32(0:grid%nm,0:grid%mm))                            ; state%cf32=0.0_kind_real
allocate(state%cf33(0:grid%nm,0:grid%mm))                            ; state%cf33=0.0_kind_real


end subroutine allocate_mg_intstate

subroutine def_mg_weights(state,grid)

implicit none

class(mgbf_state_type),intent(inout):: state
type(mgbf_grid_type),intent(in):: grid

integer(kind_int):: i,j,L

      state%p_eps(:,:)=0.0_kind_real
      state%p_del(:,:)=0.0_kind_real
      state%p_sig(:,:)=0.0_kind_real
      state%p_rho(:,:)=0.0_kind_real
          state%weig_2d(:,:,:,1)=2.0_kind_real
          state%weig_3d(:,:,:,:,1)=2.0_kind_real
      select case(grid%my_hgen)
        case(2)
          state%weig_2d(:,:,:,2)=3.0_kind_real
          state%weig_3d(:,:,:,:,2)=3.0_kind_real
        case(3)
          state%weig_2d(:,:,:,2)=4.0_kind_real
          state%weig_3d(:,:,:,:,2)=4.0_kind_real
        case default
          state%weig_2d(:,:,:,2)=5.0_kind_real
          state%weig_3d(:,:,:,:,2)=5.0_kind_real
      end select

      state%a2_diff(:,:,:,:)=0.1_kind_real
      state%a3_diff(:,:,:,:,:)=0.1_kind_real

      state%b2_diff(:,:,:,1)=1.0_kind_real
      state%b3_diff(:,:,:,:,1)=1.0_kind_real


      select case(grid%my_hgen)
        case(2)
          state%b2_diff(:,:,:,2)=20.0_kind_real
          state%b3_diff(:,:,:,:,2)=20.0_kind_real
        case(3)
          state%b2_diff(:,:,:,2)=40.0_kind_real
          state%b3_diff(:,:,:,:,2)=40.0_kind_real
        case default
          state%b2_diff(:,:,:,2)=60.0_kind_real
          state%b3_diff(:,:,:,:,2)=60.0_kind_real
      end select



          do L=1,grid%lm
            state%pasp1(1,1,L)=grid%pasp0
          enddo

          do j=0,grid%jm
          do i=0,grid%im
            state%pasp2(1,1,i,j)=grid%pasp0*(1.+state%p_del(i,j))
            state%pasp2(2,2,i,j)=grid%pasp0*(1.-state%p_del(i,j))
            state%pasp2(1,2,i,j)=grid%pasp0*state%p_eps(i,j)
            state%pasp2(2,1,i,j)=grid%pasp0*state%p_eps(i,j)
          end do
          end do
        do L=1,grid%lm
          do j=0,grid%jm
          do i=0,grid%im
            state%pasp3(1,1,i,j,l)=grid%pasp0*(1+state%p_del(i,j))
            state%pasp3(2,2,i,j,l)=grid%pasp0
            state%pasp3(3,3,i,j,l)=grid%pasp0*(1-state%p_del(i,j))
            state%pasp3(1,2,i,j,l)=grid%pasp0*state%p_eps(i,j)
            state%pasp3(2,1,i,j,l)=grid%pasp0*state%p_eps(i,j)
            state%pasp3(2,3,i,j,l)=grid%pasp0*state%p_sig(i,j)
            state%pasp3(3,2,i,j,l)=grid%pasp0*state%p_sig(i,j)
            state%pasp3(1,3,i,j,l)=grid%pasp0*state%p_rho(i,j)
            state%pasp3(3,1,i,j,l)=grid%pasp0*state%p_rho(i,j)
          end do
          end do
        end do

        call cholaspect1(1,grid%lm,state%pasp1)
        call cholaspect2(0,grid%im,0,grid%jm,state%pasp2)
        call cholaspect3(0,grid%im,0,grid%jm,1,grid%lm,state%pasp3)

        call getlinesum1(grid%hz,1,grid%lm,grid%p,state%pasp1,state%ss1)
        call getlinesum2(grid%hx,0,grid%im,grid%hy,0,grid%jm,grid%p,state%pasp2,state%ss2)
        call getlinesum3(grid%hx,0,grid%im,grid%hy,0,grid%jm,grid%hz,1,grid%lm,grid%p,state%pasp3,state%ss3)

end subroutine def_mg_weights

end module type_mgbf_state
