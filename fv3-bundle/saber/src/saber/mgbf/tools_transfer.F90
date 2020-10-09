!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        module tools_transfer
!***********************************************************************
!                                                                      !
!  Transfer data between analysis and filter grid                      !
!                                                                      !
!                                                     M. Rancic (2020) !
!***********************************************************************
use tools_kinds, only: kind_real,kind_int
use tools_interpolate,only: lsqr_adjoint_xyk,lsqr_forward_xyk
use tools_interpolate,only: lsqr_adjoint,lsqr_forward
use tools_bocos, only: bocosHn,bocosHnT,v02v
use type_mgbf_grid
use type_mgbf_state
use type_mpl
!TEST


implicit none
integer(kind_int):: n,m,l,k,i,j

public anal_to_filt1_2d
public filt_to_anal1_2d

public anal_to_filt2_2d
public filt_to_anal2_2d

public anal_to_filt2_2d_1_3d
public filt_to_anal2_2d_1_3d

public anal_to_filt
public filt_to_anal


interface compose
   module procedure compose2d_2,compose2d_3,compose2d_4                 &
                   ,compose3d_2,compose3d_3,compose3d_4                 &
                   ,compose3d_5,compose3d_6
endinterface

interface decompose
   module procedure decompose2d_2,decompose2d_3,decompose2d_4           &
                   ,decompose3d_2,decompose3d_3,decompose3d_4           &
                   ,decompose3d_5,decompose3d_6
endinterface

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine anal_to_filt1_2d &                   
!***********************************************************************
!                                                                      !
!  Transfer data from analysis to first generaton of filter grid       !
!  Only 1 2d variable                                                  !
!                                                                      !
!***********************************************************************
(mpl,grid,state)
implicit none

type(mpl_type),intent(inout) :: mpl
type(mgbf_grid_type),intent(inout) :: grid
type(mgbf_state_type),intent(inout) :: state
real(kind_real),allocatable,dimension(:,:,:):: PKA
real(kind_real),allocatable,dimension(:,:,:):: PKF

!TEST
real(kind_real),allocatable,dimension(:,:):: PF
!TEST

!----------------------------------------------------------------------
          if(grid%im==grid%nm .and. grid%jm==grid%mm) then
             state%VK2D(0:grid%im,0:grid%jm,1)=state%PA1(0:grid%nm,0:grid%mm)
             return
          endif
!***
!*** Create composite variables: state%PA1,state%PA2 -> PKA 
!***

    allocate(PKA(0:grid%nm,0:grid%mm,   grid%km2))                         

!TEST
    PKA(:,:,1)=state%PA1(:,:)
!TEST
!
!    call compose(state%PA1,state%PA2,PKA,nm,mm,km2)


!***
!*** Adjoint interpolate and build composite variables
!***
    allocate(PKF(-grid%ib:grid%im+grid%ib,-grid%jb:grid%jm+grid%jb,1))


!                                                 call btim(  aintp_tim)
!T        call lsqr_adjoint(PKA,PKF)   
        call lsqr_adjoint_xyk(grid,state,PKA,PKF)   
!                                                 call etim(  aintp_tim)

    deallocate(PKA)

!***
!***  Apply adjoint lateral bc on PKF and WKF
!***
    
         call bocosHnT(mpl,grid,PKF,grid%km2,grid%im,grid%jm,1 ,grid%ib,grid%jb,grid%Fimax,grid%Fjmax,1,1)  

!TEST
!        allocate(PF(0:im,0:jm))
!        PF(0:im,0:jm)=PKF(0:im,0:jm,1)
!        call out_2d(PF,0,0,im,jm,'f')
!        call finishMPI
!TEST

!***
!***  Form state%VK2D and state%VK3D
!***
      state%VK2D=0

      state%VK2D(0:grid%im,0:grid%jm,1:grid%km2)=PKF(0:grid%im,0:grid%jm,1:grid%km2)

    deallocate(PKF)

!----------------------------------------------------------------------
                        endsubroutine anal_to_filt1_2d

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine filt_to_anal1_2d &
!***********************************************************************
!                                                                      !
!  Transfer data from filter to analysis grid in x-direction only      !
!  2 2d variables                                                      !
!                                                                      !
!***********************************************************************
(mpl,grid,state)
implicit none

type(mpl_type),intent(inout) :: mpl
type(mgbf_grid_type),intent(inout) :: grid
type(mgbf_state_type),intent(inout) :: state

real(kind_real),allocatable,dimension(:,:,:):: PKF

real(kind_real),allocatable,dimension(:,:,:):: PKA
!TEST
real(kind_real),allocatable,dimension(:,:):: PF
!TEST

!----------------------------------------------------------------------
          if(grid%im==grid%nm .and. grid%jm==grid%mm) then
             state%PA1(0:grid%nm,0:grid%mm)=state%VK2D(0:grid%im,0:grid%jm,1 )
             return
          endif
!***
!***  Define PKF and WKF
!***
    allocate(PKF(-grid%ib:grid%im+grid%ib,-grid%jb:grid%jm+grid%jb,  1))

      PKF(0:grid%im,0:grid%jm,   1)= state%VK2D(0:grid%im,0:grid%jm,   1)


!***
!***  Supply boundary conditions for PKF, WKF
!***
         call bocosHn(mpl,grid,PKF,grid%km2,grid%im,grid%jm,1 ,grid%ib,grid%jb,grid%Fimax,grid%Fjmax,1,1)

!TEST
        allocate(PF(0:grid%im,0:grid%jm))
        PF(0:grid%im,0:grid%jm)=PKF(0:grid%im,0:grid%jm,1)
!        call out_2d(PF,0,0,im,jm,'f')
!        call finishMPI
!TEST

!***
!*** Interpolate to analysis grid composite variables
!***

    allocate(PKA(0:grid%nm,0:grid%mm,   grid%km2))

!                                                 call btim(   intp_tim)
         call lsqr_forward_xyk(grid,state,PKF,PKA)
!                                                 call etim(   intp_tim)
    deallocate(PKF)
!***
!*** Decompose composite variables: PA1,PA2 <- PKA
!***

!    call decompose(PKA,PA1,PA2,nm,mm,km2)

     state%PA1(0:grid%nm,0:grid%mm)=PKA(0:grid%nm,0:grid%mm,1)

    deallocate(PKA)


!----------------------------------------------------------------------
                        endsubroutine filt_to_anal1_2d
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine anal_to_filt2_2d &
!***********************************************************************
!                                                                      !
!  Transfer data from analysis to first generaton of filter grid       !
!  Only 2 2d variables                                                 !
!                                                                      !
!**********************************************************************
(mpl,grid,state)
implicit none

type(mpl_type),intent(inout) :: mpl
type(mgbf_grid_type),intent(inout) :: grid
type(mgbf_state_type),intent(inout) :: state

real(kind_real),allocatable,dimension(:,:,:):: PKA

real(kind_real),allocatable,dimension(:,:,:):: PKF

!----------------------------------------------------------------------
          if(grid%im==grid%nm .and. grid%jm==grid%mm) then
             state%VK2D(0:grid%im,0:grid%jm,1 )=state%PA1(0:grid%nm,0:grid%mm)
             state%VK2D(0:grid%im,0:grid%jm,2 )=state%PA2(0:grid%nm,0:grid%mm)
             return
          endif
!***
!*** Create composite variables: state%PA1,state%PA2 -> PKA
!***

    allocate(PKA(0:grid%nm,0:grid%mm,   grid%km2))

    call compose(state%PA1,state%PA2,PKA,grid%nm,grid%mm,grid%km2)


!***
!*** Adjoint interpolate and build composite variables
!***
    allocate(PKF(-grid%ib:grid%im+grid%ib,-grid%jb:grid%jm+grid%jb,   grid%km2))


!                                                 call btim(  aintp_tim)
        call lsqr_adjoint_xyk(grid,state,PKA,PKF)
!                                                 call etim(  aintp_tim)

    deallocate(PKA)

!***
!***  Apply adjoint lateral bc on PKF and WKF
!***

         call bocosHnT(mpl,grid,PKF,grid%km2,grid%im,grid%jm,1 ,grid%ib,grid%jb,grid%Fimax,grid%Fjmax,1,1)

!***
!***  Form state%VK2D and state%VK3D
!***
      state%VK2D=0

      state%VK2D(0:grid%im,0:grid%jm,1:grid%km2)=PKF(0:grid%im,0:grid%jm,1:grid%km2)

    deallocate(PKF)

!----------------------------------------------------------------------
                        endsubroutine anal_to_filt2_2d

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine filt_to_anal2_2d &
!***********************************************************************
!                                                                      !
!  Transfer data from filter to analysis grid in x-direction only      !
!  2 2d variables                                                      !
!                                                                      !
!***********************************************************************
(mpl,grid,state)
implicit none

type(mpl_type),intent(inout) :: mpl
type(mgbf_grid_type),intent(inout) :: grid
type(mgbf_state_type),intent(inout) :: state

real(kind_real),allocatable,dimension(:,:,:):: PKF

real(kind_real),allocatable,dimension(:,:,:):: PKA

!----------------------------------------------------------------------
          if(grid%im==grid%nm .and. grid%jm==grid%mm) then
             state%PA1(0:grid%nm,0:grid%mm)=state%VK2D(0:grid%im,0:grid%jm,1 )
             state%PA2(0:grid%nm,0:grid%mm)=state%VK2D(0:grid%im,0:grid%jm,2 )
             return
          endif
!***
!***  Define PKF and WKF
!***
    allocate(PKF(-grid%ib:grid%im+grid%ib,-grid%jb:grid%jm+grid%jb,1:grid%km2))

      PKF(0:grid%im,0:grid%jm,1:grid%km2)=state%VK2D(0:grid%im,0:grid%jm,1:grid%km2)


!***
!***  Supply boundary conditions for PKF, WKF
!***
         call bocosHn(mpl,grid,PKF,grid%km2,grid%im,grid%jm,1,grid%ib,grid%jb,grid%Fimax,grid%Fjmax,1,1)

!***
!*** Interpolate to analysis grid composite variables
!***

    allocate(PKA(0:grid%nm,0:grid%mm,1:grid%km2))

!                                                 call btim(   intp_tim)
         call lsqr_forward_xyk(grid,state,PKF,PKA)
!                                                 call etim(   intp_tim)
    deallocate(PKF)
!***
!*** Decompose composite variables: state%PA1,state%PA2 <- PKA
!***

    call decompose(PKA,state%PA1,state%PA2,grid%nm,grid%mm,grid%km2)


    deallocate(PKA)


!----------------------------------------------------------------------
                        endsubroutine filt_to_anal2_2d

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine anal_to_filt2_2d_1_3d &
!***********************************************************************
!                                                                      !
!  Transfer data from analysis to first generaton of filter grid       !
!  (2 2D variables and 1 3D variable)
!                                                                      !
!***********************************************************************
(mpl,grid,state)
implicit none

type(mpl_type),intent(inout) :: mpl
type(mgbf_grid_type),intent(inout) :: grid
type(mgbf_state_type),intent(inout) :: state

real(kind_real),allocatable,dimension(:,:,:):: PKA
real(kind_real),allocatable,dimension(:,:,:,:):: WKA

real(kind_real),allocatable,dimension(:,:,:):: PKF
real(kind_real),allocatable,dimension(:,:,:,:):: WKF

!----------------------------------------------------------------------
!TEST
      if(grid%km2==2 .and. grid%km3==1 .and. grid%nm==grid%im .and. grid%mm==grid%jm) then
        state%VK2D(0:grid%im,0:grid%jm,1)=state%PA1(0:grid%im,0:grid%jm)
        state%VK2D(0:grid%im,0:grid%jm,2)=state%PA2(0:grid%im,0:grid%jm)
        state%VK3D(0:grid%im,0:grid%jm,1:grid%lm,1)=state%WA1(0:grid%im,0:grid%jm,1:grid%lm)
          return
      endif
!TEST
!
! Make sure that subdomains of all variables share cogrid%mmon edges
!(We assume those are state%PA30 and state%WA30. Adjust that in GSI! )

!            state%PA3(1:grid%nm,1:grid%mm)=state%PA30(1:grid%nm,1:grid%mm)
!            state%WA2(1:grid%nm,1:grid%mm,1:grid%lm)=state%WA20(1:grid%nm,1:grid%mm,1:grid%lm)
!
!          call v02v(state%PA3,grid%nm,grid%mm,1)
!          call v02v(state%WA2,grid%nm,grid%mm,grid%lm)

!***
!*** Create composite variables: state%PA1, ..., PA4, state%WA1, ..., state%WA6
!-> PKA and WKA
!***

    allocate(PKA(0:grid%nm,0:grid%mm,   grid%km2))
    allocate(WKA(0:grid%nm,0:grid%mm,grid%lm,grid%km3))

    call compose(state%PA1,state%PA2,PKA,grid%nm,grid%mm,grid%km2)
!    call
!    compose(state%WA1,state%WA2,state%WA3,state%WA4,state%WA5,state%WA6,WKA,grid%nm,grid%mm,grid%lm,grid%km3)

          WKA(0:grid%nm,0:grid%mm,1:grid%lm,1)=state%WA1(0:grid%nm,0:grid%mm,1:grid%lm)


!***
!*** Adjoint interpolate and build composite variables
!***
    allocate(PKF(-grid%ib:grid%im+grid%ib,-grid%jb:grid%jm+grid%jb,   grid%km2))
    allocate(WKF(-grid%ib:grid%im+grid%ib,-grid%jb:grid%jm+grid%jb,grid%lm,grid%km3))


!                                                 call btim(  aintp_tgrid%im)
        call lsqr_adjoint_xyk(grid,state,PKA,PKF)
        call lsqr_adjoint_xyk(grid,state,WKA,WKF)
!                                                call etim(  aintp_tgrid%im)
!TEST
!        call out_3d2(WKF,-grid%ib,-grid%jb,grid%im+grid%ib,grid%jm+grid%jb,grid%lm,grid%lm/2,'w')
!        call finishMPI
!TEST

    deallocate(PKA)
    deallocate(WKA)

!***
!***  Apply adjoint lateral bc on PKF and WKF
!***

         call bocosHnT(mpl,grid,PKF,grid%km2,grid%im,grid%jm,1 ,grid%ib,grid%jb,grid%Fimax,grid%Fjmax,1,1)
         call bocosHnT(mpl,grid,WKF,grid%km3,grid%im,grid%jm,grid%Lm,grid%ib,grid%jb,grid%Fimax,grid%Fjmax,1,1)
!TEST
!        call out_3d2(WKF,-grid%ib,-grid%jb,grid%im+grid%ib,grid%jm+grid%jb,grid%lm,grid%lm/2,'w')
!        call finishMPI
!TEST

!***
!***  Form state%VK2D and state%VK3D
!***
      state%VK2D=0
      state%VK3D=0

      state%VK2D(0:grid%im,0:grid%jm,     1:grid%km2)=PKF(0:grid%im,0:grid%jm,     1:grid%km2)
      state%VK3D(0:grid%im,0:grid%jm,1:grid%lm,1:grid%km3)=WKF(0:grid%im,0:grid%jm,1:grid%lm,1:grid%km3)

    deallocate(PKF)
    deallocate(WKF)


!----------------------------------------------------------------------
                        endsubroutine anal_to_filt2_2d_1_3d

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine filt_to_anal2_2d_1_3d &
!***********************************************************************
!                                                                      !
!  Transfer data from filter to analysis grid                          !
!  (2 2D variables and 1 3D variable)                                  !
!                                                                      !
!***********************************************************************
(mpl,grid,state)
implicit none

type(mpl_type),intent(inout) :: mpl
type(mgbf_grid_type),intent(inout) :: grid
type(mgbf_state_type),intent(inout) :: state

real(kind_real),allocatable,dimension(:,:,:):: PKF
real(kind_real),allocatable,dimension(:,:,:,:):: WKF

real(kind_real),allocatable,dimension(:,:,:):: PKA
real(kind_real),allocatable,dimension(:,:,:,:):: WKA

!----------------------------------------------------------------------
!TEST
      if(grid%km2==1 .and. grid%nm==grid%im .and. grid%mm==grid%jm) then
        state%PA1(0:grid%im,0:grid%jm)= state%VK2D(0:grid%im,0:grid%jm,   1)
        state%PA2(0:grid%im,0:grid%jm)= state%VK2D(0:grid%im,0:grid%jm,   2)
        state%WA1(0:grid%im,0:grid%jm,1:grid%im)= state%VK3D(0:grid%im,0:grid%jm,1:grid%lm,1)
          return
      endif
!TEST
!***
!***  Define PKF and WKF
!***
    allocate(PKF(-grid%ib:grid%im+grid%ib,-grid%jb:grid%jm+grid%jb,   grid%km2))
    allocate(WKF(-grid%ib:grid%im+grid%ib,-grid%jb:grid%jm+grid%jb,grid%lm,grid%km3))

      PKF(0:grid%im,0:grid%jm,     1:grid%km2)= state%VK2D(0:grid%im,0:grid%jm,     1:grid%km2)
      WKF(0:grid%im,0:grid%jm,1:grid%lm,1:grid%km3)= state%VK3D(0:grid%im,0:grid%jm,1:grid%lm,1:grid%km3)


!***
!***  Supply boundary conditions for PKF, WKF
!***
         call bocosHn(mpl,grid,PKF,grid%km2,grid%im,grid%jm,1 ,grid%ib,grid%jb,grid%Fimax,grid%Fjmax,1,1)
         call bocosHn(mpl,grid,WKF,grid%km3,grid%im,grid%jm,grid%Lm,grid%ib,grid%jb,grid%Fimax,grid%Fjmax,1,1)

!***
!*** Interpolate to analysis grid composite variables
!***

    allocate(PKA(0:grid%nm,0:grid%mm,   grid%km2))
    allocate(WKA(0:grid%nm,0:grid%mm,grid%lm,grid%km3))

!                                                 call btim(   intp_tim)
         call lsqr_forward_xyk(grid,state,PKF,PKA)
         call lsqr_forward_xyk(grid,state,WKF,WKA)
!                                                 call etim(   intp_tim)
    deallocate(PKF)
    deallocate(WKF)
!***
!*** Decompose composite variables: state%PA1, ..., PA4, state%WA1, ...,
!state%WA6 <- PKA and WKA
!***

    call decompose(PKA,state%PA1,state%PA2,grid%nm,grid%mm,grid%km2)
!    call
!    decompose(WKA,state%WA1,state%WA2,state%WA3,state%WA4,state%WA5,state%WA6,nm,mm,lm,km)

       state%WA1(0:grid%nm,0:grid%mm,1:grid%lm)=WKA(0:grid%nm,0:grid%jm,1:grid%lm,1)


    deallocate(PKA)
    deallocate(WKA)


! Recover original variables that do not share common edges
!(We assume those are state%PA30 and state%WA20 . Adjust that in GSI! )

!            state%PA30(1:nm,1:mm)=state%PA3(1:nm,1:mm)
!            state%WA20(1:nm,1:mm,1:lm)=state%WA2(1:nm,1:mm,1:lm)


!----------------------------------------------------------------------
                        endsubroutine filt_to_anal2_2d_1_3d
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine anal_to_filt &
!***********************************************************************
!                                                                      !
!  Transfer data from analysis to first generaton of filter grid       !
!                                                                      !
!***********************************************************************
(mpl,grid,state)
implicit none

type(mpl_type),intent(inout) :: mpl
type(mgbf_grid_type),intent(inout) :: grid
type(mgbf_state_type),intent(inout) :: state


real(kind_real),allocatable,dimension(:,:,:):: PKA
real(kind_real),allocatable,dimension(:,:,:,:):: WKA

real(kind_real),allocatable,dimension(:,:,:):: PKF
real(kind_real),allocatable,dimension(:,:,:,:):: WKF

!----------------------------------------------------------------------
!TEST
      if(grid%km2==1 .and. grid%nm==grid%im .and. grid%mm==grid%jm) then
        state%VK2D(0:grid%im,0:grid%jm,   grid%km2)=state%PA1(0:grid%im,0:grid%jm)
        state%VK3D(0:grid%im,0:grid%jm,grid%lm,grid%km3)=0.
          return
      endif
!TEST
!
! Make sure that subdomains of all variables share cogrid%mmon edges
!(We assume those are state%PA30 and state%WA30. Adjust that in GSI! )

            state%PA3(1:grid%nm,1:grid%mm)=state%PA30(1:grid%nm,1:grid%mm)
            state%WA2(1:grid%nm,1:grid%mm,1:grid%lm)=state%WA20(1:grid%nm,1:grid%mm,1:grid%lm)

          call v02v(mpl,grid,state%PA3,grid%nm,grid%mm,1)
          call v02v(mpl,grid,state%WA2,grid%nm,grid%mm,grid%lm)

!***
!*** Create composite variables: state%PA1, ..., PA4, state%WA1, ..., state%WA6
!-> PKA and WKA
!***

    allocate(PKA(0:grid%nm,0:grid%mm,   grid%km2))
    allocate(WKA(0:grid%nm,0:grid%mm,grid%lm,grid%km3))

    call compose(state%PA1,state%PA2,state%PA3,state%PA4,PKA,grid%nm,grid%mm,grid%km2)
    call compose(state%WA1,state%WA2,state%WA3,state%WA4,state%WA5,state%WA6,WKA,grid%nm,grid%mm,grid%lm,grid%km3)


!***
!*** Adjoint interpolate and build composite variables
!***
    allocate(PKF(-grid%ib:grid%im+grid%ib,-grid%jb:grid%jm+grid%jb,   grid%km2))
    allocate(WKF(-grid%ib:grid%im+grid%ib,-grid%jb:grid%jm+grid%jb,grid%lm,grid%km3))


!                                                 call btim(  aintp_tgrid%im)
        call lsqr_adjoint_xyk(grid,state,PKA,PKF)
        call lsqr_adjoint_xyk(grid,state,WKA,WKF)
!                                                 call etim(  aintp_tgrid%im)

    deallocate(PKA)
    deallocate(WKA)

!***
!***  Apply adjoint lateral bc on PKF and WKF
!***

         call bocosHnT(mpl,grid,PKF,grid%km2,grid%im,grid%jm,1 ,grid%ib,grid%jb,grid%Fimax,grid%Fjmax,1,1)
         call bocosHnT(mpl,grid,WKF,grid%km3,grid%im,grid%jm,grid%Lm,grid%ib,grid%jb,grid%Fimax,grid%Fjmax,1,1)

!***
!***  Form state%VK2D and state%VK3D
!***
      state%VK2D=0
      state%VK3D=0

      state%VK2D(0:grid%im,0:grid%jm,   1:grid%km2)=PKF(0:grid%im,0:grid%jm,   1:grid%km2)
      state%VK3D(0:grid%im,0:grid%jm,1:grid%lm,1:grid%km3)=WKF(0:grid%im,0:grid%jm,1:grid%lm,1:grid%km3)

    deallocate(PKF)
    deallocate(WKF)


!----------------------------------------------------------------------
                        endsubroutine anal_to_filt
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine filt_to_anal &
!***********************************************************************
!                                                                      !
!  Transfer data from filter to analysis grid                          !
!                                                                      !
!***********************************************************************
(mpl,grid,state)
implicit none

type(mpl_type),intent(inout) :: mpl
type(mgbf_grid_type),intent(inout) :: grid
type(mgbf_state_type),intent(inout) :: state

real(kind_real),allocatable,dimension(:,:,:):: PKF
real(kind_real),allocatable,dimension(:,:,:,:):: WKF

real(kind_real),allocatable,dimension(:,:,:):: PKA
real(kind_real),allocatable,dimension(:,:,:,:):: WKA

!----------------------------------------------------------------------
!TEST
      if(grid%km2==1 .and. grid%nm==grid%im .and. grid%mm==grid%jm) then
        state%PA1(0:grid%im,0:grid%jm)= state%VK2D(0:grid%im,0:grid%jm,   grid%km2)
!        state%VK3D(0:grid%im,0:grid%jm,grid%lm,grid%km3)=0.
          return
      endif
!TEST
!***
!***  Define PKF and WKF
!***
    allocate(PKF(-grid%ib:grid%im+grid%ib,-grid%jb:grid%jm+grid%jb,   grid%km2))
    allocate(WKF(-grid%ib:grid%im+grid%ib,-grid%jb:grid%jm+grid%jb,grid%lm,grid%km3))

      PKF(0:grid%im,0:grid%jm,   1:grid%km2)= state%VK2D(0:grid%im,0:grid%jm,   1:grid%km2)
      WKF(0:grid%im,0:grid%jm,grid%lm,1:grid%km3)= state%VK3D(0:grid%im,0:grid%jm,grid%lm,1:grid%km3)


!***
!***  Supply boundary conditions for PKF, WKF
!***
         call bocosHn(mpl,grid,PKF,grid%km2,grid%im,grid%jm,1 ,grid%ib,grid%jb,grid%Fimax,grid%Fjmax,1,1)
         call bocosHn(mpl,grid,WKF,grid%km3,grid%im,grid%jm,grid%Lm,grid%ib,grid%jb,grid%Fimax,grid%Fjmax,1,1)

!***
!*** Interpolate to analysis grid composite variables
!***

    allocate(PKA(0:grid%nm,0:grid%mm,   grid%km2))
    allocate(WKA(0:grid%nm,0:grid%mm,grid%lm,grid%km3))

!                                                 call btim(   intp_tgrid%im)
         call lsqr_forward_xyk(grid,state,PKF,PKA)
         call lsqr_forward_xyk(grid,state,WKF,WKA)
!                                                 call etim(   intp_tgrid%im)
    deallocate(PKF)
    deallocate(WKF)
!***
!*** Decompose composite variables: state%PA1, ..., PA4, state%WA1, ...,
!state%WA6 <- PKA and WKA
!***

    call decompose(PKA,state%PA1,state%PA2,state%PA3,state%PA4,grid%nm,grid%mm,grid%km2)
    call decompose(WKA,state%WA1,state%WA2,state%WA3,state%WA4,state%WA5,state%WA6,grid%nm,grid%mm,grid%lm,grid%km)


    deallocate(PKA)
    deallocate(WKA)


! Recover original variables that do not share cogrid%mmon edges
!(We assume those are state%PA30 and state%WA20 . Adjust that in GSI! )

            state%PA30(1:grid%nm,1:grid%mm)=state%PA3(1:grid%nm,1:grid%mm)
            state%WA20(1:grid%nm,1:grid%mm,1:grid%lm)=state%WA2(1:grid%nm,1:grid%mm,1:grid%lm)


!----------------------------------------------------------------------
                        endsubroutine filt_to_anal
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&



!***
!***  Follows set of compose and decompose subroutines
!***

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine compose2d_2                          &
!***********************************************************************
!                                                                      !
!   Create a composite 2d array from 2 2d arrays                       !
!                                                                      !
!***********************************************************************
(V_1,V_2,V,imax,jmax,kmax)
!-----------------------------------------------------------------------
implicit none

integer(kind_int),intent(in):: imax,jmax,kmax
real(kind_real),dimension(0:imax,0:jmax),intent(in):: V_1,V_2
real(kind_real),dimension(0:imax,0:jmax,kmax),intent(inout):: V
!-----------------------------------------------------------------------

   V(:,:,1)=V_1(:,:)
   V(:,:,2)=V_2(:,:)

   if(kmax /= 2) then
     stop 'number of arrays 2 does not match'
   endif

!-----------------------------------------------------------------------
                        endsubroutine compose2d_2

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine compose2d_3                          &
!***********************************************************************
!                                                                      !
!   Create a composite 2d array from 3 2d arrays                       !
!                                                                      !
!***********************************************************************
(V_1,V_2,V_3,V,imax,jmax,kmax)
!-----------------------------------------------------------------------
implicit none

integer(kind_int),intent(in):: imax,jmax,kmax
real(kind_real),dimension(0:imax,0:jmax),intent(in):: V_1,V_2,V_3
real(kind_real),dimension(0:imax,0:jmax,kmax),intent(inout):: V
!-----------------------------------------------------------------------

   V(:,:,1)=V_1(:,:)
   V(:,:,2)=V_2(:,:)
   V(:,:,3)=V_3(:,:)

   if(kmax /= 3) then
     stop 'number of arrays 3 does not match'
   endif
!-----------------------------------------------------------------------
                        endsubroutine compose2d_3

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine compose2d_4                          &
!***********************************************************************
!                                                                      !
!   Create a composite 2d array from 4 2d arrays                       !
!                                                                      !
!***********************************************************************
(V_1,V_2,V_3,V_4,V,imax,jmax,kmax)
!-----------------------------------------------------------------------
implicit none

integer(kind_int),intent(in):: imax,jmax,kmax
real(kind_real),dimension(0:imax,0:jmax),intent(in):: V_1,V_2,V_3,V_4
real(kind_real),dimension(0:imax,0:jmax,kmax),intent(inout):: V
!-----------------------------------------------------------------------

   V(:,:,1)=V_1(:,:)
   V(:,:,2)=V_2(:,:)
   V(:,:,3)=V_3(:,:)
   V(:,:,4)=V_4(:,:)

   if(kmax /= 4) then
     stop 'number of arrays 4 does not match'
   endif
!-----------------------------------------------------------------------
                        endsubroutine compose2d_4

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine compose3d_2                          &
!***********************************************************************
!                                                                      !
!   Create a composite 3d array from 2 3d arrays                       !
!                                                                      !
!***********************************************************************
(V_1,V_2,V,imax,jmax,lmax,kmax)
!-----------------------------------------------------------------------
implicit none

integer(kind_int),intent(in):: imax,jmax,lmax,kmax
real(kind_real),dimension(0:imax,0:jmax,lmax),intent(in):: V_1,V_2
real(kind_real),dimension(0:imax,0:jmax,lmax,kmax),intent(inout):: V
!-----------------------------------------------------------------------

   V(:,:,:,1)=V_1(:,:,:)
   V(:,:,:,2)=V_2(:,:,:)

   if(kmax /= 2) then
     stop 'number of 2 arrays does not match'
   endif

!-----------------------------------------------------------------------
                        endsubroutine compose3d_2

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine compose3d_3                          &
!***********************************************************************
!                                                                      !
!   Create a composite 3d array from 3 3d arrays                       !
!                                                                      !
!***********************************************************************
(V_1,V_2,V_3,V,imax,jmax,lmax,kmax)
!-----------------------------------------------------------------------
implicit none

integer(kind_int),intent(in):: imax,jmax,lmax,kmax
real(kind_real),dimension(0:imax,0:jmax,lmax),intent(in):: V_1,V_2,V_3
real(kind_real),dimension(0:imax,0:jmax,lmax,kmax),intent(inout):: V
!-----------------------------------------------------------------------

   V(:,:,:,1)=V_1(:,:,:)
   V(:,:,:,2)=V_2(:,:,:)
   V(:,:,:,3)=V_3(:,:,:)

   if(kmax /= 3) then
     stop 'number of arrays 3 does not match'
   endif
!-----------------------------------------------------------------------
                        endsubroutine compose3d_3

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine compose3d_4                          &
!***********************************************************************
!                                                                      !
!   Create a composite 3d array from 4 3d arrays                       !
!                                                                      !
!***********************************************************************
(V_1,V_2,V_3,V_4,V,imax,jmax,lmax,kmax)
!-----------------------------------------------------------------------
implicit none

integer(kind_int),intent(in):: imax,jmax,lmax,kmax
real(kind_real),dimension(0:imax,0:jmax,lmax),intent(in):: V_1,V_2,V_3,V_4
real(kind_real),dimension(0:imax,0:jmax,lmax,kmax),intent(inout):: V
!-----------------------------------------------------------------------

   V(:,:,:,1)=V_1(:,:,:)
   V(:,:,:,2)=V_2(:,:,:)
   V(:,:,:,3)=V_3(:,:,:)
   V(:,:,:,4)=V_4(:,:,:)

   if(kmax /= 4) then
     stop 'number of arrays 4 does not match'
   endif
!-----------------------------------------------------------------------
                        endsubroutine compose3d_4

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine compose3d_5                          &
!***********************************************************************
!                                                                      !
!   Create a composite 3d array from 5 3d arrays                       !
!                                                                      !
!***********************************************************************
(V_1,V_2,V_3,V_4,V_5,V,imax,jmax,lmax,kmax)
!-----------------------------------------------------------------------
implicit none

integer(kind_int),intent(in):: imax,jmax,lmax,kmax
real(kind_real),dimension(0:imax,0:jmax,lmax),intent(in):: V_1,V_2,V_3,V_4 &
                                                       ,V_5
real(kind_real),dimension(0:imax,0:jmax,lmax,kmax),intent(inout):: V
!-----------------------------------------------------------------------

   V(:,:,:,1)=V_1(:,:,:)
   V(:,:,:,2)=V_2(:,:,:)
   V(:,:,:,3)=V_3(:,:,:)
   V(:,:,:,4)=V_4(:,:,:)
   V(:,:,:,5)=V_5(:,:,:)

   if(kmax /= 5) then
     stop 'number of arrays 5 does not match'
   endif
!-----------------------------------------------------------------------
                        endsubroutine compose3d_5

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine compose3d_6                          &
!***********************************************************************
!                                                                      !
!   Create a composite 3d array from 6 3d arrays                       !
!                                                                      !
!***********************************************************************
(V_1,V_2,V_3,V_4,V_5,V_6,V,imax,jmax,lmax,kmax)
!-----------------------------------------------------------------------
implicit none

integer(kind_int),intent(in):: imax,jmax,lmax,kmax
real(kind_real),dimension(0:imax,0:jmax,lmax),intent(in):: V_1,V_2,V_3,V_4 &
                                                       ,V_5,V_6
real(kind_real),dimension(0:imax,0:jmax,lmax,kmax),intent(inout):: V
!-----------------------------------------------------------------------

   V(:,:,:,1)=V_1(:,:,:)
   V(:,:,:,2)=V_2(:,:,:)
   V(:,:,:,3)=V_3(:,:,:)
   V(:,:,:,4)=V_4(:,:,:)
   V(:,:,:,5)=V_5(:,:,:)
   V(:,:,:,6)=V_6(:,:,:)

   if(kmax /= 6) then
     stop 'number of arrays 6 does not match'
   endif
!-----------------------------------------------------------------------
                        endsubroutine compose3d_6

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine decompose2d_2                        &
!***********************************************************************
!                                                                      !
!   Restore 2 2d arrays from a composite 2d array                      !
!                                                                      !
!***********************************************************************
(V,V_1,V_2,imax,jmax,kmax)
!-----------------------------------------------------------------------
implicit none

integer(kind_int),intent(in):: imax,jmax,kmax
real(kind_real),dimension(0:imax,0:jmax,kmax),intent(in):: V
real(kind_real),dimension(0:imax,0:jmax),intent(out):: V_1,V_2
!-----------------------------------------------------------------------

   V_1(:,:)=V(:,:,1)
   V_2(:,:)=V(:,:,2)

   if(kmax /= 2) then
     stop 'number of arrays does not match'
   endif

!-----------------------------------------------------------------------
                        endsubroutine decompose2d_2

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine decompose2d_3                        &
!***********************************************************************
!                                                                      !
!   Restode 3 2d arrays from a composite 2d array                      !
!                                                                      !
!***********************************************************************
(V,V_1,V_2,V_3,imax,jmax,kmax)
!-----------------------------------------------------------------------
implicit none

integer(kind_int),intent(in):: imax,jmax,kmax
real(kind_real),dimension(0:imax,0:jmax,kmax),intent(in):: V
real(kind_real),dimension(0:imax,0:jmax),intent(out):: V_1,V_2,V_3
!-----------------------------------------------------------------------

   V_1(:,:)=V(:,:,1)
   V_2(:,:)=V(:,:,2)
   V_3(:,:)=V(:,:,3)

   if(kmax /= 3) then
     stop 'number of arrays does not match'
   endif

!-----------------------------------------------------------------------
                        endsubroutine decompose2d_3

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine decompose2d_4                        &
!***********************************************************************
!                                                                      !
!   Restore 4 2d arrays from a composite 2d array                      !
!                                                                      !
!***********************************************************************
(V,V_1,V_2,V_3,V_4,imax,jmax,kmax)
!-----------------------------------------------------------------------
implicit none

integer(kind_int),intent(in):: imax,jmax,kmax
real(kind_real),dimension(0:imax,0:jmax,kmax),intent(in):: V
real(kind_real),dimension(0:imax,0:jmax),intent(out):: V_1,V_2,V_3,V_4
!-----------------------------------------------------------------------

   V_1(:,:)=V(:,:,1)
   V_2(:,:)=V(:,:,2)
   V_3(:,:)=V(:,:,3)
   V_4(:,:)=V(:,:,4)

   if(kmax /= 4) then
     stop 'number of arrays does not match'
   endif
!-----------------------------------------------------------------------
                        endsubroutine decompose2d_4

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine decompose3d_2                        &
!***********************************************************************
!                                                                      !
!   Restore 2 3d arrays from a composite 3d array                      !
!                                                                      !
!***********************************************************************
(V,V_1,V_2,imax,jmax,lmax,kmax)
!-----------------------------------------------------------------------
implicit none

integer(kind_int),intent(in):: imax,jmax,lmax,kmax
real(kind_real),dimension(0:imax,0:jmax,lmax,kmax),intent(in):: V
real(kind_real),dimension(0:imax,0:jmax,lmax),intent(out):: V_1,V_2
!-----------------------------------------------------------------------

   V_1(:,:,:)=V(:,:,:,1)
   V_2(:,:,:)=V(:,:,:,2)

   if(kmax /= 2) then
     stop 'number of arrays does not match'
   endif

!-----------------------------------------------------------------------
                        endsubroutine decompose3d_2

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine decompose3d_3                        &
!***********************************************************************
!                                                                      !
!   Restore 3 3d arrays from a composite 3d array                      !
!                                                                      !
!***********************************************************************
(V,V_1,V_2,V_3,imax,jmax,lmax,kmax)
!-----------------------------------------------------------------------
implicit none

integer(kind_int),intent(in):: imax,jmax,lmax,kmax
real(kind_real),dimension(0:imax,0:jmax,lmax,kmax),intent(in):: V
real(kind_real),dimension(0:imax,0:jmax,lmax),intent(out):: V_1,V_2,V_3
!-----------------------------------------------------------------------

   V_1(:,:,:)=V(:,:,:,1)
   V_2(:,:,:)=V(:,:,:,2)
   V_3(:,:,:)=V(:,:,:,3)

   if(kmax /= 3) then
     stop 'number of arrays does not match'
   endif
!-----------------------------------------------------------------------
                        endsubroutine decompose3d_3

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine decompose3d_4                        &
!***********************************************************************
!                                                                      !
!   Restore 4 3d arrays from a composite 3d array                      !
!                                                                      !
!***********************************************************************
(V,V_1,V_2,V_3,V_4,imax,jmax,lmax,kmax)
!-----------------------------------------------------------------------
implicit none

integer(kind_int),intent(in):: imax,jmax,lmax,kmax
real(kind_real),dimension(0:imax,0:jmax,lmax,kmax),intent(in):: V
real(kind_real),dimension(0:imax,0:jmax,lmax),intent(out):: V_1,V_2,V_3,V_4
!-----------------------------------------------------------------------

   V_1(:,:,:)=V(:,:,:,1)
   V_2(:,:,:)=V(:,:,:,2)
   V_3(:,:,:)=V(:,:,:,3)
   V_4(:,:,:)=V(:,:,:,4)

   if(kmax /= 4) then
     stop 'number of arrays does not match'
   endif
!-----------------------------------------------------------------------
                        endsubroutine decompose3d_4

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine decompose3d_5                        &
!***********************************************************************
!                                                                      !
!   Restore 5 3d arrays from a composite 3d array                      !
!                                                                      !
!***********************************************************************
(V,V_1,V_2,V_3,V_4,V_5,imax,jmax,lmax,kmax)
!-----------------------------------------------------------------------
implicit none

integer(kind_int),intent(in):: imax,jmax,lmax,kmax
real(kind_real),dimension(0:imax,0:jmax,lmax,kmax),intent(in):: V
real(kind_real),dimension(0:imax,0:jmax,lmax),intent(out):: V_1,V_2,V_3,V_4 &
                                                        ,V_5
!-----------------------------------------------------------------------

   V_1(:,:,:)=V(:,:,:,1)
   V_2(:,:,:)=V(:,:,:,2)
   V_3(:,:,:)=V(:,:,:,3)
   V_4(:,:,:)=V(:,:,:,4)
   V_5(:,:,:)=V(:,:,:,5)

   if(kmax /= 5) then
     stop 'number of arrays does not match'
   endif
!-----------------------------------------------------------------------
                        endsubroutine decompose3d_5

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine decompose3d_6                        &
!***********************************************************************
!                                                                      !
!   Restore 6 3d arrays from a composite 3d array                      !
!                                                                      !
!***********************************************************************
(V,V_1,V_2,V_3,V_4,V_5,V_6,imax,jmax,lmax,kmax)
!-----------------------------------------------------------------------
implicit none

integer(kind_int),intent(in):: imax,jmax,lmax,kmax
real(kind_real),dimension(0:imax,0:jmax,lmax,kmax),intent(in):: V
real(kind_real),dimension(0:imax,0:jmax,lmax),intent(out):: V_1,V_2,V_3,V_4 &
                                                         ,V_5,V_6
!-----------------------------------------------------------------------

   V_1(:,:,:)=V(:,:,:,1)
   V_2(:,:,:)=V(:,:,:,2)
   V_3(:,:,:)=V(:,:,:,3)
   V_4(:,:,:)=V(:,:,:,4)
   V_5(:,:,:)=V(:,:,:,5)
   V_6(:,:,:)=V(:,:,:,6)

   if(kmax /= 6) then
     stop 'number of arrays does not match'
   endif
!-----------------------------------------------------------------------
                        endsubroutine decompose3d_6

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        endmodule tools_transfer
