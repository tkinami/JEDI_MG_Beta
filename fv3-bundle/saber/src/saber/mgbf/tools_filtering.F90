!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        module tools_filtering
!***********************************************************************
!                                                                      !
! Contains all multigrid filtering prodecures                          ! 
!                                                                      ! 
!                                                     M. Rancic (2020) !
!***********************************************************************
use tools_kinds, only: kind_real,kind_int
use tools_bocos, only: bocosHn,bocosHnT
use tools_generations, only: weighting1,upsending1,downsending1,differencing1
use tools_betafilter, only: vrbeta1T,vrbeta1
use tools_betafilter, only: vrbeta2T,vrbeta2
use tools_betafilter, only: vrbeta3T,vrbeta3
use type_mgbf_grid
use type_mgbf_state
use type_mpl

public mg_filtering_procedure 

private mg_filtering_procedure1
private mg_filtering_procedure2
private mg_filtering_procedure3
private mg_filtering_procedure5

private sup_vrbeta1
private sup_vrbeta1T
private sup_vrbeta3
private sup_vrbeta3T

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine mg_filtering_procedure(mpl,grid,state)                       
!***********************************************************************
!                                                                      !
! Driver for Multigrid filtering procedures:                           !
!                                                                      !
!   1, 2, 3: Original filtering with two upsendings and downsendings   !
!         1: 2d radial filter for all variables                        !
!              (comparable with recursive filter)                      !
!    -->  2: 2d radial filter with 1d in vertical for 3d variables     !
!         3: 3d radial filter for 3d variables                         !
!                                                                      !
!   4, 5, 6: Filtering with Helmholtz operator and one up/downsending  !
!         4: 2d radial filter for all variables                        !
!         5: 2d radial filter with 1d in vertical for 3d variables     !
!         6: 3d radial filter for 3d variables                         !
!                                                                      !
!  (In future, line filter will be used for all 2d and 3d filtering)   !
!                                                                      !
!***********************************************************************
IMPLICIT NONE 

type(mpl_type),intent(inout) :: mpl
type(mgbf_grid_type),intent(inout) :: grid
type(mgbf_state_type),intent(inout) :: state
integer(kind_int) mg_filt
!-----------------------------------------------------------------------
!  mg_filt=5
  mg_filt=1

      select case(mg_filt)
        case(1)
          call mg_filtering_procedure1(mpl,grid,state)
        case(2)
          call mg_filtering_procedure2(mpl,grid,state)
        case(3)
          call mg_filtering_procedure3(mpl,grid,state)
        case(5)
          call mg_filtering_procedure5(mpl,grid,state)
       end select

    


!-----------------------------------------------------------------------
                        endsubroutine mg_filtering_procedure    

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine mg_filtering_procedure1(mpl,grid,state)                        
!***********************************************************************
!                                                                      !
! Multigrid filtering procedure 1:                                     !
!                                                                      !
!     - Multiple of 2D and 3D variables                                !
!     - 2 upsendings and downsendings                                  !
!     - 2d radial filter                                               !
!                                                                      !
!***********************************************************************
!-----------------------------------------------------------------------
!use mg_intstate1, only: state%VK2D,state%VK3D,state%HK2D,state%HK3D,state%pasp1,state%pasp2,state%ss1,state%ss2
IMPLICIT NONE

type(mpl_type),intent(inout) :: mpl
type(mgbf_grid_type),intent(inout) :: grid
type(mgbf_state_type),intent(inout) :: state
integer(kind_int) L

!==================== Adjoint (Conservative step) ======================

!***
!*** Adjoint interpolate and upsend (Step 1)
!***
!                                                 call btim( upsend_tim)
       call upsending1(mpl,grid,state%VK2D,state%HK2D)
       call upsending1(mpl,grid,state%VK3D,state%HK3D)
!                                                 call etim( upsend_tim)

!***
!*** Apply weights at all generations
!***

!                                                 call btim( weight_tim)
        call weighting1(grid,state,state%VK2D,state%HK2D)
        call weighting1(grid,state,state%VK3D,state%HK3D)
!                                                 call etim( weight_tim)
!***
!*** Apply adjoint of Beta filter at all generations (Step 3)
!***

!                                                 call btim(    bfiltT_tim)


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!
! Adjoint filtering
!
      call vrbeta2T(grid%km2,grid%hx,0,grid%im,grid%hy,0,grid%jm,grid%p,state%pasp2,state%ss2,state%VK2D(:,:,:))
    do L=1,grid%lm
      call vrbeta2T(grid%km3,grid%hx,0,grid%im,grid%hy,0,grid%jm,grid%p,state%pasp2,state%ss2,state%VK3D(:,:,l,:))
    end do

  if(grid%l_hgen)  then
      call vrbeta2T(grid%km2,grid%hx,0,grid%im,grid%hy,0,grid%jm,grid%p,state%pasp2,state%ss2,state%HK2D(:,:,:))
    do L=1,grid%lm
      call vrbeta2T(grid%km3,grid%hx,0,grid%im,grid%hy,0,grid%jm,grid%p,state%pasp2,state%ss2,state%HK3D(:,:,l,:))
    end do
   endif 

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

      call bocosHnT(mpl,grid,state%VK2D,grid%km2,grid%im,grid%jm,1 ,grid%hx,grid%hy,grid%Fimax,grid%Fjmax,1,1)
      call bocosHnT(mpl,grid,state%VK3D,grid%km3,grid%im,grid%jm,grid%lm,grid%hx,grid%hy,grid%Fimax,grid%Fjmax,1,1)

      call bocosHnT(mpl,grid,state%HK2D,grid%km2,grid%im,grid%jm,1 ,grid%hx,grid%hy,grid%Fimax,grid%Fjmax,2,grid%gm)
      call bocosHnT(mpl,grid,state%HK3D,grid%km3,grid%im,grid%jm,grid%lm,grid%hx,grid%hy,grid%Fimax,grid%Fjmax,2,grid%gm)

!                                                 call etim(    bfiltT_tim)

!***
!***  Downsend, interpolate and add (Step 4)
!***  Then set to zero high generations (Step 5)
!***


!TEST
!       write(20,*)'Before downsending in fitering1'
!TEST
!                                                 call btim(   dnsend_tim)
       call downsending1(mpl,grid,state%HK2D,state%VK2D)
       call downsending1(mpl,grid,state%HK3D,state%VK3D)
!                                                 call etim(   dnsend_tim)

!ADJOINT_TEST
      
!!      call multiply_add(b0,a,res_glob)
!!    
!!     if(grid%mype==0) then
!!           print *,'sum1=',res_glob  
!!      endif
   
!ADJOINT_TEST

!==================== Forward (Smoothing step) =========================

!***
!*** Adjoint interpolate and upsend (Step 6 - same as Step 1)
!***
!                                                 call btim( upsend_tim)
         call upsending1(mpl,grid,state%VK2D,state%HK2D)
         call upsending1(mpl,grid,state%VK3D,state%HK3D)
!                                                 call etim( upsend_tim)

!***
!*** Apply Beta filter at all generations (Step 7)
!***

!                                                 call btim(    bfilt_tim)
!TEST
!       write(20,*)'Before bocosHn in fitering1'
!TEST

      call bocosHn(mpl,grid,state%VK2D,grid%km2,grid%im,grid%jm,1 ,grid%hx,grid%hy,grid%Fimax,grid%Fjmax,1,1)
      call bocosHn(mpl,grid,state%VK3D,grid%km3,grid%im,grid%jm,grid%lm,grid%hx,grid%hy,grid%Fimax,grid%Fjmax,1,1)

      call bocosHn(mpl,grid,state%HK2D,grid%km2,grid%im,grid%jm,1 ,grid%hx,grid%hy,grid%Fimax,grid%Fjmax,2,grid%gm)
      call bocosHn(mpl,grid,state%HK3D,grid%km3,grid%im,grid%jm,grid%lm,grid%hx,grid%hy,grid%Fimax,grid%Fjmax,2,grid%gm)


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!
! Filtering
!

      call vrbeta2(grid%km2,grid%hx,0,grid%im,grid%hy,0,grid%jm,grid%p,state%pasp2,state%ss2,state%VK2D(:,:,:))
   do L=1,grid%lm
      call vrbeta2(grid%km3,grid%hx,0,grid%im,grid%hy,0,grid%jm,grid%p,state%pasp2,state%ss2,state%VK3D(:,:,l,:))
   enddo

  if(grid%l_hgen)  then
      call vrbeta2(grid%km2,grid%hx,0,grid%im,grid%hy,0,grid%jm,grid%p,state%pasp2,state%ss2,state%HK2D(:,:,:))
   do L=1,grid%lm
      call vrbeta2(grid%km3,grid%hx,0,grid%im,grid%hy,0,grid%jm,grid%p,state%pasp2,state%ss2,state%HK3D(:,:,l,:))
   enddo
  endif

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff


!                                                 call etim(    bfilt_tim)
!***
!*** Apply weights at all generations again
!***

!                                                 call btim( weight_tim)
        call weighting1(grid,state,state%VK2D,state%HK2D)
        call weighting1(grid,state,state%VK3D,state%HK3D)
!                                                 call etim( weight_tim)

!***
!***  Repeat:
!***  Downsend, interpolate and add (Step 4)
!***  Then zero high generations (Step 5)
!***

!                                                 call btim(   dnsend_tim)
       call downsending1(mpl,grid,state%HK2D,state%VK2D)
       call downsending1(mpl,grid,state%HK3D,state%VK3D)
!                                                 call etim(   dnsend_tim)



!-----------------------------------------------------------------------
                        endsubroutine mg_filtering_procedure1    

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine mg_filtering_procedure2(mpl,grid,state)                        
!***********************************************************************
!                                                                      !
! Multigrid filtering procedure 2:                                     !
!                                                                      !
!     - Multiple of 2D and 3D variables                                !
!     - 2 upsendings and downsendings                                  !
!     - 2d radial filter + 1d vertical filter                          !
!                                                                      !
!***********************************************************************
!-----------------------------------------------------------------------
!use mg_intstate1, only: state%VK2D,state%VK3D,state%HK2D,state%HK3D,state%pasp1,state%pasp2,state%ss1,state%ss2
IMPLICIT NONE

type(mpl_type),intent(inout) :: mpl
type(mgbf_grid_type),intent(inout) :: grid
type(mgbf_state_type),intent(inout) :: state
integer(kind_int) L

!==================== Adjoint (Conservative step) ======================

!***
!*** Adjoint interpolate and upsend (Step 1)
!***
     
!                                                 call btim( upsend_tim)
       call upsending1(mpl,grid,state%VK2D,state%HK2D)
       call upsending1(mpl,grid,state%VK3D,state%HK3D)
!                                                 call etim( upsend_tim)


!***
!*** Apply weights at all generations
!***

!                                                 call btim( weight_tim)
        call weighting1(grid,state,state%VK2D,state%HK2D)
        call weighting1(grid,state,state%VK3D,state%HK3D)

!                                                 call etim( weight_tim)

!***
!*** Apply adjoint of Beta filter at all generations (Step 3)
!***

!                                                 call btim(    bfiltT_tim)


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!
! Adjoint filtering
!
      call vrbeta2T(grid%km2,grid%hx,0,grid%im,grid%hy,0,grid%jm,grid%p,state%pasp2,state%ss2,state%VK2D(:,:,:))
    do L=1,grid%lm
      call vrbeta2T(grid%km3,grid%hx,0,grid%im,grid%hy,0,grid%jm,grid%p,state%pasp2,state%ss2,state%VK3D(:,:,l,:))
    end do

  if(grid%l_hgen)  then
      call vrbeta2T(grid%km2,grid%hx,0,grid%im,grid%hy,0,grid%jm,grid%p,state%pasp2,state%ss2,state%HK2D(:,:,:))
    do L=1,grid%lm
      call vrbeta2T(grid%km3,grid%hx,0,grid%im,grid%hy,0,grid%jm,grid%p,state%pasp2,state%ss2,state%HK3D(:,:,l,:))
    end do
   endif

      call sup_vrbeta1T(state%VK3D,grid%km3,grid%hx,grid%hy,grid%hz,grid%im,grid%jm,grid%lm,grid%p,state%pasp1,state%ss1)
  if(grid%l_hgen)  then
      call sup_vrbeta1T(state%HK3D,grid%km3,grid%hx,grid%hy,grid%hz,grid%im,grid%jm,grid%lm,grid%p,state%pasp1,state%ss1)
   endif

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff


      call bocosHnT(mpl,grid,state%VK2D,grid%km2,grid%im,grid%jm,1 ,grid%hx,grid%hy,grid%Fimax,grid%Fjmax,1,1)
      call bocosHnT(mpl,grid,state%VK3D,grid%km3,grid%im,grid%jm,grid%lm,grid%hx,grid%hy,grid%Fimax,grid%Fjmax,1,1)

      call bocosHnT(mpl,grid,state%HK2D,grid%km2,grid%im,grid%jm,1 ,grid%hx,grid%hy,grid%Fimax,grid%Fjmax,2,grid%gm)
      call bocosHnT(mpl,grid,state%HK3D,grid%km3,grid%im,grid%jm,grid%lm,grid%hx,grid%hy,grid%Fimax,grid%Fjmax,2,grid%gm)

!                                                 call etim(    bfiltT_tim)

!***
!***  Downsend, interpolate and add (Step 4)
!***  Then set to zero high generations (Step 5)
!***


 !                                                call btim(   dnsend_tim)
       call downsending1(mpl,grid,state%HK2D,state%VK2D)
       call downsending1(mpl,grid,state%HK3D,state%VK3D)
!                                                 call etim(   dnsend_tim)

!ADJOINT_TEST
      
!!      call multiply_add(b0,a,res_glob)
!!    
!!     if(grid%mype==0) then
!!           print *,'sum1=',res_glob  
!!      endif
   
!ADJOINT_TEST
!==================== Forward (Smoothing step) =========================

!***
!*** Adjoint interpolate and upsend (Step 6 - same as Step 1)
!***
!                                                 call btim( upsend_tim)
         call upsending1(mpl,grid,state%VK2D,state%HK2D)
         call upsending1(mpl,grid,state%VK3D,state%HK3D)
!                                                 call etim( upsend_tim)

!***
!*** Apply Beta filter at all generations (Step 7)
!***

!                                                 call btim(    bfilt_tim)

      call bocosHn(mpl,grid,state%VK2D,grid%km2,grid%im,grid%jm,1 ,grid%hx,grid%hy,grid%Fimax,grid%Fjmax,1,1)
      call bocosHn(mpl,grid,state%VK3D,grid%km3,grid%im,grid%jm,grid%lm,grid%hx,grid%hy,grid%Fimax,grid%Fjmax,1,1)

      call bocosHn(mpl,grid,state%HK2D,grid%km2,grid%im,grid%jm,1 ,grid%hx,grid%hy,grid%Fimax,grid%Fjmax,2,grid%gm)
      call bocosHn(mpl,grid,state%HK3D,grid%km3,grid%im,grid%jm,grid%lm,grid%hx,grid%hy,grid%Fimax,grid%Fjmax,2,grid%gm)

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!
! Filtering
!

      call vrbeta2(grid%km2,grid%hx,0,grid%im,grid%hy,0,grid%jm,grid%p,state%pasp2,state%ss2,state%VK2D(:,:,:))
   do L=1,grid%lm
      call vrbeta2(grid%km3,grid%hx,0,grid%im,grid%hy,0,grid%jm,grid%p,state%pasp2,state%ss2,state%VK3D(:,:,l,:))
   enddo

  if(grid%l_hgen)  then
      call vrbeta2(grid%km2,grid%hx,0,grid%im,grid%hy,0,grid%jm,grid%p,state%pasp2,state%ss2,state%HK2D(:,:,:))
   do L=1,grid%lm
      call vrbeta2(grid%km3,grid%hx,0,grid%im,grid%hy,0,grid%jm,grid%p,state%pasp2,state%ss2,state%HK3D(:,:,l,:))
   enddo
  endif

      call sup_vrbeta1(state%VK3D,grid%km3,grid%hx,grid%hy,grid%hz,grid%im,grid%jm,grid%lm,grid%p,state%pasp1,state%ss1)
  if(grid%l_hgen)  then
      call sup_vrbeta1(state%HK3D,grid%km3,grid%hx,grid%hy,grid%hz,grid%im,grid%jm,grid%lm,grid%p,state%pasp1,state%ss1)
   endif 

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!                                                 call etim(    bfilt_tim)
!***
!*** Apply weights at all generations again
!***

!                                                 call btim( weight_tim)
      call weighting1(grid,state,state%VK2D,state%HK2D)
      call weighting1(grid,state,state%VK3D,state%HK3D)
!                                                 call etim( weight_tim)

!***
!***  Repeat:
!***  Downsend, interpolate and add (Step 4)
!***  Then zero high generations (Step 5)
!***

!                                                 call btim(   dnsend_tim)
       call downsending1(mpl,grid,state%HK2D,state%VK2D)
       call downsending1(mpl,grid,state%HK3D,state%VK3D)
!                                                 call etim(   dnsend_tim)



!-----------------------------------------------------------------------
                        endsubroutine mg_filtering_procedure2    

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine mg_filtering_procedure3(mpl,grid,state)                        
!***********************************************************************
!                                                                      !
! Multigrid filtering procedure 3:                                     !
!                                                                      !
!     - Multiple of 3D variables                                       !
!     - 2 upsendings and downsendings                                  !
!     - 3d radial filter                                               !
!                                                                      !
!***********************************************************************
!-----------------------------------------------------------------------
!use mg_intstate1, only: state%VK2D,state%VK3D,state%HK2D,state%HK3D,state%pasp1,state%pasp2,state%pasp3,state%ss1,state%ss2,state%ss3
IMPLICIT NONE

type(mpl_type),intent(inout) :: mpl
type(mgbf_grid_type),intent(inout) :: grid
type(mgbf_state_type),intent(inout) :: state
integer(kind_int) L

!==================== Adjoint (Conservative step) ======================

!***
!*** Adjoint interpolate and upsend (Step 1)
!***
!                                                 call btim( upsend_tim)
        call upsending1(mpl,grid,state%VK2D,state%HK2D)
        call upsending1(mpl,grid,state%VK3D,state%HK3D)
!                                                 call etim( upsend_tim)

!***
!*** Apply weights at all generations
!***

!                                                 call btim( weight_tim)
        call weighting1(grid,state,state%VK2D,state%HK2D)
        call weighting1(grid,state,state%VK3D,state%HK3D)
!                                                 call etim( weight_tim)

!***
!*** Apply adjoint of Beta filter at all generations (Step 3)
!***

!                                                 call btim(    bfiltT_tim)


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!
! Adjoint filtering
!
      call vrbeta2T(grid%km2,grid%hx,0,grid%im,grid%hy,0,grid%jm,grid%p,state%pasp2,state%ss2,state%VK2D(:,:,:))
      call sup_vrbeta3T(state%VK3D,grid%km3,grid%hx,grid%hy,grid%hz,grid%im,grid%jm,grid%lm,grid%p,state%pasp3,state%ss3)

 if(grid%l_hgen) then
       call vrbeta2T(grid%km2,grid%hx,0,grid%im,grid%hy,0,grid%jm,grid%p,state%pasp2,state%ss2,state%HK2D(:,:,:))
       call sup_vrbeta3T(state%HK3D,grid%km3,grid%hx,grid%hy,grid%hz,grid%im,grid%jm,grid%lm,grid%p,state%pasp3,state%ss3)
  endif
!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

      call bocosHnT(mpl,grid,state%VK2D,grid%km2,grid%im,grid%jm,1 ,grid%hx,grid%hy,grid%Fimax,grid%Fjmax,1,1)
      call bocosHnT(mpl,grid,state%VK3D,grid%km3,grid%im,grid%jm,grid%lm,grid%hx,grid%hy,grid%Fimax,grid%Fjmax,1,1)

      call bocosHnT(mpl,grid,state%HK2D,grid%km2,grid%im,grid%jm,1 ,grid%hx,grid%hy,grid%Fimax,grid%Fjmax,2,grid%gm)
      call bocosHnT(mpl,grid,state%HK3D,grid%km3,grid%im,grid%jm,grid%lm,grid%hx,grid%hy,grid%Fimax,grid%Fjmax,2,grid%gm)

!                                                 call etim(    bfiltT_tim)

!***
!***  Downsend, interpolate and add (Step 4)
!***  Then set to zero high generations (Step 5)
!***

!                                                 call btim(   dnsend_tim)
       call downsending1(mpl,grid,state%HK2D,state%VK2D)
       call downsending1(mpl,grid,state%HK3D,state%VK3D)

!                                                 call etim(   dnsend_tim)

!ADJOINT_TEST
      
!!      call multiply_add(b0,a,res_glob)
!!    
!!     if(grid%mype==0) then
!!           print *,'sum1=',res_glob  
!!      endif
   
!ADJOINT_TEST
!==================== Forward (Smoothing step) =========================

!***
!*** Adjoint interpolate and upsend (Step 6 - same as Step 1)
!***
!                                                 call btim( upsend_tim)
         call upsending1(mpl,grid,state%VK2D,state%HK2D)
         call upsending1(mpl,grid,state%VK3D,state%HK3D)

!                                                 call etim( upsend_tim)

!***
!*** Apply Beta filter at all generations (Step 7)
!***

!                                                 call btim(    bfilt_tim)

      call bocosHn(mpl,grid,state%VK2D,grid%km2,grid%im,grid%jm,1 ,grid%hx,grid%hy,grid%Fimax,grid%Fjmax,1,1)
      call bocosHn(mpl,grid,state%VK3D,grid%km3,grid%im,grid%jm,grid%lm,grid%hx,grid%hy,grid%Fimax,grid%Fjmax,1,1)

      call bocosHn(mpl,grid,state%HK2D,grid%km2,grid%im,grid%jm,1 ,grid%hx,grid%hy,grid%Fimax,grid%Fjmax,2,grid%gm)
      call bocosHn(mpl,grid,state%HK3D,grid%km3,grid%im,grid%jm,grid%lm,grid%hx,grid%hy,grid%Fimax,grid%Fjmax,2,grid%gm)

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!
! Filtering
!

      call vrbeta2(grid%km2,grid%hx,0,grid%im,grid%hy,0,grid%jm,grid%p,state%pasp2,state%ss2,state%VK2D(:,:,:))
      call sup_vrbeta3(state%VK3D,grid%km3,grid%hx,grid%hy,grid%hz,grid%im,grid%jm,grid%lm,grid%p,state%pasp3,state%ss3)

  if(grid%l_hgen)  then
      call vrbeta2(grid%km2,grid%hx,0,grid%im,grid%hy,0,grid%jm,grid%p,state%pasp2,state%ss2,state%HK2D(:,:,:))
      call sup_vrbeta3(state%HK3D,grid%km3,grid%hx,grid%hy,grid%hz,grid%im,grid%jm,grid%lm,grid%p,state%pasp3,state%ss3)
  endif

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff


!                                                 call etim(    bfilt_tim)
!***
!*** Apply weights at all generations again
!***

!                                                 call btim( weight_tim)
        call weighting1(grid,state,state%VK2D,state%HK2D)
        call weighting1(grid,state,state%VK3D,state%HK3D)

!                                                 call etim( weight_tim)

!***
!***  Repeat:
!***  Downsend, interpolate and add (Step 4)
!***  Then zero high generations (Step 5)
!***

!                                                 call btim(   dnsend_tim)
       call downsending1(mpl,grid,state%HK2D,state%VK2D)
       call downsending1(mpl,grid,state%HK3D,state%VK3D)

!                                                 call etim(   dnsend_tim)



!-----------------------------------------------------------------------
                        endsubroutine mg_filtering_procedure3    

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine mg_filtering_procedure5(mpl,grid,state)                        
!***********************************************************************
!                                                                      !
! Multigrid filtering procedure 2:                                     !
!                                                                      !
!     - Multiple of 2D and 3D variables                                !
!     - 1 upsending and downsending                                    !
!     - Applicaton of Hegrid%lmholtz differential operator                  !
!     - 2d radial filter + 1d vertical filter                          !
!                                                                      !
!***********************************************************************
!-----------------------------------------------------------------------
!use mg_intstate1, only: state%VK2D,state%VK3D,state%HK2D,state%HK3D,state%pasp1,state%pasp2,state%ss1,state%ss2
IMPLICIT NONE

type(mpl_type),intent(inout) :: mpl
type(mgbf_grid_type),intent(inout) :: grid
type(mgbf_state_type),intent(inout) :: state
integer(kind_int) L

!==================== Adjoint (Conservative step) ======================

!***
!*** Adjoint interpolate and upsend (Step 1)
!***
     
!                                                 call btim( upsend_tim)
       call upsending1(mpl,grid,state%VK2D,state%HK2D)
       call upsending1(mpl,grid,state%VK3D,state%HK3D)
!                                                 call etim( upsend_tim)


!***
!*** Apply adjoint of Beta filter at all generations (Step 3)
!***
!                                                 call btim(    bfiltT_tim)



!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

      call vrbeta2T(grid%km2,grid%hx,0,grid%im,grid%hy,0,grid%jm,grid%p,state%pasp2,state%ss2,state%VK2D(:,:,:))
    do L=1,grid%lm
      call vrbeta2T(grid%km3,grid%hx,0,grid%im,grid%hy,0,grid%jm,grid%p,state%pasp2,state%ss2,state%VK3D(:,:,l,:))
    end do


  if(grid%l_hgen)  then
      call vrbeta2T(grid%km2,grid%hx,0,grid%im,grid%hy,0,grid%jm,grid%p,state%pasp2,state%ss2,state%HK2D(:,:,:))
    do L=1,grid%lm
      call vrbeta2T(grid%km3,grid%hx,0,grid%im,grid%hy,0,grid%jm,grid%p,state%pasp2,state%ss2,state%HK3D(:,:,l,:))
    end do
   endif

      call sup_vrbeta1T(state%VK3D,grid%km3,grid%hx,grid%hy,grid%hz,grid%im,grid%jm,grid%lm,grid%p,state%pasp1,state%ss1)
  if(grid%l_hgen)  then
      call sup_vrbeta1T(state%HK3D,grid%km3,grid%hx,grid%hy,grid%hz,grid%im,grid%jm,grid%lm,grid%p,state%pasp1,state%ss1)
   endif


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff


      call bocosHnT(mpl,grid,state%VK2D,grid%km2,grid%im,grid%jm, 1,grid%hx,grid%hy,grid%Fimax,grid%Fjmax,1,1)
      call bocosHnT(mpl,grid,state%VK3D,grid%km3,grid%im,grid%jm,grid%lm,grid%hx,grid%hy,grid%Fimax,grid%Fjmax,1,1)

      call bocosHnT(mpl,grid,state%HK2D,grid%km2,grid%im,grid%jm, 1,grid%hx,grid%hy,grid%Fimax,grid%Fjmax,2,grid%gm)
      call bocosHnT(mpl,grid,state%HK3D,grid%km3,grid%im,grid%jm,grid%lm,grid%hx,grid%hy,grid%Fimax,grid%Fjmax,2,grid%gm)



!                                                 call etim(    bfiltT_tim)
!***
!*** Apply (a-b\nabla^2)
!***

      call bocosHn(mpl,grid,state%VK2D,grid%km2,grid%im,grid%jm,1 ,grid%hx,grid%hy,grid%Fimax,grid%Fjmax,1,1)
      call bocosHn(mpl,grid,state%VK3D,grid%km3,grid%im,grid%jm,grid%lm,grid%hx,grid%hy,grid%Fimax,grid%Fjmax,1,1)

      call bocosHn(mpl,grid,state%HK2D,grid%km2,grid%im,grid%jm,1 ,grid%hx,grid%hy,grid%Fimax,grid%Fjmax,2,grid%gm)
      call bocosHn(mpl,grid,state%HK3D,grid%km3,grid%im,grid%jm,grid%lm,grid%hx,grid%hy,grid%Fimax,grid%Fjmax,2,grid%gm)

!                                                call btim( weight_tim)

      call differencing1(grid,state,state%VK2D,state%HK2D)
      call differencing1(grid,state,state%VK3D,state%HK3D)
!                                                call etim( weight_tim)


!***
!*** Apply Beta filter at all generations (Step 7)
!***


!                                                 call btim(    bfilt_tim)

      call bocosHn(mpl,grid,state%VK2D,grid%km2,grid%im,grid%jm,1 ,grid%hx,grid%hy,grid%Fimax,grid%Fjmax,1,1)
      call bocosHn(mpl,grid,state%VK3D,grid%km3,grid%im,grid%jm,grid%lm,grid%hx,grid%hy,grid%Fimax,grid%Fjmax,1,1)

      call bocosHn(mpl,grid,state%HK2D,grid%km2,grid%im,grid%jm,1 ,grid%hx,grid%hy,grid%Fimax,grid%Fjmax,2,grid%gm)
      call bocosHn(mpl,grid,state%HK3D,grid%km3,grid%im,grid%jm,grid%lm,grid%hx,grid%hy,grid%Fimax,grid%Fjmax,2,grid%gm)

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!
! Filtering
!

      call vrbeta2(grid%km2,grid%hx,0,grid%im,grid%hy,0,grid%jm,grid%p,state%pasp2,state%ss2,state%VK2D(:,:,:))
   do L=1,grid%lm
      call vrbeta2(grid%km3,grid%hx,0,grid%im,grid%hy,0,grid%jm,grid%p,state%pasp2,state%ss2,state%VK3D(:,:,l,:))
   enddo

  if(grid%l_hgen)  then
      call vrbeta2(grid%km2,grid%hx,0,grid%im,grid%hy,0,grid%jm,grid%p,state%pasp2,state%ss2,state%HK2D(:,:,:))
   do L=1,grid%lm
      call vrbeta2(grid%km3,grid%hx,0,grid%im,grid%hy,0,grid%jm,grid%p,state%pasp2,state%ss2,state%HK3D(:,:,l,:))
   enddo
  endif

      call sup_vrbeta1(state%VK3D,grid%km3,grid%hx,grid%hy,grid%hz,grid%im,grid%jm,grid%lm,grid%p,state%pasp1,state%ss1)
  if(grid%l_hgen)  then
      call sup_vrbeta1(state%HK3D,grid%km3,grid%hx,grid%hy,grid%hz,grid%im,grid%jm,grid%lm,grid%p,state%pasp1,state%ss1)
   endif

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

!                                                 call etim(    bfilt_tim)

!***
!***  Downsend, interpolate and add (Step 4)
!***  Then zero high generations (Step 5)
!***

!                                                 call btim(   dnsend_tim)
       call downsending1(mpl,grid,state%HK2D,state%VK2D)
       call downsending1(mpl,grid,state%HK3D,state%VK3D)

!                                                 call etim(   dnsend_tim)



!-----------------------------------------------------------------------
                        endsubroutine mg_filtering_procedure5    

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine sup_vrbeta1                         &
!**********************************************************************
!                                                                     *
!     conversion of vrbeta1                                           *
!                                                                     *
!**********************************************************************
(V,kmax,hx,hy,hz,im,jm,lm, p, pasp,ss)
!----------------------------------------------------------------------
IMPLICIT NONE

integer(kind_int),intent(in):: kmax,hx,hy,hz,im,jm,lm, p
real(kind_real),DIMENSION(-hx:im+hx,-hy:jm+hy,1:lm,1:kmax),intent(inout):: V
real(kind_real),DIMENSION(1,1,1:lm), intent(in):: pasp
real(kind_real),DIMENSION(1-hz,lm+hz), intent(in):: ss

real(kind_real),DIMENSION(1:kmax,1-hz:lm+hz):: W

integer(kind_int):: i,j,L

!----------------------------------------------------------------------

        do j=0,jm
        do i=0,im
          do L=1,Lm
            W(:,L)=V(i,j,L,:)
          end do
          do L=1,hz
            W(:,1-L)=W(:,1+L)
            W(:,LM+L)=W(:,LM-L)
          end do
             call vrbeta1(kmax,hz,1,lm, p, pasp,ss,W)
          do l=1,Lm
            V(i,j,L,:)=W(:,L)
          end do
        end do
        end do

  
!----------------------------------------------------------------------
                        endsubroutine sup_vrbeta1

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine sup_vrbeta1T                        &
!**********************************************************************
!                                                                     *
!     conversion of vrbeta1T                                          *
!                                                                     *
!**********************************************************************
(V,kmax,hx,hy,hz,im,jm,lm, p, pasp,ss)
!----------------------------------------------------------------------
IMPLICIT NONE

integer(kind_int),intent(in):: kmax,hx,hy,hz,im,jm,lm, p
real(kind_real),DIMENSION(-hx:im+hx,-hy:jm+hy,1:lm,1:kmax),intent(inout):: V
real(kind_real),DIMENSION(1,1,1:lm), intent(in):: pasp
real(kind_real),DIMENSION(1-hz,lm+hz), intent(in):: ss

real(kind_real),DIMENSION(1:kmax,1-hz:lm+hz):: W

integer(kind_int):: i,j,L

!----------------------------------------------------------------------

        do j=0,jm
        do i=0,im
          do L=1,Lm
            W(:,L)=V(i,j,L,:)
          end do
          do L=1,hz
            W(:,1-L)=W(:,1+L)
            W(:,LM+L)=W(:,LM-L)
          end do
             call vrbeta1T(kmax,hz,1,lm, p, pasp,ss,W)
          do l=1,Lm
            V(i,j,L,:)=W(:,L)
          end do
        end do
        end do

!----------------------------------------------------------------------
                        endsubroutine sup_vrbeta1T

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine sup_vrbeta3                         &
!**********************************************************************
!                                                                     *
!     conversion of vrbeta3                                           *
!                                                                     *
!**********************************************************************
(V,kmax,hx,hy,hz,im,jm,lm, p, pasp,ss)
!----------------------------------------------------------------------
IMPLICIT NONE

integer(kind_int),intent(in):: kmax,hx,hy,hz,im,jm,lm, p
real(kind_real),DIMENSION(-hx:im+hx,-hy:jm+hy,1:lm,1:kmax),intent(inout):: V
real(kind_real),DIMENSION(3,3,0:im,0:jm,1:lm), intent(in):: pasp
real(kind_real),DIMENSION(0:im,0:jm,1:lm), intent(in):: ss

real(kind_real),DIMENSION(1:kmax,-hx:im+hx,-hy:jm+hy,1-hz:lm+hz):: W

integer(kind_int):: i,j,L

!----------------------------------------------------------------------

          do L=1,Lm
          do j=-hy,jm+hy
          do i=-hx,im+hx
            W(:,i,j,L)=V(i,j,L,:)
          end do
          end do
          end do

          do j=-hy,jm+hy
          do i=-hx,im+hx
            do L=1,hz
              W(:,i,j,1-L)=W(:,i,j,1+L)
              W(:,i,j,LM+L)=W(:,i,j,LM-L)
            end do
          end do
          end do
    
    
           call vrbeta3(kmax,hx,0,im, hy,0,jm, hz,1,lm, p, pasp,ss,W)

  
          do l=1,Lm
          do j=0,jm
          do i=0,im
            V(i,j,L,:)=W(:,i,j,L)
          end do
          end do
          end do

!----------------------------------------------------------------------
                        endsubroutine sup_vrbeta3
 
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine sup_vrbeta3T                        &
!**********************************************************************
!                                                                     *
!     conversion of vrbeta3                                           *
!                                                                     *
!**********************************************************************
(V,kmax,hx,hy,hz,im,jm,lm, p, pasp,ss)
!----------------------------------------------------------------------
IMPLICIT NONE

integer(kind_int),intent(in):: kmax,hx,hy,hz,im,jm,lm, p
real(kind_real),DIMENSION(-hx:im+hx,-hy:jm+hy,1:lm,1:kmax),intent(inout):: V
real(kind_real),DIMENSION(3,3,0:im,0:jm,1:lm), intent(in):: pasp
real(kind_real),DIMENSION(0:im,0:jm,1:lm), intent(in):: ss

real(kind_real),DIMENSION(1:kmax,-hx:im+hx,-hy:jm+hy,1-hz:lm+hz):: W

integer(kind_int):: i,j,l

!----------------------------------------------------------------------

          do L=1,Lm
          do j=-hy,jm+hy
          do i=-hx,im+hx
            W(:,i,j,L)=V(i,j,L,:)
          end do
          end do
          end do

          do j=-hy,jm+hy
          do i=-hx,im+hx
            do L=1,hz
              W(:,i,j,1-L)=W(:,i,j,1+L)
              W(:,i,j,LM+L)=W(:,i,j,LM-L)
            end do
          end do
          end do
    
    
           call vrbeta3T(kmax,hx,0,im, hy,0,jm, hz,1,lm, p, pasp,ss,W)

  
          do l=1,lm
          do j=0,jm
          do i=0,im
            V(i,j,l,:)=W(:,i,j,l)
          end do
          end do
          end do

!----------------------------------------------------------------------
                        endsubroutine sup_vrbeta3T
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        endmodule tools_filtering
