! (C) Copyright 2018-2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module type_mgbf_interface

use atlas_module
use fckit_configuration_module, only: fckit_configuration
use fckit_mpi_module, only: fckit_mpi_comm
use iso_c_binding
use type_fieldset, only: fieldset_type
use type_mgbf, only: mgbf_type,mgbf_registry

implicit none

private

contains
!----------------------------------------------------------------------
! Subroutine: mgbf_create_c
! Purpose: create
!----------------------------------------------------------------------
subroutine mgbf_create_c(key_mgbf,c_comm,c_afunctionspace,c_afieldset,c_conf,c_grid) bind(c,name='mgbf_create_f90')

implicit none

! Passed variables
integer(c_int),intent(inout) :: key_mgbf         ! GSIRF
type(c_ptr),intent(in),value :: c_comm           ! FCKIT MPI communicator wrapper
type(c_ptr),intent(in),value :: c_afunctionspace ! ATLAS function space
type(c_ptr),intent(in),value :: c_afieldset      ! ATLAS fieldset  (containing geometry features: area, vunit, gmask, smask)
type(c_ptr),intent(in),value :: c_conf           ! FCKIT configuration
type(c_ptr),intent(in),value :: c_grid           ! FCKIT grid configuration

! Local variables
type(mgbf_type),pointer :: mgbf
type(fckit_mpi_comm) :: f_comm
type(atlas_functionspace) :: f_afunctionspace
type(fieldset_type) :: f_fieldset
type(fckit_configuration) :: f_conf
type(fckit_configuration) :: f_grid

! Interface
call mgbf_registry%init()
call mgbf_registry%add(key_mgbf)
call mgbf_registry%get(key_mgbf,mgbf)
f_comm = fckit_mpi_comm(c_comm)
f_afunctionspace = atlas_functionspace(c_afunctionspace)
f_fieldset = atlas_fieldset(c_afieldset)
f_conf = fckit_configuration(c_conf)
f_grid = fckit_configuration(c_grid)

! Call Fortran
call mgbf%create(f_comm,f_afunctionspace,f_fieldset,f_conf,f_grid)

end subroutine mgbf_create_c

!----------------------------------------------------------------------
! Subroutine: mgbf_run_drivers_c
! Purpose: run drivers
!----------------------------------------------------------------------
subroutine mgbf_run_drivers_c(key_mgbf) bind(c,name='mgbf_run_drivers_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: key_mgbf ! MGBF

! Local variables
type(mgbf_type),pointer :: mgbf

! Interface
call mgbf_registry%get(key_mgbf,mgbf)

! Call Fortran
!call mgbf%run_drivers

end subroutine mgbf_run_drivers_c

!----------------------------------------------------------------------
! Subroutine: mgbf_apply_vbal_c
! Purpose: vertical balance application
!----------------------------------------------------------------------
subroutine mgbf_apply_vbal_c(key_mgbf,c_afieldset) bind(c,name='mgbf_apply_vbal_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: key_mgbf       ! MGBF
type(c_ptr),intent(in),value :: c_afieldset ! ATLAS fieldset pointer

! Local variables
type(mgbf_type),pointer :: mgbf
type(fieldset_type) :: f_fieldset

! Interface
call mgbf_registry%get(key_mgbf,mgbf)
f_fieldset = atlas_fieldset(c_afieldset)

! Call Fortran
!call mgbf%apply_vbal(f_fieldset)

end subroutine mgbf_apply_vbal_c

!----------------------------------------------------------------------
! Subroutine: mgbf_apply_vbal_inv_c
! Purpose: vertical balance application, inverse
!----------------------------------------------------------------------
subroutine mgbf_apply_vbal_inv_c(key_mgbf,c_afieldset) bind(c,name='mgbf_apply_vbal_inv_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: key_mgbf       ! MGBF
type(c_ptr),intent(in),value :: c_afieldset ! ATLAS fieldset pointer

! Local variables
type(mgbf_type),pointer :: mgbf
type(fieldset_type) :: f_fieldset

! Interface
call mgbf_registry%get(key_mgbf,mgbf)
f_fieldset = atlas_fieldset(c_afieldset)

! Call Fortran
!call mgbf%apply_vbal_inv(f_fieldset)

end subroutine mgbf_apply_vbal_inv_c

!----------------------------------------------------------------------
! Subroutine: mgbf_apply_vbal_ad_c
! Purpose: vertical balance application, adjoint
!----------------------------------------------------------------------
subroutine mgbf_apply_vbal_ad_c(key_mgbf,c_afieldset) bind(c,name='mgbf_apply_vbal_ad_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: key_mgbf       ! MGBF
type(c_ptr),intent(in),value :: c_afieldset ! ATLAS fieldset pointer

! Local variables
type(mgbf_type),pointer :: mgbf
type(fieldset_type) :: f_fieldset

! Interface
call mgbf_registry%get(key_mgbf,mgbf)
f_fieldset = atlas_fieldset(c_afieldset)

! Call Fortran
!call mgbf%apply_vbal_ad(f_fieldset)

end subroutine mgbf_apply_vbal_ad_c

!----------------------------------------------------------------------
! Subroutine: mgbf_apply_vbal_inv_ad_c
! Purpose: vertical balance application, inverse adjoint
!----------------------------------------------------------------------
subroutine mgbf_apply_vbal_inv_ad_c(key_mgbf,c_afieldset) bind(c,name='mgbf_apply_vbal_inv_ad_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: key_mgbf       ! MGBF
type(c_ptr),intent(in),value :: c_afieldset ! ATLAS fieldset pointer

! Local variables
type(mgbf_type),pointer :: mgbf
type(fieldset_type) :: f_fieldset

! Interface
call mgbf_registry%get(key_mgbf,mgbf)
f_fieldset = atlas_fieldset(c_afieldset)

! Call Fortran
!call mgbf%apply_vbal_inv_ad(f_fieldset)

end subroutine mgbf_apply_vbal_inv_ad_c

!----------------------------------------------------------------------
! Subroutine: mgbf_apply_stddev_c
! Purpose: standard-deviation application
!----------------------------------------------------------------------
subroutine mgbf_apply_stddev_c(key_mgbf,c_afieldset) bind(c,name='mgbf_apply_stddev_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: key_mgbf       ! MGBF
type(c_ptr),intent(in),value :: c_afieldset ! ATLAS fieldset pointer

! Local variables
type(mgbf_type),pointer :: mgbf
type(fieldset_type) :: f_fieldset

! Interface
call mgbf_registry%get(key_mgbf,mgbf)
f_fieldset = atlas_fieldset(c_afieldset)

! Call Fortran
!call mgbf%apply_stddev(f_fieldset)

end subroutine mgbf_apply_stddev_c

!----------------------------------------------------------------------
! Subroutine: mgbf_apply_stddev_inv_c
! Purpose: standard-deviation application, inverse
!----------------------------------------------------------------------
subroutine mgbf_apply_stddev_inv_c(key_mgbf,c_afieldset) bind(c,name='mgbf_apply_stddev_inv_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: key_mgbf       ! MGBF
type(c_ptr),intent(in),value :: c_afieldset ! ATLAS fieldset pointer

! Local variables
type(mgbf_type),pointer :: mgbf
type(fieldset_type) :: f_fieldset

! Interface
call mgbf_registry%get(key_mgbf,mgbf)
f_fieldset = atlas_fieldset(c_afieldset)

! Call Fortran
!call mgbf%apply_stddev_inv(f_fieldset)

end subroutine mgbf_apply_stddev_inv_c

!----------------------------------------------------------------------
! Subroutine: mgbf_apply_filter_c
! Purpose: NICAS application
!----------------------------------------------------------------------
subroutine mgbf_apply_filter_c(key_mgbf,c_afieldset) bind(c,name='mgbf_apply_filter_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: key_mgbf       ! MGBF
type(c_ptr),intent(in),value :: c_afieldset ! ATLAS fieldset pointer

! Local variables
type(mgbf_type),pointer :: mgbf
type(fieldset_type) :: f_fieldset

! Interface
call mgbf_registry%get(key_mgbf,mgbf)
f_fieldset = atlas_fieldset(c_afieldset)

! Call Fortran
call mgbf%apply_filter(f_fieldset)

end subroutine mgbf_apply_filter_c

!----------------------------------------------------------------------
! Subroutine: mgbf_get_cv_size_c
! Purpose: get control variable size
!----------------------------------------------------------------------
subroutine mgbf_get_cv_size_c(key_mgbf,n) bind(c,name='mgbf_get_cv_size_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: key_mgbf ! MGBF
integer(c_int),intent(out) :: n       ! Control variable size

! Local variables
type(mgbf_type),pointer :: mgbf

! Interface
call mgbf_registry%get(key_mgbf,mgbf)

! Call Fortran
!call mgbf%get_cv_size(n)

end subroutine mgbf_get_cv_size_c

!----------------------------------------------------------------------
! Subroutine: mgbf_apply_filter_sqrt
! Purpose: NICAS square-root application
!----------------------------------------------------------------------
subroutine mgbf_apply_filter_sqrt_c(key_mgbf,cv,c_afieldset) bind(c,name='mgbf_apply_filter_sqrt_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: key_mgbf       ! MGBF
real(c_double),intent(in) :: cv(:)          ! Control variable
type(c_ptr),intent(in),value :: c_afieldset ! ATLAS fieldset pointer

! Local variables
type(mgbf_type),pointer :: mgbf
type(fieldset_type) :: f_fieldset

! Interface
call mgbf_registry%get(key_mgbf,mgbf)
f_fieldset = atlas_fieldset(c_afieldset)

! Call Fortran
!call mgbf%apply_filter_sqrt(cv,f_fieldset)

end subroutine mgbf_apply_filter_sqrt_c

!----------------------------------------------------------------------
! Subroutine: mgbf_apply_filter_sqrt_ad_c
! Purpose: NICAS square-root adjoint application
!----------------------------------------------------------------------
subroutine mgbf_apply_filter_sqrt_ad_c(key_mgbf,c_afieldset,cv) bind(c,name='mgbf_apply_filter_sqrt_ad_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: key_mgbf       ! MGBF
type(c_ptr),intent(in),value :: c_afieldset ! ATLAS fieldset pointer
real(c_double),intent(inout) :: cv(:)       ! Control variable

! Local variables
type(mgbf_type),pointer :: mgbf
type(fieldset_type) :: f_fieldset

! Interface
call mgbf_registry%get(key_mgbf,mgbf)
f_fieldset = atlas_fieldset(c_afieldset)

! Call Fortran
!call mgbf%apply_filter_sqrt_ad(f_fieldset,cv)

end subroutine mgbf_apply_filter_sqrt_ad_c

!----------------------------------------------------------------------
! Subroutine: mgbf_randomize_c
! Purpose: NICAS randomization
!----------------------------------------------------------------------
subroutine mgbf_randomize_c(key_mgbf,c_afieldset) bind(c,name='mgbf_randomize_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: key_mgbf       ! MGBF
type(c_ptr),intent(in),value :: c_afieldset ! ATLAS fieldset pointer

! Local variables
type(mgbf_type),pointer :: mgbf
type(fieldset_type) :: f_fieldset

! Interface
call mgbf_registry%get(key_mgbf,mgbf)
f_fieldset = atlas_fieldset(c_afieldset)

! Call Fortran
!call mgbf%randomize(f_fieldset)

end subroutine mgbf_randomize_c

!----------------------------------------------------------------------
! Subroutine: mgbf_get_parameter_c
! Purpose: get a parameter
!----------------------------------------------------------------------
subroutine mgbf_get_parameter_c(key_mgbf,nstr,cstr,c_afieldset) bind(c,name='mgbf_get_parameter_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: key_mgbf           ! MGBF
integer(c_int),intent(in) :: nstr               ! Parameter name size
character(kind=c_char),intent(in) :: cstr(nstr) ! Parameter name
type(c_ptr),intent(in),value :: c_afieldset     ! ATLAS fieldset pointer

! Local variables
type(mgbf_type),pointer :: mgbf
integer :: istr
character(len=nstr) :: param
type(fieldset_type) :: f_fieldset

! Interface
call mgbf_registry%get(key_mgbf,mgbf)
param = ''
do istr=1,nstr
  param = trim(param)//cstr(istr)
end do
f_fieldset = atlas_fieldset(c_afieldset)

! Call Fortran
!call mgbf%get_parameter(param,f_fieldset)

end subroutine mgbf_get_parameter_c

!----------------------------------------------------------------------
! Subroutine: mgbf_set_parameter_c
! Purpose: set a parameter
!----------------------------------------------------------------------
subroutine mgbf_set_parameter_c(key_mgbf,nstr,cstr,c_afieldset) bind(c,name='mgbf_set_parameter_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: key_mgbf           ! MGBF
integer(c_int),intent(in) :: nstr               ! Parameter name size
character(kind=c_char),intent(in) :: cstr(nstr) ! Parameter name
type(c_ptr),intent(in),value :: c_afieldset     ! ATLAS fieldset pointer

! Local variables
type(mgbf_type),pointer :: mgbf
integer :: istr
character(len=nstr) :: param
type(fieldset_type) :: f_fieldset

! Interface
call mgbf_registry%get(key_mgbf,mgbf)
param = ''
do istr=1,nstr
  param = trim(param)//cstr(istr)
end do
f_fieldset = atlas_fieldset(c_afieldset)

! Call Fortran
!call mgbf%set_parameter(param,f_fieldset)

end subroutine mgbf_set_parameter_c

!----------------------------------------------------------------------
! Subroutine: mgbf_dealloc_c
! Purpose: deallocation
!----------------------------------------------------------------------
subroutine mgbf_dealloc_c(key_mgbf) bind(c,name='mgbf_dealloc_f90')

implicit none

! Passed variables
integer(c_int),intent(inout) :: key_mgbf ! MGBF

! Local variables
type(mgbf_type),pointer :: mgbf

! Interface
call mgbf_registry%get(key_mgbf,mgbf)

! Deallocate MGBF
call mgbf%dealloc

! Clean interface
call mgbf_registry%remove(key_mgbf)

end subroutine mgbf_dealloc_c

end module type_mgbf_interface
