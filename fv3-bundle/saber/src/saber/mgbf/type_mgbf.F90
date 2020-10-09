! (C) Copyright 2018-2019 UCAR
!
! This software is licensed under the terms of the Apache Licence
! Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module type_mgbf

use atlas_module
use fckit_configuration_module, only: fckit_configuration
use fckit_mpi_module, only: fckit_mpi_comm,fckit_mpi_sum
use tools_atlas, only: create_atlas_function_space
use tools_const, only: req,deg2rad,rad2deg
use tools_kinds,only: kind_int,kind_real
use type_fieldset, only: fieldset_type
use type_geom, only: geom_type

use tools_interpolate
use tools_transfer
use tools_filtering
use tools_sp

use type_mpl, only: mpl_type
use type_bump, only: bump_type
use type_mgbf_nam, only: mgbf_nam_type
use type_mgbf_grid, only: mgbf_grid_type
use type_mgbf_state, only: mgbf_state_type

implicit none

! MGBF derived type
type mgbf_type
  type(mpl_type) :: mpl
  type(mgbf_nam_type) :: nam
  type(mgbf_grid_type) :: grid
  type(mgbf_state_type) :: state
  type(bump_type) :: cs2gg
  type(atlas_functionspace) :: afunctionspace_cs2gg
  integer :: ngrid
  integer :: idrt
  type(bump_type) :: gg2cs
  type(atlas_functionspace) :: afunctionspace_gg2cs
contains
  procedure :: create => mgbf_create
  procedure :: setup => mgbf_setup
  procedure :: apply_filter => mgbf_apply_filter
  procedure :: dealloc => mgbf_dealloc
end type mgbf_type

integer,parameter :: dmsvali = -999           ! Default missing value for integers
real(kind_real),parameter :: dmsvalr = -999.0 ! Default missing value for reals

private
public :: mgbf_type
public :: mgbf_registry

! MGBF registry
#define LISTED_TYPE mgbf_type

! Linked list interface - defines registry_t type
#include "saber/util/linkedList_i.f"

! Global registry
type(registry_t) :: mgbf_registry

contains

!----------------------------------------------------------------------
! Linked list implementation
!----------------------------------------------------------------------
#include "saber/util/linkedList_c.f"

!----------------------------------------------------------------------
! Subroutine: mgbf_create
! Purpose: create
!----------------------------------------------------------------------
subroutine mgbf_create(mgbf,comm,afunctionspace,fieldset,conf,grid)

implicit none

! Passed variables
class(mgbf_type),intent(inout) :: mgbf               ! MGBF
type(fckit_mpi_comm),intent(in) :: comm                ! FCKIT MPI communicator wrapper
type(atlas_functionspace),intent(in) :: afunctionspace ! ATLAS function space
type(fieldset_type),intent(in) :: fieldset             ! Fieldset containing geometry elements
type(fckit_configuration),intent(in) :: conf           ! FCKIT configuration
type(fckit_configuration),intent(in) :: grid           ! FCKIT grid configuration

! Local variables
integer :: lmsvali, llunit
real(kind_real) :: lmsvalr

call mgbf%nam%init(comm%size())
call mgbf%grid%init(comm%size(),comm%rank())

! Read configuration
call mgbf%nam%from_conf(mgbf%grid,conf)

! Read grid configuration
call mgbf%nam%from_conf(mgbf%grid,grid)

! Set missing values
lmsvali = dmsvali
lmsvalr = dmsvalr
if (conf%has('msvali')) call conf%get_or_die('msvali',lmsvali)
if (conf%has('msvalr')) call conf%get_or_die('msvalr',lmsvalr)

! Set log unit
llunit = lmsvali
if (conf%has('lunit')) call conf%get_or_die('lunit',llunit)

! Setup MGBF
call mgbf%setup(comm,afunctionspace,fieldset, &
&               lunit=llunit,msvali=lmsvali,msvalr=lmsvalr)

end subroutine mgbf_create

!----------------------------------------------------------------------
! Subroutine: mgbf_setup
! Purpose: setup
!----------------------------------------------------------------------
subroutine mgbf_setup(mgbf,f_comm,afunctionspace,fieldset, &
                     & nobs,lonobs,latobs,lunit,msvali,msvalr)

implicit none

! Passed variables
class(mgbf_type),intent(inout) :: mgbf               ! MGBF
type(fckit_mpi_comm),intent(in) :: f_comm              ! FCKIT MPI communicator wrapper
type(atlas_functionspace),intent(in) :: afunctionspace ! ATLAS functionspace
type(fieldset_type),intent(in),optional :: fieldset    ! Fieldset containing geometry elements
integer,intent(in),optional :: nobs                    ! Number of observations
real(kind_real),intent(in),optional :: lonobs(:)       ! Observations longitude (in degrees)
real(kind_real),intent(in),optional :: latobs(:)       ! Observations latitude (in degrees)
integer,intent(in),optional :: lunit                   ! Listing unit
integer,intent(in),optional :: msvali                  ! Missing value for integers
real(kind_real),intent(in),optional :: msvalr          ! Missing value for reals

! Local variables
integer :: lmsvali,iv,its
real(kind_real) :: lmsvalr
real(kind_real),pointer :: real_ptr(:,:)
character(len=1024) :: fieldname
character(len=1024),parameter :: subr = 'mgbf_setup'
type(atlas_field) :: afield

! Initialize MPL
call mgbf%mpl%init(f_comm)

! Set missing values
mgbf%mpl%msv%vali = dmsvali
mgbf%mpl%msv%valr = dmsvalr
if (present(msvali)) mgbf%mpl%msv%vali = msvali
if (present(msvalr)) mgbf%mpl%msv%valr = msvalr

! Initialize listing
mgbf%mpl%lunit = mgbf%mpl%msv%vali
if (present(lunit)) mgbf%mpl%lunit = lunit
mgbf%mpl%verbosity = 'all'
mgbf%mpl%black = ' '
mgbf%mpl%green = ' '
mgbf%mpl%peach = ' '
mgbf%mpl%aqua = ' '
mgbf%mpl%purple = ' '
mgbf%mpl%err = ' '
mgbf%mpl%wng = ' '

! Check namelist parameters
!write(mgbf%mpl%info,'(a)') '-------------------------------------------------------------------'
!call mgbf%mpl%flush
!write(mgbf%mpl%info,'(a)') '--- Check namelist parameters'
!call mgbf%mpl%flush
!call mgbf%nam%check(mgbf%mpl,mgbf%grid,mgbf%berr)
!if(mgbf%mpl%main) call mgbf%nam%write(mgbf%mpl,mgbf%grid)

call mgbf%grid%setup(mgbf%nam%nv,mgbf%nam%variables)

call mgbf%state%setup(mgbf%grid)

call lsqr_mg_coef(mgbf%grid,mgbf%state)

call mgbf_bump_init(mgbf,f_comm,afunctionspace)

end subroutine mgbf_setup

!----------------------------------------------------------------------
! Subroutine: mgbf_apply_filter
! Purpose: apply recursive filter
!----------------------------------------------------------------------
subroutine mgbf_apply_filter(mgbf,fieldset)

implicit none

class(mgbf_type), intent(inout) :: mgbf
type(fieldset_type),intent(inout) :: fieldset             ! Fieldset containing geometry elements

! Local variable
integer :: its,iv
real(kind_real) :: fld_mga(mgbf%cs2gg%geom%nmga,mgbf%cs2gg%geom%nl0,mgbf%nam%nv)
integer :: i,j,k,ij

real(kind_real) :: fld_ngrid(mgbf%ngrid,mgbf%cs2gg%geom%nl0)

integer :: n,t

! Initialize fieldset
call fieldset%init(mgbf%mpl,mgbf%cs2gg%geom%nmga,mgbf%cs2gg%geom%nl0, &
 & mgbf%cs2gg%geom%gmask_mga,mgbf%nam%variables(1:mgbf%nam%nv),     &
 & mgbf%nam%lev2d)

! Fieldset to Fortran array on model grid
call fieldset%to_array(mgbf%mpl,fld_mga)

do n=1,mgbf%nam%nv
  select case(trim(mgbf%nam%variables(n)))
    case('ua','sf','psi')
     call fld_interp_cs2gg(mgbf,mgbf%cs2gg%geom%nl0,fld_mga(:,:,n),mgbf%state%WA1) 
    case('va','vp','chi')
     call fld_interp_cs2gg(mgbf,mgbf%cs2gg%geom%nl0,fld_mga(:,:,n),mgbf%state%WA20)
    case('t','tv')
     call fld_interp_cs2gg(mgbf,mgbf%cs2gg%geom%nl0,fld_mga(:,:,n),mgbf%state%WA3)
    case('q','sphum','rh')
     call fld_interp_cs2gg(mgbf,mgbf%cs2gg%geom%nl0,fld_mga(:,:,n),mgbf%state%WA4)
    case("cw","liq_wat")
     call fld_interp_cs2gg(mgbf,mgbf%cs2gg%geom%nl0,fld_mga(:,:,n),mgbf%state%WA5)
    case("ice_wat")
     call fld_interp_cs2gg(mgbf,mgbf%cs2gg%geom%nl0,fld_mga(:,:,n),mgbf%state%WA6)
    case('ps')
     call fld_interp_cs2gg(mgbf,1,fld_mga(:,:,n),mgbf%state%PA1)
  end select
end do

call anal_to_filt(mgbf%mpl,mgbf%grid,mgbf%state)

call mg_filtering_procedure(mgbf%mpl,mgbf%grid,mgbf%state)

call filt_to_anal(mgbf%mpl,mgbf%grid,mgbf%state)

do n=1,mgbf%nam%nv
  select case(trim(mgbf%nam%variables(n)))
    case('ua','sf','psi')
     call fld_interp_gg2cs(mgbf,mgbf%cs2gg%geom%nl0,mgbf%state%WA1,fld_mga(:,:,n))
    case('va','vp','chi')
     call fld_interp_gg2cs(mgbf,mgbf%cs2gg%geom%nl0,mgbf%state%WA20,fld_mga(:,:,n))
    case('t','tv')
     call fld_interp_gg2cs(mgbf,mgbf%cs2gg%geom%nl0,mgbf%state%WA3,fld_mga(:,:,n))
    case('q','sphum','rh')
     call fld_interp_gg2cs(mgbf,mgbf%cs2gg%geom%nl0,mgbf%state%WA4,fld_mga(:,:,n))
    case("cw","liq_wat")
     call fld_interp_gg2cs(mgbf,mgbf%cs2gg%geom%nl0,mgbf%state%WA5,fld_mga(:,:,n))
    case("ice_wat")
     call fld_interp_gg2cs(mgbf,mgbf%cs2gg%geom%nl0,mgbf%state%WA6,fld_mga(:,:,n))
    case('ps')
     call fld_interp_gg2cs(mgbf,1,mgbf%state%PA1,fld_mga(:,:,n))
  end select
end do

! Fortran array on model grid to fieldset
call fieldset%from_array(mgbf%mpl,fld_mga)

end subroutine mgbf_apply_filter

!----------------------------------------------------------------------
! Subroutine: mgbf_dealloc
! Purpose: release memory (full)
!----------------------------------------------------------------------
subroutine mgbf_dealloc(mgbf)

implicit none

! Passed variables
class(mgbf_type),intent(inout) :: mgbf ! MGBF

! Release memory
call mgbf%cs2gg%dealloc
call mgbf%gg2cs%dealloc
call mgbf%mpl%final

end subroutine mgbf_dealloc

subroutine mgbf_bump_init(mgbf,f_comm,afunctionspace)

implicit none

type(mgbf_type), intent(inout) :: mgbf
type(fckit_mpi_comm),intent(in) :: f_comm                     ! FCKIT MPI communicator wrapper
type(atlas_functionspace),intent(in) :: afunctionspace ! ATLAS functionspace

character(len=1024),parameter :: subr = 'mgbf_bump_init'
real(kind=kind_real), pointer :: real_ptr(:,:)
type(atlas_fieldset) :: afieldset
type(atlas_field) :: afield

type(geom_type) :: geom ! dummy to copy latlon

real(kind_real),allocatable :: lons_f(:),lats_f(:)
real(kind_real),allocatable :: lons_m(:),lats_m(:)
real(kind_real),allocatable :: rlons(:),rlats(:)
real(kind_real),allocatable :: slat(:),wlat(:)
real(kind_real) :: dlat,dlon
integer :: IMAX,JMAX
integer :: i,j,ij,ii,jj

! analysis to filter
call mgbf%cs2gg%nam%init(f_comm%size())

mgbf%cs2gg%nam%prefix = trim(subr)
mgbf%cs2gg%nam%default_seed = .true.
mgbf%cs2gg%nam%new_obsop = .true.
mgbf%cs2gg%nam%write_obsop = .false.
mgbf%cs2gg%nam%verbosity = "none"
mgbf%cs2gg%nam%nl = mgbf%grid%lm
mgbf%cs2gg%nam%nv = 1
mgbf%cs2gg%nam%variables(1) = "var"

afieldset = atlas_fieldset()

mgbf%ngrid = mgbf%grid%nm*mgbf%grid%mm
IMAX = mgbf%grid%nm*mgbf%grid%nxm
JMAX = mgbf%grid%mm*mgbf%grid%mym
allocate(rlons(1:IMAX))
allocate(rlats(1:JMAX))
allocate(SLAT(1:JMAX))
allocate(WLAT(1:JMAX))
allocate(lons_f(mgbf%ngrid))
allocate(lats_f(mgbf%ngrid))

! Create global longitude
dlon=360./IMAX
do i=1,IMAX
  rlons(i) = (i-1)*dlon
end do
! Create global latitude
mgbf%IDRT = 256 !   0: Equarlly-spaced grid with the pole (can't be used in bump) 
                !   4: Gaussian grid
                ! 256: Equarlly-spaced grid without the pole
call splat(mgbf%IDRT,JMAX,SLAT,WLAT)
rlats=0.0_kind_real
do j=1,JMAX/2
  rlats(JMAX+1-j) = asin(SLAT(j))*rad2deg
  rlats(j) = -asin(SLAT(j))*rad2deg
end do
deallocate(SLAT)
deallocate(WLAT)
ij=0
do j=1,mgbf%grid%mm
  jj=(mgbf%grid%mype/mgbf%grid%nxm)*mgbf%grid%mm+j
  if(jj<1) jj=1
  if(jj>mgbf%grid%mm*mgbf%grid%mym) jj=mgbf%grid%mm*mgbf%grid%mym
  do i=1,mgbf%grid%nm
    ii=mod(mgbf%grid%mype,mgbf%grid%nxm)*mgbf%grid%nm+i
    if(ii<1) ii=ii+mgbf%grid%nm*mgbf%grid%nxm
    if(ii>mgbf%grid%nm*mgbf%grid%nxm) ii=ii-mgbf%grid%nm*mgbf%grid%nxm
    ij = ij+1
    lons_f(ij)=rlons(ii)
    lats_f(ij)=rlats(jj)
  end do
end do

call mgbf%cs2gg%setup(f_comm,afunctionspace, &
& nobs=mgbf%ngrid,lonobs=lons_f,latobs=lats_f)
call mgbf%cs2gg%run_drivers()

! Copy model lat-lon
call geom%from_atlas(mgbf%mpl,afunctionspace)
allocate(lons_m(mgbf%cs2gg%geom%nmga))
allocate(lats_m(mgbf%cs2gg%geom%nmga))
do i=1,mgbf%cs2gg%geom%nmga
  lons_m(i)=geom%lon_mga(i)*rad2deg
  lats_m(i)=geom%lat_mga(i)*rad2deg
end do
call geom%dealloc
afield = atlas_field(name="lonlat", kind=atlas_real(kind_real),shape=(/2,mgbf%cs2gg%geom%nmga/))
call afield%data(real_ptr)
real_ptr(1,:) = lons_m(:)
real_ptr(2,:) = lats_m(:)
mgbf%afunctionspace_cs2gg = atlas_functionspace_pointcloud(afield)
call afield%final()

call afieldset%final()

! filter to analysis
call mgbf%gg2cs%nam%init(f_comm%size())

mgbf%gg2cs%nam%prefix = trim(subr)
mgbf%gg2cs%nam%default_seed = .true.
mgbf%gg2cs%nam%new_obsop = .true.
mgbf%gg2cs%nam%write_obsop = .false.
mgbf%gg2cs%nam%verbosity = "none"
mgbf%gg2cs%nam%nl = mgbf%nam%nl
mgbf%gg2cs%nam%nv = 1
mgbf%gg2cs%nam%variables(1) = "var"

afield = atlas_field(name="lonlat", kind=atlas_real(kind_real),shape=(/2,mgbf%ngrid/))
call afield%data(real_ptr)
real_ptr(1,:) = lons_f(:)
real_ptr(2,:) = lats_f(:)
mgbf%afunctionspace_gg2cs = atlas_functionspace_pointcloud(afield)
call afield%final()

call mgbf%gg2cs%setup(                  &
&     f_comm,mgbf%afunctionspace_gg2cs, &
&     nobs=mgbf%cs2gg%geom%nmga,        &
&     lonobs=lons_m,latobs=lats_m)
call mgbf%gg2cs%run_drivers()

call mgbf%gg2cs%partial_dealloc()

deallocate(lons_m)
deallocate(lats_m)
deallocate(lons_f)
deallocate(lats_f)

end subroutine mgbf_bump_init 

!----------------------------------------------------------------------
! Subroutine: fld_interp_cs2gg
! Purpose: fields interpolation
!----------------------------------------------------------------------
subroutine fld_interp_cs2gg(mgbf,nz,fld_mga,WA)

implicit none
type(mgbf_type),intent(inout) :: mgbf
integer,intent(in) :: nz
real(kind_real),intent(in) :: fld_mga(mgbf%cs2gg%geom%nmga,mgbf%cs2gg%geom%nl0)
real(kind_real),intent(out) :: WA(mgbf%grid%nm,mgbf%grid%mm,nz)
real(kind_real) :: fld_ngrid(mgbf%ngrid,mgbf%cs2gg%geom%nl0)

integer :: i,j,ij 

fld_ngrid=0.0
WA=0.0
call bump_apply(mgbf%cs2gg,mgbf%afunctionspace_cs2gg, &
&               mgbf%cs2gg%geom%nl0,                &
&               fld_mga(:,:),                     &
&               mgbf%ngrid,fld_ngrid)
ij=0
do j=1,mgbf%grid%mm
  do i=1,mgbf%grid%nm
    ij=ij+1
    WA(i,j,1:nz)=fld_ngrid(ij,1:nz)
  end do
end do

end subroutine fld_interp_cs2gg

!----------------------------------------------------------------------
! Subroutine: fld_interp_gg2cs
! Purpose: fields interpolation
!----------------------------------------------------------------------
subroutine fld_interp_gg2cs(mgbf,nz,WA,fld_mga)

implicit none
type(mgbf_type),intent(inout) :: mgbf
integer,intent(in) :: nz
real(kind_real),intent(in) :: WA(mgbf%grid%nm,mgbf%grid%mm,nz)
real(kind_real),intent(out) :: fld_mga(mgbf%cs2gg%geom%nmga,mgbf%cs2gg%geom%nl0)
real(kind_real) :: fld_ngrid(mgbf%ngrid,mgbf%cs2gg%geom%nl0)

integer :: i,j,ij

fld_ngrid=0.0
fld_mga=0.0
ij=0
do j=1,mgbf%grid%mm
  do i=1,mgbf%grid%nm
    ij=ij+1
    fld_ngrid(ij,1:nz)=WA(i,j,1:nz)
  end do
end do
call bump_apply(mgbf%gg2cs,mgbf%afunctionspace_gg2cs, &
&               mgbf%cs2gg%geom%nl0,                &
&               fld_ngrid,                        &
&               mgbf%cs2gg%geom%nmga,fld_mga(:,:))

end subroutine fld_interp_gg2cs

!----------------------------------------------------------------------
! Subroutine: bump_apply
! Purpose: bump_apply_obsop wrapper
!----------------------------------------------------------------------
subroutine bump_apply(bump,afunctionspace,nz,fld_in,ngrid_ou,fld_ou)

implicit none

type(bump_type),intent(inout) :: bump
type(atlas_functionspace),intent(in) :: afunctionspace ! ATLAS functionspace
integer,intent(in) :: nz,ngrid_ou
real(kind_real),intent(in ) :: fld_in(bump%geom%nmga,nz)
real(kind_real),intent(out) :: fld_ou(ngrid_ou,nz)

real(kind_real), pointer :: real_ptr_2(:,:)
type(atlas_field) :: afield
type(atlas_fieldset) :: afieldset

integer :: k

bump%geom%nl0 = nz

afieldset = atlas_fieldset()
afield = afunctionspace%create_field(name='var', kind=atlas_real(kind_real),levels=nz)
call afieldset%add(afield)

call afield%data(real_ptr_2)
do k=1,nz
  real_ptr_2(k,:) = fld_in(:,k)
end do

call bump%apply_obsop(afieldset,fld_ou)

call afieldset%final()
call afield%final()

end subroutine bump_apply

end module type_mgbf
