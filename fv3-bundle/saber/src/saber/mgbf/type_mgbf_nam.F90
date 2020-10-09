! (C) Copyright 2018-2019 UCAR
!
! This software is licensed under the terms of the Apache Licence
! Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
module type_mgbf_nam

use fckit_configuration_module, only: fckit_configuration,fckit_yamlconfiguration
use fckit_pathname_module, only : fckit_pathname
use iso_c_binding
use tools_const, only: pi,req,deg2rad,rad2deg
use tools_kinds,only: kind_real
use type_mpl, only: mpl_type
use type_mgbf_grid, only: mgbf_grid_type

implicit none

integer,parameter :: nvmax = 20      ! Maximum number of variables
integer,parameter :: nlmax = 200     ! Maximum number of levels
integer,parameter :: nc3max = 1000   ! Maximum number of classes
integer,parameter :: nscalesmax = 5  ! Maximum number of variables
integer,parameter :: ndirmax = 300   ! Maximum number of diracs
integer,parameter :: nldwvmax = 99   ! Maximum number of local diagnostic profiles
integer,parameter :: nprociomax = 20 ! Maximum number of I/O tasks

type mgbf_nam_type
  ! general_param
  integer :: nprocio                                   ! Number of IO processors
  ! model_param
  integer :: nl                                        ! Number of levels
  integer :: levs(nlmax)                               ! Levels
  character(len=1024) :: lev2d                         ! Level for 2D variables ('first' or 'last')
  logical :: logpres                                   ! Use pressure logarithm as vertical coordinate (model level if .false.)
  integer :: nv                                        ! Number of variables
  character(len=1024),dimension(nvmax) :: variables      ! Variables names
  logical :: nomask                                    ! Do not use geometry mask
contains
  procedure :: init => mgbf_nam_init
  procedure :: read => mgbf_nam_read
  procedure :: read_yaml => mgbf_nam_read_yaml
  procedure :: bcast => mgbf_nam_bcast
  procedure :: from_conf => mgbf_nam_from_conf
  procedure :: check => mgbf_nam_check
  procedure :: write => mgbf_nam_write
end type mgbf_nam_type

private
public :: nvmax,nlmax,nc3max,nscalesmax,ndirmax,nldwvmax
public :: mgbf_nam_type

contains

!----------------------------------------------------------------------
! Subroutine: mgbf_nam_init
! Purpose: intialize
!----------------------------------------------------------------------
subroutine mgbf_nam_init(nam,nproc)

implicit none

! Passed variable
class(mgbf_nam_type),intent(inout) :: nam ! Namelist
integer,intent(in) :: nproc        ! Number of MPI task

! Local variable
integer :: il,iv,ildwv

! general_param default
nam%nprocio = min(nproc,nprociomax)

! model_param default
nam%nl = 0
do il=1,nlmax
   nam%levs(il) = il
end do
nam%lev2d = 'first'
nam%logpres = .false.
nam%nv = 0
do iv=1,nvmax
   nam%variables(iv) = ''
end do
nam%nomask = .false.

end subroutine mgbf_nam_init

!----------------------------------------------------------------------
! Subroutine: mgbf_nam_read
! Purpose: read
!----------------------------------------------------------------------
subroutine mgbf_nam_read(nam,mpl,grid,namelname)

implicit none

! Passed variable
class(mgbf_nam_type),intent(inout) :: nam     ! Namelist
type(mpl_type),intent(inout) :: mpl      ! MPI data
type(mgbf_grid_type),intent(inout) :: grid
character(len=*),intent(in) :: namelname ! Namelist name

! Local variables
integer :: iv
character(len=1024),parameter :: subr = 'nam_read'
integer :: lunit,i

! Namelist variables
integer :: nl
integer :: levs(nlmax)
character(len=1024) :: lev2d
integer :: nv
character(len=1024),dimension(nvmax) :: variables

integer :: nA_max0
integer :: mA_max0
integer :: nm0
integer :: mm0
integer :: im00
integer :: jm00
integer :: nxm
integer :: mym
integer :: lm

namelist/model_param/nl,      &
                   & levs,    &
                   & lev2d,   &
                   & nv,      &
                   & variables

namelist/grid_param/  nA_max0, &
                    & nA_max0, &
                    & nm0,     &
                    & mm0,     &
                    & im00,    &
                    & jm00,    &
                    & nxm,     &
                    & mym,     &
                    & lm


if(mpl%main) then
  ! Open namelist
  call mpl%newunit(lunit)
  open(unit=lunit,file=trim(namelname),status='old',action='read')

  ! model_param default
  nl = 0
  levs = 0
  lev2d = 'first'
  nv = 0
  do iv=1,nvmax
    variables(iv) = ''
  end do
  ! model_param
  read(lunit,nml=model_param)
  if (nl>nlmax) call mpl%abort(subr,'nl is too large')
  if (nv>nvmax) call mpl%abort(subr,'nv is too large')
  nam%nl = nl
  if (nl>0) nam%levs(1:nl) = levs(1:nl)
  nam%lev2d = lev2d
  nam%nv = nv
  if (nv>0) nam%variables(1:nv) = variables(1:nv)

  ! grid_param default
  nA_max0 = 1800
  mA_max0 = 1060

  nm0 = 1804
  mm0 = 1072
  im00=1760
  jm00=960

  nxm = 8
  mym = 6

  lm = 64

  ! grid_param
  read(lunit,nml=grid_param)
  grid%nA_max0 = nA_max0
  grid%mA_max0 = mA_max0
  grid%nm0 = nm0
  grid%mm0 = mm0
  grid%im00= im00
  grid%jm00= jm00
  grid%nxm = nxm
  grid%mym = mym
  grid%lm = lm

  ! Close namelist
  close(unit=lunit)
end if

end subroutine mgbf_nam_read

!----------------------------------------------------------------------
! Subroutine: mgbf_nam_read_yaml
! Purpose: read YAML file
!----------------------------------------------------------------------
subroutine mgbf_nam_read_yaml(nam,mpl,grid,yamlname)

implicit none

! Passed variable
class(mgbf_nam_type),intent(inout) :: nam       ! Namelist
type(mpl_type),intent(inout) :: mpl        ! MPI data
type(mgbf_grid_type),intent(inout) :: grid
character(len=*),intent(inout) :: yamlname ! YAML name

! Local variables
type(fckit_configuration) :: conf

if (mpl%main) then
   ! Set fckit configuration from yamlname
   conf = fckit_yamlconfiguration(fckit_pathname(trim(yamlname)))

   ! Convert fckit configuration to namelist
   call nam%from_conf(grid,conf)
end if

end subroutine mgbf_nam_read_yaml

!----------------------------------------------------------------------
! Subroutine: mgbf_nam_bcast
! Purpose: broadcast
!----------------------------------------------------------------------
subroutine mgbf_nam_bcast(nam,mpl,grid)

implicit none

! Passed variable
class(mgbf_nam_type),intent(inout) :: nam ! Namelist
type(mpl_type),intent(inout) :: mpl  ! MPI data
type(mgbf_grid_type),intent(inout) :: grid

! model_param
call mpl%f_comm%broadcast(nam%nl,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%levs,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%lev2d,mpl%rootproc-1)
call mpl%f_comm%broadcast(nam%nv,mpl%rootproc-1)
call mpl%broadcast(nam%variables,mpl%rootproc-1)

! grid_param
call mpl%f_comm%broadcast(grid%nA_max0,mpl%rootproc-1)
call mpl%f_comm%broadcast(grid%mA_max0,mpl%rootproc-1)
call mpl%f_comm%broadcast(grid%nm0,mpl%rootproc-1)
call mpl%f_comm%broadcast(grid%mm0,mpl%rootproc-1)
call mpl%f_comm%broadcast(grid%im00,mpl%rootproc-1)
call mpl%f_comm%broadcast(grid%jm00,mpl%rootproc-1)
call mpl%f_comm%broadcast(grid%nxm,mpl%rootproc-1)
call mpl%f_comm%broadcast(grid%mym,mpl%rootproc-1)
call mpl%f_comm%broadcast(grid%lm,mpl%rootproc-1)

end subroutine mgbf_nam_bcast

!----------------------------------------------------------------------
! Subroutine: mgbf_nam_from_conf
! Purpose: intialize from configuration
!----------------------------------------------------------------------
subroutine mgbf_nam_from_conf(nam,grid,conf)

implicit none

! Passed variable
class(mgbf_nam_type),intent(inout) :: nam         ! Namelist
type(mgbf_grid_type),intent(inout) :: grid
type(fckit_configuration),intent(in) :: conf ! Configuration

! Local variables
integer,allocatable :: integer_array(:)
real(kind_real),allocatable :: real_array(:)
logical,allocatable :: logical_array(:)
character(len=:),allocatable :: str
character(len=:),allocatable :: str_array(:)

! model_param
if (conf%has("nl")) call conf%get_or_die("nl",nam%nl)
if (conf%has("levs")) then
   call conf%get_or_die("levs",integer_array)
   nam%levs(1:nam%nl) = integer_array(1:nam%nl)
end if
if (conf%has("lev2d")) then
   call conf%get_or_die("lev2d",str)
   nam%lev2d = str
end if
if (conf%has("logpres")) call conf%get_or_die("logpres",nam%logpres)
if (conf%has("nv")) call conf%get_or_die("nv",nam%nv)
if (conf%has("variables")) then
   call conf%get_or_die("variables",str_array)
   nam%variables(1:nam%nv) = str_array(1:nam%nv)
end if

! mgbf_grid_param
if (conf%has("nA_max0")) call conf%get_or_die("nA_max0",grid%nA_max0)
if (conf%has("mA_max0")) call conf%get_or_die("mA_max0",grid%mA_max0)
if (conf%has("nm0"))     call conf%get_or_die("nm0",grid%nm0)
if (conf%has("mm0"))     call conf%get_or_die("mm0",grid%mm0)
if (conf%has("im00"))    call conf%get_or_die("im00",grid%im00)
if (conf%has("jm00"))    call conf%get_or_die("jm00",grid%jm00)
if (conf%has("nxm"))     call conf%get_or_die("nxm",grid%nxm)
if (conf%has("mym"))     call conf%get_or_die("mym",grid%mym)
if (conf%has("lm"))      call conf%get_or_die("lm",grid%lm)

end subroutine mgbf_nam_from_conf

!----------------------------------------------------------------------
! Subroutine: mgbf_nam_check
! Purpose: check namelist parameters
!----------------------------------------------------------------------
subroutine mgbf_nam_check(nam,mpl,grid)

implicit none

! Passed variable
class(mgbf_nam_type),intent(inout) :: nam ! Namelist
type(mpl_type),intent(inout) :: mpl  ! MPI data
type(mgbf_grid_type),intent(inout) :: grid

! Local variables
integer :: iv,its,il,idir,ildwv
character(len=2) :: ivchar,itschar,ildwvchar
character(len=1024),parameter :: subr = 'nam_check'

end subroutine mgbf_nam_check

!----------------------------------------------------------------------
! Subroutine: mgbf_nam_write
! Purpose: write namelist parameters into a log file or into a NetCDF file
!----------------------------------------------------------------------
subroutine mgbf_nam_write(nam,mpl,grid,ncid)

implicit none

! Passed variable
class(mgbf_nam_type),intent(in) :: nam   ! Namelist
type(mpl_type),intent(inout) :: mpl ! MPI data
type(mgbf_grid_type),intent(inout) :: grid
integer,intent(in),optional :: ncid ! NetCDF file ID

! Local variables
integer :: lncid
real(kind_real),allocatable :: londir(:),latdir(:),lon_ldwv(:),lat_ldwv(:)

! Set ncid
lncid = mpl%msv%vali
if (present(ncid)) lncid = ncid

! model_param
if (mpl%msv%is(lncid)) then
  write(mpl%info,'(a7,a)') '','Model parameters'
  call mpl%flush
end if
call mpl%write(lncid,'mgbf_nam','nl',nam%nl)
call mpl%write(lncid,'mgbf_nam','levs',nam%nl,nam%levs(1:nam%nl))
call mpl%write(lncid,'mgbf_nam','lev2d',nam%lev2d)
call mpl%write(lncid,'mgbf_nam','nv',nam%nv)
call mpl%write(lncid,'mgbf_nam','variables',nam%nv,nam%variables(1:nam%nv))

! grid_param
if (mpl%msv%is(lncid)) then
  write(mpl%info,'(a7,a)') '','Grid parameters'
  call mpl%flush
end if
call mpl%write(lncid,'mgbf_nam','nA_max0',grid%nA_max0)
call mpl%write(lncid,'mgbf_nam','mA_max0',grid%mA_max0)
call mpl%write(lncid,'mgbf_nam','nm0',grid%nm0)
call mpl%write(lncid,'mgbf_nam','mm0',grid%mm0)
call mpl%write(lncid,'mgbf_nam','im00',grid%im00)
call mpl%write(lncid,'mgbf_nam','jm00',grid%jm00)
call mpl%write(lncid,'mgbf_nam','nxm',grid%nxm)
call mpl%write(lncid,'mgbf_nam','mym',grid%mym)
call mpl%write(lncid,'mgbf_nam','lm',grid%lm)

end subroutine mgbf_nam_write

end module type_mgbf_nam
