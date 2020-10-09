! (C) Copyright 2018-2019 UCAR
!
! This software is licensed under the terms of the Apache Licence
! Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
module type_mgbf_grid

use tools_kinds, only: kind_int,kind_real
use type_mpl, only: mpl_type

implicit none

! subdomain derived type
type mgbf_grid_type
  integer :: npe
  integer :: mype
  integer :: nc3d
  integer :: nc2d
!
  integer :: gm   !
  integer :: nA_max0
  integer :: mA_max0
  integer :: nm0  !
  integer :: mm0  !
  integer :: nxm  !
  integer :: mym  !
  integer :: nm   !
  integer :: mm   !
  integer :: im00
  integer :: jm00
  integer :: im
  integer :: jm
  integer :: ib
  integer :: jb
  integer :: nb
  integer :: mb
  integer:: hx,hy,hz
  integer:: p
  integer:: pasp0
  integer, allocatable, dimension(:):: maxpe_fgen
  integer, allocatable, dimension(:):: ixm,jym,nxy
  integer, allocatable, dimension(:):: im0,jm0
  integer, allocatable, dimension(:):: Fimax,Fjmax
  integer, allocatable, dimension(:):: FimaxL,FjmaxL
  integer:: maxpe_filt
  integer:: imL,jmL
  integer:: imH,jmH
  integer:: lm            ! number of vertical layers
  integer:: km            ! number of 3d variables
  integer:: km3           ! number of 3d variables
  integer:: km2           ! number of 2d variables
  integer:: lm05          ! half of vertical levels
  real(kind_real):: lengthx,lengthy,x0,y0
  real(kind_real):: dxf,dyf,dxa,dya
  integer:: npadx         ! x padding on analysis grid
  integer:: mpady         ! y padding on analysis grid
  integer:: ipadx         ! x padding on filter decomposition
  integer:: jpady         ! y padding on filter deocmposition
  integer:: nx
  integer:: my
  integer:: my_hgen
  integer:: mype_hgen
  logical:: l_hgen

  logical,dimension(2):: Flwest,Fleast,Flnorth,Flsouth
  logical,dimension(2):: Flcorner_sw,Flcorner_nw,Flcorner_se,Flcorner_ne
  integer(kind_int),dimension(2):: Fitarg_n,Fitarg_e,Fitarg_s,Fitarg_w
  integer(kind_int),dimension(2):: Fitarg_sw,Fitarg_se,Fitarg_ne,Fitarg_nw
  logical,dimension(2):: Flsendup_sw,Flsendup_se,Flsendup_nw,Flsendup_ne
  integer(kind_int),dimension(2):: Fitarg_up
  integer(kind_int):: itargdn_sw,itargdn_se,itargdn_ne,itargdn_nw
  integer(kind_int):: itarg_wA,itarg_eA,itarg_sA,itarg_nA
  logical:: lwestA,leastA,lsouthA,lnorthA
  integer(kind_int) ix,jy
  integer(kind_int),dimension(2):: mype_filt

contains
  procedure :: dealloc => mgbf_grid_dealloc
  procedure :: init => mgbf_grid_init
  procedure :: setup => mgbf_grid_setup
end type mgbf_grid_type

private
public :: mgbf_grid_type

contains

!----------------------------------------------------------------------
! Subroutine: mgbf_grid_dealloc
! Purpose: dealloc
!----------------------------------------------------------------------
subroutine mgbf_grid_dealloc(grid)

implicit none

! Passed variables
class(mgbf_grid_type),intent(inout) :: grid               ! GRID

! Local variables

if(allocated(grid%maxpe_fgen)) deallocate(grid%maxpe_fgen)
if(allocated(grid%ixm))        deallocate(grid%ixm)
if(allocated(grid%jym))        deallocate(grid%jym)
if(allocated(grid%nxy))        deallocate(grid%nxy)
if(allocated(grid%im0))        deallocate(grid%im0)
if(allocated(grid%jm0))        deallocate(grid%jm0)
if(allocated(grid%Fimax))      deallocate(grid%Fimax)
if(allocated(grid%Fjmax))      deallocate(grid%Fjmax)
if(allocated(grid%FimaxL))     deallocate(grid%FimaxL)
if(allocated(grid%FjmaxL))     deallocate(grid%FjmaxL)

end subroutine mgbf_grid_dealloc

!----------------------------------------------------------------------
! Subroutine: mgbf_grid_init
! Purpose: init
!----------------------------------------------------------------------
subroutine mgbf_grid_init(grid,nproc,myproc)

implicit none

! Passed variables
class(mgbf_grid_type),intent(inout) :: grid               ! GRID
integer,intent(in) :: nproc ! Number of MPI task
integer,intent(in) :: myproc ! My MPI task

grid%npe = nproc
grid%mype = myproc

grid%nA_max0 = 1800
grid%mA_max0 = 1060

grid%nm0 = 1804
grid%mm0 = 1072
grid%im00=1760
grid%jm00=960

grid%lm = 64

grid%nxm = 8
grid%mym = 6


end subroutine mgbf_grid_init

!----------------------------------------------------------------------
! Subroutine: mgbf_grid_setup
! Purpose: setup
!----------------------------------------------------------------------
subroutine mgbf_grid_setup(grid,nv,variables)

implicit none

! Passed variables
class(mgbf_grid_type),intent(inout) :: grid               ! GRID
integer,intent(in) :: nv
character(len=*),intent(in) :: variables(nv)

! Local variables
integer :: n,nn,k,g

grid%nc3d=0
grid%nc2d=0
do n=1,nv
  select case (trim(variables(n)))
  case("ua","va","sf","psi","vp","chi","t","tv",  &
 &     "q","sphum","rh","cw","liq_wat","ice_wat", &
 &     "oz","o3","o3mr")
    grid%nc3d = grid%nc3d + 1
  case("ps","sst")
    grid%nc2d = grid%nc2d + 1
  case default
    grid%nc3d = grid%nc3d + 1
  end select
end do

call init_mg_parameter(grid)

call init_mg_MPI(grid)

call init_mg_domain(grid)
call init_topology_2d(grid)

end subroutine mgbf_grid_setup

subroutine init_mg_parameter(grid)
implicit none
type(mgbf_grid_type),intent(inout) :: grid
integer :: g

call def_maxgen(grid)
if(grid%gm>4) then
  grid%gm=4
end if

grid%nm = grid%nm0/grid%nxm
grid%mm = grid%mm0/grid%mym

grid%im = grid%im00/grid%nxm
grid%jm = grid%jm00/grid%mym

if(grid%im*grid%nxm /= grid%im00 ) then
  write(6,*) 'im,nxm,im00=',grid%im,grid%nxm,grid%im00
  stop 'im00 is not divisible by nxm'
endif

if(grid%jm*grid%mym /= grid%jm00 ) then
  write(6,*) 'jm,mym,jm00=',grid%jm,grid%mym,grid%jm00
  stop 'jm00 is not divisible by mym'
endif

!
! Set number of processors in higher generations
!

allocate(grid%ixm(grid%gm))
allocate(grid%jym(grid%gm))
allocate(grid%nxy(grid%gm))
allocate(grid%maxpe_fgen(0:grid%gm))
allocate(grid%im0(grid%gm))
allocate(grid%jm0(grid%gm))
allocate(grid%Fimax(grid%gm))
allocate(grid%Fjmax(grid%gm))
allocate(grid%FimaxL(grid%gm))
allocate(grid%FjmaxL(grid%gm))

call def_ngens(grid%ixm,grid%gm,grid%nxm)
call def_ngens(grid%jym,grid%gm,grid%mym)

do g=1,grid%gm
  grid%nxy(g)=grid%ixm(g)*grid%jym(g)
end do

grid%maxpe_fgen(0)= 0
do g=1,grid%gm
  grid%maxpe_fgen(g)=grid%maxpe_fgen(g-1)+grid%nxy(g)
end do

grid%maxpe_filt=grid%maxpe_fgen(grid%gm)

grid%im0(1)=grid%im00
do g=2,grid%gm
  grid%im0(g)=(grid%im0(g-1)+1)/2
end do

grid%jm0(1)=grid%jm00
do g=2,grid%gm
  grid%jm0(g)=(grid%jm0(g-1)+1)/2
end do

do g=1,grid%gm
  grid%Fimax(g)=grid%im0(g)-grid%im*(grid%ixm(g)-1)
  grid%Fjmax(g)=grid%jm0(g)-grid%jm*(grid%jym(g)-1)
!TEST
!      write(15,*)'Fimax(',g,')=',Fimax(g)
!      write(15,*)'Fjmax(',g,')=',Fjmax(g)
!TEST
end do

!
!  Double check this - for now should be fine !!!!!
!
do g=1,grid%gm
  grid%FimaxL(g)=grid%Fimax(g)/2
  grid%FjmaxL(g)=grid%Fjmax(g)/2
end do

!***
!***     Vertical distribution
!***

! lm = 1
grid%lm05 = grid%lm/2

grid%km = 6
grid%km2= 4
grid%km3= grid%km

!grid%km = grid%nc3d
!grid%km2= grid%nc2d
!  km2= 1
!grid%km3= grid%km


!***
!*** Filter related parameters
!**
grid%lengthx = 6.      ! arbitrary chosen scale of the domain
grid%lengthy = 6.      ! arbitrary chosen scale of the domain

grid%x0 = -3.
grid%y0 = -3.

grid%ib=4
grid%jb=4

grid%dxa = grid%lengthx/grid%nm
grid%dxf = grid%lengthx/grid%jm
grid%nb = 2*grid%dxf/grid%dxa

grid%dya = grid%lengthy/grid%mm
grid%dyf = grid%lengthy/grid%jm
grid%mb = 2*grid%dyf/grid%dya

grid%imL=grid%im/2
grid%jmL=grid%jm/2

grid%imH=2*grid%im
grid%jmH=2*grid%jm

!  pasp0=1.
!  pasp0=9.
grid%pasp0 = 10.
grid%pasp0 = 18.


grid%hx=8
grid%hx=9
grid%hx=12
grid%hy=grid%hx
grid%hz=grid%hx+2

grid%p = 2                !  Exponent of Beta function

end subroutine init_mg_parameter

subroutine init_mg_domain(grid)
implicit none

class(mgbf_grid_type),intent(inout) :: grid

call init_domain(grid)

end subroutine init_mg_domain

subroutine def_maxgen(grid)

implicit none

type(mgbf_grid_type),intent(inout) :: grid

integer:: npx,npy,gx,gy

npx = grid%nxm;  gx=1
Do
  npx = (npx + 1)/2
  gx = gx + 1
  if(npx == 1) exit
end do

npy = grid%mym;  gy=1
Do
  npy = (npy + 1)/2
  gy = gy + 1
  if(npy == 1) exit
end do

grid%gm = Min(gx,gy)

end subroutine def_maxgen

subroutine def_ngens(nsm,gm,nsm0)
implicit none
integer, intent(in):: gm,nsm0
integer, dimension(gm), intent(out):: nsm
integer:: g

nsm(1)=nsm0
do g=2,gm
  nsm(g) = (nsm(g-1) + 1)/2
end do

end subroutine def_ngens

subroutine init_mg_MPI(grid)
implicit none

type(mgbf_grid_type),intent(inout) :: grid
integer :: g

grid%nx = mod(grid%mype,grid%nxm)+1
grid%my = (grid%mype/grid%nxm)+1

grid%mype_hgen=-1
grid%my_hgen=-1

if( grid%mype < grid%maxpe_filt-grid%nxy(1)) then
  grid%mype_hgen=grid%mype+grid%nxy(1)
endif
do g=1,grid%gm
  if(grid%maxpe_fgen(g-1)<= grid%mype_hgen .and. grid%mype_hgen< grid%maxpe_fgen(g)) then
    grid%my_hgen=g
  endif
enddo
grid%l_hgen = grid%mype_hgen >-1


end subroutine init_mg_MPI

subroutine init_domain(grid)
implicit none

class(mgbf_grid_type),intent(inout) :: grid

integer :: n,nstrd,i,j
logical :: F=.false., T=.true.

integer :: loc_pe,g

grid%Flwest(1)=grid%nx.eq.1
grid%Fleast(1)=grid%nx.eq.grid%nxm
grid%Flsouth(1)=grid%my.eq.1
grid%Flnorth(1)=grid%my.eq.grid%mym

if(grid%l_hgen) then

  loc_pe=grid%mype_hgen-grid%maxpe_fgen(grid%my_hgen-1)
  grid%jy=loc_pe/grid%ixm(grid%my_hgen)+1
  grid%ix=mod(loc_pe,grid%ixm(grid%my_hgen))+1

  grid%Flwest(2)=grid%ix.eq.1
  grid%Fleast(2)=grid%ix.eq.grid%ixm(grid%my_hgen)
  grid%Flsouth(2)=grid%jy.eq.1
  grid%Flnorth(2)=grid%jy.eq.grid%jym(grid%my_hgen)

else

  grid%jy = -1
  grid%ix = -1

  grid%Flwest(2)=F
  grid%Fleast(2)=F
  grid%Flsouth(2)=F
  grid%Flnorth(2)=F

endif

do g=1,2

  grid%Flcorner_sw(g)=F
  grid%Flcorner_se(g)=F
  grid%Flcorner_nw(g)=F
  grid%Flcorner_ne(g)=F

  if(grid%Flsouth(g).and.grid%Flwest(g))then
    grid%Flcorner_sw(g)=T
  endif
  if(grid%Flsouth(g).and.grid%Fleast(g))then
    grid%Flcorner_se(g)=T
  endif
  if(grid%Flnorth(g).and.grid%Flwest(g))then
    grid%Flcorner_nw(g)=T
  endif
  if(grid%Flnorth(g).and.grid%Fleast(g))then
    grid%Flcorner_ne(g)=T
  endif

enddo

grid%mype_filt(1)=grid%mype
grid%mype_filt(2)=grid%mype_hgen

!
! Communication params for analysis grid
!
if(grid%nx==1) then
  grid%itarg_wA=-1
else
  grid%itarg_wA=grid%mype-1
endif

if(grid%nx==grid%nxm) then
  grid%itarg_eA=-1
else
  grid%itarg_eA=grid%mype+1
endif

if(grid%my==1) then
  grid%itarg_sA=-1
else
  grid%itarg_sA=grid%mype-grid%nxm
endif

if(grid%my==grid%mym) then
  grid%itarg_nA=-1
else
  grid%itarg_nA=grid%mype+grid%nxm
endif

grid%lwestA=grid%nx.eq.1
grid%leastA=grid%nx.eq.grid%nxm
grid%lsouthA=grid%my.eq.1
grid%lnorthA=grid%my.eq.grid%mym

end subroutine init_domain

subroutine init_topology_2d(grid)

implicit none

type(mgbf_grid_type),intent(inout) :: grid

logical:: F=.false., T=.true.

integer(kind_int) :: mx2,my2,ix_up,jy_up,ix_dn,jy_dn
integer(kind_int) :: g,naux,nx_up,my_up

do g = 1,2
!***
!*** Send WEST
!***
  if(grid%Flwest(g)) then
    grid%Fitarg_w(g) = -1
  else
    if(g==1.or.grid%l_hgen) then
      grid%Fitarg_w(g) = grid%mype_filt(g)-1
    else
      grid%Fitarg_w(g) = -1
    end if
  end if
!***
!*** Send EAST
!***
  if(grid%Fleast(g)) then
    grid%Fitarg_e(g) = -1
  else
    if(g==1.or.grid%l_hgen) then
      grid%Fitarg_e(g) = grid%mype_filt(g)+1
    else
      grid%Fitarg_e(g) = -1
    end if
  end if

!***
!*** Send SOUTH
!***

  if(grid%Flsouth(g)) then
    grid%Fitarg_s(g)=-1
  else
    select case(g)
      case(1)
        naux = grid%nxm
      case(2)
        if(grid%l_hgen) then
          naux = grid%ixm(grid%my_hgen)
        end if
    end select
    if(g==1.or.grid%l_hgen) then
      grid%Fitarg_s(g)=grid%mype_filt(g)-naux
    else
      grid%Fitarg_s(g)=-1
    end if
  end if

!***
!*** Send NORTH
!***
  if(grid%Flnorth(g)) then
    grid%Fitarg_n(g)=-1
  else
    select case(g)
      case(1)
        naux = grid%nxm
      case(2)
        if(grid%l_hgen) then
          naux = grid%ixm(grid%my_hgen)
        end if
    end select
    if(g==1.or.grid%l_hgen) then
      grid%Fitarg_n(g)=grid%mype_filt(g)+naux
    else
      grid%Fitarg_n(g)=-1
    end if
  end if

!***
!*** Send SOUTH-WEST
!***

  if(grid%Flsouth(g).and.grid%Flwest(g)) then
    grid%Fitarg_sw(g)=-1
  else if(grid%Flsouth(g)) then
    grid%Fitarg_sw(g)=grid%Fitarg_w(g)
  else if(grid%Flwest(g)) then
    grid%Fitarg_sw(g)=grid%Fitarg_s(g)
  else
    grid%Fitarg_sw(g)=grid%Fitarg_s(g)-1
  end if
  if(g>1 .and. .not.grid%l_hgen) then
    grid%Fitarg_sw(g)=-1
  end if

!***
!*** Send SOUTH-EAST
!***

  if(grid%Flsouth(g).and.grid%Fleast(g)) then
    grid%Fitarg_se(g)=-1
  else if(grid%Flsouth(g)) then
    grid%Fitarg_se(g)=grid%Fitarg_e(g)
  else if(grid%Fleast(g)) then
    grid%Fitarg_se(g)=grid%Fitarg_s(g)
  else
    grid%Fitarg_se(g)=grid%Fitarg_s(g)+1
  end if
  if(g>1 .and. .not.grid%l_hgen) then
    grid%Fitarg_se(g)=-1
  end if

!***
!*** Send NORTH-WEST
!***
  if(grid%Flnorth(g).and.grid%Flwest(g)) then
    grid%Fitarg_nw(g)=-1
  else if(grid%Flnorth(g)) then
    grid%Fitarg_nw(g)=grid%Fitarg_w(g)
  else if(grid%Flwest(g)) then
    grid%Fitarg_nw(g)=grid%Fitarg_n(g)
  else
    grid%Fitarg_nw(g)=grid%Fitarg_n(g)-1
  end if
  if(g>1 .and. .not.grid%l_hgen) then
    grid%Fitarg_nw(g)=-1
  end if


!***
!*** Send NORTH-EAST
!***

  if(grid%Flnorth(g).and.grid%Fleast(g)) then
    grid%Fitarg_ne(g)=-1
  else if(grid%Flnorth(g)) then
    grid%Fitarg_ne(g)=grid%Fitarg_e(g)
  else if(grid%Fleast(g)) then
    grid%Fitarg_ne(g)=grid%Fitarg_n(g)
  else
    grid%Fitarg_ne(g)=grid%Fitarg_n(g)+1
  end if
  if(g>1 .and. .not.grid%l_hgen) then
    grid%Fitarg_ne(g)=-1
  end if


end do

!-----------------------------------------------------------------------
!
! Upsending flags
!

mx2=mod(grid%nx,2)
my2=mod(grid%my,2)

if(mx2==1.and.my2==1) then
  grid%Flsendup_sw(1)=T
else if(mx2==0.and.my2==1) then
  grid%Flsendup_se(1)=T
else if(mx2==1.and.my2==0) then
  grid%Flsendup_nw(1)=T
else
  grid%Flsendup_ne(1)=T
end if

nx_up=(grid%nx-1)/2   !+1
my_up=(grid%my-1)/2   !+1


grid%Fitarg_up(1)=grid%maxpe_fgen(1)+my_up*grid%ixm(2)+nx_up

if(grid%l_hgen.and.grid%my_hgen < grid%gm) then

  mx2=mod(grid%ix,2)
  my2=mod(grid%jy,2)

  if(mx2==1.and.my2==1) then
    grid%Flsendup_sw(2)=T
  else if(mx2==0.and.my2==1) then
    grid%Flsendup_se(2)=T
  else if(mx2==1.and.my2==0) then
    grid%Flsendup_nw(2)=T
  else
    grid%Flsendup_ne(2)=T
  end if

  ix_up=(grid%ix-1)/2   !+1
  jy_up=(grid%jy-1)/2   !+1

  grid%Fitarg_up(2)=grid%maxpe_fgen(grid%my_hgen)+jy_up*grid%ixm(grid%my_hgen+1)+ix_up

else
  grid%Flsendup_sw(2)=F
  grid%Flsendup_se(2)=F
  grid%Flsendup_nw(2)=F
  grid%Flsendup_ne(2)=F

  grid%Fitarg_up(2)=-1

end if

if(grid%my_hgen > 1) then

  ix_dn = 2*grid%ix-1
  jy_dn = 2*grid%jy-1

  grid%itargdn_sw=grid%maxpe_fgen(grid%my_hgen-2)+(jy_dn-1)*grid%ixm(grid%my_hgen-1)+(ix_dn-1)
  grid%itargdn_nw=grid%itargdn_sw+grid%ixm(grid%my_hgen-1)
  grid%itargdn_se=grid%itargdn_sw+1
  grid%itargdn_ne=grid%itargdn_nw+1

  if(grid%Fimax(grid%my_hgen) <= grid%imL .and. grid%Fleast(2)) then
    grid%itargdn_se=-1
    grid%itargdn_ne=-1
  end if
  if(grid%Fjmax(grid%my_hgen) <= grid%jmL .and. grid%Flnorth(2)) then
    grid%itargdn_nw=-1
    grid%itargdn_ne=-1
  end if

else

  grid%itargdn_sw=-1
  grid%itargdn_se=-1
  grid%itargdn_nw=-1
  grid%itargdn_ne=-1

end if
call real_itarg(grid%Fitarg_w(2),grid%nxy(1))
call real_itarg(grid%Fitarg_e(2),grid%nxy(1))
call real_itarg(grid%Fitarg_s(2),grid%nxy(1))
call real_itarg(grid%Fitarg_n(2),grid%nxy(1))

call real_itarg(grid%Fitarg_sw(2),grid%nxy(1))
call real_itarg(grid%Fitarg_se(2),grid%nxy(1))
call real_itarg(grid%Fitarg_nw(2),grid%nxy(1))
call real_itarg(grid%Fitarg_ne(2),grid%nxy(1))

if(grid%itargdn_sw .ge. grid%maxpe_fgen(1)) call real_itarg(grid%itargdn_sw,grid%nxy(1))
if(grid%itargdn_se .ge. grid%maxpe_fgen(1)) call real_itarg(grid%itargdn_se,grid%nxy(1))
if(grid%itargdn_nw .ge. grid%maxpe_fgen(1)) call real_itarg(grid%itargdn_nw,grid%nxy(1))
if(grid%itargdn_ne .ge. grid%maxpe_fgen(1)) call real_itarg(grid%itargdn_ne,grid%nxy(1))

call real_itarg(grid%Fitarg_up(1),grid%nxy(1))
call real_itarg(grid%Fitarg_up(2),grid%nxy(1))

end subroutine init_topology_2d

subroutine real_itarg(itarg,nxy1)
implicit none
integer(kind_int), intent(inout):: itarg
integer(kind_int), intent(in   ):: nxy1
if(itarg>-1) then
  itarg = itarg-nxy1
end if

end subroutine real_itarg

end module type_mgbf_grid

