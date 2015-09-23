module nmlvars
implicit none
!---------written by Felix Schueller (FSS)-----------------
! -INPUT:
! -OUTPUT:
! -DESCRIPTION: Used for namelist reading of MMBM
! -TODO:
! -Last modified:  Fri Dec 12, 2014  10:59
!!@author Felix Schueller
!----------------------------------------------------------

! Variable which are defined from the NAMELIST
!FSS---IMPORTANT: Every time you add a variable please also modifiy 
! subroutine wr_setup_nc@nc_io and the namelist file 
type :: nc_variable
    character(len=32) :: name = '---'
    character(len=32) :: units = '---'
    integer :: id
    logical :: output = .false.
end type nc_variable

type :: nc_data
    real, dimension(:,:), allocatable :: data_2d
    real, dimension(:,:,:), allocatable :: data_3d
end type nc_data

type(nc_variable),dimension(20),save :: nc_var
type(nc_data),dimension(20) :: nc_var_data

integer :: nx
integer :: ny
integer :: ncx,ncy
character(len=240) :: regions_string
integer,dimension(:),allocatable :: regions
integer :: n_closest_glac
integer :: sim_year_end, sim_year_start, precomp_year_start, precomp_year_end
character(len=120)  :: nc_outfile,met_type,anom_type
character(len=240)  :: ncf_crosval,ncf_clim,ncf_anom
integer :: nc_flag
integer :: length_of_climatology
real :: A_V_scaling_error, L_V_scaling_error, area_error_rgi
real :: tau_A_error, tau_L_error
real :: gamma_glacier, c_a_glacier, q_glacier, c_l_glacier
real :: gamma_icecap, c_a_icecap, q_icecap, c_l_icecap
real :: t_precip_solid, t_melt
real :: clim_precip_factor

!FSS---IMPORTANT: Every time you add a variable please also modifiy 
! subroutine wr_setup_nc@nc_io and the namelist file 
namelist /CONTROL/  regions_string, & 
                    nc_outfile, &
                    sim_year_start, &
                    sim_year_end, &
                    precomp_year_start, &
                    precomp_year_end, &
                    ncf_clim, &
                    met_type, &
                    anom_type, &
                    ncf_anom, &
                    ncf_crosval, &
                    A_V_scaling_error, &
                    L_V_scaling_error, &
                    area_error_rgi, &
                    tau_L_error, &
                    tau_A_error, &
                    gamma_glacier, &
                    c_a_glacier, &
                    q_glacier, &
                    c_l_glacier, &
                    gamma_icecap, &
                    c_a_icecap,  &
                    q_icecap,  &
                    c_l_icecap, &
                    n_closest_glac, &
                    length_of_climatology, &
                    t_melt, &
                    clim_precip_factor, &
                    t_precip_solid

!FSS---IMPORTANT: Every time you add a variable please also modifiy 
! subroutine wr_setup_nc@nc_io and the namelist file 

public :: get_namelist

contains

subroutine get_namelist(infile)
!---------written by Felix Schueller (FSS)-----------------
! INPUT:
!   infile - filename with namelist variables
! OUTPUT:
!   set variables  which are globally available
! DESCRIPTION:
! TODO:
! last modified: 
!----------------------------------------------------------
    character(len=*),intent(in) :: infile
    character(len=240) :: regions_tmp
    integer :: inunit
    integer :: iostat
    integer :: pos,i


    ! Open namelist file
    inunit = 21

    open(inunit,file=infile,status="OLD",iostat=iostat)

    if (iostat /= 0) then
        print*,"ERROR with INFILE: ",infile
        stop
    end if


    read(inunit, nml=CONTROL)

    close(inunit)

    ! FSS--- convert region string to integer array
    regions_tmp = regions_string
    i = 0
    do while (index(regions_tmp,',') /= 0) 
        i = i + 1
        pos = index(regions_tmp, ",")
        regions_tmp = regions_tmp(pos+1:)
    end do 

    allocate(regions(i+1))
    read(regions_string,*) regions

end subroutine get_namelist


end module nmlvars
