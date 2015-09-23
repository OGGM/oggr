!---------written by Felix Oesterle (FSO)-----------------
! -DESCRIPTION: precompute file for Fortran MMBM
! -TODO:
! -Last modified:  Wed Sep 23, 2015  15:58
! @author Felix Oesterle
! !--------------------------------------------------------
program mmbm_precomp
!use ieee_arithmetic !provides ieee_is_nan, if not available use if(x /= x) as
!                    !NaN test
use nmlvars
use txt_io
use helper 
!use mmbm
!use omp_lib ! openMP
implicit none

character(len=120) :: argu
character(len=120) :: infile
character(len=120) :: gi_file
character(len=120) :: ncf 


type(rgi_region),dimension(:), allocatable :: rgi

integer :: j,i,year


real, dimension(:), allocatable :: clim_lat,clim_lon,clim_time
real, dimension(:), allocatable :: anom_lat,anom_lon,anom_time
real, dimension(:,:,:), allocatable :: clim_precip,clim_temp
real, dimension(:,:,:), allocatable :: anom_precip,anom_temp
real, dimension(:,:), allocatable :: clim_elev

integer, dimension(:), allocatable :: t_star_rgi, meteo_years
integer, dimension(:), allocatable :: glacier_id
real, dimension(:), allocatable :: mu_rgi, turnover
real, dimension(:), allocatable :: model_bias_rgi
real, dimension(:), allocatable :: mse_mu_star_cross
real, dimension(:,:), allocatable :: precip_per_glac, temp_per_glac_tongue
real, dimension(:), allocatable :: meteo_months, lapse_temp

integer :: g_dimid,ncid, mu_id, t_id
 
real :: start
character(len=240)  :: datadir, outdir


real :: nan

integer:: t1, t2
integer ::  clock_rate, clock_max 
call system_clock ( t1, clock_rate, clock_max )

!FSO---directory setup 
datadir = 'data/'
outdir = 'results/'

! FSO--- get a NaN variable
nan = 0
nan = nan/nan

! FSO---get_command_argument is Fortran 2003 standard, so not all compilers
! might support it. If not supported: use f2kcli (google/yahoo or whatever it)....
! Enables model to be run with -f INFILE argument
call get_command_argument(1,argu)                                                      
if (trim(argu) == "-f") then
    call get_command_argument(2,infile)
    write(*,*)
    write(*,fmt="(a,a)") "Using specified infile: ", infile
else 
    write(*,"(a,a)") "-----------------------------------------------------------------------"
    write(*,"(a,a)")  "No infile specified, using default: MMBM_Settings.nml"
    write(*,"(a,a)")  "Give command line option -f INFILE if you want to use a different one!"
    write(*,"(a,a)") "-----------------------------------------------------------------------"
    infile = "MMBM_Settings.nml"
end if

call cpu_time(start) 

! FSO---Open and read setup from namelist INFILE
call get_namelist(infile)

! FSO--- LOAD METEOROLOGICAL INPUT ClIMATOLOGY ----------
call read_clim(ncf_clim,datadir,met_type,clim_lat,clim_lon,clim_time, &
    clim_precip,clim_temp,clim_elev)

! FSO---arbitrary factor
clim_precip=clim_precip*clim_precip_factor


! FSO--- LOAD RGI------------
allocate(rgi(size(regions)))
do i = 1,size(regions)
    write(*,"(a,a)") "--------------------------------------------------------------------"
    write(*,"(a,i5)")  "Loading region",regions(i)
    ! rgi type contains 5 arrays:
    ! lat and lon
    ! info(glacier_id,ice_sheet)
    ! props(area,zmin,zmax,zmean,zstd)
    ! nan_idx
    ! year
    ! this makes it less readable but easier to work on
    write(gi_file,'(a,i2.2,a)') 'rgi/region',regions(i),'.txt'
    call get_rgi(gi_file,regions(i),rgi(i)) 
end do



do year = precomp_year_start,precomp_year_end
    !write(*,*) "Precomputing year  ", year
    ! FSO--- LOAD METEOROLOGICAL INPUT Anomalies----------
    ! FSO--- check type of anomaly and read accordingly
    ! FSO--- not via openmp parallelizable, as netcdf is read. parallel netcdf
    ! is needed...
    if (anom_type == 'timeseries') then
        call read_anom_timeseries(year,datadir,met_type,anom_lat,anom_lon, &
            anom_time,anom_precip,anom_temp,glacier_id) !anom[time,glacier,1]
    else if (anom_type == 'grid') then
        call read_anom_year(year,datadir,met_type,anom_lat,anom_lon,anom_time, &
            anom_precip,anom_temp) !anom[lat,lon,time]
    else 
        write(*,*) "Unknown anom_type, exiting"
        stop
    end if

    do i = 1,size(regions)
        write(*,"(a,a)") "--------------------------------------------------------------------"
        write(*,"(a,i5,a,i5)")  "Precomputing region",regions(i),' Year ',year

        ! Do precomputation here
        ! i.e. extract all data for each glacier
        call precompute(met_type,anom_type,rgi(i),clim_lat,clim_lon,clim_time, &
            anom_lat,anom_lon, anom_time,clim_precip,clim_temp,anom_precip, &
            anom_temp,clim_elev,glacier_id)
        !call system_clock ( t2, clock_rate, clock_max )
        !write(*,'("Wall time precompute = ",f9.3," seconds.")')real ( t2 - t1 ) / real ( clock_rate )
    end do 
    deallocate(anom_lat,anom_lon,anom_time)
    deallocate(anom_precip,anom_temp)
end do

!

do i = 1,size(regions)
    write(*,"(a,a)") "--------------------------------------------------------------------"
    write(*,"(a,i5)")  "Precomputing crossval region",regions(i)

    call precompute_crosval(ncf_crosval,datadir,rgi(i))

    !call system_clock ( t2, clock_rate, clock_max )
    !write(*,'("Wall time precompute = ",f9.3," seconds.")')real ( t2 - t1 ) / real ( clock_rate )
end do

! FSO--- if met_type is CRU, precompute mu_rgi as well 
if (met_type == 'CRU') then
    do i = 1,size(regions)
        write(*,"(a,a)") "--------------------------------------------------------------------"
        write(*,"(a,i5)")  "Precomputing mu_rgi region",regions(i)

        ! FSO--- LOAD METEOROLOGICAL INPUT ----------
        ! read precomputed values
        ! data contains NaN's, so make sure to treat them correctly
        write(ncf,'(a,i2.2,a,i2.2,a)') 'precomputed/region',rgi(i)%region,'/mb_dist_',rgi(i)%region,'.nc'
        !write(ncf,'(a,i2.2,a)') 'precomputed/mb_dist_',rgi(i)%region,'.nc'
        call rd_1D_var_I(ncf,'t_star_rgi',t_star_rgi)
        call rd_1D_var(ncf,'model_bias_rgi',model_bias_rgi)
        call rd_1D_var(ncf,'mse_mu_star_cross',mse_mu_star_cross)

        call read_precomputed_met(rgi(i)%region,met_type,1901,2009, &
            meteo_months,lapse_temp,precip_per_glac,temp_per_glac_tongue,meteo_years)
        
        allocate(mu_rgi(size(t_star_rgi)))
        allocate(turnover(size(t_star_rgi)))
        ! FSO--- calc model parameter mu_rgi and turnover
        do j = 1, size(t_star_rgi)
            call calc_mu_rgi(length_of_climatology,t_star_rgi(j), meteo_years, & 
                temp_per_glac_tongue(j,:), &
                precip_per_glac(j,:),rgi(i)%nan_idx(j), &
                mu_rgi(j),turnover(j))

        end do 
        
        write(ncf,'(a,i2.2,a,i2.2,a)') 'precomputed/region',rgi(i)%region,'/mu_rgi_',rgi(i)%region,'.nc'
        !write(ncf,'(a,i2.2,a)') 'precomputed/mu_rgi_',rgi(i)%region,'.nc'

        !FSO---create new file, overwrite if exists 
        call check( nf90_create(ncf, NF90_NETCDF4, ncid) )
        
        !FSO---Define the dimensions. NetCDF will hand back an ID for each.
        call check( nf90_def_dim(ncid, "glacier", size(t_star_rgi,1), g_dimid) )
        
        !FSO---Define the variable. 
        call check( nf90_def_var(ncid, "mu_rgi", nf90_double, g_dimid, mu_id) )
        call check( nf90_def_var_deflate(ncid, mu_id, shuffle = 1, deflate =1 , &
            deflate_level = 5))
        
        call check( nf90_def_var(ncid, "turnover", nf90_double, g_dimid, t_id) )
        call check( nf90_def_var_deflate(ncid, t_id, shuffle = 1, deflate =1 , &
            deflate_level = 5))
        
        
        ! FSO---put attributes 
        call check( nf90_put_att(ncid, nf90_global,"region",rgi(i)%region))
        !FSO---End define mode. This tells netCDF we are done defining metadata
        call check( nf90_enddef(ncid) )

        !FSO---write variables 
        call check( nf90_put_var(ncid, mu_id,mu_rgi)) 
        call check( nf90_put_var(ncid, t_id,turnover)) 

        call finalize_nc(ncid)

        deallocate(mu_rgi,turnover)
     end do 
end if 

end program mmbm_precomp
