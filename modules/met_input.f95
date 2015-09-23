module helper 
!---------written by Felix Oesterle (FSO)-----------------
! -DESCRIPTION: This module contains various subroutines 
! regarding cru data, netcdf input output, precompute and various
! -TODO: split netcdf stuff to own module
! -Last modified:  Tue Aug 25, 2015  10:49
!@author Felix Oesterle
!----------------------------------------------------------
use netcdf
use nmlvars
use txt_io
use m_mrgrnk
use mmbm
!use ieee_arithmetic !provides ieee_is_nan, if not available use if(x /= x) as
implicit none

include 'geodesic.inc'

public :: read_precomputed_met
public :: precompute
public :: precompute_crosval
public :: read_clim
public :: read_anom_year
public :: read_anom_timeseries
public :: write_results_to_nc
public :: set_nan_by_idx
public :: rd_lat_lon
public :: rd_time
public :: rd_1D_var
public :: rd_1D_var_I
public :: rd_2D_var
public :: rd_3D_var
public :: create_out_nc
public :: wr_setup_nc
!public :: get_dims_nc
public :: finalize_nc
!public :: prep_ideal_nc
public :: ls_poly


public :: check

contains

subroutine read_precomputed_met(region,met_type,sim_year_start,sim_year_end, &
    meteo_months,lapse_temp,precip_per_glac,temp_per_glac_tongue,meteo_years)
!---------written by Felix Oesterle (FSO)-----------------
! -INPUT:
! -OUTPUT:
! -DESCRIPTION:
! -TODO:
! @author Felix Oesterle
! !--------------------------------------------------------
integer, intent(in) :: region
character(*), intent(in) :: met_type
integer, intent(in) :: sim_year_start,sim_year_end
real, dimension(:,:), allocatable,intent(out) :: precip_per_glac, temp_per_glac_tongue
real, dimension(:), allocatable,intent(out) ::  meteo_months, lapse_temp
integer, dimension(:), allocatable, intent(out) :: meteo_years
real, dimension(:), allocatable :: tmp1d
real, dimension(:,:), allocatable :: tmp2d
character(len=240) :: ncf
integer :: n_years,i,year,nglac,n_months

n_years = (sim_year_end - sim_year_start) + 1


! FSO--- read last file to get dimensions and check for partial year
write(ncf,'(a,i2.2,a,a,a,i2.2,a,i4,a)') 'precomputed/region',region,'/', &
        adjustl(trim(met_type)),'_precomp_',region,'_',int(sim_year_end),'.nc'

call rd_1D_var(ncf,'lapse_temp',lapse_temp)
nglac = size(lapse_temp,1)

! FSO--- check for partial year
call rd_1D_var(ncf,'months',tmp1d)
if (size(tmp1d)< 12) then
    n_months = (n_years - 1) * 12 + size(tmp1d)
    print*, "ENCOUNTERED PARTIAL YEAR"
else
    n_months = n_years * 12
end if
deallocate(tmp1d)

allocate(meteo_years(n_years))
allocate(meteo_months(n_months))
allocate(precip_per_glac(nglac,n_months))
allocate(temp_per_glac_tongue(nglac,n_months))

! FSO--- read all yearly files
i = 0
do year = sim_year_start,sim_year_end

    meteo_years(i+1) = year
    write(ncf,'(a,i2.2,a,a,a,i2.2,a,i4,a)') 'precomputed/region',region,'/', &
        adjustl(trim(met_type)),'_precomp_',region,'_',int(year),'.nc'
    call rd_1D_var(ncf,'months',tmp1d)
    meteo_months(i*12+1:(i+1)*12) = tmp1d 
    deallocate(tmp1d)
    
    call rd_2D_var(ncf,'precip_per_glac',tmp2d)
    precip_per_glac(:,i*12+1:(i+1)*12) = tmp2d 
    deallocate(tmp2d)
    
    call rd_2D_var(ncf,'temp_per_glac_tongue',tmp2d)
    temp_per_glac_tongue(:,i*12+1:(i+1)*12) = tmp2d 
    deallocate(tmp2d)
    i = i+1 
end do 

end subroutine read_precomputed_met

subroutine precompute_crosval(ncf,datadir,rgi)
!---------written by Felix Oesterle (FSO)-----------------
! -INPUT:
! -OUTPUT:
! -DESCRIPTION:
! -TODO:
! @author Felix Oesterle
! !--------------------------------------------------------
    character(*), intent(in) :: ncf
    character(*), intent(in) :: datadir 
    type(rgi_region), intent(in) :: rgi

    character(len=240) :: nc_out
    character(len=240) :: ncf_cv
    integer :: t_star_rgi_idx
    real, dimension(:), allocatable :: lon_mb, lat_mb, mean_dist
    real, dimension(:), allocatable :: weight_dist, model_bias_rgi
    real, dimension(:,:), allocatable :: dist
    real, dimension(:,:), allocatable :: model_bias 
    real, dimension(:), allocatable :: mse_mu_star_cross

    integer, dimension(:), allocatable :: rank_idx,cru_years
    integer, dimension(:), allocatable :: t_star_rgi, t_star

    integer :: g_dimid,ncid, t_star_id, mb_id, t_dimid, t_id,mse_id
    integer :: mse_dimid,glac_id
    
    real :: a, f, lat1, lon1, azi1, lat2, lon2, azi2, s12, &
         dummy
    integer :: omask
    integer :: i, j ! WGS84 values
    a = 6378137.
    f = 1.0/298.257223563
    omask = 0

    lat1=30.
    lon1=0.0
    lat2=29.5
    lon2=179.5
    call invers(a, f, lat1, lon1, lat2, lon2, & 
        s12, azi1, azi2, omask, dummy, dummy, dummy, dummy, dummy)
    
    write(ncf_cv,'(a,a)') adjustl(trim(datadir)),adjustl(trim(ncf))
    print*,ncf_cv
    !call rd_time(ncf_cv,'years',cru_years) 
    call rd_1D_var(ncf_cv,'lon',lon_mb)
    call rd_1D_var(ncf_cv,'lat',lat_mb)
    call rd_2D_var(ncf_cv,'model_bias',model_bias)
    print*,'Converting input model bias from mm to SI m'
    model_bias = model_bias / 1000.
    call rd_1D_var(ncf_cv,'mse_mu_star_cross',mse_mu_star_cross)
    call rd_1D_var_I(ncf_cv,'t_star',t_star)
    call rd_1D_var_I(ncf_cv,'time',cru_years)
    
    !(N glaciers, M mb measurements) Iceland (320,255)
    allocate(dist(size(rgi%info,1),size(lon_mb,1)))
    allocate(rank_idx(size(lon_mb,1)))
    allocate(mean_dist(size(rgi%info,1)))
    allocate(t_star_rgi(size(rgi%info,1)))
    allocate(model_bias_rgi(size(rgi%info,1)))
    allocate(weight_dist(n_closest_glac))
    
    ! FSO--- loop over all mass balance measurements
    do j = 1, size(lon_mb,1)
        ! FSO--- loop over all glaciers
        do i = 1, size(rgi%info,1)
            call invers(a, f, rgi%lat(i), rgi%lon(i), lat_mb(j), lon_mb(j), & 
                s12, azi1, azi2, omask, dummy, dummy, dummy, dummy, dummy)
            dist(i,j) = s12 / 1000.
        end do
    end do 
    !stop

    ! FSO--- loop over all glaciers
    do i = 1, size(rgi%info,1)
        ! FSO--- rank values and get index
        call r_mrgrnk(dist(i,:),rank_idx)
        mean_dist(i) = sum(dist(i,rank_idx(1:n_closest_glac)))/n_closest_glac

        ! FSO--- this is done to avoid division by 0 in next step 
        if (dist(i,rank_idx(1)) < 0.00001) then
            dist(i,rank_idx(1)) = 0.0001
        end if
        
        ! FSO--- weigh by distance
        weight_dist= 1./dist(i,rank_idx(1:n_closest_glac))  &
            *(1./sum(1./dist(i,rank_idx(1:n_closest_glac))));
!print*,weight_dist
!print*
!print*,dist(i,rank_idx(1:n_closest_glac))
        ! FSO--- determine t_star
        !if (dist(i,rank_idx(1)) < 0.00001) then
        !    t_star_rgi(i) = t_star(rank_idx(1))
        !else
        t_star_rgi(i)=nint(sum(t_star(rank_idx(1:n_closest_glac)) &
            * weight_dist));
        !end if
        
        t_star_rgi_idx = -99
        do j = 1,size(cru_years,1)
            if (cru_years(j) == t_star_rgi(i)) then
                t_star_rgi_idx = j
                exit
            end if
        end do 
        
        !if (dist(i,rank_idx(1)) < 0.00001) then
        !    model_bias_rgi(i) = model_bias(t_star_rgi_idx,rank_idx(1))
        !else
        model_bias_rgi(i)=sum(model_bias(t_star_rgi_idx,rank_idx(1:n_closest_glac)) &
            *weight_dist)
        !end if 
!print*,'T star idx',t_star_rgi_idx
!print*
!print*,model_bias(t_star_rgi_idx,rank_idx(1:n_closest_glac))
!print*
!print*,'model bias',model_bias_rgi(i)

    end do
    
    write(nc_out,'(a,i2.2,a,i2.2,a)') 'precomputed/region',rgi%region, &
        '/mb_dist_',rgi%region,'.nc'
    !write(nc_out,'(a,i2.2,a)') 'precomputed/mb_dist_',rgi%region,'.nc'
   
    !FSO---create new file, overwrite if exists 
    call check( nf90_create(nc_out, NF90_NETCDF4, ncid) )
    
    !FSO---Define the dimensions. NetCDF will hand back an ID for each.
    call check( nf90_def_dim(ncid, "glacier", size(model_bias_rgi,1), g_dimid) )
    call check( nf90_def_dim(ncid, "time", size(cru_years,1), t_dimid) )
    call check( nf90_def_dim(ncid, "mb_measurements", size(mse_mu_star_cross), mse_dimid) )
    
    !FSO---Define the variable. 
    call check( nf90_def_var(ncid, "t_star_rgi", nf90_int, g_dimid, t_star_id) )
    call check( nf90_def_var_deflate(ncid, t_star_id, shuffle = 1, deflate =1 , &
        deflate_level = 5))
    
    call check( nf90_def_var(ncid, "model_bias_rgi", nf90_double, g_dimid, mb_id) )
    call check( nf90_def_var_deflate(ncid, mb_id, shuffle = 1, deflate =1 , &
        deflate_level = 5))
    
    call check( nf90_def_var(ncid, "mse_mu_star_cross", nf90_double, mse_dimid, mse_id) )
    call check( nf90_def_var_deflate(ncid, mse_id, shuffle = 1, deflate =1 , &
        deflate_level = 5))
    
    call check( nf90_def_var(ncid, "time", nf90_double, t_dimid, t_id) )
    call check( nf90_def_var_deflate(ncid, t_id, shuffle = 1, deflate =1 , &
        deflate_level = 5))
    
    call check( nf90_def_var(ncid, "glacier", nf90_int, g_dimid, glac_id) )
    call check( nf90_def_var_deflate(ncid, glac_id, shuffle = 1, deflate =1 , &
        deflate_level = 5))
    
    ! FSO---put attributes 
    call check( nf90_put_att(ncid, nf90_global,"region",rgi%region))
    !FSO---End define mode. This tells netCDF we are done defining metadata
    call check( nf90_enddef(ncid) )

    !FSO---write variables 
    call check( nf90_put_var(ncid, t_star_id,t_star_rgi)) 
    call check( nf90_put_var(ncid, mb_id,model_bias_rgi)) 
    call check( nf90_put_var(ncid, mse_id,mse_mu_star_cross)) 
    call check( nf90_put_var(ncid, t_id,cru_years)) 
    call check( nf90_put_var(ncid, glac_id,rgi%info(:,1))) 
    !call check( nf90_sync(ncid))

    call finalize_nc(ncid)

end subroutine precompute_crosval

subroutine read_clim(ncf_clim,datadir,met_type,clim_lat,clim_lon,clim_time, &
    clim_precip,clim_temp,clim_elev)
    
    character(*), intent(in) :: ncf_clim,datadir
    character(*), intent(in) :: met_type
    real, dimension(:), allocatable, intent(out) :: clim_lat,clim_lon,clim_time
    real, dimension(:,:,:), allocatable, intent(out) :: clim_precip,clim_temp
    real, dimension(:,:), allocatable, intent(out) :: clim_elev

    character(len=240) ::ncf
    
    ! FSO--- read climatology data from nc files
    write(ncf,'(a,a)') adjustl(trim(datadir)),adjustl(trim(ncf_clim))
    print*,trim(ncf)
    call rd_lat_lon(ncf,clim_lat,clim_lon) 
    call rd_time(ncf,'month',clim_time) 
    call rd_3D_var(ncf,'tp',clim_precip) !clim_precip(lon,lat,months)
    call rd_3D_var(ncf,'t2m',clim_temp) !clim_precip(lon,lat,months)
    call rd_2D_var(ncf,'elev',clim_elev) ! elev(lon,lat)
   
end subroutine read_clim

subroutine read_anom_year(year,datadir,met_type, &
    anom_lat,anom_lon,anom_time,anom_precip,anom_temp)
    
    integer, intent(in) :: year
    character(*), intent(in) :: datadir
    character(*), intent(in) :: met_type
    real, dimension(:), allocatable, intent(out) :: anom_lat,anom_lon,anom_time
    real, dimension(:,:,:), allocatable, intent(out) :: anom_precip,anom_temp

    character(len=240) :: ncf_anom
        
    ! FSO--- read yearly data
    write(ncf_anom,'(a,a,a,i4,a)') adjustl(trim(datadir)),adjustl(trim(met_type)), &
        '_anom_interp_',year,'.nc'
    print*,trim(ncf_anom)
    call rd_lat_lon(ncf_anom,anom_lat,anom_lon) 
    call rd_time(ncf_anom,'time',anom_time) 
    call rd_3D_var(ncf_anom,'tp',anom_precip) 
    call rd_3D_var(ncf_anom,'t2m',anom_temp) 

end subroutine read_anom_year

subroutine read_anom_timeseries(year,datadir,met_type, &
    anom_lat,anom_lon,anom_time,anom_precip,anom_temp,glacier_id)
    
    integer, intent(in) :: year
    character(*), intent(in) :: datadir
    character(*), intent(in) :: met_type
    real, dimension(:), allocatable, intent(out) :: anom_lat,anom_lon,anom_time
    real, dimension(:,:,:), allocatable, intent(out) :: anom_precip,anom_temp
    integer, dimension(:), allocatable, intent(out) :: glacier_id 
    real, dimension(:,:), allocatable :: p_ts,t_ts

    character(len=240) :: ncf_anom
        
    ! FSO--- read data
    write(ncf_anom,'(a,a,a,i4,a)') adjustl(trim(datadir)),adjustl(trim(met_type)), &
        '_anom_timeseries_',year,'.nc'
    call rd_1D_var(ncf_anom,'latitude',anom_lat) 
    call rd_1D_var(ncf_anom,'longitude',anom_lon) 
    call rd_1D_var_I(ncf_anom,'glacier',glacier_id) 
    call rd_time(ncf_anom,'time',anom_time) 
    call rd_2D_var(ncf_anom,'tp',p_ts)  ! anom[time,glacier]
    call rd_2D_var(ncf_anom,'t2m',t_ts)  !anom[time,glacier]

    allocate(anom_precip(size(p_ts,1),size(p_ts,2),1))
    allocate(anom_temp(size(t_ts,1),size(t_ts,2),1))
    anom_precip(:,:,1) = p_ts
    anom_temp(:,:,1) = t_ts

end subroutine read_anom_timeseries


subroutine precompute(met_type,anom_type,rgi,clim_lat,clim_lon,clim_time, &
    anom_lat,anom_lon,anom_time,clim_precip,clim_temp,anom_precip,anom_temp,clim_elev, &
    glacier_id)
!---------written by Felix Oesterle (FSO)-----------------
! -INPUT: Anomalies and Climatologies, Glacier inventory data
! -OUTPUT:
! -DESCRIPTION:
! if anom_type is time series anom_precip and anom_temp have dims (time, glacier,1)
! @author Felix Oesterle
! !--------------------------------------------------------
    character(*), intent(in) :: met_type, anom_type
    type(rgi_region), intent(in) :: rgi
    real, dimension(:), intent(in) :: clim_lat,clim_lon,clim_time
    real, dimension(:), intent(in) :: anom_lat,anom_lon,anom_time
    real, dimension(:,:,:), intent(in) :: clim_precip,clim_temp
    real, dimension(:,:,:), intent(in) :: anom_precip,anom_temp
    real, dimension(:,:), intent(in) :: clim_elev
    integer, dimension(:),intent(in), optional :: glacier_id
    
    integer :: i,j,n_glac,n_month
    integer :: l,m
    integer :: grid_dist
    integer :: nmonths
    integer,dimension(1) :: clim_lon_idx, clim_lat_idx
    integer,dimension(1) :: anom_lon_idx, anom_lat_idx
    
    real :: dd
    real :: nan
    real :: lapse_precip
    real, dimension(9) :: elev_surround,temp_clim_surround
    
    real, dimension(:), allocatable :: height_cru_clim_glacier, lapse_temp,ev_tmp,t_tmp
    real, dimension(:), allocatable :: temp_range, factor
    real, dimension(:), allocatable :: c
    real, dimension(:,:), allocatable :: precip_per_glac
    real, dimension(:,:), allocatable :: temp_per_glac, temp_per_glac_tongue
    real, dimension(:,:), allocatable :: clim_temp_mean


    logical, dimension(9) :: nan_mask
    logical, dimension(:,:), allocatable :: clim_is_nan
    
    character(len=240) :: nc_out
    integer :: ncid, g_dimid, t_dimid, precip_id, lt_id, months_id,glac_id
    integer :: temp_id
    integer, dimension(2) :: dimids

    integer :: lonub,lonlb,latub,latlb
    
    ! FSO--- get a NaN variable
    nan = 0
    nan = nan/nan
    !write(*,*) "Precomputing CRU" 
    
    ! FSO--- if anom_type = timeseries, glacier_id has to be set
    if (anom_type == 'timeseries' .and. .not. present(glacier_id)) then
        print*,'ERROR: Glacier IDs have to be set for timeseries. Aborting'
        stop
    end if 

    ! FSO--- number of glaciers
    n_glac = size(rgi%info,1)
    
    allocate(clim_temp_mean(size(clim_temp,1),size(clim_temp,2)))
    allocate(clim_is_nan(size(clim_temp,1),size(clim_temp,2)))
    clim_is_nan = .false.
    ! FSO--- generate a mask to indicate where clim values are NaN
    do j = 1,size(clim_is_nan,2)
        do i = 1, size(clim_is_nan,1)
            if (any(clim_temp(i,j,:) /= clim_temp(i,j,:))) clim_is_nan(i,j) = .true.
        end do 
    end do
    
    ! FSO--- clim_temp == clim_temp generates a mask where nan's are false
    clim_temp_mean = sum(clim_temp,3,clim_temp == clim_temp) / 12.0
    where (clim_is_nan) clim_temp_mean = nan

    ! FSO--- holds data per glacier
    allocate(precip_per_glac(n_glac,size(anom_time,1)))
    allocate(temp_per_glac(n_glac,size(anom_time,1)))
    allocate(temp_per_glac_tongue(n_glac,size(anom_time,1)))
    allocate(height_cru_clim_glacier(n_glac))
    allocate(lapse_temp(n_glac))
    allocate(temp_range(n_glac))
    allocate(factor(size(anom_time,1)))

    n_glac = size(rgi%info,1)

    ! FSO--- i => index of glaciers
    do i = 1,n_glac
        !print*,i
        clim_lon_idx = minloc(abs(clim_lon - rgi%lon(i)))
        clim_lat_idx = minloc(abs(clim_lat - rgi%lat(i)))
        anom_lon_idx = minloc(abs(anom_lon - rgi%lon(i)))
        anom_lat_idx = minloc(abs(anom_lat - rgi%lat(i)))
        
        height_cru_clim_glacier(i)=clim_elev(clim_lon_idx(1),clim_lat_idx(1))

        if (anom_type == 'timeseries' .and. rgi%info(i,1) /= glacier_id(i) ) then
            print*, 'ERROR: GI glacier ID and timeseries glacier ID dont match!'
            stop
        end if


        ! FSO--- loop over all times
        do j = 1,size(anom_time)
            ! FSO--- get number of months
            n_month = NINT((anom_time(j)-int(anom_time(j)) +1./24.)*12.)
            if (anom_type == 'timeseries') then
                precip_per_glac(i,j) = anom_precip(j,i,1) + &
                    clim_precip(clim_lon_idx(1),clim_lat_idx(1),n_month)
                temp_per_glac(i,j) = anom_temp(j,i,1) + &
                    clim_temp(clim_lon_idx(1),clim_lat_idx(1),n_month)
            else if (anom_type == 'grid')  then
                precip_per_glac(i,j) = anom_precip(anom_lon_idx(1),anom_lat_idx(1),j) + &
                    clim_precip(clim_lon_idx(1),clim_lat_idx(1),n_month)
                temp_per_glac(i,j) = anom_temp(anom_lon_idx(1),anom_lat_idx(1),j) + &
                    clim_temp(clim_lon_idx(1),clim_lat_idx(1),n_month)
            end if
        end do 
    
        ! FSO--- calculate lapse rates from surrounding points
        grid_dist=1;! FSO---set as input parameter 
        lonlb = clim_lon_idx(1)-grid_dist
        lonub = clim_lon_idx(1)+grid_dist
        latlb = clim_lat_idx(1)-grid_dist
        latub = clim_lat_idx(1)+grid_dist
        if (lonlb < 1) then
            lonlb = 1
            lonub = lonub+1
        end if
        if (latlb < 1) then
            latlb = 1
            latub = latub+1
        end if
        
        if (lonub > size(clim_elev,1)) then
            lonub = size(clim_elev,1) 
            lonlb = lonlb-1
        end if
        if (latub > size(clim_elev,2)) then
            latub = size(clim_elev,2) 
            latlb = latlb-1
        end if
        
        elev_surround=reshape(clim_elev(lonlb:lonub,latlb:latub),(/9/))
        temp_clim_surround=reshape(clim_temp_mean(lonlb:lonub,latlb:latub),(/9/))

        nan_mask = elev_surround == elev_surround ! true for real, false for nan
        
        ! FSO--- count NaN's and if only 1 or 0 heights are available, extend
        ! region
        if (count(nan_mask) < 2) then
            grid_dist = grid_dist + 1
            lonlb = clim_lon_idx(1)-grid_dist
            lonub = clim_lon_idx(1)+grid_dist
            latlb = clim_lat_idx(1)-grid_dist
            latub = clim_lat_idx(1)+grid_dist
            if (lonlb < 1) then
                lonlb = 1
                lonub = lonub+1
            end if
            if (latlb < 1) then
                latlb = 1
                latub = latub+1
            end if
            if (lonub > size(clim_elev,1)) then
                lonub = size(clim_elev,1) 
                lonlb = lonlb-1
            end if
            if (latub > size(clim_elev,2)) then
                latub = size(clim_elev,2) 
                latlb = latlb-1
            end if
            elev_surround=reshape(clim_elev(lonlb:lonub,latlb:latub),(/9/))
            temp_clim_surround=reshape(clim_temp_mean(lonlb:lonub,latlb:latub),(/9/))
            grid_dist = grid_dist - 1
        end if
    
    
        ! FSO--- if still less than 2 elevations surrounding, set fixed lapse
        ! rate
        if (count(nan_mask) < 2) then
            !print*, 'Setting fixed lapse rate'
            lapse_temp(i)=0.0065 ! K/m
        else
            call remove_nan(elev_surround,ev_tmp)
            call remove_nan(temp_clim_surround,t_tmp)
            ! FSO--- least square polynomial fitting
            m = 1 ! order of fit 
            call LS_POLY(m,0.0,size(ev_tmp),l,ev_tmp,t_tmp,c,dd)
            !write(*,'(f10.4)') -1.0 * c(1)
            lapse_temp(i)=-1.0*c(1);
        end if

        !write(*,*) 'LapseTemp', lapse_temp(i)

        ! calculate temp over glacier tongue by assuming lapse-rate
        !% TODO FSO -- what happens if tongue height changes?
        ! FSO--- loop over all times
        do j = 1,size(anom_time)
            temp_per_glac_tongue(i,j) = temp_per_glac(i,j)+ &
                lapse_temp(i)*(height_cru_clim_glacier(i)-rgi%props(i,2)) !rgi%props(i,2) is zmin
            ! FSO--- height correction of precipitation with lapse rate 
            lapse_precip = -0.0003 !fixed lapse rate of 3% / 100m 
            precip_per_glac(i,j)=precip_per_glac(i,j)*(1.0+lapse_precip* &
                (height_cru_clim_glacier(i)-(rgi%props(i,2)+rgi%props(i,3))/2.))
                                            !   zmin           zmax
        end do 
 
        ! remove non-solid precipitation
                       !zmax              zmin
        temp_range(i)=(rgi%props(i,3)-rgi%props(i,2))*lapse_temp(i)
        factor=-1.0/temp_range(i)*(temp_per_glac_tongue(i,:)) + 1.0 + &
            t_precip_solid / temp_range(i)

        factor=min(1.0,factor)
        factor=max(0.0,factor)
        precip_per_glac(i,:) = precip_per_glac(i,:) * factor
        if (allocated(ev_tmp)) deallocate(ev_tmp)
        if (allocated(t_tmp)) deallocate(t_tmp)
        if (allocated(c)) deallocate(c)

    
    end do ! end of per galcier loop

    !! FSO--- alert NaN's  TODO convert to function/subroutine
    !write(*,*) " Precip per glacier NaN's"
    !do i = 1 , size(precip_per_glac,1)
    !    if (any(precip_per_glac(i,:) /= precip_per_glac(i,:))) then
    !        print*,i, "NaN"
    !    end if
    !end do 

    write(nc_out,'(a,i2.2,a,a,a,i2.2,a,i4,a)') 'precomputed/region',rgi%region,'/', &
        adjustl(trim(met_type)),'_precomp_',rgi%region,'_',int(anom_time(1)),'.nc'
    !print*,nc_out
   
    !FSO---create new file, overwrite if exists 
    call check( nf90_create(nc_out, NF90_NETCDF4, ncid) )
    
    !FSO---Define the dimensions. NetCDF will hand back an ID for each.
    call check( nf90_def_dim(ncid, "glacier", size(precip_per_glac,1), g_dimid) )
    call check( nf90_def_dim(ncid, "months", NF90_UNLIMITED, t_dimid) )
    dimids =  (/ g_dimid, t_dimid /) 
    
    !FSO---Define the variable. 
    call check( nf90_def_var(ncid, "precip_per_glac", nf90_double, dimids, precip_id) )
    call check( nf90_def_var_deflate(ncid, precip_id, shuffle = 1, deflate =1 , &
        deflate_level = 3))
    
    call check( nf90_def_var(ncid, "temp_per_glac_tongue", nf90_double, dimids, temp_id) )
    call check( nf90_def_var_deflate(ncid, temp_id, shuffle = 1, deflate =1 , &
        deflate_level = 3))
    
    call check( nf90_def_var(ncid, "lapse_temp", nf90_double, g_dimid, lt_id) )
    call check( nf90_def_var_deflate(ncid, lt_id, shuffle = 1, deflate =1 , &
        deflate_level = 3))

    call check( nf90_def_var(ncid, "months", nf90_double, t_dimid, months_id) )
    call check( nf90_def_var_deflate(ncid, months_id, shuffle = 1, deflate =1 , &
        deflate_level = 3))
    
    call check( nf90_def_var(ncid, "glacier", nf90_int, g_dimid, glac_id) )
    call check( nf90_def_var_deflate(ncid, glac_id, shuffle = 1, deflate =1 , &
        deflate_level = 3))
    !
    ! FSO---put attributes 
    call check( nf90_put_att(ncid, nf90_global,"region",rgi%region))
    
    !FSO---End define mode. This tells netCDF we are done defining metadata
    call check( nf90_enddef(ncid) )

    !FSO---write variables 
    call check( nf90_put_var(ncid, precip_id,precip_per_glac)) 
    call check( nf90_put_var(ncid, temp_id,temp_per_glac_tongue)) 
    call check( nf90_put_var(ncid, lt_id,lapse_temp)) 
    call check( nf90_put_var(ncid, months_id,anom_time)) 
    call check( nf90_put_var(ncid, glac_id,rgi%info(:,1))) 
    
    call finalize_nc(ncid)
    deallocate(height_cru_clim_glacier, lapse_temp)
    deallocate(temp_range, factor)
    deallocate(precip_per_glac)
    deallocate(temp_per_glac, temp_per_glac_tongue)
    deallocate(clim_temp_mean)
    deallocate(clim_is_nan)



end subroutine precompute

subroutine write_errors_to_nc(nc_out,years,rgi,V,dV,dL, L, A, mb_modeled)
!---------written by Felix Oesterle (FSO)-----------------
! -INPUT:
! -OUTPUT:
! -DESCRIPTION: all inputs are errors
! -TODO:
! -Last modified:  Fri Nov 07, 2014  15:39
! @author Felix Oesterle
!--------------------------------------------------------
    character(len=*), intent(in) :: nc_out
    type(rgi_region), intent(in) :: rgi
    real, dimension(:,:), intent(inout) :: V,dV, dL, L, A, mb_modeled
    integer, dimension(:), intent(in) :: years

    integer :: ncid, g_dimid, t_dimid
    integer :: V_id, dV_id, dL_id, L_id, mb_id,A_id
    integer :: y_id,g_id,i
    integer, dimension(2) :: dimids

    print*, 'Writing results to ', nc_out

    ! FSO--- set NaN's for glaciers with nan_idx = 1
    call set_nan_by_idx(V,rgi%nan_idx)
    call set_nan_by_idx(A,rgi%nan_idx)
    call set_nan_by_idx(dV,rgi%nan_idx)
    call set_nan_by_idx(L,rgi%nan_idx)
    call set_nan_by_idx(dL,rgi%nan_idx)
    call set_nan_by_idx(mb_modeled,rgi%nan_idx)

    !FSO---create new file, overwrite if exists 
    call check( nf90_create(nc_out, NF90_NETCDF4, ncid) )
    
    !FSO---Define the dimensions. NetCDF will hand back an ID for each.
    call check( nf90_def_dim(ncid, 'glacier', size(A,1), g_dimid) )
    call check( nf90_def_dim(ncid, 'years', size(years,1), t_dimid) )
    dimids =  (/ g_dimid, t_dimid /) 
    
    !FSO---Define the variable. 
    call check( nf90_def_var(ncid, "A_error", nf90_double, dimids, a_id) )
    call check( nf90_def_var_deflate(ncid, a_id, shuffle = 1, deflate =1 , &
        deflate_level = 3))
    
    call check( nf90_def_var(ncid, "V_error", nf90_double, dimids, V_id) )
    call check( nf90_def_var_deflate(ncid, V_id, shuffle = 1, deflate =1 , &
        deflate_level = 3))
    
    call check( nf90_def_var(ncid, "dV_cum_error", nf90_double, dimids, dV_id) )
    call check( nf90_def_var_deflate(ncid, dV_id, shuffle = 1, deflate =1 , &
        deflate_level = 3))

    call check( nf90_def_var(ncid, "L_error", nf90_double, dimids, L_id) )
    call check( nf90_def_var_deflate(ncid, L_id, shuffle = 1, deflate =1 , &
        deflate_level = 3))
    
    call check( nf90_def_var(ncid, "dL_cum_error", nf90_double, dimids, dL_id) )
    call check( nf90_def_var_deflate(ncid, dL_id, shuffle = 1, deflate =1 , &
        deflate_level = 3))
    
    call check( nf90_def_var(ncid, "mb_error", nf90_double, dimids, mb_id) )
    call check( nf90_def_var_deflate(ncid, mb_id, shuffle = 1, deflate =1 , &
        deflate_level = 3))

    call check( nf90_def_var(ncid, "years", nf90_double, t_dimid, y_id) )
    call check( nf90_def_var_deflate(ncid, y_id, shuffle = 1, deflate =1 , &
        deflate_level = 3))
    
    call check( nf90_def_var(ncid, "glacier", nf90_double, g_dimid, g_id) )
    call check( nf90_def_var_deflate(ncid, g_id, shuffle = 1, deflate =1 , &
        deflate_level = 3))
    
    ! FSO---put attributes 
    call check( nf90_put_att(ncid, nf90_global,"region",rgi%region))
    
    !FSO---End define mode. This tells netCDF we are done defining metadata
    call check( nf90_enddef(ncid) )

    !FSO---write variables 
    call check( nf90_put_var(ncid, a_id,A)) 
    call check( nf90_put_var(ncid, V_id,V)) 
    call check( nf90_put_var(ncid, dV_id,dV)) 
    call check( nf90_put_var(ncid, L_id,L)) 
    call check( nf90_put_var(ncid, dL_id,dL)) 
    call check( nf90_put_var(ncid, mb_id,mb_modeled)) 
    call check( nf90_put_var(ncid, y_id,years)) 
    call check( nf90_put_var(ncid, g_id, [ (i,i=1,size(A,1)) ]     )) 

    call finalize_nc(ncid)


end subroutine write_errors_to_nc

subroutine write_results_to_nc(nc_out,years,months,rgi,V,dV,dL, L, A, mb_modeled,mb_mod_mly,zterm)
!---------written by Felix Oesterle (FSO)-----------------
! -INPUT:
! -OUTPUT:
! -DESCRIPTION:
! -TODO:
! -Last modified:  Thu Aug 27, 2015  11:45
! @author Felix Oesterle
!--------------------------------------------------------
    character(len=*), intent(in) :: nc_out
    type(rgi_region), intent(in) :: rgi
    real, dimension(:,:), intent(inout) :: V,dV, dL, L, A, mb_modeled,mb_mod_mly,zterm
    integer, dimension(:), intent(in) :: years
    real, dimension(:), intent(in) :: months

    integer :: ncid, g_dimid, t_dimid, m_dimid
    integer :: V_id, dV_id, dL_id, L_id, mb_id,A_id,mmb_id,m_id,zterm_id
    integer :: y_id,g_id,i
    integer, dimension(2) :: dimids

    print*, 'Writing results to ', nc_out

    ! FSO--- set NaN's for glaciers with nan_idx = 1
    call set_nan_by_idx(V,rgi%nan_idx)
    call set_nan_by_idx(A,rgi%nan_idx)
    call set_nan_by_idx(dV,rgi%nan_idx)
    call set_nan_by_idx(L,rgi%nan_idx)
    call set_nan_by_idx(dL,rgi%nan_idx)
    call set_nan_by_idx(mb_modeled,rgi%nan_idx)
    call set_nan_by_idx(zterm,rgi%nan_idx)
    call set_nan_by_idx(mb_mod_mly,rgi%nan_idx)

    !FSO---create new file, overwrite if exists 
    call check( nf90_create(nc_out, NF90_NETCDF4, ncid) )
    
    !FSO---Define the dimensions. NetCDF will hand back an ID for each.
    call check( nf90_def_dim(ncid, 'glacier', size(A,1), g_dimid) )
    call check( nf90_def_dim(ncid, 'years', size(years,1), t_dimid) )
    call check( nf90_def_dim(ncid, 'months', size(months,1), m_dimid) )
    dimids =  (/ g_dimid, t_dimid /) 
    
    !FSO---Define the variable. 
    call check( nf90_def_var(ncid, "A", nf90_double, dimids, a_id) )
    call check( nf90_def_var_deflate(ncid, a_id, shuffle = 1, deflate =1 , &
        deflate_level = 3))
    call check( nf90_put_att(ncid,a_id,'units','m^2'))
    call check( nf90_put_att(ncid,a_id,'longname','area'))
    
    call check( nf90_def_var(ncid, "V", nf90_double, dimids, V_id) )
    call check( nf90_def_var_deflate(ncid, V_id, shuffle = 1, deflate =1 , &
        deflate_level = 3))
    call check( nf90_put_att(ncid,V_id,'units','m^3'))
    call check( nf90_put_att(ncid,V_id,'longname','volume'))
    
    call check( nf90_def_var(ncid, "dV", nf90_double, dimids, dV_id) )
    call check( nf90_def_var_deflate(ncid, dV_id, shuffle = 1, deflate =1 , &
        deflate_level = 3))
    call check( nf90_put_att(ncid,dV_id,'units','m^3'))
    call check( nf90_put_att(ncid,dV_id,'longname','volume change'))

    call check( nf90_def_var(ncid, "L", nf90_double, dimids, L_id) )
    call check( nf90_def_var_deflate(ncid, L_id, shuffle = 1, deflate =1 , &
        deflate_level = 3))
    call check( nf90_put_att(ncid,L_id,'units','m'))
    call check( nf90_put_att(ncid,L_id,'longname','length'))
    
    call check( nf90_def_var(ncid, "dL", nf90_double, dimids, dL_id) )
    call check( nf90_def_var_deflate(ncid, dL_id, shuffle = 1, deflate =1 , &
        deflate_level = 3))
    call check( nf90_put_att(ncid,dL_id,'units','m'))
    call check( nf90_put_att(ncid,dL_id,'longname','length change'))
    
    call check( nf90_def_var(ncid, "massbalance_modeled", nf90_double, dimids, mb_id) )
    call check( nf90_def_var_deflate(ncid, mb_id, shuffle = 1, deflate =1 , &
        deflate_level = 3))
    call check( nf90_put_att(ncid,mb_id,'units','m'))
    call check( nf90_put_att(ncid,mb_id,'longname','mass balance yearly'))
    
    call check( nf90_def_var(ncid, "z_term", nf90_double, dimids, zterm_id) )
    call check( nf90_def_var_deflate(ncid, zterm_id, shuffle = 1, deflate =1 , &
        deflate_level = 3))
    call check( nf90_put_att(ncid,zterm_id,'units','m'))
    call check( nf90_put_att(ncid,zterm_id,'longname','terminus height'))
    
    call check( nf90_def_var(ncid, "massbalance_modeled_monthly", nf90_double, &
        (/ g_dimid, m_dimid /), mmb_id) )
    call check( nf90_def_var_deflate(ncid, mmb_id, shuffle = 1, deflate =1 , &
        deflate_level = 3))
    call check( nf90_put_att(ncid,mmb_id,'units','m/m^2'))
    call check( nf90_put_att(ncid,mmb_id,'longname','monthly mass balance'))

    call check( nf90_def_var(ncid, "years", nf90_double, t_dimid, y_id) )
    call check( nf90_def_var_deflate(ncid, y_id, shuffle = 1, deflate =1 , &
        deflate_level = 3))
    call check( nf90_put_att(ncid,y_id,'units','decimal year'))
    call check( nf90_put_att(ncid,y_id,'longname','years'))
    
    call check( nf90_def_var(ncid, "months", nf90_double, m_dimid, m_id) )
    call check( nf90_def_var_deflate(ncid, m_id, shuffle = 1, deflate =1 , &
        deflate_level = 3))
    call check( nf90_put_att(ncid,m_id,'units','decimal year'))
    
    call check( nf90_def_var(ncid, "glacier", nf90_int, g_dimid, g_id) )
    call check( nf90_def_var_deflate(ncid, g_id, shuffle = 1, deflate =1 , &
        deflate_level = 3))
    
    ! FSO---put attributes 
    call check( nf90_put_att(ncid, nf90_global,"region",rgi%region))
    
    print*, "BIG TODO: add proper description to variables in netcdf output"
    
    !FSO---End define mode. This tells netCDF we are done defining metadata
    call check( nf90_enddef(ncid) )

    !FSO---write variables 
    call check( nf90_put_var(ncid, a_id,A)) 
    call check( nf90_put_var(ncid, V_id,V)) 
    call check( nf90_put_var(ncid, dV_id,dV)) 
    call check( nf90_put_var(ncid, L_id,L)) 
    call check( nf90_put_var(ncid, dL_id,dL)) 
    call check( nf90_put_var(ncid, mb_id,mb_modeled)) 
    call check( nf90_put_var(ncid, mmb_id,mb_mod_mly)) 
    call check( nf90_put_var(ncid, y_id,years)) 
    call check( nf90_put_var(ncid, m_id,months)) 
    call check( nf90_put_var(ncid, zterm_id,zterm)) 
    !call check( nf90_put_var(ncid, g_id, [ (i,i=1,size(A,1)) ]     )) 
    call check( nf90_put_var(ncid, g_id, rgi%info(:,1)  )) 

    call finalize_nc(ncid)


end subroutine write_results_to_nc

subroutine set_nan_by_idx(var,idx)
!---------written by Felix Oesterle (FSO)-----------------
! -INPUT:
! -OUTPUT:
! -DESCRIPTION:
! -TODO:
! @author Felix Oesterle
!--------------------------------------------------------
    real, dimension(:,:), intent(inout) :: var
    integer, dimension(:), intent(in) :: idx

    integer :: i
    real :: nan
    ! FSO--- get a NaN variable
    nan = 0
    nan = nan/nan

    do i = 1, size(var,1)
        if (idx(i) == 1) then
            var(i,:) = nan
        end if
    end do 

end subroutine set_nan_by_idx
    
subroutine rd_1D_var_I(ncf,var_name,var)
!---------written by Felix Oesterle (FSO)-----------------
! -INPUT: ncf: name of netcdf file
!   var_name : guess what
! -OUTPUT:
!   var : values of var_name variable
! -DESCRIPTION:
! -TODO:
! @author Felix Oesterle
! !--------------------------------------------------------
    character(*), intent(in) :: ncf
    character(*), intent(in) :: var_name
    integer, dimension(:), allocatable,intent(out) :: var
    
    character(len=120) :: vn
    integer :: ncid
    integer :: varid, nd,n1
    integer,dimension(2):: dids
    character(len=90) :: lnm
    call check (nf90_open(ncf, nf90_write,ncid))
    call check (nf90_inq_varid(ncid, var_name, varid) )
    ! FSO--- get info about variable and allocate
    call check (nf90_inquire_variable(ncid,varid,name=vn,ndims=nd,dimids=dids))
    !print*,trim(vn),nd,dids
    call check (nf90_inquire_dimension(ncid, dids(1), lnm, n1))
    allocate(var(n1))
    call check (nf90_get_var(ncid, varid, var)) 
    call check (nf90_close(ncid))
end subroutine rd_1D_var_I


subroutine rd_1D_var(ncf,var_name,var)
!---------written by Felix Oesterle (FSO)-----------------
! -INPUT:
! -OUTPUT:
! -DESCRIPTION:
! -TODO:
! @author Felix Oesterle
! !--------------------------------------------------------
    character(*), intent(in) :: ncf
    character(*), intent(in) :: var_name
    character(len=120) :: vn
    real, dimension(:), allocatable,intent(out) :: var
    integer :: ncid
    integer :: varid, nd,n1
    integer,dimension(2):: dids
    character(len=90) :: lnm
    call check (nf90_open(ncf, nf90_write,ncid))
    call check (nf90_inq_varid(ncid, var_name, varid) )
    ! FSO--- get info about variable and allocate
    call check (nf90_inquire_variable(ncid,varid,name=vn,ndims=nd,dimids=dids))
    !print*,trim(vn),nd,dids
    call check (nf90_inquire_dimension(ncid, dids(1), lnm, n1))
    allocate(var(n1))
    call check (nf90_get_var(ncid, varid, var)) 
    call check (nf90_close(ncid))
end subroutine rd_1D_var


subroutine rd_2D_var(ncf,var_name,var)
!---------written by Felix Oesterle (FSO)-----------------
! -INPUT:
! -OUTPUT:
! -DESCRIPTION:
! -TODO:
! @author Felix Oesterle
! !--------------------------------------------------------
    character(*), intent(in) :: ncf
    character(*), intent(in) :: var_name
    character(len=120) :: vn
    real, dimension(:,:), allocatable,intent(out) :: var
    integer :: ncid
    integer :: varid, nd,n1,n2
    integer,dimension(2) :: dids
    character(len=90) :: lnm
    call check (nf90_open(ncf, nf90_write,ncid))
    call check (nf90_inq_varid(ncid, var_name, varid) )
    ! FSO--- get info about variable and allocate
    call check (nf90_inquire_variable(ncid,varid,name=vn,ndims=nd,dimids=dids))
    !print*,trim(vn),nd,dids
    call check (nf90_inquire_dimension(ncid, dids(1), lnm, n1))
    call check (nf90_inquire_dimension(ncid, dids(2), lnm, n2))
    allocate(var(n1,n2))
    call check (nf90_get_var(ncid, varid, var)) 
    call check (nf90_close(ncid))
end subroutine rd_2D_var


subroutine rd_3D_var(ncf,var_name,var)
!---------written by Felix Oesterle (FSO)-----------------
! -INPUT:
! -OUTPUT:
! -DESCRIPTION:
! -TODO:
! @author Felix Oesterle
! !--------------------------------------------------------
    character(*), intent(in) :: ncf
    character(*), intent(in) :: var_name
    character(len=120) :: vn
    real, dimension(:,:,:), allocatable,intent(out) :: var
    integer :: ncid
    integer :: varid, nd,n1,n2,n3
    integer,dimension(3) :: dids
    character(len=90) :: lnm
    call check (nf90_open(ncf, nf90_write,ncid))
    call check (nf90_inq_varid(ncid, var_name, varid) )
    ! FSO--- get info about variable and allocate
    call check (nf90_inquire_variable(ncid,varid,name=vn,ndims=nd,dimids=dids))
    call check (nf90_inquire_dimension(ncid, dids(1), lnm, n1))
    call check (nf90_inquire_dimension(ncid, dids(2), lnm, n2))
    call check (nf90_inquire_dimension(ncid, dids(3), lnm, n3))
    allocate(var(n1,n2,n3))
    call check (nf90_get_var(ncid, varid, var)) 
    call check (nf90_close(ncid))
end subroutine rd_3D_var



subroutine rd_lat_lon(ncf,lat,lon)
!---------written by Felix Oesterle (FSO)-----------------
! -INPUT:
! -OUTPUT:
! -DESCRIPTION:
! -TODO:
! @author Felix Oesterle
! !--------------------------------------------------------
    character(*), intent(in) :: ncf
    real, dimension(:), allocatable,intent(out) :: lat,lon
    integer :: ncid
    integer :: ny,nx
    integer :: t
    integer :: y_dimid,x_dimid
    character(len=90) :: lnm
    call check (nf90_open(ncf, NF90_NOWRITE,ncid))
    call check (nf90_inq_dimid(ncid, "latitude", y_dimid) )
    call check (nf90_inquire_dimension(ncid, y_dimid, lnm, ny))
    call check (nf90_inq_dimid(ncid, "longitude", x_dimid) )
    call check (nf90_inquire_dimension(ncid, x_dimid, lnm, nx))
    allocate(lat(ny),lon(nx))
    call check (nf90_get_var(ncid, y_dimid, lat)) 
    call check (nf90_get_var(ncid, x_dimid, lon)) 
    call check (nf90_close(ncid))
end subroutine rd_lat_lon

subroutine rd_time(ncf,time_name,time)
!---------written by Felix Oesterle (FSO)-----------------
! -INPUT:
! -OUTPUT:
! -DESCRIPTION:
! -TODO:
! @author Felix Oesterle
! !--------------------------------------------------------
    character(*), intent(in) :: ncf
    character(*), intent(in) :: time_name
    real, dimension(:), allocatable,intent(out) :: time
    integer :: ncid
    integer :: nt
    integer :: time_dimid
    character(len=90) :: lnm
    call check (nf90_open(ncf, nf90_write,ncid))
    call check (nf90_inq_dimid(ncid, time_name, time_dimid) )
    call check (nf90_inquire_dimension(ncid,time_dimid, lnm, nt))
    !print*,nt
    allocate(time(nt))
    call check (nf90_get_var(ncid, time_dimid, time)) 
    call check (nf90_close(ncid))
end subroutine rd_time


subroutine create_out_nc(ncf,nx,ny,ncid,dimids)
!---------written by Felix Oesterle (FSO)-----------------
!! -INPUT:
!! -OUTPUT:
!! -DESCRIPTION: creates a nc_file if nc_diff_out is true in namelistfile
!! -TODO:
!!@author Felix Oesterle
!----------------------------------------------------------
    character(*), intent(in) :: ncf
    integer,intent(inout) :: ncid
    integer, intent(in) :: nx, ny
    integer, dimension(2), intent(inout) :: dimids
    integer :: x_dimid, y_dimid,lat_varid,lon_varid 
    integer :: x_varid,y_varid 
    real, dimension(:,:), allocatable :: lat,lon
    real, dimension(:), allocatable :: x,y 
    real :: dx,dy
    character(len=90) :: yu,xu
    character(len=90) :: latu,lonu
    character(len=90) :: x_lnm,y_lnm
    character(len=90) :: lat_lnm,lon_lnm
print*,"in create out"
    ! FSO---sync topography ncfile  
    call check( nf90_sync(ncid))
   

     ! FSO---get dimension: lat 
    call check (nf90_inq_varid(ncid, "lat", y_dimid) )
    allocate(y(ny))
    call check (nf90_get_var(ncid,y_dimid,y))
    call check (nf90_get_att(ncid,y_dimid,"units",yu))
    call check (nf90_get_att(ncid,y_dimid,"long_name",y_lnm))
    call check (nf90_get_att(ncid,y_dimid,"gridspace_lat",dy))

    ! FSO---get dimension: lon 
    call check (nf90_inq_varid(ncid, "lon", x_dimid) )
    allocate(x(nx))
    call check (nf90_get_var(ncid,x_dimid,x))
    call check (nf90_get_att(ncid,x_dimid,"units",xu))
    call check (nf90_get_att(ncid,x_dimid,"long_name",x_lnm))
    call check (nf90_get_att(ncid,x_dimid,"gridspace_lon",dx))
    
    ! FSO---get variable: lat
    call check (nf90_inq_varid(ncid, "lat", lat_varid) )
    allocate(lat(nx,ny))
    call check (nf90_get_var(ncid,lat_varid,lat))
    call check (nf90_get_att(ncid,lat_varid,"units",latu))
    call check (nf90_get_att(ncid,lat_varid,"long_name",lat_lnm))

    ! FSO---get variable: lon
    call check (nf90_inq_varid(ncid, "lon", lon_varid) )
    allocate(lon(nx,ny))
    call check (nf90_get_var(ncid,lon_varid,lon))
    call check (nf90_get_att(ncid,lon_varid,"units",lonu))
    call check (nf90_get_att(ncid,lon_varid,"long_name",lon_lnm))
    
    ! FSO---close topography ncfile  
    call check (nf90_close(ncid))
   

    print*, "creating new nc file: ", ncf
    !FSO---create new file, overwrite if exists 
    call check( nf90_create(ncf, NF90_CLOBBER, ncid) )
    
    !FSO---Define the dimensions. NetCDF will hand back an ID for each.
    call check( nf90_def_dim(ncid, "lon", nx, x_dimid) )
    call check( nf90_def_dim(ncid, "lat", ny, y_dimid) )
    
    ! Define the coordinate variables. They will hold the coordinate
    ! information, that is, the lon and lat. A varid is
    ! returned for each.
    call check( nf90_def_var(ncid, "lat", NF90_REAL, y_dimid, y_varid) )
    call check( nf90_def_var(ncid, "lon", NF90_REAL, x_dimid, x_varid) )

    ! assign units attributes to coordinate var data. this attaches a
    ! text attribute to each of the coordinate variables, containing the
    ! units.
    call check( nf90_put_att(ncid, y_varid, "units", yu) )
    call check( nf90_put_att(ncid, y_varid, "long_name", y_lnm) )
    call check( nf90_put_att(ncid, y_varid, "gridspace_lat", dy) )
    call check( nf90_put_att(ncid, x_varid, "units", xu) )
    call check( nf90_put_att(ncid, x_varid, "long_name", x_lnm) )
    call check( nf90_put_att(ncid, x_varid, "gridspace_lon", dx) )


    !FSO---The dimids array is used to pass the IDs of the dimensions of
    ! the variables. Note that in fortran arrays are stored in
    ! column-major format.
    dimids =  (/ x_dimid, y_dimid /) !FSO: do not trust ncbrowse order of dims, there seems to be
    !a bug in it...

    ! Define secondary coord variables. They will hold the coordinate
    ! information, that is, the latitudes and longitudes. A varid is
    ! returned for each.
    call check( nf90_def_var(ncid, "lat", NF90_REAL, dimids, lat_varid) )
    call check( nf90_def_var(ncid, "lon", NF90_REAL, dimids, lon_varid) )

    ! Assign units attributes to secondary coordinate var data. This attaches a
    ! text attribute to each of the coordinate variables, containing the
    ! units.
    call check( nf90_put_att(ncid, lat_varid, "units", latu) )
    call check( nf90_put_att(ncid, lat_varid, "long_name", lat_lnm) )
    call check( nf90_put_att(ncid, lon_varid, "long_name", lon_lnm) )
    call check( nf90_put_att(ncid, lon_varid, "units", lonu) )


    !FSO---End define mode. This tells netCDF we are done defining metadata
    call check( nf90_enddef(ncid) )
    
    ! FSO---put values from topography file
    call check( nf90_put_var(ncid, x_varid,x))
    call check( nf90_put_var(ncid, y_varid,y))
    
    ! FSO---put values from topography file
    call check( nf90_put_var(ncid, lat_varid,lat))
    call check( nf90_put_var(ncid, lon_varid,lon))
    
    call check( nf90_sync(ncid))
    print*,"FSO: created new nc file"

end subroutine create_out_nc

subroutine wr_setup_nc(ncid)
!---------written by Felix Oesterle (FSO)-----------------
!! -INPUT: ncid: id of netcdf file
!! -OUTPUT:
!! -DESCRIPTION:
!! -TODO:
!! -Last modified:  Wed Jun 25, 2014  13:44
!!@author Felix Oesterle
!----------------------------------------------------------
    use nmlvars
    integer,intent(in) :: ncid

    
    
    !call check( nf90_redef(ncid))
    !
    !call check( nf90_put_att(ncid, nf90_global,"Nlon",nx))
    !call check( nf90_put_att(ncid, nf90_global,"Nlat",ny))
    !call check( nf90_put_att(ncid, nf90_global,"Undef",undef))
    !if(use_fall) then
    !    call check( nf90_put_att(ncid, nf90_global,"use_fall","true"))
    !else
    !    call check( nf90_put_att(ncid, nf90_global,"use_fall","false"))
    !end if
    !
    call check( nf90_enddef(ncid) )
end subroutine wr_setup_nc


!subroutine get_dims_nc(ncf,nx,ny,delx,dely)
!!---------written by Felix Oesterle (FSO)-----------------
!!! -INPUT: ncf : name of netcdf file
!!! -OUTPUT: nx,ny,delx,dely
!!! -DESCRIPTION:
!!! -TODO:
!!! -Last modified:  Wed Sep 23, 2015  15:57
!!!@author Felix Oesterle
!!----------------------------------------------------------
!    character(*), intent(in) :: ncf
!    integer, intent(inout) :: nx, ny
!    real, intent(inout) :: delx, dely
!    integer :: y_dimid,x_dimid,ncid
!    character(len=90) :: lnm
!
!    ! FSO--- not sur what is going on here, possibly a bug in netcdf ?!
!    ! There seems to be a mixup in the dimensions
!    ! The structure of an example file looks like this (ncdump -h ... )    
!    !dimensions:
!        !lon = 55 ;
!        !lat = 78 ;
!    !variables:
!        !double lat(lat) ;
!                !lat:long_name = "northward distance from southwest corner of domain in projection coordinates" ;
!                !lat:axis = "lat" ;
!                !lat:gridspace_lat = 500. ;
!        !double lon(lon) ;
!                !lon:long_name = "eastward distance from southwest corner of domain in projection coordinates" ;
!                !lon:axis = "lon" ;
!                !lon:gridspace_lon = 333. ;
!        !double topo(lat, lon) ;
!                !topo:coordinates = "lat lon" ;
!        !double lat(lat, lon) ;
!                !lat:units = "degrees_north" ;
!                !lat:long_name = "Latitude" ;
!        !double lon(lat, lon) ;
!                !lon:units = "degrees_east" ;
!                !lon:long_name = "Longitude" ;
! 
!    
!    ! NOTE that lat has attribute gridspace_lat !!!
!    ! In the code below I get y_dimid from the dimension 'lat', however this one
!    ! has the attribute gridspace_lon !!! Don't know why....
!    
!    
!    
!    
!    print*, "in get_dims"
!    call check (nf90_open(ncf, nf90_write,ncid))
!    call check (nf90_inq_dimid(ncid, "lat", y_dimid) )
!    call check (nf90_inquire_dimension(ncid, y_dimid, lnm, ny))
!    call check (nf90_get_att(ncid,y_dimid,"gridspace_lon",delx)) !BUg IN NETCDF ?!?
!    
!    call check (nf90_inq_dimid(ncid, "lon", x_dimid) )
!    call check (nf90_inquire_dimension(ncid, x_dimid, lnm, nx))
!    call check (nf90_get_att(ncid,x_dimid,"gridspace_lat",dely)) !BUG IN NETCDF ?!?!
!
!
!    print*, "In get dims"
!    print*, "dely",dely,ny
!    print*, "delx",delx,nx
!end subroutine get_dims_nc

subroutine finalize_nc(ncid)
    integer, intent(in) :: ncid
    call check( nf90_sync(ncid))
    call check (nf90_close(ncid))
end subroutine finalize_nc

  
!


!------------------------------------------------
! FSO: subroutine in conjunction with netcdf stuff. Taken from netcdf examples
!------------------------------------------------
subroutine check(status)
integer, intent ( in) :: status

    if(status /= nf90_noerr) then
        print *, trim(nf90_strerror(status))
        stop "Stopped because of netcdf call"
    end if

end subroutine check

subroutine check_exi(status,ncid,name,id)
integer, intent(in) :: status
integer, intent(in) :: ncid
character(len=*) :: name
integer, intent(inout) :: id

    if (status .eq. nf90_enameinuse) then 
        print*, "WARNING: Got an existing variable in netcdf file, will overwrite..."
        call check(nf90_inq_varid(ncid,name,id))
    else if(status /= nf90_noerr) then
        print*, status 
        print*, trim(nf90_strerror(status))
        stop "Stopped because of netcdf call"
     end if

end subroutine check_exi


!*****************************************************************
!*         LEAST SQUARES POLYNOMIAL FITTING PROCEDURE            *
!* ------------------------------------------------------------- *
!* This program least squares fits a polynomial to input data.   *
!* forsythe orthogonal polynomials are used in the fitting.      *
!* The number of data points is n.                               *
!* The data is input to the subroutine in x[i], y[i] pairs.      *
!* The coefficients are returned in c[i],                        *
!* the smoothed data is returned in v[i],                        *
!* the order of the fit is specified by m.                       *
!* The standard deviation of the fit is returned in d.           *
!* There are two options available by use of the parameter e:    *
!*  1. if e = 0, the fit is to order m,                          *
!*  2. if e > 0, the order of fit increases towards m, but will  *
!*     stop if the relative standard deviation does not decrease *
!*     by more than e between successive fits.                   *
!* The order of the fit then obtained is l.                      *
!*****************************************************************
Subroutine LS_POLY(m,e1,n,l,x,y,c,dd)
  integer :: si 
  integer, intent(inout) :: m
  !Labels: 10,15,20,30,50
  real,dimension(:),intent(in) :: x,y
  real,dimension(:), allocatable, intent(out) :: c
  real,dimension(:), allocatable :: v, a, b
  real, dimension(:), allocatable :: d,c2,e,f
  integer :: i,l,l2,n,n1
  real :: a1,a2,b1,b2,c1,dd,d1,e1,f1,f2,v1,v2,w,vv
 
  si = size(x)
  
  allocate(c(0:si))
  allocate(v(si),a(si),b(si)) 
  allocate(d(si),c2(si),e(si),f(si))
  
  n1 = m + 1; l=0
  v1 = 1.7
  ! Initialize the arrays
  do i=1, n1
    a(i) = 0.0; b(i) = 0.0; f(i) = 0.0
  end do
  do i=1, n
    v(i) = 0.0; d(i) = 0.0
  end do
  d1 = sqrt(float(n)); w = d1;
  do i=1, n
    e(i) = 1.0 / w
  end do
  f1 = d1; a1 = 0.0
  do i=1, n
    a1 = a1 + x(i) * e(i) * e(i)
  end do
  c1 = 0.0
  do i=1, n
    c1 = c1 + y(i) * e(i)
  end do
  b(1) = 1.0 / f1; f(1) = b(1) * c1
  do i=1, n
    v(i) = v(i) + e(i) * c1
  end do
  m = 1
! Save latest results
10 do i=1, l
    c2(i) = c(i)
  end do
  l2 = l; v2 = v1; f2 = f1; a2 = a1; f1 = 0.0
  do i=1, n
    b1 = e(i)
    e(i) = (x(i) - a2) * e(i) - f2 * d(i)
    d(i) = b1
    f1 = f1 + e(i) * e(i)
  end do
  f1 = sqrt(f1)
  do i=1, n
    e(i) = e(i) / f1
  end do
  a1 = 0.0
  do i=1, n  
    a1 = a1 + x(i) * e(i) * e(i)
  end do
  c1 = 0.0
  do i=1, n  
    c1 = c1 + e(i) * y(i)
  end do
  m = m + 1; i = 0
15 l = m - i; b2 = b(l); d1 = 0.0
  if (l > 1)  d1 = b(l - 1)
  d1 = d1 - a2 * b(l) - f2 * a(l)
  b(l) = d1 / f1; a(l) = b2; i = i + 1
  if (i.ne.m) goto 15
  do i=1, n
    v(i) = v(i) + e(i) * c1 
  end do
  do i=1, n1
    f(i) = f(i) + b(i) * c1
    c(i) = f(i)
  end do
  vv = 0.0
  do i=1, n
    vv = vv + (v(i) - y(i)) * (v(i) - y(i))
  end do
  !Note the division is by the number of degrees of freedom
  !vv = dsqrt(vv / dfloat(n - l - 1)); l = m
  vv = sqrt(vv / float(n - l - 1)); l = m
  if (e1.eq.0.0) goto 20
  !Test for minimal improvement
  if (abs(v1 - vv) / vv < e1) goto 50
  !if error is larger, quit
  if (e1 * vv > e1 * v1) goto 50
  v1 = vv
20 if (m.eq.n1) goto 30
  goto 10
!Shift the c[i] down, so c(0) is the constant term
30 do i=1, l  
    c(i - 1) = c(i)
  end do
  c(l) = 0.0
  ! l is the order of the polynomial fitted
  l = l - 1; dd = vv
  return
! Aborted sequence, recover last values
50 l = l2; vv = v2
  do i=1, l  
    c(i) = c2(i)
  end do
  goto 30
end subroutine



end module helper 
