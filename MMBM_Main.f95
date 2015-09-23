!---------written by Felix Oesterle (FSO)-----------------
! -DESCRIPTION: main file for Fortran MMBM
! -TODO:
! -Last modified:  Wed Sep 23, 2015  15:57
! @author Felix Oesterle
! !--------------------------------------------------------
program mmbm_main
!use ieee_arithmetic !provides ieee_is_nan, if not available use if(x /= x) as
!                    !NaN test
use nmlvars
use txt_io
use helper 
use mmbm
use omp_lib ! openMP
implicit none

character(len=120) :: argu
character(len=120) :: infile,ncf
character(len=120) :: result_ncf,gi_file

real, dimension(:), allocatable :: model_bias_rgi
real, dimension(12) :: mb_months_base

type(rgi_region) :: rgi

integer, dimension(:), allocatable :: t_star_rgi, meteo_years
integer, dimension(:), allocatable :: n_months_above_freezing
integer, dimension(:,:), allocatable :: n_months_above_freezing_per_glac
!integer, dimension(12) :: months_idx
integer :: i,j,k,idx
integer :: n_years, n_glac, n_months
integer :: sim_years

real :: gamma, c_a, q, c_l
real :: A_test
real :: A_begin
real :: start, start_loop, finish
real, dimension(:), allocatable :: V, dL, L, A, mb_modeled,mb_mod_mly
real, dimension(:), allocatable :: z_terminus
real, dimension(:), allocatable :: V_error, dL_cum_error, L_error, A_error
real, dimension(:), allocatable :: mb_error, dV_cum_error
real, dimension(:), allocatable :: A_1, meteo_months, lapse_temp
real, dimension(:), allocatable :: mu_rgi, turnover
real, dimension(:), allocatable :: mse_mu_star_cross 
real, dimension(:,:), allocatable :: precip_per_glac, temp_per_glac_tongue
real, dimension(:,:), allocatable :: V_per_glac, L_per_glac
real, dimension(:,:), allocatable :: dV_per_glac, dL_per_glac
real, dimension(:,:), allocatable :: dA_per_glac, A_per_glac
real, dimension(:,:), allocatable :: mb_modeled_per_glac 
real, dimension(:,:), allocatable :: z_terminus_per_glac 
real, dimension(:,:), allocatable :: mb_mod_mly_per_glac 

real, dimension(:,:), allocatable :: V_error_per_glac, L_error_per_glac
real, dimension(:,:), allocatable :: A_error_per_glac, mb_error_per_glac
real, dimension(:,:), allocatable :: dV_cum_error_per_glac, dL_cum_error_per_glac

character(len=240) :: nc_out
character(len=240)  :: datadir, outdir
integer :: ncid, g_dimid, t_dimid, precip_id, lt_id, months_id
integer :: temp_id
integer, dimension(2) :: dimids


integer :: tnr,t

real :: nan

integer:: t1, t2
integer ::  clock_rate, clock_max 
call system_clock ( t1, clock_rate, clock_max )

! FSO--- get a NaN variable
nan = 0
nan = nan/nan

!FSO---directory setup 
datadir = 'data'
outdir = 'results'

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

write(*,"(a,a)") "-----------------------------------------------------------------------"
write(*,"(a,i5)")  "Simulating region",regions(1)

! FSO--- LOAD RGI------------
! rgi type contains 5 arrays:
! lat and lon
! info(glacier_id,ice_sheet)
! props(area,zmin,zmax,zmean,zstd)
! nan_idx
! year
! this makes it less readable but easier to work on
write(gi_file,'(a,i2.2,a)') 'rgi/region',regions(1),'.txt'
call get_rgi(gi_file,regions(1),rgi) 

! FSO--- LOAD METEOROLOGICAL INPUT ----------
! read precomputed values
! data contains NaN's, so make sure to treat them correctly
write(ncf,'(a,i2.2,a,i2.2,a)') 'precomputed/region',rgi%region,'/mb_dist_',rgi%region,'.nc'
!write(ncf,'(a,i2.2,a)') 'precomputed/mb_dist_',rgi%region,'.nc'
call rd_1D_var_I(ncf,'t_star_rgi',t_star_rgi)
call rd_1D_var(ncf,'model_bias_rgi',model_bias_rgi)
call rd_1D_var(ncf,'mse_mu_star_cross',mse_mu_star_cross)

call read_precomputed_met(rgi%region,met_type,sim_year_start,sim_year_end, &
    meteo_months,lapse_temp,precip_per_glac,temp_per_glac_tongue,meteo_years)

! FSO--- read precomputed mu_rgi
! mu_rgi contains only CRU data for now
write(ncf,'(a,i2.2,a,i2.2,a)') 'precomputed/region',rgi%region,'/mu_rgi_',rgi%region,'.nc'
write(*,*) ncf
call rd_1D_var(ncf,'mu_rgi',mu_rgi)
call rd_1D_var(ncf,'turnover',turnover)

n_years = size(meteo_years,1)
n_months = size(meteo_months,1)

write(*,*) "Number of years: ", n_years
write(*,*) "Number of months: ", n_months 

! FSO--- allocate additional parameters that are computed at runtime
n_glac = size(precip_per_glac,1)

! FSO--- 
allocate(V_per_glac(n_glac,n_years), L_per_glac(n_glac,n_years))
allocate(dL_per_glac(n_glac,n_years))
allocate(dA_per_glac(n_glac,n_years), A_per_glac(n_glac,n_years))
! FSO--- the next will have a NaN at the beginning 
allocate(mb_modeled_per_glac(n_glac,n_years))
allocate(z_terminus_per_glac(n_glac,n_years))
allocate(mb_mod_mly_per_glac(n_glac,n_months))
allocate(dV_per_glac(n_glac,n_years))
allocate(n_months_above_freezing_per_glac(n_glac,n_years))

allocate(V_error_per_glac(n_glac,n_years), L_error_per_glac(n_glac,n_years))
allocate(A_error_per_glac(n_glac,n_years), mb_error_per_glac(n_glac,n_years))
allocate(dV_cum_error_per_glac(n_glac,n_years), dL_cum_error_per_glac(n_glac,n_years))


! FSO--- set environment variable:
! set -x OMP_NUM_THREADS 4
! or
!call omp_set_num_threads( 20 )
t = omp_get_max_threads()
print*, 'Number of openmp threads: ', t


! FSO--- START FIRST MODEL ITERATION
! loop over all glaciers

! FSO--- parallelization with openMP
!$omp parallel do private(i,k,gamma,c_a,q,c_l,mb_months_base,A_begin,A_1,start_loop), &
!$omp& private(sim_years,V,dL,L,A,mb_modeled,mb_mod_mly,n_months_above_freezing,finish), &
!$omp& private(idx,j,A_test,V_error,A_error,L_error,mb_error,dV_cum_error,dL_cum_error,z_terminus)
do i=1,size(rgi%info,1)
    call cpu_time(start_loop) 
    !tnr = omp_get_thread_num()
    !write( *, * ) 'Thread', tnr, ':',  i
    ! FSO--- set parameters according to icecap status
    if (rgi%info(i,2)==1) then
        gamma=gamma_icecap
        c_a=c_a_icecap
        q=q_icecap
        c_l=c_l_icecap
    else
        gamma=gamma_glacier
        c_a=c_a_glacier
        q=q_glacier
        c_l=c_l_glacier
    end if

    ! FSO---mass balance years (NH: oct - sept, SH: apr - mar)
    if (rgi%lat(i) >= 0) then
        mb_months_base = [((-3./12.+1./24.) + real(i) * 1./12.,i=0,11)]
    else 
        mb_months_base = [((-9./12.+1./24.) + real(i) * 1./12.,i=0,11)]
    end if

    ! FSO--- use inventory area as first guess for starting area
    A_begin=rgi%props(i,1) !rgi area

    !rgi_year: year of inventory 
    !min(meteo_years): first year of driving dataset
    sim_years = rgi%year(i)-minval(meteo_years)
    
    
    if (sim_years < 0) then
        print*,rgi%year(i),minval(meteo_years)
        print*,'RGI year is not within meteo_years, STOPPING'
        print*,'Need to implement filling with CRU data'
        stop
    else if (sim_years == 0) then
        print*, 'Got area right a begin, not iterating'
    else 
        ! FSO---iterativly compute area at beginning 
        if (.not. allocated(A_1)) allocate(A_1(50))

        do k = 1,50 ! start loop iteration


            call mmb_model(A_begin,meteo_years(1:sim_years), meteo_months, &
                precip_per_glac(i,:), &
                temp_per_glac_tongue(i,:), mb_months_base, &
                c_a, c_l, q, gamma, rgi%props(i,3), rgi%props(i,2), &
                lapse_temp(i), mu_rgi(i), turnover(i), model_bias_rgi(i), &
                V, dL, L, A, mb_modeled,mb_mod_mly,z_terminus,n_months_above_freezing)
            
            ! FSO--- A_1 contains history of already tried starting areas
            A_1(k)=A(1) !could also be A_begin...

            idx = 0
            do j = 1,size(meteo_years,1)
                if (meteo_years(j) == rgi%year(i)) then
                    idx = j
                    exit
                end if
            end do 
            
            ! FSO--- compare inventory area with the simulated area of the same year
            A_test=A(idx)

            ! FSO ---  props(i,1) => area
            if (abs(A_test/rgi%props(i,1)-1.0)<0.01) then
                exit ! exit iteration loop
            end if
            
            ! FSO--- set new beginning area
            if (A_test<rgi%props(i,1) .and. k==1) then
                A_begin = A_1(k)*2.0
            else if (A_test>rgi%props(i,1) .and. k==1) then
                A_begin = A_1(k)*1./2.
            else if (A_test<rgi%props(i,1) .and. k/=1 ) then
                if (A_1(k) == maxval(A_1(1:k))) then
                    A_begin=A_1(k)+rgi%props(i,1)
                else
                    A_begin=A_1(k)+abs(A_1(k)-A_1(k-1))/2.0
                end if 
            else if (A_test>rgi%props(i,1) .and. k/=1) then
                A_begin=A_1(k)-abs(A_1(k)-A_1(k-1))/2.0
            end if 

        end do ! loop iteration

        ! FSO--- if 50 iterations are not enough disregard glacier
        if (k-1==50) then
            rgi%nan_idx(i)=1
        end if

        ! FSO--- set starting area 
        A_begin = A(1)
    end if 
    ! FSO--- model the glacier for the full period 
    call mmb_model(A_begin,meteo_years(1:size(meteo_years,1)-1), &
            meteo_months, precip_per_glac(i,:), &
            temp_per_glac_tongue(i,:), mb_months_base, &
            c_a, c_l, q, gamma, rgi%props(i,3), rgi%props(i,2), &
            lapse_temp(i), mu_rgi(i), turnover(i), model_bias_rgi(i), &
            V, dL, L, A, mb_modeled,mb_mod_mly,z_terminus,n_months_above_freezing)

    V_per_glac(i,:) = V
    dL_per_glac(i,:) = dL 
    L_per_glac(i,:) = L
    A_per_glac(i,:) = A 
    mb_modeled_per_glac(i,:) = mb_modeled 
    z_terminus_per_glac(i,:) = z_terminus
    mb_mod_mly_per_glac(i,:) = mb_mod_mly 
    dV_per_glac(i,1) = NaN
    dV_per_glac(i,2:) = mb_modeled(2:)*A(1:size(A,1)-1)*1.0e-6
    
    ! FSO--- this line is base on Ben's code, although I don't get it yet. TODO.
    ! Why is last year just doubled?
    !n_months_above_freezing(size(n_months_above_freezing,1))= &
        !n_months_above_freezing(size(n_months_above_freezing,1)-1)
    n_months_above_freezing_per_glac(i,:) = n_months_above_freezing
   
    !! props(area,zmin,zmax,zmean,zstd)
    !print*, "TODO CHECK ERRROR MODEL FOR units. (especially mse_mu_star_cross)"
    !call error_model(meteo_years,rgi%year(i),rgi%props(i,1),rgi%props(i,2), &
    !    rgi%props(i,3), A,V,L,dL,mb_modeled,mse_mu_star_cross, &
    !    n_months_above_freezing,mu_rgi(i),c_l, c_a, &
    !    turnover(i), gamma, q, lapse_temp(i), V_error, A_error, L_error, &
    !    mb_error, dV_cum_error, dL_cum_error)

    !V_error_per_glac(i,:) = V_error
    !A_error_per_glac(i,:) = A_error
    !L_error_per_glac(i,:) = L_error
    !mb_error_per_glac(i,:) = mb_error
    !dV_cum_error_per_glac(i,:) = dV_cum_error 
    !dL_cum_error_per_glac(i,:) = dL_cum_error 


    deallocate(V)
    deallocate(dL)
    deallocate(L)
    deallocate(A)
    deallocate(mb_modeled)
    deallocate(mb_mod_mly)
    deallocate(z_terminus)
    !deallocate(V_error)
    !deallocate(L_error)
    !deallocate(A_error)
    !deallocate(mb_error)
    !deallocate(dL_cum_error)
    !deallocate(dV_cum_error)

    if (i/1000.==nint(i/1000.)) then
        call cpu_time(finish)
        write(*,'("Glacier ",i5," loop time = ",f9.3,"s.")') &
            i, finish-start_loop
    end if
end do ! loop over all glacier
!$omp end parallel do

! FSO--- write to netcdf file
! FSO: Dimensions conversion as not all vars are SI yet
! BIG TODO!!
write(result_ncf,'(a,a,a,a,i2.2,a)') trim(outdir),'/',adjustl(trim(met_type)),'_results_',rgi%region,'.nc'
A_per_glac = A_per_glac * 1000. * 1000. !km^2 -> m^2
dL_per_glac = dL_per_glac * 1000. !km -> m
dV_per_glac = dV_per_glac * 1000. * 1000. *1000. !km^3 -> m^3
L_per_glac = L_per_glac * 1000. !km -> m
V_per_glac = V_per_glac * 1000. * 1000. *1000. !km^3 -> m^3

call write_results_to_nc(result_ncf,meteo_years,meteo_months,rgi,V_per_glac, &
    dV_per_glac, dL_per_glac,L_per_glac, &
    A_per_glac,mb_modeled_per_glac,mb_mod_mly_per_glac,z_terminus_per_glac)

!write(result_ncf,'(a,a,i2.2,a)') trim(outdir),'/mmbm_errors_',rgi%region,'.nc'
!call write_errors_to_nc(result_ncf,meteo_years,rgi,V_error_per_glac, &
!    dV_cum_error_per_glac, dL_cum_error_per_glac,L_error_per_glac, &
!    A_error_per_glac,mb_error_per_glac)


call cpu_time(finish) 
call system_clock ( t2, clock_rate, clock_max )
write(*,'("Wall time(system clock) = ",f9.3," seconds.")')real ( t2 - t1 ) / real ( clock_rate )
write(*,'("CPU time = ",f9.3," seconds.")')finish-start

end program mmbm_main
