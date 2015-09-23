module mmbm
!---------written by Felix Schueller (FSS)-----------------
! -DESCRIPTION: contains the main model and general helper routines
! -TODO:
! -Last modified:  Fri Nov 07, 2014  15:35
! @author Felix Schueller
!--------------------------------------------------------
use nmlvars
implicit none
public :: mmb_model
public :: error_model
public :: calc_mu_rgi
public :: remove_nan

contains

subroutine error_model(meteo_years,gi_year,gi_area,gi_zmin,gi_zmax, &
    A,V,L,dL,mb_modeled,mse_mu_star_cross,n_months_above_freezing,mu_rgi,c_l, c_a, &
    turnover, gam, q, lapse_temp, V_error, A_error, L_error, mb_error, dV_cum_error, &
    dL_cum_error)
!---------written by Felix Schueller (FSS)-----------------
! -INPUT: 
!   meteo_years 
!   gi_year rgi_year(i)
!   gi_area  rgi_area(i)
!   gi_zmax rgi_zmax(i)
!   V(i,:)
!   L(i,:)
!   dL(i,:)
!   mb_modeled(i,:)
!   mse_mu_star_cross
!   n_months_above_freezing(i,:)
!   mu_rgi(i)
!   c_l,q
!   lapse_temp : lapse_temp(i)
!   turnover, gam
! -OUTPUT:
! -DESCRIPTION:
!  calculate errors
! -TODO:
! @author Felix Schueller
!--------------------------------------------------------
    integer, dimension(:), intent(in) :: meteo_years
    integer, intent(in) :: gi_year
    real, intent(in) ::  gi_area, gi_zmin, gi_zmax
    real, dimension(:), intent(in) :: A, V, L, dL, mb_modeled
    real, dimension(:), intent(in) :: mse_mu_star_cross
    integer, dimension(:), intent(in) :: n_months_above_freezing
    real, intent(in) :: mu_rgi, c_l, c_a, turnover, gam
    real, intent(in) :: lapse_temp, q
    
    real, dimension(:), allocatable, intent(out) :: V_error, L_error
    real, dimension(:), allocatable, intent(out) :: A_error, mb_error
    real, dimension(:), allocatable, intent(out) :: dV_cum_error, dL_cum_error
    
    integer :: sim_years
    integer :: start_j_idx
    integer :: j
    real :: A_relax, dA
    real :: tau_A, tau_L
    real :: dA_error, dL_error
    real :: dTz_error, dV_error
    real :: L_relax
    real :: dz, mean_mse
    real :: nan
    ! FSS--- get a NaN variable
    nan = 0
    nan = nan/nan

    sim_years = size(meteo_years)

    ! FSS--- allocate arrays 
    allocate(A_error(sim_years))
    allocate(V_error(sim_years))
    allocate(L_error(sim_years))
    allocate(mb_error(sim_years)) ! will have a NaN at the beginning 
    allocate(dV_cum_error(sim_years)) !     "
    allocate(dL_cum_error(sim_years)) !     "
    mb_error = nan
    dV_cum_error = nan
    dL_cum_error = nan
    
    
    ! FSS--- get index of rgi year in meteo_years
    start_j_idx = -99
    do j = 1,size(meteo_years,1)
        if (meteo_years(j) == gi_year) then 
            start_j_idx = j
            exit
        end if
    end do 

    mean_mse = sum(mse_mu_star_cross)/size(mse_mu_star_cross,1)
    
    A_error(start_j_idx)=gi_area * area_error_rgi/100.0
    V_error(start_j_idx)=V(start_j_idx) * A_V_scaling_error
    L_error(start_j_idx)=L(start_j_idx) * L_V_scaling_error
    mb_error = NaN
    mb_error(start_j_idx)=mean_mse ** (1./2.)


    !FSS TODO: need +1 ?
    dV_cum_error(start_j_idx)= min(V(start_j_idx), &
        abs(mb_modeled(start_j_idx+1)*A(start_j_idx)*1.0e-6))* &
        ((mb_error(start_j_idx)/ &
        max(10.0,abs(mb_modeled(start_j_idx+1))))**2.0+ &
        (A_error(start_j_idx)/max(1.0e-6,A(start_j_idx)))**2.0)**(1./2.)
    dL_cum_error(start_j_idx)=abs(dL(start_j_idx)* &
        ((tau_L_error)**2.0+(L_V_scaling_error)**2.0)**(1./2.))

!write(*,'(2f15.6)') dV_cum_error(start_j_idx), dL_cum_error(start_j_idx) 

    do j=start_j_idx,sim_years-1
!        print*, 
!print*, meteo_years(j)
!delta_z=((zmax-zmin) /  ((V(1)/c_l)**(1./q)))*(((V(1)/c_l)**(1.0/q))-L(j))
!delta_T_z= -delta_z*lapse_temp

        dz=((gi_zmax-gi_zmin)/((V(1)/c_l)**(1.0/q)))*(L(start_j_idx)-L(j+1))
        dTz_error=lapse_temp*dz*L_error(j)/max(1.0e-3,L(j))
        
        mb_error(j+1)=(mean_mse+ &
                n_months_above_freezing(j+1)*(mu_rgi*dTz_error)**2.)**(1./2.)
!write(*,'(a,i4,f15.6)') ' nm(j+1)', n_months_above_freezing(j+1),mb_modeled(j+1)
!write(*,'(a,f15.6,a,f15.6)') ' mb_error(j)', mb_error(j),' mb_error(j+1)',mb_error(j+1)


!write(*,'(3f15.6)') V(j), mb_modeled(j+1), A(j)
        dV_error=min(V(j),abs(mb_modeled(j+1)*A(j)*1.e-6))* &
                ((mb_error(j)/max(10.0,abs(mb_modeled(j+1))))**2+ &
                (A_error(j)/max(1.e-6,A(j)))**2.)**(1./2.)
        V_error(j+1)=(dV_error**2.+V_error(j)**2)**(1./2.)

        L_relax=(V(j+1)/c_l)**(1./q)

        dL_error=abs(dL(j)*((tau_L_error)**2+(L_V_scaling_error)**2)**(1./2.))

        dL_cum_error(j+1)=((dL_cum_error(j))**2+dL_error**2)**(1./2.)
        L_error(j+1)=((L_error(j))**2+(dL_error)**2)**(1./2.)

        A_relax=(V(j+1)/c_a)**(1.0/gam)
        tau_L=max(1.0,V(j)/(turnover*A(j)*1.0e-6))
!write(*,'(a,f9.6)') 'tau_L', tau_L
        tau_A=max(1.0,tau_L*(c_l**(2.0/q))/(c_a**(1.0/gam))*(V(j)**(1.0/gam-2.0/q)))
!write(*,'(a,f9.6)') 'tau_A', tau_A
        
        dA=1./tau_A*(A_relax-A(j))
        dA_error=abs(dA*((tau_A_error)**2+(A_V_scaling_error)**2)**(1./2.))
        A_error(j+1)=((A_error(j))**2+(dA_error)**2)**(1./2.)
!write(*,'(a,f9.6)') 'dV_error', dV_error
!write(*,'(a,f9.6)') 'V_error', V_error(j+1)
!write(*,'(a,f9.6)') 'L_relax', L_relax 
!write(*,'(a,f9.6)') 'dL_error', dL_error 
!write(*,'(a,f9.6)') 'dL_cum_error', dL_cum_error(j+1)
!write(*,'(a,f9.6)') 'L_error', L_error(j+1)
!write(*,'(a,f9.6)') 'A_rlax', A_relax 
!write(*,'(a,f9.6)') 'dA', dA
!write(*,'(a,f9.6)') 'dA_error', dA_error
!write(*,'(a,f9.6)') 'A_error', A_error(j+1)

    !    if V(i,j+1)==0;
        if (abs(V(j+1)) < 1.0e-12) then
            L_error(j+1)=L_error(j)
            A_error(j+1)=A_error(j)
        end if
    end do
!stop

    do j=start_j_idx,size(meteo_years,1)-2
        dV_cum_error(j+1)=((dV_cum_error(j))**2+(min(V(j+1), &
            abs(mb_modeled(j+2)*A(j+1)*1.e-6))* &
            ((mb_error(j+1)/max(10.0,abs(mb_modeled(j+2))))**2+ &
            (A_error(j+1)/max(1.0e-6,A(j+1)))**2)**(1./2.))**2)**(1./2.)
!write(*,'(a,f9.6)') 'dV_cum_error', dV_cum_error(j+1)
    end do

    do j=start_j_idx-1,1,-1
!print*, 
!print*, meteo_years(j)
        dz=((gi_zmax-gi_zmin)/((V(1)/c_l)**(1./q)))*(L(start_j_idx)-L(j))
        dTz_error=lapse_temp*dz*L_error(j+1)/max(1.e-3,L(j+1))
        mb_error(j)=(mean_mse+n_months_above_freezing(j)* &
            (mu_rgi*dTz_error)**2)**(1./2.)
!write(*,'(a,f15.6)') 'mb_error', mb_error(j)

        dV_error=min(V(j+1),abs(mb_modeled(j+1)*A(j+1)*1.e-6))* &
            ((mb_error(j+1)/max(10.0,abs(mb_modeled(j+1))))**2+ &
            (A_error(j+1)/max(1.e-6,A(j+1)))**2)**(1./2.)
!write(*,'(a,f15.6)') 'dV_error', dV_error
        V_error(j)=(dV_error**2+V_error(j+1)**2)**(1./2.)
!write(*,'(a,f15.6)') 'V_error', V_error(j)
        L_relax=(V(j)/c_l)**(1./q)
        dL_error=abs(dL(j+1)*((tau_L_error)**2+(L_V_scaling_error)**2)**(1./2.))
!write(*,'(a,f15.6)') 'dL_error', dL_error
        dL_cum_error(j)=((dL_cum_error(j+1))**2+dL_error**2)**(1./2.)
!write(*,'(a,f15.6)') 'dL_cum_error', dL_cum_error(j)
        L_error(j)=((L_error(j+1))**2+(dL_error)**2)**(1./2.)

        tau_L=max(1.0,V(j)/(turnover*A(j)*1.0e-6))
        tau_A=max(1.0,tau_L*(c_l**(2.0/q))/(c_a**(1.0/gam))*(V(j)**(1.0/gam-2.0/q)))
!write(*,'(a,f9.6)') 'tau_A', tau_A
        A_relax=(V(j)/c_a)**(1./gam)
        dA=1.0/tau_A*(A_relax-A(j+1))
        dA_error=abs(dA*((tau_A_error)**2+(A_V_scaling_error)**2)**(1./2.))
!write(*,'(a,f15.6)') 'dA_error', dA_error
        A_error(j)=((A_error(j+1))**2+(dA_error)**2)**(1./2.)
!write(*,'(a,f15.6)') 'A_error', A_error(j)

        dV_cum_error(j)=((dV_cum_error(j+1))**2+(min(V(j),&
            abs(mb_modeled(j+1)*A(j)*1.e-6))* &
            ((mb_error(j)/max(10.,abs(mb_modeled(j+1))))**2+ &
            (A_error(j)/max(1.0e-6,A(j)))**2)**(1./2.))**2)**(1./2.)
!write(*,'(a,f9.6)') 'dV_cum_error', dV_cum_error(j+1)

    end do
    ! FSS--- shifting arrays as error model code is based on Ben's code. See
    ! description of subroutine mmbm_model
    dV_cum_error(2:size(dV_cum_error,1)) = dV_cum_error(1:size(dV_cum_error,1)-1)
    dV_cum_error(1) = NaN
end subroutine error_model

subroutine mmb_model(A_begin,meteo_years, meteo_months, precip_all, &
    temp_all_tongue, mb_months_base, &
    c_a,c_l,q,gam,zmax, zmin, lapse_temp,mu_gi,turnover,model_bias_gi,V,dL,L,A,mb_modeled, &
    mb_mod_mly,z_terminus, n_months_above_freezing)
!---------written by Felix Schueller (FSS)-----------------
! -INPUT:
!   A_begin : the area to begin with
!   meteo_years : years of driving meteorological input
!   meteo_months : months of driving meteorological input, i.e meteo_months and
!       precip and temp have to have the same size!
!   precip, temp : driving meteorological data
!   mb_months_base : decimal number of months from which to generate mass
!       balance (differs between southern and northern hemisphere)
!   c_a, c_l, q , gam, : glacier parameters
!   zmax, zmin : info about vertical dimension of glacier
!   lapse_temp
!   mu_gi : 
!   model_bias_gi
!   turnover 
! -OUTPUT:
!   V, dL, L, A, mb_modeled,mb_mod_mly
!   n_months_above_freezing
! -DESCRIPTION:
!   Area, Volume, Length will have one more point than years as a beginning data
!   point is added. dV, dA, dL and mass balance are set to the same length, but
!   with a NaN at the beginning. Mass balance of 1902 means mass balance between
!   Oct 1901 and Sept 1902 (northern hemisphere). i.e
!   year: 1901 1902 1903 ...
!   Area: 0.3   0.4  0.3 ...
!   mb  :  NaN  300  210 ...
! -TODO:
! -Last modified:  Mon Nov 10, 2014  11:48
! @author Felix Schueller
!--------------------------------------------------------
    integer, dimension(:), intent(in) :: meteo_years
    real, intent(in) :: c_a,c_l,q,gam, zmax, zmin, lapse_temp, mu_gi
    real, intent(in) :: A_begin
    real, intent(in) :: turnover, model_bias_gi
    real, dimension(:), intent(in) :: meteo_months
    real, dimension(:), intent(in) :: precip_all, temp_all_tongue
    real, dimension(:), intent(in) :: mb_months_base
    real, dimension(:), allocatable, intent(out) :: V, dL, L, A, mb_modeled
    real, dimension(:), allocatable, intent(out) :: mb_mod_mly, z_terminus
    integer, dimension(:), allocatable, intent(out) :: n_months_above_freezing

    real, dimension(12) :: mb_months
    integer, dimension(12) :: tmp_idx
    real :: turnover2
    
    real, dimension(:), allocatable :: precip, temp
    integer, dimension(:), allocatable :: months_idx
    
    integer :: sim_years, sim_months
    integer :: j,i
    real :: delta_z, delta_T_z
    real :: tau_A, tau_L
    integer :: last_idx
    
    real :: nan
    ! FSS--- get a NaN variable
    nan = 0
    nan = nan/nan

    sim_years = size(meteo_years)
    sim_months= size(meteo_months)

    ! FSS--- allocate arrays 
    if (allocated(V)) deallocate(V)
    if (allocated(dL)) deallocate(dL)
    if (allocated(L)) deallocate(L)
    if (allocated(A)) deallocate(A)
    if (allocated(mb_modeled)) deallocate(mb_modeled)
    if (allocated(mb_mod_mly)) deallocate(mb_mod_mly)
    if (allocated(z_terminus)) deallocate(z_terminus)
    allocate(V(sim_years+1))
    allocate(dL(sim_years+1))
    allocate(L(sim_years+1))
    allocate(A(sim_years+1))
    allocate(z_terminus(sim_years+1))
    allocate(mb_modeled(sim_years+1)) ! will have a NaN at the beginning 
    allocate(mb_mod_mly(sim_months)) ! will have a NaN at the beginning 
    allocate(n_months_above_freezing(sim_years+1))! will have a NaN at the beginning 
    mb_modeled = nan
    mb_mod_mly = nan
    n_months_above_freezing = -99 
    A = A_begin
    z_terminus = zmin

    !FSS TODO: check if V(1) is correct
    V(1)=c_a * (A(1)**gam)
    dL=0.0
    L(1)=(V(1)/c_l)**(1.0/q)
!print*,"A(1)", A(1)
!print*,"V(1)", V(1)
!print*,'L(1)', L(1)
!print*,'lapse_temp', lapse_temp
!
    !do j = 1,gi%year(i)-minval(meteo_years)
    do j = 1,sim_years

        mb_months= mb_months_base + (meteo_years(j)+1)

        tmp_idx = 0
        last_idx = -1
        do i=1,12
            tmp_idx(i)=minloc(abs(meteo_months-mb_months(i)),1)
            if (i>1 .and. last_idx == tmp_idx(i)) then
                !print*,"uh oh, partial yaer"
                exit 
            else
                last_idx = tmp_idx(i)
            end if 
        end do
        if (allocated(precip)) deallocate(precip)
        if (allocated(temp)) deallocate(temp)
        if (allocated(months_idx)) deallocate(months_idx)
        allocate(precip(i-1),temp(i-1),months_idx(i-1))
        months_idx = tmp_idx(:i-1)
        
        !print*,nint((meteo_months(1)-int(meteo_months(1))+1./24.) *12)
        
        ! FSS--- TODO: really V(1) ?
        ! properties don't change ?
        delta_z=((zmax-zmin) /  ((V(1)/c_l)**(1./q)))*(((V(1)/c_l)**(1.0/q))-L(j))
                                !|   L0                       L0                  |
        z_terminus(j) = zmin + delta_z
        delta_T_z= -delta_z*lapse_temp
        precip=precip_all(months_idx)
        temp=temp_all_tongue(months_idx)
        !print*,temp

        ! FSS--- beware the +1 in mb_modeled! this differs from original code. I
        ! added it to be able to have arrays with same lengths
        precip=precip_all(months_idx)
        mb_modeled(j+1)=sum(precip) - mu_gi * &
            (sum(max(0.0,delta_T_z+temp-t_melt)))-model_bias_gi
!write(*,'(a,i5,a,f8.2,a,f8.2,a,f8.4,a,f8.4,a,f8.4,a,f8.4)') "year",meteo_years(j),' p', &
    !sum(precip), &
    !' mu_gi', mu_gi,' delta_T_z',delta_T_z, ' temp',(sum(max(0.0,delta_T_z+temp-t_melt))),' model_bias_gi', &
    !model_bias_gi, ' mb', mb_modeled(j+1) 
           
!write(*,'(a,i5,a,f8.2,a,f8.2,f8.2,a,f8.4,a,f8.4)') "year",meteo_years(j),' zmin',zmin, &
!    ' deltaZ', delta_z, zmax,' delta_T_z',delta_T_z, ' mb', mb_modeled(j+1) 

        ! FSS--- monthly mass balances
        do i = 1,size(months_idx)
            !mb_mod_mly(months_idx(i)) = sum(precip(1:i)) - mu_gi * &
            !    (sum(max(0.0,delta_T_z+temp(1:i)-t_melt)))- (i * model_bias_gi/12.)
            mb_mod_mly(months_idx(i)) = precip(i) - mu_gi * &
                (max(0.0,delta_T_z+temp(i)-t_melt))- (model_bias_gi/12.)
        end do
        
        ! FSS--- n_months_above_freezing is used in error computation
        n_months_above_freezing(j+1)=count(delta_T_z+temp-t_melt>0)

        ! FSS--- beware the +1 in mb_modeled! this differs from original code. I
        ! added it to be able to have arrays with same lengths
        ! BIG TODO V/A etc are still in km ...
        ! * 1000 is to convert mb [m] in mb [mm]
        V(j+1)=max(0.0,V(j)+mb_modeled(j+1)*1000.0*A(j)*1.0e-6)
        turnover2 = turnover * 1000.0
        tau_L=max(1.0,V(j)/(turnover2*A(j)*1.0e-6)) !convert A to m
        tau_A=max(1.0,tau_L*(c_l**(2.0/q))/(c_a**(1.0/gam))*(V(j)**(1.0/gam-2.0/q)))

        L(j+1)=max(0.0,L(j)+1/tau_L*(((V(j+1)/c_l)**(1.0/q))-L(j)))
        dL(j+1)=1.0/tau_L*(((V(j+1)/c_l)**(1./q))-L(j))


        A(j+1)=max(0.0,A(j)+1.0/tau_A*(((V(j+1)/c_a)**(1.0/gam))-A(j)))
        
        if (abs(V(j+1)) < 1.0e-12) then
            L(j+1)=0.
            dL(j+1)=0.
            A(j+1)=0.
        end if

    end do ! end modeling
end subroutine mmb_model

subroutine calc_mu_rgi(l_o_c,t_star_rgi,meteo_years,T_tongue,Prec,nan_idx,mu_rgi,turnover)
!---------written by Felix Schueller (FSS)-----------------
! -INPUT:
!   l_o_c : length_of_climatology
!   t_star_rgi
!   meteo_years
!   Temp at glacier tongue
!   precipation
!   nan_idx 
! -OUTPUT:
!   nan_idx 
!   mu_rgi
!   turnover
! -DESCRIPTION: calculates climatologies within length_of_climatology window 
!   if lower or upper margin are hit, less years than l_o_c are taken
! @author Felix Schueller
!--------------------------------------------------------
    integer, intent(in) :: l_o_c
    integer, intent(in) :: t_star_rgi
    integer, dimension(:), intent(in) :: meteo_years
    real, dimension(:), intent(in) :: T_tongue
    real, dimension(:), intent(in) :: Prec 
    integer, intent(inout) :: nan_idx
    real, intent(out) :: mu_rgi
    real, intent(out) :: turnover

    integer :: t_star_rgi_idx
    integer :: j,i
    integer :: lowbound, upbound
    real, dimension(:), allocatable :: precip, temp,temp_tmp, precip_tmp
    real, dimension(12) :: t_tongue_clim, p_clim
    
    ! FSS--- get index of rgi year in meteo_years
    t_star_rgi_idx = -99
    do j = 1,size(meteo_years,1)
        if (meteo_years(j) == t_star_rgi) then
            t_star_rgi_idx = j
            exit
        end if
    end do 

    ! FSS--- throw error if not match was found
    if (t_star_rgi_idx == -99) then
        print*,'T_star rgi',t_star_rgi
        print*, 'ERROR: T_star_rgi is outside of meteo_years STOPPING'
        stop
    end if
    ! FSS--- lower boundary should not go below index 1
    lowbound = max(1,((t_star_rgi_idx-1)*12+1) - 12* (l_o_c -1)/2)
    ! FSS--- upper boundary not outside of length of data array
    upbound = min(size(T_tongue,1),(t_star_rgi_idx+(l_o_c-1)/2)*12)

    allocate(temp( (upbound - lowbound) +1))
    allocate(precip( (upbound - lowbound) +1))


    ! select data from time window
    ! get months from l.._o_clima.. number of years 
    temp=T_tongue(lowbound:upbound)
    precip = Prec(lowbound:upbound) 

    ! calculate monthly climatologies
    do i=1,12
        ! FSS--- extract matching months and remove nans
        call remove_nan(temp([(i+j, j = 0,size(temp,1)-1,12)]),temp_tmp) 
        call remove_nan(precip([(i+j, j = 0,size(precip,1)-1,12)]),precip_tmp) 
        ! FSS--- mean
        t_tongue_clim(i) = sum(temp_tmp)/size(temp_tmp,1)
        p_clim(i) = sum(precip_tmp)/size(precip_tmp,1)
    end do

    ! calculate model parameter
    mu_rgi = sum(p_clim)/sum(max(0.0,t_tongue_clim-t_melt))
    if (mu_rgi /= mu_rgi) then
        nan_idx = 1
    end if

    turnover=max(0.01,sum(p_clim))

    deallocate(precip, temp)

end subroutine calc_mu_rgi

subroutine remove_nan(ain,aout)
    real, dimension(:), intent(in) :: ain
    real, dimension(:), allocatable :: aout
    logical, dimension(size(ain)) :: nan_mask
    integer :: i,j
    
    nan_mask = ain == ain ! true for real, false for nan
    if (allocated(aout)) deallocate(aout)
    allocate( aout(count(nan_mask)) )
    j = 1
    do i = 1,size(ain)
        if (nan_mask(i)) then
            aout(j) = ain(i)
            j = j + 1
        end if
    end do
    
end subroutine remove_nan


end module mmbm
