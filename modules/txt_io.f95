module txt_io
implicit none
!---------written by Felix Schueller (FSS)-----------------
! -DESCRIPTION: subroutines concerning RGI IO 
! -Last modified:  Tue Jun 17, 2014  16:22
! @author Felix Schueller
!--------------------------------------------------------
type :: rgi_region
    character(len=60) :: region_name
    integer :: region
    integer, dimension(:,:), allocatable :: info
    real, dimension(:), allocatable :: lat, lon
    real, dimension(:,:), allocatable :: props 
    integer, dimension(:), allocatable :: nan_idx,year
    
end type rgi_region

type :: rgi_in
    type(rgi_in), pointer :: next => null()
    type(rgi_in), pointer :: prev => null()
    integer :: id
    real :: lat, lon, area
    real :: zmax, zmin, zmean, zstd
    integer :: is
    integer :: year
end type rgi_in
    
public :: get_rgi

contains

subroutine get_rgi(infile,reg,rgi)
!---------written by Felix Schueller (FSS)-----------------
! -INPUT:
! -OUTPUT:
! -DESCRIPTION: get Randolph Glacier inventory from txt file
! -Last modified:  Tue Nov 25, 2014  11:19
! @author Felix Schueller
!--------------------------------------------------------
    character(len=*), intent(in) :: infile
    integer, intent(in) :: reg
    type(rgi_region), intent(out) :: rgi
    type(rgi_in),pointer :: fi 
    type(rgi_in),pointer :: cu
    character(len=40) :: line
    integer :: n,ios,i
    allocate(fi)
    cu => fi
    open(unit=42,file=infile,status='old')
    read(42,*) !header line
    n = 0 
    do 
        read(42,*,iostat=ios) line,cu%id,cu%lon,cu%lat,cu%area,cu%zmin,cu%zmax, &
            cu%zmean,cu%zstd,cu%is,cu%year
        if (ios /= 0) exit
        n = n+1
        allocate(cu%next)
        cu => cu%next
        cu%next => null()
    end do
    close(42)

    allocate(rgi%lat(n),rgi%lon(n))
    allocate(rgi%props(n,5))
    allocate(rgi%info(n,2))
    allocate(rgi%year(n))
    allocate(rgi%nan_idx(n))
    rgi%nan_idx = 0

    rgi%region = reg
    
    write(*,*) 'RGI has ', n,'glaciers'
    ! FSS--- convert all to an array (pointers would work too...)
    cu => fi
    do i = 1,n
        rgi%info(i,1) = cu%id
        rgi%info(i,2) = cu%is
        rgi%year(i) = cu%year
        rgi%lat(i) = cu%lat
        rgi%lon(i) = cu%lon
        rgi%props(i,1) = cu%area
        rgi%props(i,2) = cu%zmin
        rgi%props(i,3) = cu%zmax
        rgi%props(i,4) = cu%zmean
        rgi%props(i,5) = cu%zstd
        cu => cu%next
    end do
    
    ! FSS--- get NaN's in properties
    do i = 1 , size(rgi%lat,1)
        if (any(rgi%props(i,:) /= rgi%props(i,:))) then
            !print*,i, "NaN"
            rgi%nan_idx(i) = 1
        end if
    end do 


    ! FSS--- tidy up
    cu => fi
    do while (associated(cu%next))
        fi => cu%next
        deallocate(cu)
        cu => fi
    end do 

end subroutine get_rgi
end module txt_io










