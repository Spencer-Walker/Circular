function forward_atomic(n,l,m,lmax,Ns)
  implicit none
  integer, intent(in)  :: n, l, m, lmax, Ns
  integer :: forward_atomic

  forward_atomic = Ns*( l**2 + l + m) + n

end function forward_atomic

function forward(n,l,m,NL,Ns,mmax) 
  implicit none
  integer, intent(in)  :: n,l,m,Nl,Ns,mmax
  integer :: forward

  forward = (mmax+m)*Nl*Ns + (l-abs(m))*Ns + n

end function forward

function reverse(i,Nl,Ns,mmax) 
  implicit none
  integer, intent(in)  :: Nl,Ns,mmax,i
  integer :: reverse(3), n,l,m

  m = floor(dble(i)/dble(Nl*Ns))-mmax
  l = floor(dble(i-(mmax+m)*Nl*Ns)/Ns)+abs(m)
  n = i - (mmax+m)*Nl*Ns - (l-abs(m))*Ns
  reverse(1) = n
  reverse(2) = l
  reverse(3) = m
end function reverse

function ThreeJ(l1,l2,l3,m1,m2,m3)
  use fgsl
  use iso_c_binding
  implicit none 
  integer, intent(in) :: l1, l2, l3, m1, m2, m3
  real(fgsl_double) :: ThreeJ 
  ThreeJ = fgsl_sf_coupling_3j(2*l1,2*l2,2*l3,2*m1,2*m2,2*m3)
  return 
end function ThreeJ

program main
#include <slepc/finclude/slepceps.h>
use slepceps
use hdf5
use mpi
use fgsl
use iso_c_binding
  implicit none

interface 
  function reverse(i,Nl,Ns,mmax) 
    implicit none
    integer, intent(in)  :: Nl,Ns,mmax,i
    integer :: reverse(3), n,l,m
  end function reverse
end interface

  integer(HID_T) :: memtype, eps_group_id, param_file_id, eps_dat_id
  integer(SIZE_T), parameter :: sdim = 300 
  integer(HSIZE_T)  :: dims(1)
  integer :: forward, forward_atomic, ret(3)
  PetscInt, parameter :: dp = kind(1.d0)
  PetscInt :: Ns, nev, Nl, tmp_int, h5_err, num_proc, proc_id, Na
  PetscInt :: Istart, Iend, n, l, m, Ntot, Ntot_atomic, mmax, i, j, row(1), lmax
  PetscInt, allocatable :: col(:),ix(:), ix_insert(:)
  character(len = 15) :: strl ! file name w/o .h5
  character(len = 200)   :: file_name 
  character(len = 300)   :: tmp_character
  character(len = 6)    :: fmt  ! format descriptor
  PetscBool      :: rot_present 
  PetscReal,  parameter :: pi = 3.141592653589793238462643383279502884197169
  PetscReal :: rot_theta, start_time, end_time, kabs
  PetscReal :: angular 
  real(fgsl_double) :: ThreeJ
  PetscErrorCode :: ierr
  PetscScalar :: k, norm
  PetscScalar, allocatable :: val(:), Lx_basis(:,:), Z_basis(:,:), overlap(:,:), insert(:)
  Mat            :: Lx, Z, S
  Mat, allocatable :: S_atomic(:) 
  Vec :: V, W 
  Vec, allocatable :: U(:,:), US(:)
  PetscViewer :: viewer
  
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !     Beginning of program
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  call CPU_TIME(start_time)

  call SlepcInitialize(PETSC_NULL_CHARACTER,ierr)
  if (ierr .ne. 0) then
    call PetscPrintf(MPI_COMM_WORLD, "SlepcInitialize failed\n", ierr)
    CHKERRA(ierr)
    stop
  end if
  call MPI_Comm_rank(PETSC_COMM_WORLD,proc_id,ierr);CHKERRA(ierr)
  call MPI_Comm_size(PETSC_COMM_WORLD,num_proc,ierr);CHKERRA(ierr)

  ! Initialize hdf5 
  call h5open_f( h5_err)
  if ( h5_err .ne. 0 ) then
    call PetscPrintf(MPI_COMM_WORLD, 'h5open_f failed\n', ierr);CHKERRA(ierr)
    stop
  end if

  call h5fopen_f( trim("parameters.h5"), H5F_ACC_RDONLY_F, param_file_id, h5_err)
  call h5gopen_f(param_file_id, "EPS", eps_group_id, h5_err)

  dims(1) = 1
  call H5Tcopy_f(H5T_FORTRAN_S1, memtype, h5_err)
  call H5Tset_size_f(memtype, sdim, h5_err)

  call h5dopen_f(eps_group_id, "rot_present", eps_dat_id, h5_err)
  call h5dread_f(eps_dat_id, H5T_NATIVE_INTEGER, tmp_int, dims, h5_err)
  call h5dclose_f( eps_dat_id, h5_err)  
  
  if (tmp_int .eq. 1) then
    rot_present = .true.
    call h5dopen_f(eps_group_id, "rot_theta", eps_dat_id, h5_err)
    call h5dread_f(eps_dat_id, H5T_NATIVE_DOUBLE, rot_theta, dims, h5_err)
    call h5dclose_f( eps_dat_id, h5_err)
  else
    rot_present = .false.
  end if
  if ( rot_present .eqv. .false.) then
    rot_theta = 0.d0
  end if

  call h5dopen_f(eps_group_id, "k", eps_dat_id, h5_err)
  call h5dread_f(eps_dat_id, H5T_NATIVE_DOUBLE, kabs, dims, h5_err)
  call h5dclose_f( eps_dat_id, h5_err)

  k = zexp(-(0d0,1d0)*rot_theta)*kabs

  dims(1) = 1
  call h5dopen_f(eps_group_id, "Nl", eps_dat_id, h5_err)
  call h5dread_f(eps_dat_id, H5T_NATIVE_INTEGER, Nl, dims, h5_err)
  call h5dclose_f( eps_dat_id, h5_err)

  call h5dopen_f(eps_group_id, "Ns", eps_dat_id, h5_err)
  call h5dread_f(eps_dat_id, H5T_NATIVE_INTEGER, Ns, dims, h5_err)
  call h5dclose_f( eps_dat_id, h5_err)

  call h5dopen_f(eps_group_id, "Na", eps_dat_id, h5_err)
  call h5dread_f(eps_dat_id, H5T_NATIVE_INTEGER, Na, dims, h5_err)
  call h5dclose_f( eps_dat_id, h5_err)

  call h5dopen_f(eps_group_id, "nev", eps_dat_id, h5_err)
  call h5dread_f(eps_dat_id, H5T_NATIVE_INTEGER, nev, dims, h5_err)
  call h5dclose_f( eps_dat_id, h5_err)

  call h5dopen_f(eps_group_id, "m_max", eps_dat_id, h5_err)
  call h5dread_f(eps_dat_id, H5T_NATIVE_INTEGER, mmax, dims, h5_err)
  call h5dclose_f( eps_dat_id, h5_err)

  lmax = Nl - 1 + mmax 

  Ntot = forward(Ns-1,NL-1+abs(mmax),mmax,Nl,Ns,mmax)+1
  Ntot_atomic  = Ns

  call MatCreate(PETSC_COMM_WORLD,S,ierr);CHKERRA(ierr)
  call MatSetSizes(S,PETSC_DECIDE,PETSC_DECIDE,Ntot,Ntot,ierr);CHKERRA(ierr)
  call MatSetFromOptions(S,ierr);CHKERRA(ierr)
  call MatSetUp(S,ierr);CHKERRA(ierr)
  
  call MatCreate(PETSC_COMM_WORLD,Lx,ierr);CHKERRA(ierr)
  call MatSetSizes(Lx,PETSC_DECIDE,PETSC_DECIDE,Ntot,Ntot,ierr);CHKERRA(ierr)
  call MatSetFromOptions(Lx,ierr);CHKERRA(ierr)
  call MatSetUp(Lx,ierr);CHKERRA(ierr)

  call MatCreate(PETSC_COMM_WORLD,Z,ierr);CHKERRA(ierr)
  call MatSetSizes(Z,PETSC_DECIDE,PETSC_DECIDE,Ntot,Ntot,ierr);CHKERRA(ierr)
  call MatSetFromOptions(Z,ierr);CHKERRA(ierr)
  call MatSetUp(Z,ierr);CHKERRA(ierr)

!-----------------------------------------------------------------------------------------
! Build Overlap Matrix
!-----------------------------------------------------------------------------------------

  call MatGetOwnershipRange(S,Istart,Iend,ierr);CHKERRA(ierr)
  allocate(col(Ns),val(Ns))
  do i = Istart,Iend-1
    ret = reverse(i,Nl,Ns,mmax)
    n = ret(1)
    l = ret(2)
    m = ret(3) 
    if (n .eq. 0) then 
      row(1) = i
      col(1) = i
      col(2) = i+1
      val(1) = 1d0
      val(2) = -0.5d0*dsqrt(dble(n+2*l+2)/dble(n+l+2))*dsqrt(dble(n+1)/dble(n+l+1))
      call MatSetValues(S,1,row,2,col,val,ADD_VALUES,ierr);CHKERRA(ierr)
    else if (n .eq. Ns-1) then
      row(1) = i
      col(1) = i-1
      col(2) = i
      val(1) = -0.5d0*dsqrt(dble(n+2*l+1)/dble(n+l+1))*dsqrt(dble(n)/dble(n+l))
      val(2) = 1d0
      call MatSetValues(S,1,row,2,col,val,ADD_VALUES,ierr);CHKERRA(ierr)
    else 
      row(1) = i
      col(1) = i - 1
      col(2) = i
      col(3) = i + 1
      val(1) = -0.5d0*dsqrt(dble(n)/dble(n+l))*dsqrt(dble(n+2*l+1)/dble(n+l+1))
      val(2) = 1d0
      val(3) = -0.5d0*dsqrt(dble(n+2*l+2)/dble(n+l+2))*dsqrt(dble(n+1)/dble(n+l+1))
      call MatSetValues(S,1,row,3,col,val,ADD_VALUES,ierr);CHKERRA(ierr)
    end if
  end do 

  call MatAssemblyBegin(S,MAT_FINAL_ASSEMBLY,ierr);CHKERRA(ierr)
  call MatAssemblyEnd(S,MAT_FINAL_ASSEMBLY,ierr);CHKERRA(ierr)


!-----------------------------------------------------------------------------------------
! Build Z Matrix
!-----------------------------------------------------------------------------------------
  call MatGetOwnershipRange(Z,Istart,Iend,ierr);CHKERRA(ierr)
  do i = Istart,Iend-1
    ret = reverse(i,Nl,Ns,mmax)
    n = ret(1)
    l = ret(2)
    m = ret(3) 
    row(1) = i
    if (l .lt. abs(m)+Nl-1) then 
      angular = (-1d0)**dble(m)*dsqrt(dble(2*l+1)*dble(2*l+3))*ThreeJ(l,1,l+1,0,0,0)*ThreeJ(l,1,l+1,-m,0,m)
      if (n .eq. 0) then
        col(1) = forward(n,l+1,m,Nl,Ns,mmax)
        col(2) = forward(n+1,l+1,m,Nl,Ns,mmax)

        val(1) = 0.5d0*k**(-1d0)*dble(2*n+l+2)*dsqrt(dble((n+2*l+2)*(n+2*l+3))/dble((n+l+1)*(n+l+2)))
        val(2) = -0.25d0*k**(-1d0)*dsqrt(dble((n+1)*(n+2*l+4)*(n+2*l+3)*(n+2*l+2))/dble((n+l+1)*(n+l+3)))

        val = angular*val
        call MatSetValues(Z,1,row,2,col,val,ADD_VALUES,ierr);CHKERRA(ierr)
        call MatSetValues(Z,2,col,1,row,val,ADD_VALUES,ierr);CHKERRA(ierr)

      else if (n .eq. 1) then 
        col(1) = forward(n-1,l+1,m,Nl,Ns,mmax)
        col(2) = forward(n,l+1,m,Nl,Ns,mmax)
        col(3) = forward(n+1,l+1,m,Nl,Ns,mmax)

        val(1) = -1.5d0*k**(-1d0)*dsqrt(dble(n*(n+2*l+2)))
        val(2) = 0.5d0*k**(-1d0)*dble(2*n+l+2)*dsqrt(dble((n+2*l+2)*(n+2*l+3))/dble((n+l+1)*(n+l+2)))
        val(3) = -0.25d0*k**(-1d0)*dsqrt(dble((n+1)*(n+2*l+4)*(n+2*l+3)*(n+2*l+2))/dble((n+l+1)*(n+l+3)))
        
        val = angular*val
        call MatSetValues(Z,1,row,3,col,val,ADD_VALUES,ierr);CHKERRA(ierr)
        call MatSetValues(Z,3,col,1,row,val,ADD_VALUES,ierr);CHKERRA(ierr)

      else if (n .eq. 2) then
        col(1) = forward(n-2,l+1,m,Nl,Ns,mmax)
        col(2) = forward(n-1,l+1,m,Nl,Ns,mmax)
        col(3) = forward(n,l+1,m,Nl,Ns,mmax)
        col(4) = forward(n+1,l+1,m,Nl,Ns,mmax)

        val(1) = 0.5d0*k**(-1d0)*dble(2*n+3*l+2)*dsqrt(dble(n*(n-1))/dble((n+l+1)*(n+l)))
        val(2) = -1.5d0*k**(-1d0)*dsqrt(dble(n*(n+2*l+2)))
        val(3) = 0.5d0*k**(-1d0)*dble(2*n+l+2)*dsqrt(dble((n+2*l+2)*(n+2*l+3))/dble((n+l+1)*(n+l+2)))
        val(4) = -0.25d0*k**(-1d0)*dsqrt(dble((n+1)*(n+2*l+4)*(n+2*l+3)*(n+2*l+2))/dble((n+l+1)*(n+l+3)))

        val = angular*val
        call MatSetValues(Z,1,row,4,col,val,ADD_VALUES,ierr);CHKERRA(ierr)
        call MatSetValues(Z,4,col,1,row,val,ADD_VALUES,ierr);CHKERRA(ierr)

      else if (n .eq. Ns-1) then 
        col(1) = forward(n-3,l+1,m,Nl,Ns,mmax)
        col(2) = forward(n-2,l+1,m,Nl,Ns,mmax)
        col(3) = forward(n-1,l+1,m,Nl,Ns,mmax)
        col(4) = forward(n,l+1,m,Nl,Ns,mmax)

        val(1) = -0.25d0*k**(-1d0)*dsqrt(dble(n*(n-1)*(n-2)*(n+2*l+1))/dble((n+l+1)*(n+l-1)))
        val(2) = 0.5d0*k**(-1d0)*dble(2*n+3*l+2)*dsqrt(dble(n*(n-1))/dble((n+l+1)*(n+l)))
        val(3) = -1.5d0*k**(-1d0)*dsqrt(dble(n*(n+2*l+2)))
        val(4) = 0.5d0*k**(-1d0)*dble(2*n+l+2)*dsqrt(dble((n+2*l+2)*(n+2*l+3))/dble((n+l+1)*(n+l+2)))

        val = angular*val
        call MatSetValues(Z,1,row,4,col,val,ADD_VALUES,ierr);CHKERRA(ierr)
        call MatSetValues(Z,4,col,1,row,val,ADD_VALUES,ierr);CHKERRA(ierr)

      else 
        col(1) = forward(n-3,l+1,m,Nl,Ns,mmax)
        col(2) = forward(n-2,l+1,m,Nl,Ns,mmax)
        col(3) = forward(n-1,l+1,m,Nl,Ns,mmax)
        col(4) = forward(n,l+1,m,Nl,Ns,mmax)
        col(5) = forward(n+1,l+1,m,Nl,Ns,mmax)

        val(1) = -0.25d0*k**(-1d0)*dsqrt(dble(n*(n-1)*(n-2)*(n+2*l+1))/dble((n+l+1)*(n+l-1)))
        val(2) = 0.5d0*k**(-1d0)*dble(2*n+3*l+2)*dsqrt(dble(n*(n-1))/dble((n+l+1)*(n+l)))
        val(3) = -1.5d0*k**(-1d0)*dsqrt(dble(n*(n+2*l+2)))
        val(4) = 0.5d0*k**(-1d0)*dble(2*n+l+2)*dsqrt(dble((n+2*l+2)*(n+2*l+3))/dble((n+l+1)*(n+l+2)))
        val(5) = -0.25d0*k**(-1d0)*dsqrt(dble((n+1)*(n+2*l+4)*(n+2*l+3)*(n+2*l+2))/dble((n+l+1)*(n+l+3)))

        val = angular*val 
        call MatSetValues(Z,1,row,5,col,val,ADD_VALUES,ierr);CHKERRA(ierr)
        call MatSetValues(Z,5,col,1,row,val,ADD_VALUES,ierr);CHKERRA(ierr)
      end if 
    end if 
  end do 

  call MatAssemblyBegin(Z,MAT_FINAL_ASSEMBLY,ierr);CHKERRA(ierr)
  call MatAssemblyEnd(Z,MAT_FINAL_ASSEMBLY,ierr);CHKERRA(ierr)


!-----------------------------------------------------------------------------------------
! Build Lx Matrix
!-----------------------------------------------------------------------------------------
  call MatGetOwnershipRange(Lx,Istart,Iend,ierr);CHKERRA(ierr)
  do i = Istart,Iend-1
    ret = reverse(i,Nl,Ns,mmax)
    n = ret(1)
    l = ret(2)
    m = ret(3) 

    if (m .lt. mmax .and. m .lt. l .and. l .lt. abs(m+1)+Nl) then 
      row(1) = i
      angular = dsqrt(dble((l+m+1)*(l-m)))

      if (n .eq. 0) then 
        col(1) = forward(n,l,m+1,Nl,Ns,mmax)
        col(2) = forward(n+1,l,m+1,Nl,Ns,mmax)

        val(1) = 1d0
        val(2) = -0.5d0*dsqrt(dble(n+2*l+2)/dble(n+l+2))*dsqrt(dble(n+1)/dble(n+l+1))

        val = angular*val
        call MatSetValues(Lx,1,row,2,col,val,ADD_VALUES,ierr);CHKERRA(ierr)
        call MatSetValues(Lx,2,col,1,row,val,ADD_VALUES,ierr);CHKERRA(ierr)

      else if (n .eq. Ns-1) then 
        col(1) = forward(n-1,l,m+1,Nl,Ns,mmax)
        col(2) = forward(n,l,m+1,Nl,Ns,mmax)
 
        val(1) = -0.5d0*dsqrt(dble(n+2*l+1)/dble(n+l+1))*dsqrt(dble(n)/dble(n+l))
        val(2) = 1d0

        val = angular*val
        call MatSetValues(Lx,1,row,2,col,val,ADD_VALUES,ierr);CHKERRA(ierr)
        call MatSetValues(Lx,2,col,1,row,val,ADD_VALUES,ierr);CHKERRA(ierr)

      else 
        col(1) = forward(n-1,l,m+1,Nl,Ns,mmax)
        col(2) = forward(n,l,m+1,Nl,Ns,mmax)
        col(3) = forward(n+1,l,m+1,Nl,Ns,mmax)

        val(1) = -0.5d0*dsqrt(dble(n)/dble(n+l))*dsqrt(dble(n+2*l+1)/dble(n+l+1))
        val(2) = 1d0
        val(3) = -0.5d0*dsqrt(dble(n+2*l+2)/dble(n+l+2))*dsqrt(dble(n+1)/dble(n+l+1))

        val = angular*val
        call MatSetValues(Lx,1,row,3,col,val,ADD_VALUES,ierr);CHKERRA(ierr)
        call MatSetValues(Lx,3,col,1,row,val,ADD_VALUES,ierr);CHKERRA(ierr)
      end if 
    end if 
  end do 

  deallocate(col,val)

  call MatAssemblyBegin(Lx,MAT_FINAL_ASSEMBLY,ierr);CHKERRA(ierr)
  call MatAssemblyEnd(Lx,MAT_FINAL_ASSEMBLY,ierr);CHKERRA(ierr)

!-----------------------------------------------------------------------------------------
! Vec_Load
!-----------------------------------------------------------------------------------------
  print*, 'Normalizing Stark States'
  file_name = 'stark_eigenvectors.bin'
  allocate(US(0:nev-1))
  call PetscViewerBinaryOpen(MPI_COMM_WORLD, file_name, FILE_MODE_READ, viewer, ierr);CHKERRA(ierr)

  call VecCreate(MPI_COMM_WORLD, US(0), ierr); CHKERRA(ierr)
  call VecLoad(US(0), viewer, ierr);CHKERRA(ierr)
  call VecDuplicate(US(0), V, ierr); CHKERRA(ierr)
  call VecCopy(US(0), V, ierr); CHKERRA(ierr)
  call MatMult(S, US(0), V, ierr);CHKERRA(ierr)
  call VecTDot( V, US(0), norm, ierr);CHKERRA(ierr)
  call VecScale(US(0), 1/zsqrt(norm), ierr);CHKERRA(ierr)
 
  do i = 1, nev-1
    call VecCreate(MPI_COMM_WORLD, US(i), ierr); CHKERRA(ierr)
    call VecLoad(US(i), viewer, ierr);CHKERRA(ierr)
    call VecCopy(US(i), V, ierr); CHKERRA(ierr)
    call MatMult(S, US(i), V, ierr);CHKERRA(ierr)
    call VecTDot( V, US(i), norm, ierr);CHKERRA(ierr)
    call VecScale(US(i), 1/zsqrt(norm), ierr);CHKERRA(ierr)
  end do
  call PetscViewerDestroy(viewer,ierr);CHKERRA(ierr)

  print*, 'Normalizing atomic States'
  allocate(U(0:lmax,0:Na-1))
  allocate(S_atomic(0:lmax))
  allocate(col(Ns))
  allocate(val(Ns))
  do l = 0, lmax
    print*, 'l = ', l
    if ( l .le. 9 ) then
      fmt = '(I1.1)'
    else if ( l .le. 99 ) then
      fmt = '(I2.2)'
    else
      fmt = '(I3.3)'
    end if

    write(strl,fmt) l
    call MatCreate(PETSC_COMM_WORLD,S_atomic(l),ierr);CHKERRA(ierr)
    call MatSetSizes(S_atomic(l),PETSC_DECIDE,PETSC_DECIDE,Ntot_atomic,Ntot_atomic,ierr);CHKERRA(ierr)
    call MatSetFromOptions(S_atomic(l),ierr);CHKERRA(ierr)
    call MatSetUp(S_atomic(l),ierr);CHKERRA(ierr)
    call MatGetOwnershipRange(S_atomic(l),Istart,Iend,ierr);CHKERRA(ierr)
    do i = Istart,Iend-1
      n = i 
      if (n .eq. 0) then
        row(1) = i
        col(1) = i
        col(2) = i+1
        val(1) = 1d0
        val(2) = -0.5d0*dsqrt(dble(n+2*l+2)/dble(n+l+2))*dsqrt(dble(n+1)/dble(n+l+1))
        call MatSetValues(S_atomic(l),1,row,2,col,val,ADD_VALUES,ierr);CHKERRA(ierr)
      else if (n .eq. Ns-1) then
        row(1) = i
        col(1) = i-1
        col(2) = i
        val(1) = -0.5d0*dsqrt(dble(n+2*l+1)/dble(n+l+1))*dsqrt(dble(n)/dble(n+l))
        val(2) = 1d0
        call MatSetValues(S_atomic(l),1,row,2,col,val,ADD_VALUES,ierr);CHKERRA(ierr)
      else
        row(1) = i
        col(1) = i - 1
        col(2) = i
        col(3) = i + 1
        val(1) = -0.5d0*dsqrt(dble(n)/dble(n+l))*dsqrt(dble(n+2*l+1)/dble(n+l+1))
        val(2) = 1d0
        val(3) = -0.5d0*dsqrt(dble(n+2*l+2)/dble(n+l+2))*dsqrt(dble(n+1)/dble(n+l+1))
        call MatSetValues(S_atomic(l),1,row,3,col,val,ADD_VALUES,ierr);CHKERRA(ierr)
      end if
    end do 
    
    call MatAssemblyBegin(S_atomic(l),MAT_FINAL_ASSEMBLY,ierr);CHKERRA(ierr)
    call MatAssemblyEnd(S_atomic(l),MAT_FINAL_ASSEMBLY,ierr);CHKERRA(ierr)

    file_name = 'atomic_eigenvectors_l'//trim(strl)//'.bin'
    call PetscViewerBinaryOpen(MPI_COMM_WORLD, file_name, FILE_MODE_READ, viewer, ierr);CHKERRA(ierr)

    i = 0       

    call VecCreate(MPI_COMM_WORLD, U(l,i), ierr); CHKERRA(ierr)
    call VecLoad(U(l,i), viewer, ierr);CHKERRA(ierr)
    if (l .eq. 0) then
      call VecDuplicate(U(l,i), V, ierr); CHKERRA(ierr)
    end if 
    call VecCopy(U(l,i), V, ierr); CHKERRA(ierr)
    call MatMult(S_atomic(l), U(l,i), V, ierr);CHKERRA(ierr)
    call VecTDot( V, U(l,i), norm, ierr);CHKERRA(ierr)
    call VecScale(U(l,i), 1/zsqrt(norm), ierr);CHKERRA(ierr)    

    do i = 1, Na-1
      call VecCreate(MPI_COMM_WORLD, U(l,i), ierr); CHKERRA(ierr)
      call VecLoad(U(l,i), viewer, ierr);CHKERRA(ierr)
      call VecCopy(U(l,i), V, ierr); CHKERRA(ierr)
      call MatMult(S_atomic(l), U(l,i), V, ierr);CHKERRA(ierr)
      call VecTDot( V, U(l,i), norm, ierr);CHKERRA(ierr)
      call VecScale(U(l,i), 1/zsqrt(norm), ierr);CHKERRA(ierr)
    end do
    call PetscViewerDestroy(viewer,ierr);CHKERRA(ierr)
  end do 
  deallocate(col)
  deallocate(val)

  print*, 'Saving stark states'
  file_name = 'stark_eigenvectors.h5'
  call PetscViewerHDF5Open(MPI_COMM_WORLD, file_name, FILE_MODE_WRITE, viewer, ierr);CHKERRA(ierr)
  do i = 0, nev-1
    call VecView(US(i),viewer,ierr);CHKERRA(ierr)
  end do
  call PetscViewerDestroy(viewer,ierr);CHKERRA(ierr)

  print*, 'Saving atomic states'
  do l = 0, lmax
    print*, 'l = ', l 
    if ( l .le. 9 ) then
      fmt = '(I1.1)'
    else if ( l .le. 99 ) then
      fmt = '(I2.2)'
    else
      fmt = '(I3.3)'
    end if

    write(strl,fmt) l
    file_name = 'atomic_eigenvectors_l'//trim(strl)//'.h5'
    call PetscViewerHDF5Open(MPI_COMM_WORLD, file_name, FILE_MODE_WRITE, viewer, ierr);CHKERRA(ierr)
    do i = 0, Na-1
      call VecView(U(l,i),viewer,ierr);CHKERRA(ierr)
    end do
    call PetscViewerDestroy(viewer,ierr);CHKERRA(ierr)
  end do
!-----------------------------------------------------------------------------------------
!  Stark - atomic overlap matrix
!-----------------------------------------------------------------------------------------


  allocate(ix(0:Ns-1))
  allocate(val(0:Ns-1))
  allocate(insert(0:Ns-1))
  allocate(ix_insert(0:Ns-1))
  do i = 0,Ns-1
    ix(i) = i
  end do

  print*, 'Saving stark-atomic overlap' 
  allocate(overlap(1,0:Na*(lmax+1)**2-1))
  overlap = (0d0,0d0)

  print*, Na*(lmax+1)**2-1

  open( unit = 10, file = "stark-atomic-overlap-real.csv")
  open( unit = 20, file = "stark-atomic-overlap-imag.csv")

  call VecDuplicate(US(0), V, ierr); CHKERRA(ierr)
  call VecDuplicate(US(0), W, ierr); CHKERRA(ierr)
  do i = 0, nev-1
    print*, 'stark i = ', i, 'out of', nev-1
    do m = -mmax, mmax
      do l = abs(m), abs(m) + Nl-1 
        do j = 0, Na-1
          call VecSet(V,(0d0,0d0),ierr);CHKERRA(ierr)
          call VecGetValues(U(l,j), Ns, ix, val, ierr);CHKERRA(ierr)
          ix_insert(:) = ix(:) + forward(0,l,m,Nl,Ns,mmax) 
          call VecSetValues(V, Ns, ix_insert, val, INSERT_VALUES, ierr)
          CHKERRA(ierr)
          call VecAssemblyBegin(V,ierr);CHKERRA(ierr)
          call VecAssemblyEnd(V,ierr);CHKERRA(ierr)
          call MatMult(S, V, W, ierr);CHKERRA(ierr)
          call VecTDot(W, US(i), norm, ierr);CHKERRA(ierr) 
          overlap(1,forward_atomic(j,l,m,lmax,Na)) = norm
        end do 
      end do 
    end do 
    write(10,*) (realpart(overlap(1,j)), j = 0,Na*(lmax+1)**2-1)
    write(20,*) (imagpart(overlap(1,j)), j = 0,Na*(lmax+1)**2-1)
  end do

  deallocate(overlap)
  deallocate(ix_insert)
  deallocate(val)
  deallocate(insert)
  deallocate(ix)

  print*, 'Saving stark overlap'
  open( unit = 10, file = "stark-overlap-real.csv")
  open( unit = 20, file = "stark-overlap-imag.csv")
  allocate(overlap(0:nev-1,0:nev-1))
  call VecDuplicate(US(0), V, ierr); CHKERRA(ierr)
  call VecCopy(US(0), V, ierr); CHKERRA(ierr)
  do j = 0, nev-1
    do i = 0, nev-1
        call MatMult(S, US(j), V, ierr);CHKERRA(ierr)
        call VecTDot( V, US(i), norm, ierr);CHKERRA(ierr)
        overlap(i,j) = norm
    end do
  end do

  do i = 0, nev-1
    write(10,*) (realpart(overlap(i,j)), j = 0,nev-1)
    write(20,*) (imagpart(overlap(i,j)), j = 0,nev-1)
  end do

  deallocate(overlap)



!-----------------------------------------------------------------------------------------
!  Stark observables
!-----------------------------------------------------------------------------------------
  print*, 'Saving stark observables'
  allocate(Z_basis(0:nev-1,0:nev-1), Lx_basis(0:nev-1,0:nev-1))
  open( unit = 10, file = "Z_stark-real.csv")
  open( unit = 20, file = "Z_stark-imag.csv")
  open( unit = 30, file = "Lx_stark-real.csv")
  open( unit = 40, file = "Lx_stark-imag.csv")
  do j = 0, nev-1
    do i = 0, nev-1
        call MatMult(Z, US(j), V, ierr);CHKERRA(ierr)
        call VecTDot( V, US(i), norm, ierr);CHKERRA(ierr)
        Z_basis(i,j) = norm

        call MatMult(Lx, US(j), V, ierr);CHKERRA(ierr)
        call VecTDot( V, US(i), norm, ierr);CHKERRA(ierr)
        Lx_basis(i,j) = norm
    end do
  end do

  do i = 0, nev-1
    write(10,*) (realpart(Z_basis(i,j)), j = 0,nev-1)
    write(20,*) (imagpart(Z_basis(i,j)), j = 0,nev-1)
    write(30,*) (realpart(Lx_basis(i,j)), j = 0,nev-1)
    write(40,*) (imagpart(Lx_basis(i,j)), j = 0,nev-1)
  end do 

  deallocate(Z_basis,Lx_basis)

  ! Closes the hdf5 file
  call h5gclose_f( eps_group_id, h5_err)
  call h5fclose_f( param_file_id, h5_err)
  call CPU_TIME(end_time)
  write(tmp_character, "(ES9.2)")  end_time-start_time
  call PetscPrintf(MPI_COMM_WORLD, 'time   :'//trim(tmp_character)//"\n", ierr)
  CHKERRA(ierr)
  call SlepcFinalize(ierr)


end program main 
