function forward(n,l,m,Nl,Ns,mmax) 
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
use forpy_mod
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
  integer :: forward, ret(3)
  PetscInt, parameter :: dp = kind(1.d0)
  PetscInt :: Ns, nev, ncv, mpd, maxits, Nl, tmp_int, h5_err, num_proc, proc_id
  PetscInt :: Istart, Iend, n, l, m, Nshell, Ntot, mmax, i, row(1), j, ierror
  PetscInt :: n1
  PetscInt, allocatable :: col(:)
  character(len = 24)   :: file_name 
  character(len = 300)   :: tmp_character
  PetscBool      :: rot_present, eps_two_sided 
  PetscReal,  parameter :: pi = 3.141592653589793238462643383279502884197169
  PetscReal :: rot_theta, tol, start_time, end_time, eps_target_components(2), kabs
  PetscReal :: C0, C, Zc, F, w, angular 
  PetscReal, allocatable ::  a(:), b(:)
  real(fgsl_double) :: ThreeJ
  PetscErrorCode :: ierr
  EPSProblemType :: eps_problem 
  PetscScalar :: eps_target, k, g, q, q2, q3 
  PetscScalar, allocatable :: val(:)
  Mat            :: H, S
  EPS            :: eps
  EPSType :: eps_type   
  EPSWhich :: eps_which
  EPSBalance :: eps_balance 
  PetscViewer :: viewer
  type(module_py) :: opt
  type(tuple) :: args
  type(object) :: retval  
  type(list) :: paths
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
  call h5dopen_f(eps_group_id, "EPSSetBalance", eps_dat_id, h5_err)
  call h5dread_f(eps_dat_id, memtype, tmp_character,dims, h5_err)
  call h5dclose_f( eps_dat_id, h5_err)
  if (trim(tmp_character) .eq. "EPS_BALANCE_TWOSIDE") then
    eps_balance = EPS_BALANCE_TWOSIDE
  else if (trim(tmp_character) .eq. "EPS_BALANCE_ONESIDE") then
    eps_balance = EPS_BALANCE_ONESIDE
  else if (trim(tmp_character) .eq. "EPS_BALANCE_NONE") then
    eps_balance = EPS_BALANCE_NONE
  else 
    call PetscPrintf(MPI_COMM_WORLD, "EPSBalance not supported defaulting to EPS_BALANCE_NONE\n", ierr)
    CHKERRA(ierr)
    eps_balance = EPS_BALANCE_NONE
  end if 
  
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

  call h5dopen_f(eps_group_id, "F", eps_dat_id, h5_err)
  call h5dread_f(eps_dat_id, H5T_NATIVE_DOUBLE, F, dims, h5_err)
  call h5dclose_f( eps_dat_id, h5_err)

  call h5dopen_f(eps_group_id, "w", eps_dat_id, h5_err)
  call h5dread_f(eps_dat_id, H5T_NATIVE_DOUBLE, w, dims, h5_err)
  call h5dclose_f( eps_dat_id, h5_err)

  call h5dopen_f(eps_group_id, "k", eps_dat_id, h5_err)
  call h5dread_f(eps_dat_id, H5T_NATIVE_DOUBLE, kabs, dims, h5_err)
  call h5dclose_f( eps_dat_id, h5_err)

  k = zexp(-(0d0,1d0)*rot_theta)*kabs

  call h5dopen_f(eps_group_id, "C0", eps_dat_id, h5_err)
  call h5dread_f(eps_dat_id, H5T_NATIVE_DOUBLE, C0, dims, h5_err)
  call h5dclose_f( eps_dat_id, h5_err)
  
  call h5dopen_f(eps_group_id, "C", eps_dat_id, h5_err)
  call h5dread_f(eps_dat_id, H5T_NATIVE_DOUBLE, C, dims, h5_err)
  call h5dclose_f( eps_dat_id, h5_err)

  call h5dopen_f(eps_group_id, "Zc", eps_dat_id, h5_err)
  call h5dread_f(eps_dat_id, H5T_NATIVE_DOUBLE, Zc, dims, h5_err)
  call h5dclose_f( eps_dat_id, h5_err)

  call h5dopen_f(eps_group_id, "Nshell", eps_dat_id, h5_err)
  call h5dread_f(eps_dat_id, H5T_NATIVE_INTEGER, Nshell, dims, h5_err)
  call h5dclose_f( eps_dat_id, h5_err)  

  allocate(a(Nshell),b(Nshell))
  dims(1) = Nshell
  call h5dopen_f(eps_group_id, "a", eps_dat_id, h5_err)
  call h5dread_f(eps_dat_id, H5T_NATIVE_DOUBLE, a, dims, h5_err)
  call h5dclose_f( eps_dat_id, h5_err)

  call h5dopen_f(eps_group_id, "b", eps_dat_id, h5_err)
  call h5dread_f(eps_dat_id, H5T_NATIVE_DOUBLE, b, dims, h5_err)
  call h5dclose_f( eps_dat_id, h5_err)


  dims(1) = 1
  call h5dopen_f(eps_group_id, "Nl", eps_dat_id, h5_err)
  call h5dread_f(eps_dat_id, H5T_NATIVE_INTEGER, Nl, dims, h5_err)
  call h5dclose_f( eps_dat_id, h5_err)

  call h5dopen_f(eps_group_id, "Ns", eps_dat_id, h5_err)
  call h5dread_f(eps_dat_id, H5T_NATIVE_INTEGER, Ns, dims, h5_err)
  call h5dclose_f( eps_dat_id, h5_err)

  call h5dopen_f(eps_group_id, "m_max", eps_dat_id, h5_err)
  call h5dread_f(eps_dat_id, H5T_NATIVE_INTEGER, mmax, dims, h5_err)
  call h5dclose_f( eps_dat_id, h5_err)

  call h5dopen_f(eps_group_id, "nev", eps_dat_id, h5_err)
  call h5dread_f(eps_dat_id, H5T_NATIVE_INTEGER, nev, dims, h5_err)
  call h5dclose_f( eps_dat_id, h5_err)

  call h5dopen_f(eps_group_id, "max_its", eps_dat_id, h5_err)
  call h5dread_f(eps_dat_id, H5T_NATIVE_INTEGER, maxits, dims, h5_err)
  call h5dclose_f( eps_dat_id, h5_err)

  call h5dopen_f(eps_group_id, "abs_tol", eps_dat_id, h5_err)
  call h5dread_f(eps_dat_id, H5T_NATIVE_DOUBLE, tol, dims, h5_err)
  call h5dclose_f( eps_dat_id, h5_err)

  call h5dopen_f(eps_group_id, "ncv", eps_dat_id, h5_err)
  call h5dread_f(eps_dat_id, H5T_NATIVE_INTEGER, ncv, dims, h5_err)
  call h5dclose_f( eps_dat_id, h5_err)
  if (ncv .eq. -1) then
    ncv = PETSC_DEFAULT_INTEGER
  end if 

  call h5dopen_f(eps_group_id, "mpd", eps_dat_id, h5_err)
  call h5dread_f(eps_dat_id, H5T_NATIVE_INTEGER, mpd, dims, h5_err)
  call h5dclose_f( eps_dat_id, h5_err)
  if (mpd .eq. -1) then
    mpd = PETSC_DEFAULT_INTEGER
  end if 

  call h5dopen_f(eps_group_id, "EPSSetProblemType", eps_dat_id, h5_err)
  call h5dread_f(eps_dat_id, memtype, tmp_character, dims, h5_err)
  call h5dclose_f( eps_dat_id, h5_err)
  if (trim(tmp_character) .eq. "EPS_NHEP" ) then
    eps_problem = EPS_NHEP
  else if (trim(tmp_character) .eq. "EPS_PGNHEP") then
    eps_problem = EPS_PGNHEP
  else if (trim(tmp_character) .eq. "EPS_HEP") then
    eps_problem = EPS_HEP
  else if (trim(tmp_character) .eq. "EPS_GHEP" ) then
    eps_problem = EPS_GHEP
  else if (trim(tmp_character) .eq. "EPS_GNHEP") then
    eps_problem = EPS_GNHEP
  else 
    call PetscPrintf(MPI_COMM_WORLD, "EPSProblemType not supported defaulting to EPS_NHEP\n", ierr)
    CHKERRA(ierr)
    eps_problem = EPS_NHEP
  end if 
  

  dims(1) = 2
  call h5dopen_f(eps_group_id, "EPSSetTarget", eps_dat_id, h5_err)
  call h5dread_f(eps_dat_id, H5T_NATIVE_DOUBLE, eps_target_components, dims, h5_err)
  call h5dclose_f( eps_dat_id, h5_err)
  eps_target = dcmplx(eps_target_components(1),eps_target_components(2))

  dims(1) = 1
  call h5dopen_f(eps_group_id, "EPSSetTwoSided", eps_dat_id, h5_err)
  call h5dread_f(eps_dat_id, H5T_NATIVE_INTEGER, tmp_int, dims, h5_err)
  call h5dclose_f( eps_dat_id, h5_err)
  if (tmp_int .eq. 1) then
    eps_two_sided = PETSC_TRUE
  else
    eps_two_sided = PETSC_FALSE
  end if 

  call h5dopen_f(eps_group_id, "EPSSetType", eps_dat_id, h5_err)
  call h5dread_f(eps_dat_id, memtype, tmp_character, dims, h5_err)
  call h5dclose_f( eps_dat_id, h5_err)
  if (trim(tmp_character) .eq. "EPSPOWER") then 
    eps_type = EPSPOWER
  else if (trim(tmp_character) .eq. "EPSPOWER") then 
    eps_type = EPSPOWER
  else if (trim(tmp_character) .eq. "EPSSUBSPACE") then 
    eps_type = EPSSUBSPACE    
  else if (trim(tmp_character) .eq. "EPSARNOLDI") then 
    eps_type = EPSARNOLDI
  else if (trim(tmp_character) .eq. "EPSLANCZOS") then 
    eps_type = EPSLANCZOS
  else if (trim(tmp_character) .eq. "EPSKRYLOVSCHUR") then 
    eps_type = EPSKRYLOVSCHUR
  else if (trim(tmp_character) .eq. "EPSGD") then 
    eps_type = EPSGD
  else if (trim(tmp_character) .eq. "EPSJD") then 
    eps_type = EPSJD
  else if (trim(tmp_character) .eq. "EPSRQCG") then 
    eps_type = EPSRQCG
  else if (trim(tmp_character) .eq. "EPSLOBPCG") then 
    eps_type = EPSLOBPCG
  else if (trim(tmp_character) .eq. "EPSCISS") then 
    eps_type = EPSCISS
  else if (trim(tmp_character) .eq. "EPSLAPACK") then 
    eps_type = EPSLAPACK
  else if (trim(tmp_character) .eq. "EPSARPACK") then 
    eps_type = EPSARPACK
  else if (trim(tmp_character) .eq. "EPSBLZPACK") then 
    eps_type = EPSBLZPACK
  else if (trim(tmp_character) .eq. "EPSTRLAN") then 
    eps_type = EPSTRLAN
  else if (trim(tmp_character) .eq. "EPSBLOPEX") then 
    eps_type = EPSBLOPEX
  else if (trim(tmp_character) .eq. "EPSPRIMME") then 
    eps_type = EPSPRIMME
  else if (trim(tmp_character) .eq. "EPSFEAST") then 
    eps_type = EPSFEAST
  else 
    call PetscPrintf(MPI_COMM_WORLD, "EPSType not supported defaulting to EPSKRYLOVSCHUR\n", ierr)
    CHKERRA(ierr)
    eps_type = EPSKRYLOVSCHUR
  end if 

  call h5dopen_f(eps_group_id, "EPSSetWhichEigenpairs", eps_dat_id, h5_err)
  call h5dread_f(eps_dat_id, memtype, tmp_character, dims, h5_err)
  call h5dclose_f( eps_dat_id, h5_err)
  if (trim(tmp_character) .eq. "EPS_LARGEST_MAGNITUDE") then 
    eps_which = EPS_LARGEST_MAGNITUDE 
  else if (trim(tmp_character) .eq. "EPS_SMALLEST_MAGNITUDE") then 
    eps_which = EPS_SMALLEST_MAGNITUDE
  else if (trim(tmp_character) .eq. "EPS_LARGEST_REAL") then 
    eps_which = EPS_LARGEST_REAL   
  else if (trim(tmp_character) .eq. "EPS_SMALLEST_REAL") then 
    eps_which = EPS_SMALLEST_REAL
  else if (trim(tmp_character) .eq. "EPS_LARGEST_IMAGINARY") then 
    eps_which = EPS_LARGEST_IMAGINARY
  else if (trim(tmp_character) .eq. "EPS_SMALLEST_IMAGINARY") then 
    eps_which = EPS_SMALLEST_IMAGINARY
  else if (trim(tmp_character) .eq. "EPS_TARGET_MAGNITUDE") then 
    eps_which = EPS_TARGET_MAGNITUDE
  else if (trim(tmp_character) .eq. "EPS_TARGET_REAL") then 
    eps_which = EPS_TARGET_REAL
  else if (trim(tmp_character) .eq. "EPS_TARGET_IMAGINARY") then 
    eps_which = EPS_TARGET_IMAGINARY
  else if (trim(tmp_character) .eq. "EPS_ALL") then 
    eps_which = EPS_ALL
  else
    call PetscPrintf(MPI_COMM_WORLD, "EPSWhich not supported defaulting to EPS_SMALLEST_REAL\n", ierr)
    CHKERRA(ierr)
    eps_which = EPS_SMALLEST_REAL
  end if 

  Ntot = forward(Ns-1,NL-1+abs(mmax),mmax,Nl,Ns,mmax)+1

  call MatCreate(PETSC_COMM_WORLD,S,ierr);CHKERRA(ierr)
  call MatSetSizes(S,PETSC_DECIDE,PETSC_DECIDE,Ntot,Ntot,ierr);CHKERRA(ierr)
  call MatSetFromOptions(S,ierr);CHKERRA(ierr)
  call MatSetUp(S,ierr);CHKERRA(ierr)
  
  call MatCreate(PETSC_COMM_WORLD,H,ierr);CHKERRA(ierr)
  call MatSetSizes(H,PETSC_DECIDE,PETSC_DECIDE,Ntot,Ntot,ierr);CHKERRA(ierr)
  call MatSetFromOptions(H,ierr);CHKERRA(ierr)
  call MatSetUp(H,ierr);CHKERRA(ierr)
  
  call MatGetOwnershipRange(H,Istart,Iend,ierr);CHKERRA(ierr)

!-----------------------------------------------------------------------------------------
! Build Overlap Matrix
!-----------------------------------------------------------------------------------------
  
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
! Build Field-Free Matrix
!-----------------------------------------------------------------------------------------
  do i = Istart,Iend-1
    ret = reverse(i,Nl,Ns,mmax)
    n = ret(1)
    l = ret(2)
    m = ret(3) 
    
    if (n .eq. 0) then 
      ! Main diagonal
      row(1) = i
      col(1) = i
      col(2) = i+1
      val(1) = 0.5d0*k**2d0 - k*C0/dble(n+l+1)
      val(2) = 0.25d0*k**2d0*dsqrt(dble(n+2*l+2)/dble(n+l+2))*dsqrt(dble(n+1)/dble(n+l+1))
      call MatSetValues(H,1,row,2,col,val,ADD_VALUES,ierr);CHKERRA(ierr)
  
      
    else if (n .eq. Ns-1) then
      row(1) = i
      col(1) = i-1
      col(2) = i
      val(1) = 0.25d0*k**2d0*dsqrt(dble(n)/dble(n+l))*dsqrt(dble(n+2*l+1)/dble(n+l+1))
      val(2) = 0.5d0*k**2d0 - k*C0/dble(n+l+1)
      call MatSetValues(H,1,row,2,col,val,ADD_VALUES,ierr);CHKERRA(ierr)
    else 
      row(1) = i
      col(1) = i - 1
      col(2) = i
      col(3) = i + 1
      val(1) = 0.25d0*k**2d0*dsqrt(dble(n)/dble(n+l))*dsqrt(dble(n+2*l+1)/dble(n+l+1))
      val(2) = 0.5d0*k**2d0 - k*C0/dble(n+l+1)
      val(3) = 0.25d0*k**2d0*dsqrt(dble(n+2*l+2)/dble(n+l+2))*dsqrt(dble(n+1)/dble(n+l+1))
      call MatSetValues(H,1,row,3,col,val,ADD_VALUES,ierr);CHKERRA(ierr)
    end if
  end do
  
  !Short range potential
  g = 2*(k/C)
  ierror = forpy_initialize()
  ierror = get_sys_path(paths)
  ierror = paths%append("/home/becker/spwa4419/Documents/Circular/")
  ierror = import_py(opt, "mymodule")
  if (dabs(Zc) .gt. 1d-15) then
    g = 2*(k/C)
    do i = Istart, Iend-1
      ret = reverse(i,Nl,Ns,mmax)
      n1 = ret(1)
      l = ret(2)
      m = ret(3)
      row(1)  = i
      do n = 0, Ns-1
        ierror = tuple_create(args, 6)
        ierror = args%setitem(0,n)
        ierror = args%setitem(1,n1)
        ierror = args%setitem(2,l)
        ierror = args%setitem(3,k)
        ierror = args%setitem(4,Zc)
        ierror = args%setitem(5,g)
        ierror = call_py(retval, opt, "yukawa", args)
        ierror = cast(q, retval)
        call args%destroy
        call retval%destroy
        col(n+1) = forward(n,l,m,Nl,Ns,mmax)
        val(n+1) = q
      end do
      call MatSetValues(H,1,row,Ns,col,val,ADD_VALUES,ierr);CHKERRA(ierr)
    end do
  end if 

  if (Nshell .ne. 0) then
    do j = 1, Nshell
      g = 2d0*(k/b(j))
      do i = Istart, Iend-1
        ret = reverse(i,Nl,Ns,mmax)
        n1 = ret(1)
        l = ret(2)
        m = ret(3)
        row(1)  = i
        val = 0
        do n = 0, Ns-1
          ierror = tuple_create(args, 5)
          ierror = args%setitem(0,n)
          ierror = args%setitem(1,n1)
          ierror = args%setitem(2,l)
          ierror = args%setitem(3,g)
          ierror = args%setitem(4,a(j))
          ierror = call_py(retval, opt, "shell", args)
          ierror = cast(q, retval)
          call args%destroy
          call retval%destroy
          col(n+1) =  forward(n,l,m,Nl,Ns,mmax)
          val(n+1) =  q
        end do
        call MatSetValues(H,1,row,Ns,col,val,ADD_VALUES,ierr);CHKERRA(ierr)
      end do
    end do
  end if 

  ! DC Stark Shifts
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

        val = F*angular*val
        call MatSetValues(H,1,row,2,col,val,ADD_VALUES,ierr);CHKERRA(ierr)
        call MatSetValues(H,2,col,1,row,val,ADD_VALUES,ierr);CHKERRA(ierr)

      else if (n .eq. 1) then 
        col(1) = forward(n-1,l+1,m,Nl,Ns,mmax)
        col(2) = forward(n,l+1,m,Nl,Ns,mmax)
        col(3) = forward(n+1,l+1,m,Nl,Ns,mmax)

        val(1) = -1.5d0*k**(-1d0)*dsqrt(dble(n*(n+2*l+2)))
        val(2) = 0.5d0*k**(-1d0)*dble(2*n+l+2)*dsqrt(dble((n+2*l+2)*(n+2*l+3))/dble((n+l+1)*(n+l+2)))
        val(3) = -0.25d0*k**(-1d0)*dsqrt(dble((n+1)*(n+2*l+4)*(n+2*l+3)*(n+2*l+2))/dble((n+l+1)*(n+l+3)))
        
        val = F*angular*val
        call MatSetValues(H,1,row,3,col,val,ADD_VALUES,ierr);CHKERRA(ierr)
        call MatSetValues(H,3,col,1,row,val,ADD_VALUES,ierr);CHKERRA(ierr)

      else if (n .eq. 2) then
        col(1) = forward(n-2,l+1,m,Nl,Ns,mmax)
        col(2) = forward(n-1,l+1,m,Nl,Ns,mmax)
        col(3) = forward(n,l+1,m,Nl,Ns,mmax)
        col(4) = forward(n+1,l+1,m,Nl,Ns,mmax)

        val(1) = 0.5d0*k**(-1d0)*dble(2*n+3*l+2)*dsqrt(dble(n*(n-1))/dble((n+l+1)*(n+l)))
        val(2) = -1.5d0*k**(-1d0)*dsqrt(dble(n*(n+2*l+2)))
        val(3) = 0.5d0*k**(-1d0)*dble(2*n+l+2)*dsqrt(dble((n+2*l+2)*(n+2*l+3))/dble((n+l+1)*(n+l+2)))
        val(4) = -0.25d0*k**(-1d0)*dsqrt(dble((n+1)*(n+2*l+4)*(n+2*l+3)*(n+2*l+2))/dble((n+l+1)*(n+l+3)))

        val = F*angular*val
        call MatSetValues(H,1,row,4,col,val,ADD_VALUES,ierr);CHKERRA(ierr)
        call MatSetValues(H,4,col,1,row,val,ADD_VALUES,ierr);CHKERRA(ierr)

      else if (n .eq. Ns-1) then 
        col(1) = forward(n-3,l+1,m,Nl,Ns,mmax)
        col(2) = forward(n-2,l+1,m,Nl,Ns,mmax)
        col(3) = forward(n-1,l+1,m,Nl,Ns,mmax)
        col(4) = forward(n,l+1,m,Nl,Ns,mmax)

        val(1) = -0.25d0*k**(-1d0)*dsqrt(dble(n*(n-1)*(n-2)*(n+2*l+1))/dble((n+l+1)*(n+l-1)))
        val(2) = 0.5d0*k**(-1d0)*dble(2*n+3*l+2)*dsqrt(dble(n*(n-1))/dble((n+l+1)*(n+l)))
        val(3) = -1.5d0*k**(-1d0)*dsqrt(dble(n*(n+2*l+2)))
        val(4) = 0.5d0*k**(-1d0)*dble(2*n+l+2)*dsqrt(dble((n+2*l+2)*(n+2*l+3))/dble((n+l+1)*(n+l+2)))

        val = F*angular*val
        call MatSetValues(H,1,row,4,col,val,ADD_VALUES,ierr);CHKERRA(ierr)
        call MatSetValues(H,4,col,1,row,val,ADD_VALUES,ierr);CHKERRA(ierr)

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

        val = F*angular*val 
        call MatSetValues(H,1,row,5,col,val,ADD_VALUES,ierr);CHKERRA(ierr)
        call MatSetValues(H,5,col,1,row,val,ADD_VALUES,ierr);CHKERRA(ierr)
      end if 
    end if 
  end do 

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

        val = -0.5d0*w*angular*val
        call MatSetValues(H,1,row,2,col,val,ADD_VALUES,ierr);CHKERRA(ierr)
        call MatSetValues(H,2,col,1,row,val,ADD_VALUES,ierr);CHKERRA(ierr)

      else if (n .eq. Ns-1) then 
        col(1) = forward(n-1,l,m+1,Nl,Ns,mmax)
        col(2) = forward(n,l,m+1,Nl,Ns,mmax)

        val(1) = -0.5d0*dsqrt(dble(n+2*l+1)/dble(n+l+1))*dsqrt(dble(n)/dble(n+l))
        val(2) = 1d0

        val = -0.5d0*w*angular*val
        call MatSetValues(H,1,row,2,col,val,ADD_VALUES,ierr);CHKERRA(ierr)
        call MatSetValues(H,2,col,1,row,val,ADD_VALUES,ierr);CHKERRA(ierr)

      else 
        col(1) = forward(n-1,l,m+1,Nl,Ns,mmax)
        col(2) = forward(n,l,m+1,Nl,Ns,mmax)
        col(3) = forward(n+1,l,m+1,Nl,Ns,mmax)

        val(1) = -0.5d0*dsqrt(dble(n)/dble(n+l))*dsqrt(dble(n+2*l+1)/dble(n+l+1))
        val(2) = 1d0
        val(3) = -0.5d0*dsqrt(dble(n+2*l+2)/dble(n+l+2))*dsqrt(dble(n+1)/dble(n+l+1))

        val = -0.5d0*w*angular*val
        call MatSetValues(H,1,row,3,col,val,ADD_VALUES,ierr);CHKERRA(ierr)
        call MatSetValues(H,3,col,1,row,val,ADD_VALUES,ierr);CHKERRA(ierr)
      end if 
    end if 
  end do 

  deallocate(a,b)
  deallocate(col,val)

  call MatAssemblyBegin(H,MAT_FINAL_ASSEMBLY,ierr);CHKERRA(ierr)
  call MatAssemblyEnd(H,MAT_FINAL_ASSEMBLY,ierr);CHKERRA(ierr)
!-----------------------------------------------------------------------------------------
! Create Eigensolver
!-----------------------------------------------------------------------------------------
  call EPSCreate(PETSC_COMM_WORLD,eps,ierr);CHKERRA(ierr)

  call EPSSetOperators(eps,H,S,ierr);CHKERRA(ierr)

  call EPSSetWhichEigenpairs(eps,eps_which,ierr);CHKERRA(ierr)
  call EPSSetType(eps,eps_type,ierr);CHKERRA(ierr)

  call EPSSetTolerances(eps,tol,maxits,ierr);CHKERRA(ierr)
  call EPSSetProblemType(eps,eps_problem,ierr);CHKERRA(ierr)
  call EPSSetDimensions(eps,nev,ncv,mpd,ierr);CHKERRA(ierr)

  !     ** Set solver parameters at runtime
  call EPSSetFromOptions(eps,ierr);CHKERRA(ierr)

  call EPSSolve(eps,ierr);CHKERRA(ierr)


  file_name = 'stark_eigenvectors.bin'
  call PetscViewerBinaryOpen(MPI_COMM_WORLD, file_name, FILE_MODE_WRITE, viewer, ierr)
  CHKERRA(ierr)
  call EPSVectorsView(eps, viewer, ierr);CHKERRA(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRA(ierr)


  ! Closes the hdf5 file
  call h5gclose_f( eps_group_id, h5_err)
  call h5fclose_f( param_file_id, h5_err)
  call CPU_TIME(end_time)
  write(tmp_character, "(ES9.2)")  end_time-start_time
  call PetscPrintf(MPI_COMM_WORLD, 'time   :'//trim(tmp_character)//"\n", ierr)
  CHKERRA(ierr)
  call SlepcFinalize(ierr)
  call forpy_finalize

end program main 
