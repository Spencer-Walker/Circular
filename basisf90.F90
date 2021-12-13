
program main
#include <slepc/finclude/slepceps.h>
use slepceps
use hdf5
use mpi
use fgsl
use iso_c_binding
  implicit none
  PetscInt,   parameter :: dp = kind(1.d0)
  PetscInt  :: nmax, INFO
  PetscReal :: Rmax
  character(len = 15) :: label ! file name w/o .h5
  integer(HID_T), allocatable :: psi_space_right(:), psi_dset_right(:)
  integer(HID_T), allocatable :: ener_space(:), ener_dset(:)
  integer(HID_T), allocatable :: psi_space_left(:), psi_dset_left(:)
  integer(HID_T) :: file_id, eps_group_id, param_file_id, eps_dat_id
  PetscInt :: h5_err, num_proc, proc_id
  integer(HSIZE_T)  :: ener_dims(1:2)
  integer(HSIZE_T)  :: psi_dims(1:3) 
  integer(HID_T)    :: h5_kind, plist_id, memtype
  character(len = 24)   :: file_name 
  character(len = 3)    :: strl ! file number
  character(len = 6)    :: fmt  ! format descriptor 
  integer(SIZE_T), parameter :: sdim      = 300 
  character(len = 300)   :: tmp_character
  integer(HSIZE_T)  :: dims(1), dims2(2)
  Mat            :: H
  Vec            :: Vi
  EPS            :: eps
  EPSType        :: tname
  PetscBool      :: cap_present, rot_present
  PetscInt       :: num_points, i, j, Istart, Iend, l_start, l_stop, l_stride
  PetscInt       :: nev, ncv, mpd, ncon, maxits, lmax, l
  PetscInt       :: row(1), col(3), tmp_int
  PetscInt, allocatable :: ix(:), IPIV(:)
  PetscReal :: cap_eta, rot_theta
  PetscErrorCode :: ierr
  EPSProblemType :: eps_problem
  PetscScalar    :: value(3), rint 
  PetscScalar, allocatable :: u_right(:,:), u_left(:,:), E_right(:)
  PetscScalar, allocatable ::V(:),S(:,:),UL(:,:),UR(:,:)
  PetscReal      :: dr, gobbler, tol, start_time, end_time
  PetscReal, allocatable :: EE_right(:,:), uu_right(:,:,:), VV(:,:),r(:)
  PetscReal, allocatable :: uu_left(:,:,:)
  EPSBalance :: eps_balance 
  PetscReal,  parameter :: pi = 3.141592653589793238462643383279502884197169
  PetscBool,  parameter :: eps_world = PETSC_FALSE
  PetscScalar :: eps_target
  PetscReal :: eps_target_components(2)
  PetscBool :: eps_two_sided
  EPSType :: eps_type     
  EPSWhich :: eps_which
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
  
  call h5dopen_f(eps_group_id, "cap_present", eps_dat_id, h5_err)
  call h5dread_f(eps_dat_id, H5T_NATIVE_INTEGER, tmp_int, dims, h5_err)
  call h5dclose_f( eps_dat_id, h5_err)
  if (tmp_int .eq. 1) then
    cap_present = .true.
  else 
    cap_present = .false.
  end if 
  if ( cap_present .eqv. .false.) then
    gobbler = 1.0d0
    cap_eta = 0.0d0
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

  call h5dopen_f(eps_group_id, "cap_eta", eps_dat_id, h5_err)
  call h5dread_f(eps_dat_id, H5T_NATIVE_DOUBLE, cap_eta, dims, h5_err)
  call h5dclose_f( eps_dat_id, h5_err)

  call h5dopen_f(eps_group_id, "gobbler", eps_dat_id, h5_err)
  call h5dread_f(eps_dat_id, H5T_NATIVE_DOUBLE, gobbler, dims, h5_err)
  call h5dclose_f( eps_dat_id, h5_err)

  call h5dopen_f(eps_group_id, "l_max", eps_dat_id, h5_err)
  call h5dread_f(eps_dat_id, H5T_NATIVE_INTEGER, lmax, dims, h5_err)
  call h5dclose_f( eps_dat_id, h5_err)

  call h5dopen_f(eps_group_id, "n_max", eps_dat_id, h5_err)
  call h5dread_f(eps_dat_id, H5T_NATIVE_INTEGER, nmax, dims, h5_err)
  call h5dclose_f( eps_dat_id, h5_err)

  call h5dopen_f(eps_group_id, "R_max", eps_dat_id, h5_err)
  call h5dread_f(eps_dat_id, H5T_NATIVE_DOUBLE, Rmax, dims, h5_err)
  call h5dclose_f( eps_dat_id, h5_err)

  call h5dopen_f(eps_group_id, "delta_x", eps_dat_id, h5_err)
  call h5dread_f(eps_dat_id, H5T_NATIVE_DOUBLE, dr, dims, h5_err)
  call h5dclose_f( eps_dat_id, h5_err)

  call h5dopen_f(eps_group_id, "max_its", eps_dat_id, h5_err)
  call h5dread_f(eps_dat_id, H5T_NATIVE_INTEGER, maxits, dims, h5_err)
  call h5dclose_f( eps_dat_id, h5_err)

  call h5dopen_f(eps_group_id, "abs_tol", eps_dat_id, h5_err)
  call h5dread_f(eps_dat_id, H5T_NATIVE_DOUBLE, tol, dims, h5_err)
  call h5dclose_f( eps_dat_id, h5_err)

  call h5dopen_f(eps_group_id, "label", eps_dat_id, h5_err)
  call h5dread_f(eps_dat_id, memtype, label, dims, h5_err)
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

  ! Add the .h5 extension to the file label
  file_name = trim(label)//'.h5'
  allocate(psi_space_right(0:lmax))
  allocate(psi_dset_right(0:lmax))
  allocate(ener_space(0:lmax))
  allocate(ener_dset(0:lmax))
  allocate(psi_space_left(0:lmax))
  allocate(psi_dset_left(0:lmax))

  num_points = nint(Rmax/dr)
  Rmax = num_points*dr
 
  allocate(r(0:num_points-1))
  dims(1) = num_points-1
  call h5dopen_f(eps_group_id, "r", eps_dat_id, h5_err)
  call h5dread_f(eps_dat_id, H5T_NATIVE_DOUBLE, r, dims, h5_err)
  call h5dclose_f( eps_dat_id, h5_err)

  call h5pcreate_f( H5P_FILE_ACCESS_F, plist_id, h5_err)
  call h5pset_fapl_mpio_f( plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, h5_err)
  call h5fcreate_f( file_name, H5F_ACC_TRUNC_F, file_id, h5_err, &
  & access_prp = plist_id)
  h5_kind = h5kind_to_type( dp, H5_REAL_KIND)

  ! Iterate over all l and create a file for each energy and wfn
  do l = 0, lmax
    ! Figures out how many characters are needed to represent l
    ! (if l = 0 1 character if l = 10 2 characters, ... etc. )
    if ( l .le. 9 ) then
      fmt = '(I1.1)'
    else if ( l .le. 99 ) then
      fmt = '(I2.2)'
    else
      fmt = '(I3.3)'
    end if
    
    ! Create string strl for number l
    write(strl,fmt) l

    ! For each l there will be an array of nmax-l energies
    ener_dims(1) = nmax-l
    ener_dims(2) = 2
    ! There will be nmax-l wfns, each with dim int(Rmax/dr)
    psi_dims(1) = int(Rmax/dr)
    psi_dims(2) = nmax-l
    psi_dims(3) = 2
    ! Creates an hdf5 interface that knows information about the data
    call h5screate_simple_f( 2, ener_dims, ener_space(l), h5_err)
    
    ! Connects the datasapce to the hdf5 file and gives it a name 
    call h5dcreate_f( file_id, "Energy_l"//trim(strl), h5_kind, ener_space(l), &
    & ener_dset(l), h5_err)

    ! Does the same as before for the energy but with the wfns 
    call h5screate_simple_f( 3, psi_dims, psi_space_right(l), h5_err)
    call h5dcreate_f( file_id, "Psi_r_l"//trim(strl), h5_kind, psi_space_right(l), &
    & psi_dset_right(l), h5_err)

    ! Does the same as before for the energy but with the wfns 
    call h5screate_simple_f( 3, psi_dims, psi_space_left(l), h5_err)
    call h5dcreate_f( file_id, "Psi_l_l"//trim(strl), h5_kind, psi_space_left(l), &
    & psi_dset_left(l), h5_err)

  end do 

  if (eps_world .eqv. PETSC_FALSE) then
    l_start = proc_id
    l_stop = lmax
    l_stride = num_proc
  else
    l_start = 0
    l_stop = lmax
    l_stride = 1 
  end if 

  do l = l_start,l_stop,l_stride
    write(tmp_character, "(A2,I4)")  'lB', l
    call PetscPrintf(MPI_COMM_SELF, trim(tmp_character)//"\n", ierr)
    CHKERRA(ierr)
    nev = nmax - l
    allocate(V(0:num_points-1))
    allocate(VV(0:num_points-1,2))
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !     Compute the operator matrix that defines the eigensystem, Ax=kx
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (eps_world .eqv. PETSC_TRUE) then
      call MatCreate(PETSC_COMM_WORLD,H,ierr);CHKERRA(ierr)
    else
      call MatCreate(PETSC_COMM_SELF,H,ierr);CHKERRA(ierr)
    end if
    call MatSetSizes(H,PETSC_DECIDE,PETSC_DECIDE,num_points,num_points,ierr);CHKERRA(ierr)
    call MatSetFromOptions(H,ierr);CHKERRA(ierr)
    call MatSetUp(H,ierr);CHKERRA(ierr)
    call MatGetOwnershipRange(H,Istart,Iend,ierr);CHKERRA(ierr)

    dims2(1) = num_points
    dims2(2) = 2
    call h5dopen_f(eps_group_id, "V", eps_dat_id, h5_err)
    call h5dread_f(eps_dat_id, H5T_NATIVE_DOUBLE, VV, dims2, h5_err)
    call h5dclose_f( eps_dat_id, h5_err)
    V = dcmplx(VV(:,1),VV(:,2))
    V = V + zexp((0d0,-2d0)*rot_theta)* dble(l*(l+1))/(2d0*r**2d0)

    if (Istart .eq. 0) then
      row(1) = 0
      col(1) = 0
      col(2) = 1
      value(1) =  zexp((0d0,-2d0)*rot_theta)*1.0d0/(dr**2.0d0) + V(Istart)
      value(2) = -zexp((0d0,-2d0)*rot_theta)*0.5d0/(dr**2.0d0)
      call MatSetValues(H,1,row,2,col,value,INSERT_VALUES,ierr);CHKERRA(ierr)
      Istart = Istart+1
    end if

    if (Iend .eq. num_points) then
      row(1) = num_points-1
      col(1) = num_points-2
      col(2) = num_points-1
      value(1) = -zexp((0d0,-2d0)*rot_theta)*0.5d0/(dr**2.0d0)
      value(2) =  zexp((0d0,-2d0)*rot_theta)*1.0d0/(dr**2.0d0) + V(num_points-1) 
      call MatSetValues(H,1,row,2,col,value,INSERT_VALUES,ierr);CHKERRA(ierr)
      Iend = Iend-1
    end if

    value(1) = -zexp((0d0,-2d0)*rot_theta)*0.5d0/(dr**2.0d0)
    value(2) =  zexp((0d0,-2d0)*rot_theta)*1.0d0/(dr**2.0d0) 
    value(3) = -zexp((0d0,-2d0)*rot_theta)*0.5d0/(dr**2.0d0)
    do i=Istart,Iend-1
      row(1) = i
      col(1) = i-1
      col(2) = i
      col(3) = i+1
      value(2) =  zexp((0d0,-2d0)*rot_theta)*1.0d0/(dr**2.0d0) + V(i)
      call MatSetValues(H,1,row,3,col,value,INSERT_VALUES,ierr);CHKERRA(ierr)
    end do

    call MatAssemblyBegin(H,MAT_FINAL_ASSEMBLY,ierr);CHKERRA(ierr)
    call MatAssemblyEnd(H,MAT_FINAL_ASSEMBLY,ierr);CHKERRA(ierr)
    call MatCreateVecs(H,Vi,PETSC_NULL_VEC,ierr);CHKERRA(ierr)

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !     Create the eigensolver and display info
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    !     ** Create eigensolver context
    if (eps_world .eqv. PETSC_TRUE) then
      call EPSCreate(PETSC_COMM_WORLD,eps,ierr);CHKERRA(ierr)
    else
      call EPSCreate(PETSC_COMM_SELF,eps,ierr);CHKERRA(ierr)
    end if 

    call EPSSetTwoSided(eps,eps_two_sided,ierr);CHKERRA(ierr)

    !     ** Set operators. In this case, it is a standard eigenvalue problem
    call EPSSetOperators(eps,H,PETSC_NULL_MAT,ierr);CHKERRA(ierr)
    call EPSSetBalance(eps,eps_balance,PETSC_DEFAULT_INTEGER,PETSC_DEFAULT_REAL,ierr);CHKERRA(ierr)

    call EPSSetWhichEigenpairs(eps,eps_which,ierr);CHKERRA(ierr)
    call EPSSetType(eps,eps_type,ierr);CHKERRA(ierr)

    if (eps_which .eq. EPS_TARGET_MAGNITUDE ) then
      call EPSSetTarget(eps,eps_target,ierr);CHKERRA(ierr)
    end if 

    call EPSSetTolerances(eps,tol,maxits,ierr);CHKERRA(ierr)
    call EPSSetProblemType(eps,eps_problem,ierr);CHKERRA(ierr)
    call EPSSetDimensions(eps,nev,ncv,mpd,ierr);CHKERRA(ierr)

    !     ** Set solver parameters at runtime
    call EPSSetFromOptions(eps,ierr);CHKERRA(ierr)

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !     Solve the eigensystem
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    call EPSSolve(eps,ierr);CHKERRA(ierr)

    !     ** Optional: Get some information from the solver and display it
    call EPSGetType(eps,tname,ierr);CHKERRA(ierr)

    call EPSGetDimensions(eps,nev,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,ierr);CHKERRA(ierr)

    call EPSGetConverged(eps,ncon,ierr);CHKERRA(ierr)

    if (ncon < nev) then
      call PetscPrintf(MPI_COMM_SELF, "ncon < nev\n", ierr)
      CHKERRA(ierr)
      write(tmp_character, "(I4,I4)")  ncon, nev
      call PetscPrintf(MPI_COMM_SELF,"ncon, nev = "//tmp_character//"\n", ierr)
      CHKERRA(ierr)
      stop
    end if

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !     Store solution
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  

    allocate(u_right(num_points,nev))
    allocate(E_right(nev))
    allocate(ix(num_points))
    allocate(uu_right(num_points,nev,2))
    allocate(EE_right(nev,2))
    allocate(u_left(num_points,nev))
    allocate(uu_left(num_points,nev,2))
    allocate(S(nev,nev))
    allocate(UL(nev,nev))
    allocate(UR(nev,nev))
    allocate(IPIV(nev))

    do i = 1, num_points
      ix(i) = i
    end do 

    do i = 0, nev-1
      call EPSGetEigenpair(eps,i,E_right(i+1),PETSC_NULL_SCALAR,Vi,PETSC_NULL_VEC,ierr);CHKERRA(ierr)
      call VecGetValues(Vi,num_points,ix,u_right(1:num_points,i+1),ierr);CHKERRA(ierr) 

      u_right(num_points,i+1) = (0d0,0d0)
      rint = sum(u_right(:,i+1)*u_right(:,i+1))
      u_right(:,i+1) = u_right(:,i+1)/zsqrt(rint)

      if (eps_two_sided .eqv. PETSC_TRUE) then
        call EPSGetLeftEigenvector(eps,i,Vi,PETSC_NULL_VEC,ierr);CHKERRA(ierr)
        call VecGetValues(Vi,num_points,ix,u_left(1:num_points,i+1),ierr);CHKERRA(ierr) 
        u_left(num_points,i+1) = (0d0,0d0)
        rint  = sum(u_left(:,i+1)*u_left(:,i+1))
        u_left(:,i+1) = dconjg(u_left(:,i+1)/zsqrt(rint))

        rint  = sum(u_left(:,i+1)*u_right(:,i+1))
        u_left(:,i+1) = u_left(:,i+1)/zsqrt(rint)
        u_right(:,i+1) = u_right(:,i+1)/zsqrt(rint)
      else 
        do j = 0,i-1
          u_right(:,i+1) = u_right(:,i+1) - sum(u_right(:,i+1)*u_right(:,j+1))*u_right(:,j+1)
        end do
        rint = sum(u_right(:,i+1)*u_right(:,i+1))
        u_right(:,i+1) = u_right(:,i+1)/zsqrt(rint)
      end if
    end do 

#if 0
    if(eps_two_sided .eqv. PETSC_TRUE ) then
      S = matmul(transpose(u_left),u_right)

      call zgetrf(nev,nev,S,nev,IPIV,INFO)

      if (INFO .ne. 0) then
        call PetscPrintf(MPI_COMM_SELF, "'LU' factorization failed\n", ierr)
        CHKERRA(ierr)
        write(tmp_character, "(A10,I4)")  "exit code ", INFO
        call PetscPrintf(MPI_COMM_SELF, trim(tmp_character)//"\n", ierr)
        CHKERRA(ierr)
      end if 

      do i = 1, nev
        if (IPIV(i) .ne. i) then
          call PetscPrintf(MPI_COMM_SELF, "Pivits are hapening\n", ierr)
          CHKERRA(ierr)
          stop
        end if 
      end do 

      UR = 0.d0
      UL = 0.d0
      do i = 1,nev
        do j = i,nev
          UR(i,j) = S(i,j)
        end do
      end do 
      
      do j = 1, nev
        do i = j+1,nev
          UL(i,j) = S(i,j)
        end do
      end do 

      do i = 1,nev
        UL(i,i) = 1.d0
      end do 

      call ztrtri('U','N',nev,UR,nev,INFO)
      call ztrtri('L','N',nev,UL,nev,INFO)
      
      if (INFO .ne. 0) then
        call PetscPrintf(MPI_COMM_SELF, 'Inversion failed\n', ierr)
        CHKERRA(ierr)
        stop
      end if 


      u_left = matmul(u_left,transpose(UL))

      u_right = matmul(u_right,UR)
    else
      u_left = u_right
    end if 
#endif 
    u_left = u_right
    ! There will be nmax-l eigenstates for each l 
    ener_dims(1) = nmax-l
    ener_dims(2) = 2
    ! Each of the nmax-l eigenstates will have a wfn of length int(Rmax/h)
    psi_dims(1) = int(Rmax/dr)
    psi_dims(2) = nmax-l
    psi_dims(3) = 2
    
    uu_right(:,:,1) = real(real(u_right))
    uu_right(:,:,2) = real(aimag(u_right))

    EE_right(:,1) = real(real(E_right))
    EE_right(:,2) = real(aimag(E_right))

    ! Writes E_right to Energy_l#l (l#l can be l0, l1, ... etc ) 
    call h5dwrite_f( ener_dset(l), h5_kind, EE_right, ener_dims, h5_err)
    
    ! Writes u_right to Psi_l#l
    call h5dwrite_f( psi_dset_right(l), h5_kind, uu_right, psi_dims, h5_err)

    uu_left(:,:,1) = real(real(u_left))
    uu_left(:,:,2) = real(aimag(u_left))

    ! Writes u_left to Psi_l_l#l
    call h5dwrite_f( psi_dset_left(l), h5_kind, uu_left, psi_dims, h5_err)

    
    deallocate(u_right)
    deallocate(uu_right)
    deallocate(EE_right)
    deallocate(ix)
    deallocate(u_left)
    deallocate(E_right)
    deallocate(uu_left)
    deallocate(S)
    deallocate(UL)
    deallocate(UR)
    deallocate(IPIV)
    
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !     Display solution and clean up
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    call EPSDestroy(eps,ierr);CHKERRA(ierr)
    call MatDestroy(H,ierr);CHKERRA(ierr)
    call MatDestroy(H,ierr);CHKERRA(ierr)

    deallocate(V)
    deallocate(VV)
  end do 

    ! For each l close the hdf5 objects that were created in initialize
  do l = 0, lmax
    ! Closes the energy dataset
    call h5dclose_f( ener_dset(l), h5_err)
    ! Closes the energy dataspace
    call h5sclose_f( ener_space(l), h5_err)
    ! Closes the psi dataset
    call h5dclose_f( psi_dset_right(l), h5_err)
    ! Closes the psi dataspace
    call h5sclose_f( psi_space_right(l), h5_err)
    ! Closes the psi dataset
    call h5dclose_f( psi_dset_left(l), h5_err)
    ! Closes the psi dataspace
    call h5sclose_f( psi_space_left(l), h5_err)
  end do 
  
  deallocate(r)
  deallocate(psi_space_right)
  deallocate(psi_dset_right)
  deallocate(ener_space)
  deallocate(ener_dset)
  deallocate(psi_space_left)
  deallocate(psi_dset_left)

  ! Closes the hdf5 file
  call h5fclose_f( file_id, h5_err)
  call h5gclose_f( eps_group_id, h5_err)
  call h5fclose_f( param_file_id, h5_err)
  call CPU_TIME(end_time)
  write(tmp_character, "(ES9.2)")  end_time-start_time
  call PetscPrintf(MPI_COMM_WORLD, 'time   :'//trim(tmp_character)//"\n", ierr)
  CHKERRA(ierr)
  call SlepcFinalize(ierr)


end program main 
