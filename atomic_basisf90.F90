program main
#include <slepc/finclude/slepceps.h>
use slepceps
use hdf5
use mpi
use fgsl
use forpy_mod
use iso_c_binding
use forpy_mod
  implicit none
  integer(HID_T) :: memtype, eps_group_id, param_file_id, eps_dat_id
  integer(SIZE_T), parameter :: sdim = 300 
  integer(HSIZE_T)  :: dims(1)
  PetscInt, parameter :: dp = kind(1.d0)
  PetscInt :: Ns, nev, ncv, mpd, maxits, tmp_int, h5_err, num_proc, proc_id
  PetscInt :: n, l, Nshell, Ntot, mmax, i, j, row(1), lmax, Nl, ierror
  PetscInt, allocatable :: col(:)
  character(len = 6)    :: fmt,strl  ! format descriptor
  character(len = 300)   :: tmp_character, file_name
  PetscBool      :: rot_present, eps_two_sided 
  PetscReal,  parameter :: pi = 3.141592653589793238462643383279502884197169
  PetscReal :: rot_theta, tol, start_time, end_time, eps_target_components(2), kabs
  PetscReal :: C0, C, Zc 
  PetscReal, allocatable ::  a(:), b(:)
  PetscErrorCode :: ierr
  EPSProblemType :: eps_problem 
  PetscScalar :: eps_target, k, eigr, eigi, q, q2, q3, g
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
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !     Beginning of program
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  call CPU_TIME(start_time)

  ierror = forpy_initialize()
  ierror = import_py(opt, "mpmath")

  call SlepcInitialize(PETSC_NULL_CHARACTER,ierr)
  if (ierr .ne. 0) then
    call PetscPrintf(MPI_COMM_WORLD, "SlepcInitialize failed\n", ierr)
    CHKERRA(ierr)
    stop
  end if
  call MPI_Comm_rank(PETSC_COMM_SELF,proc_id,ierr);CHKERRA(ierr)
  call MPI_Comm_size(PETSC_COMM_SELF,num_proc,ierr);CHKERRA(ierr)

  ! Initialize hdf5 
  call h5open_f( h5_err)
  if ( h5_err .ne. 0 ) then
    call PetscPrintf(MPI_COMM_SELF, 'h5open_f failed\n', ierr);CHKERRA(ierr)
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

  print*, 'a =', a
  print*, 'b =', b

  dims(1) = 1
  call h5dopen_f(eps_group_id, "Nl", eps_dat_id, h5_err)
  call h5dread_f(eps_dat_id, H5T_NATIVE_INTEGER, Nl, dims, h5_err)
  call h5dclose_f( eps_dat_id, h5_err)

  call h5dopen_f(eps_group_id, "m_max", eps_dat_id, h5_err)
  call h5dread_f(eps_dat_id, H5T_NATIVE_INTEGER, mmax, dims, h5_err)
  call h5dclose_f( eps_dat_id, h5_err)

  lmax = Nl - 1 + mmax 
   
  call h5dopen_f(eps_group_id, "Ns", eps_dat_id, h5_err)
  call h5dread_f(eps_dat_id, H5T_NATIVE_INTEGER, Ns, dims, h5_err)
  call h5dclose_f( eps_dat_id, h5_err)
  
  call h5dopen_f(eps_group_id, "Na", eps_dat_id, h5_err)
  call h5dread_f(eps_dat_id, H5T_NATIVE_INTEGER, nev, dims, h5_err)
  call h5dclose_f( eps_dat_id, h5_err)

  mmax = 0 ! Spherically symmetric problem

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

  tmp_character = "EPS_SMALLEST_REAL"

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


          Ntot = Ns


        do l = 0, lmax
          print*, 'l =', l
          if ( l .le. 9 ) then
            fmt = '(I1.1)'
          else if ( l .le. 99 ) then
            fmt = '(I2.2)'
          else
            fmt = '(I3.3)'
          end if

          write(strl,fmt) l


          call MatCreate(PETSC_COMM_SELF,S,ierr);CHKERRA(ierr)
          call MatSetSizes(S,PETSC_DECIDE,PETSC_DECIDE,Ntot,Ntot,ierr);CHKERRA(ierr)
          call MatSetFromOptions(S,ierr);CHKERRA(ierr)
          call MatSetUp(S,ierr);CHKERRA(ierr)
          
          call MatCreate(PETSC_COMM_SELF,H,ierr);CHKERRA(ierr)
          call MatSetSizes(H,PETSC_DECIDE,PETSC_DECIDE,Ntot,Ntot,ierr);CHKERRA(ierr)
          call MatSetFromOptions(H,ierr);CHKERRA(ierr)
          call MatSetUp(H,ierr);CHKERRA(ierr)
          

        !-----------------------------------------------------------------------------------------
        ! Build Overlap Matrix
        !-----------------------------------------------------------------------------------------
          allocate(col(Ns),val(Ns))

           
          do i = 0,Ns-1
            n = i  
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
          do i = 0,Ns-1  
            n = i 
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

          if (dabs(Zc) .gt. 1d-15) then
            g = 2*(k/C) 
            do i = 0, Ns-1
              row(1)  = i
              do n = 0, Ns-1
                ierror = tuple_create(args, 4)
                ierror = args%setitem(0, -i)
                ierror = args%setitem(1, -n)
                ierror = args%setitem(2, -n-i-2*l-1)
                ierror = args%setitem(3, 1-g**2)
                ierror = call_py(retval, opt, "hyp2f1", args)
                ierror = cast_nonstrict(q, retval)
                col(n+1) = n
                val(n+1) = -k*dsqrt(dgamma(dble(n+1))*dgamma(dble(i+1))/(dble((n+l+1)*(i+l+1)) &
                &*dgamma(dble(n+2*l+2))*dgamma(dble(i+2*l+2))))*Zc*dgamma(dble(n+i+2*l+2)) &
                &/(dgamma(dble(n+1))*dgamma(dble(i+1)))*(g**dble(2*l+2)/(g+1)**dble(n+i+2*l+2))*q
              end do
              call MatSetValues(H,1,row,Ns,col,val,ADD_VALUES,ierr);CHKERRA(ierr)  
            end do 
          end if 

          if (Nshell .ne. 0) then
            do j = 1, Nshell 
              g = 2d0*(k/b(j))
              do i = 0, Ns-1
                row(1)  = i
                do n = 0, Ns-1
                  ierror = tuple_create(args, 4)
                  ierror = args%setitem(0, -i)
                  ierror = args%setitem(1, -n)
                  ierror = args%setitem(2, -n-i-2*l-1)
                  ierror = args%setitem(3, 1-g**2)
                  ierror = call_py(retval, opt, "hyp2f1", args)
                  ierror = cast_nonstrict(q, retval)

                  ierror = tuple_create(args, 4)
                  ierror = args%setitem(0, -i)
                  ierror = args%setitem(1, -n-1)
                  ierror = args%setitem(2, -n-i-2*l-2)
                  ierror = args%setitem(3, 1-g**2)
                  ierror = call_py(retval, opt, "hyp2f1", args)
                  ierror = cast_nonstrict(q2, retval)

          ierror = tuple_create(args, 4)
          ierror = args%setitem(0, -i)
          ierror = args%setitem(1, -n+1)
          ierror = args%setitem(2, -n-i-2*l)
          ierror = args%setitem(3, 1-g**2)
          ierror = call_py(retval, opt, "hyp2f1", args)
          ierror = cast_nonstrict(q3, retval)
        
          col(n+1) = n
          val(n+1) = -dsqrt(dgamma(dble(n+1))*dgamma(dble((i+1))) &
                & /(dble(n+l+1)*dble(i+l+1)*dgamma(dble(n+2*l+2))*dgamma(dble(i+2*l+2))))&
                & *(0.5d0*a(j))*( 2d0*dble(n+l+1)*(dgamma(dble(n+i+2*l+2)) &
                & /(dgamma(dble(n+1))*dgamma(dble(i+1))))*(g**dble(2*l+2)/(g+1d0)**dble(n+i+2*l+2))*q& 
                & -dble(n+1)*(dgamma(dble(n+i+2*l+3))/(dgamma(dble(n+2))*dgamma(dble(i+1))))&
                & *(g**dble(2*l+2)/(g+1d0)**dble(n+i+2*l+3))*q2&
                & -dble(n+2*l+1)*(dgamma(dble(n+i+2*l+1))/(dgamma(dble(n))*dgamma(dble(i+1))))&
                & *(g**dble(2*l+2)/(g+1)**dble(n+i+2*l+1)) *q3)
        end do
        call MatSetValues(H,1,row,Ns,col,val,ADD_VALUES,ierr);CHKERRA(ierr)
      end do
    end do 
  end if
  deallocate(col,val)
  call MatAssemblyBegin(H,MAT_FINAL_ASSEMBLY,ierr);CHKERRA(ierr)
  call MatAssemblyEnd(H,MAT_FINAL_ASSEMBLY,ierr);CHKERRA(ierr)

!-----------------------------------------------------------------------------------------
! Create Eigensolver
!-----------------------------------------------------------------------------------------
  call EPSCreate(PETSC_COMM_SELF,eps,ierr);CHKERRA(ierr)
  call EPSSetOperators(eps,H,S,ierr);CHKERRA(ierr)

  call EPSSetWhichEigenpairs(eps,eps_which,ierr);CHKERRA(ierr)
  call EPSSetType(eps,eps_type,ierr);CHKERRA(ierr)


  call EPSSetTolerances(eps,tol,maxits,ierr);CHKERRA(ierr)
  call EPSSetProblemType(eps,eps_problem,ierr);CHKERRA(ierr)
  call EPSSetDimensions(eps,nev,ncv,mpd,ierr);CHKERRA(ierr)


  !     ** Set solver parameters at runtime
  call EPSSetFromOptions(eps,ierr);CHKERRA(ierr)

  call EPSSolve(eps,ierr);CHKERRA(ierr)

  file_name = 'atomic_eigenvectors_l'//trim(strl)//'.bin'
  call PetscViewerBinaryOpen(MPI_COMM_WORLD, file_name, FILE_MODE_WRITE, viewer, ierr)
  CHKERRA(ierr)
  call EPSVectorsView(eps, viewer, ierr);CHKERRA(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRA(ierr)


  if (proc_id .eq. 0 ) then 
    open( unit = 10, file = 'atomic_eigenvalues_l'//trim(strl)//'-real.csv')
    open( unit = 20, file = 'atomic_eigenvalues_l'//trim(strl)//'-imag.csv') 
    do i = 0, nev-1
      call EPSGetEigenvalue(eps, i, eigr, eigi, ierr);CHKERRA(ierr)
      write(10,*) realpart(eigr)
      write(20,*) imagpart(eigr)
    end do 
  end if 

  call EPSDestroy(eps,ierr);CHKERRA(ierr) 
  call MatDestroy(S,ierr);CHKERRA(ierr)
  call MatDestroy(H,ierr);CHKERRA(ierr)

end do 

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
