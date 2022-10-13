module hamil
  use prec
  use fileio
  use couplings
  use constants
  implicit none

  type TDKS
    integer :: ndim
    ! _[c,p,n] means current, previous, next
    complex(kind=q), allocatable, dimension(:) :: psi_c
    complex(kind=q), allocatable, dimension(:) :: psi_p
    complex(kind=q), allocatable, dimension(:) :: psi_n
    complex(kind=q), allocatable, dimension(:,:) :: psi_a
    ! the result of hamiltonian acting on a vector
    complex(kind=q), allocatable, dimension(:) :: hpsi
    ! population
    real(kind=q), allocatable, dimension(:,:) :: pop_a
    real(kind=q), allocatable, dimension(:) :: norm

    complex(kind=q), allocatable, dimension(:,:) :: ham_c
    ! complex(kind=q), allocatable, dimension(:,:) :: ham_p
    ! complex(kind=q), allocatable, dimension(:,:) :: ham_n

    ! KS eigenvalues
    real(kind=q), allocatable, dimension(:,:) :: eigKs
    ! Non-adiabatic couplings
    complex(kind=q), allocatable, dimension(:,:,:) :: NAcoup
    ! Spin-orbital couplings
    complex(kind=q), allocatable, dimension(:,:,:) :: SOcoup

    ! surface hopping related

    ! Bkm = REAL(CONJG(Akm) * Ckm)
    real(kind=q), allocatable, dimension(:) :: Bkm
    real(kind=q), allocatable, dimension(:,:) :: sh_pops
    real(kind=q), allocatable, dimension(:,:) :: sh_prop

    ! whether the memory has been allocated
    logical :: LALLO = .FALSE.

  end type

  contains

  subroutine initTDKS(ks, inp, olap)
    implicit none

    type(TDKS), intent(inout)  :: ks
    type(overlap), intent(in)  :: olap
    type(namdInfo), intent(in) :: inp

    integer :: i, j, N, istat
    integer, allocatable, dimension(:) :: inibs

    ! memory allocation
    ks%ndim = inp%NBASIS
    N = inp%NBASIS

    if (.NOT. ks%LALLO) then
      allocate(ks%psi_c(N))
      allocate(ks%psi_p(N))
      allocate(ks%psi_n(N))
      allocate(ks%hpsi(N))
      allocate(ks%psi_a(N, inp%NAMDTIME))
      allocate(ks%pop_a(N, inp%NAMDTIME))
      allocate(ks%norm(inp%NAMDTIME))

      allocate(ks%ham_c(N,N))

      allocate(ks%eigKs(N, inp%NAMDTIME))
      allocate(ks%NAcoup(N,N, inp%NAMDTIME))
      if (inp%SOCTYPE==2) then
        allocate(ks%SOcoup(N,N, inp%NAMDTIME))
      end if

      allocate(ks%sh_pops(N, inp%NAMDTIME))
      allocate(ks%sh_prop(N,N))
      allocate(ks%Bkm(N))
      ! allocate(ks%ham_p(N,N))
      ! allocate(ks%ham_n(N,N))
      ks%LALLO = .TRUE.
    end if

    ! cero = (0, 0)
    ! uno = (1, 0)
    ks%psi_c = cero
    ks%psi_p = cero
    ks%psi_n = cero
    ! ks%ham_c = cero
    ! ks%ham_p = cero
    ! ks%ham_n = cero

    allocate(inibs(inp%NINIBS))
    do i=1, inp%NINIBS
      if (inp%SOCTYPE==1) then
        istat = inp%INIBAND(i) - inp%BMIN + 1
      else if (inp%SOCTYPE==2) then
        if (inp%INISPIN(i) == 1) then
          istat = inp%INIBAND(i) - inp%BMINU + 1
        else
          istat = inp%INIBAND(i) - inp%BMIND + inp%BMAXU - inp%BMINU + 2
        end if
      end if
      inibs(i) = istat
    end do

    ks%psi_c(inibs) = uno
    ks%psi_c = ks%psi_c / SQRT(REAL(inp%NINIBS))

    do i=1, inp%NAMDTIME
      j = int( mod(inp%NAMDTINI+i-1, inp%NSW-1) )
      if (j==0) j = inp%NSW - 1
      ! We don't need all the information, only a section of it
      ks%eigKs(:,i) = olap%Eig(:, j)
      ! Divide by 2 * POTIM here, because we didn't do this in the calculation
      ! of couplings
      ks%NAcoup(:,:,i) = olap%Dij(:,:, j) / (2*inp%POTIM)
      if (inp%SOCTYPE==2) then
        ks%SOcoup(:,:,i) = olap%Sij(:,:, j)
      end if

    end do

  end subroutine


  ! constructing the hamiltonian
  subroutine make_hamil(TION, TELE, ks, inp)
    implicit none

    type(TDKS), intent(inout) :: ks
    type(namdInfo), intent(in) :: inp
    integer, intent(in) :: TION, TELE

    integer :: i

    ! the hamiltonian contains two parts, which are obtained by interpolation
    ! method between two ionic tims step

    ! The non-adiabatic coupling part
    ks%ham_c(:,:) = ks%NAcoup(:,:,TION) + &
        (ks%NAcoup(:,:,TION+1) - ks%NAcoup(:,:,TION)) * TELE / inp%NELM

    ! multiply by -i * hbar
    ks%ham_c = -imgUnit * hbar * ks%ham_c

    if (inp%SOCTYPE==2) then
      ks%ham_c(:,:) = ks%ham_c(:,:) + ks%SOcoup(:,:,TION) + &
        (ks%SOcoup(:,:,TION+1) - ks%SOcoup(:,:,TION)) * TELE / inp%NELM
    end if

    ! the energy eigenvalue part
    do i=1, ks%ndim
      ks%ham_c(i,i) = ks%ham_c(i,i) + ks%eigKs(i,TION) + &
        (ks%eigKs(i,TION+1) - ks%eigKs(i,TION)) * TELE / inp%NELM
    end do
  end subroutine

  ! Acting the hamiltonian on the state vector
  subroutine hamil_act(ks)
    implicit none
    type(TDKS), intent(inout) :: ks
    integer :: i, j, N
    complex(kind=q) :: tmp

    N = ks%ndim
    do i=1, N
      tmp = cero
      do j=1, N
        tmp = tmp + ks%ham_c(i,j) * ks%psi_c(j)
      end do
      ks%hpsi(i) = tmp
    end do

  end subroutine

end module
