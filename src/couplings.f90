module couplings
  use prec
  use constants
  use lattice
  use wavecar
  use fileio

  implicit none

  type overlap
    integer :: ISPIN
    integer :: NBANDS
    integer :: TSTEPS
    real(kind=q) :: dt
    complex(kind=q), allocatable, dimension(:,:,:) :: Dij
    complex(kind=q), allocatable, dimension(:,:,:) :: Sij
    real(kind=q), allocatable, dimension(:,:) :: Eig
  end type

  contains

  subroutine CoupFromFile(olap)
    implicit none
    type(overlap), intent(inout) :: olap
    integer :: i, j, k, ierr, irec
    integer :: irecordL
    real(kind=q) :: recordL, rnbands, rnsw, rdt
    ! to find out the record length
    real(kind=q), allocatable, dimension(:)  :: values

    open(unit=20, file='COUPCAR', access='direct', form='unformatted', &
         status = 'unknown', recl=256, iostat=ierr)
    if(ierr /= 0) then
      write(*,*) "File I/O error with COUPCAR"
      stop
    end if

    read(unit=20,rec=1) recordL, rnbands, rnsw, rdt
    ! write(*,*) recordL, rnbands, rnsw, rdt

    if (olap%NBANDS /= NINT(rnbands) .or. &
        olap%TSTEPS /= NINT(rnsw)) then
      ! write(*,*) olap%NBANDS, NINT(rnbands), olap%TSTEPS, NINT(rnsw)
      write(*,*) "The COUPCAR seems to be wrong..."
      stop
    end if

    close(20)
    irecordL = NINT(recordL)
    open(unit=20, file='COUPCAR', access='direct', form='unformatted', &
         status = 'unknown', recl=irecordL, iostat=ierr)

    write(*,*)
    write(*,'(A)') "------------------------------------------------------------"
    write(*,*) "Reading couplings from COUPCAR..."
    irec = 1
    do k=1, olap%TSTEPS-1
      irec = irec + 1
      read(unit=20, rec=irec) ((olap%Dij(i,j,k), i=1, olap%NBANDS), j=1, olap%NBANDS), &
                               (olap%Eig(i,k), i=1, olap%NBANDS)
    end do
    write(*,*) "Done..."
    write(*,'(A)') "------------------------------------------------------------"
    write(*,*)
    close(20)
  end subroutine

  subroutine CoupToFile(olap)
    implicit none
    type(overlap), intent(in) :: olap

    ! Couplings are save to a binary file
    integer :: recordLA, recordLB, recordL
    integer :: i, j, k, ierr, irec
    ! to find out the record length
    real(kind=q), allocatable, dimension(:)  :: valuesA
    complex(kind=q), allocatable, dimension(:)  :: valuesB

    allocate(valuesA(olap%NBANDS))
    allocate(valuesB(olap%NBANDS * olap%NBANDS))
    inquire (iolength=recordLA) valuesA
    inquire (iolength=recordLB) valuesB
    recordL = recordLA + recordLB
    deallocate(valuesA)
    deallocate(valuesB)

    open(unit=20, file='COUPCAR', access='direct', form='unformatted', &
         status='unknown', recl=recordL, iostat=ierr)
    if(ierr /= 0) then
        write(*,*) "File I/O error with COUPCAR"
        stop
    end if

    irec = 1
    write(unit=20, rec=irec) real(recordL, kind=q),     &
                             real(olap%NBANDS, kind=q), &
                             real(olap%TSTEPS, kind=q), &
                             real(olap%dt, kind=q)
    do k=1, olap%TSTEPS-1
      irec = irec + 1
      write(unit=20, rec=irec) ((olap%Dij(i,j,k), i=1, olap%NBANDS), j=1, olap%NBANDS), &
                               (olap%Eig(i,k), i=1, olap%NBANDS)
    end do
    close(unit=20)

  end subroutine

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Calculate Nonadiabatic Couplings From Two WAVECARs
  ! <psi_i(t)| d/dt |(psi_j(t))> ~=~
  !                                 (<psi_i(t)|psi_j(t+dt)> -
  !                                  <psi_j(t)|psi_i(t+dt)>) / (2dt)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine CoupIJ(waveA, waveB, Cij, is)
    implicit none

    type(waveinfo), intent(in) :: waveA, waveB
    integer, intent(in) :: is

    integer :: s, i, j, dbs
    type(psi) :: ket
    ! <psi_i(t)| d/dt |(psi_j(t))>
    ! The coupling as defined above is a real number
    complex(kind=q), dimension(:,:), intent(inout)  :: Cij
    ! stores plane wave coefficients of the wavefunctions
    complex(kind=qs), allocatable, dimension(:,:) :: crA
    complex(kind=qs), allocatable, dimension(:,:) :: crB
    complex(kind=q) :: pij, pji

    allocate(crA(waveA%NPLWS(1), waveA%NBANDS))
    allocate(crB(waveB%NPLWS(1), waveB%NBANDS))

    ! read in all the wavefunctions
    ! the Gamma point WAVECAR has only ONE kpoint
    do i=1, waveA%NBANDS
        ! i-th band, firtst kpoint, first spin
        call setKet(ket, i, 1, is)
        ! the coefficients are normalized in the LOADWAVE subroutine
        ! here, we don't have to worry about normalization problem
        call LOADWAVE(crA(:,i), ket, waveA)
        call LOADWAVE(crB(:,i), ket, waveB)
    end do

    ! Initialization
    Cij = 0

    ! <psi_i(t)| d/dt |(psi_j(t))> ~=~
    !                             (<psi_i(t)|psi_j(t+dt)> - <psi_j(t)|psi_i(t+dt)>) / (2dt)
    if (is==1) write(*,*) "#", trim(waveA%WAVECAR)
    do i=1, waveA%NBANDS
      ! write(*,*) "C_ZERO: ", REAL(crA(1,i)), AIMAG(crA(1,i))
      do j=i, waveA%NBANDS
        ! <psi_i(t)|psi_j(t+dt)>
        pij = SUM(CONJG(crA(:,i)) * crB(:,j))
        ! <psi_j(t)|psi_i(t+dt)>
        pji = SUM(CONJG(crB(:,i)) * crA(:,j))
        ! Not devided by 2 * dt
        Cij(i, j) = pij - pji
        if ( i /= j ) then
          Cij(j, i) = -CONJG(Cij(i, j))
        end if
      end do
    end do

    ! Cij = -0.658218 * Cij / 2.
    ! do i = 1, waveA%NBANDS
    !   write(*,*) (Cij(i,j), j=1, waveA%NBANDS)
    ! end do

    deallocate(crA, crB)

  end subroutine CoupIJ

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Open WAVECARs and do some checking
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine initAB(FileA, FileB, waveA, waveB)
    implicit none

    character(len=*), intent(in) :: FileA, FileB
    type(waveinfo), intent(inout)  :: waveA, waveB

    integer, parameter :: IUA = 12
    integer, parameter :: IUB = 13

    waveA%IU = IUA
    waveA%WAVECAR = FileA
    waveB%IU = IUB
    waveB%WAVECAR = FileB

    call sysinfo(waveA)
    call sysinfo(waveB)

    if (waveA%ISPIN /= waveB%ISPIN) stop
    if (waveA%NKPTS /= waveB%NKPTS) stop
    if (waveA%NBANDS /= waveB%NBANDS) stop
    ! only one K-point
    if (waveA%NPLWS(1) /= waveB%NPLWS(1)) stop

  end subroutine


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! close WAVECARs and free memories
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine finishAB(waveA, waveB)
    implicit none

    type(waveinfo), intent(inout)  :: waveA, waveB

    call freemem(waveA)
    call freemem(waveB)

    call closewav(waveA%IU)
    call closewav(waveB%IU)

  end subroutine

  subroutine TDCoupIJ(rundir, olap, olap_sec, inp)
    implicit none
    character(len=*), intent(in) :: rundir
    type(overlap), intent(inout) :: olap
    type(overlap), intent(inout) :: olap_sec
    type(namdInfo), intent(in) :: inp

    logical :: lcoup
    integer :: i, j, nsw, ndigit, nb, is, nspin
    character(len=256) :: fileA, fileB, buf, tmp
    type(waveinfo) :: waveA, waveB

    nb = inp%NBANDS
    nspin = inp%SOCTYPE

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Initialization
    olap%ISPIN = inp%SOCTYPE
    olap%NBANDS = inp%NBANDS
    olap%TSTEPS = inp%NSW
    olap%dt = inp%POTIM
    allocate(olap%Dij(olap%NBANDS, olap%NBANDS*olap%ISPIN, olap%TSTEPS-1))
    allocate(olap%Eig(olap%NBANDS*olap%ISPIN, olap%TSTEPS-1))

    olap_sec%NBANDS = inp%NBASIS
    olap_sec%TSTEPS = inp%NSW
    olap_sec%dt = inp%POTIM
    allocate(olap_sec%Dij(olap_sec%NBANDS, olap_sec%NBANDS, olap_sec%TSTEPS-1))
    allocate(olap_sec%Eig(olap_sec%NBANDS, olap_sec%TSTEPS-1))
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    nsw = olap%TSTEPS
    write(buf,*) nsw
    ndigit = len_trim(adjustl(buf))
    write(buf,*) '(I0.', ndigit, ')'

    inquire(file='COUPCAR', exist=lcoup)
    if (lcoup) then
      ! file containing couplings exists, then read it
      if (inp%LCPTXT) then
        call readNaEig(olap_sec, inp)
      else
        call CoupFromFile(olap)
        call copyToSec(olap, olap_sec, inp)
        call writeNaEig(olap_sec, inp)
      end if
    else
      ! create the couplings from the wavefunctions
      write(*,*)
      write(*,'(A)') "------------------------------------------------------------"
      do i=1, nsw-1
        ! wavefunction at t
        write(tmp, buf) i
        fileA = trim(rundir) // '/' // trim(adjustl(tmp)) // '/WAVECAR'
        ! wavefunction at t + dt
        write(tmp, buf) i + 1
        fileB = trim(rundir) // '/' // trim(adjustl(tmp)) // '/WAVECAR'

        call initAB(trim(fileA), trim(fileB), waveA, waveB)

        if (olap%NBANDS /= waveA%NBANDS) then
          write(*,*) "No. of bands does NOT match! WAVECAR: ", trim(adjustl(tmp))
          stop
        end if
        if (olap%ISPIN /= waveA%ISPIN) then
          write(*,*) "No. of spin components does NOT match! WAVECAR: ", &
                      trim(adjustl(tmp))
          stop
        end if
        do is=1, nspin
          call CoupIJ(waveA, waveB, olap%Dij(:,1+nb*(is-1):nb*is,i), is)
          olap%Eig(1+nb*(is-1):nb*is,i) = waveA%BANDS(:,1,is)
        end do

        call finishAB(waveA, waveB)
      end do
      ! After reading, write the couplings to disk
      call CoupToFile(olap)
      call copyToSec(olap, olap_sec, inp)
      call writeNaEig(olap_sec, inp)
    end if

    deallocate(olap%Dij, olap%Eig)
  end subroutine

  subroutine copyToSec(olap, olap_sec, inp)
    implicit none

    type(overlap), intent(in) :: olap
    type(overlap), intent(inout) :: olap_sec
    type(namdInfo), intent(in) :: inp

    integer, allocatable, dimension(:) :: basisU, basisD, basis
    integer :: i, nbasU, nbasD, nbas, nb

    if (inp%SOCTYPE==1) then

      olap_sec%Eig = olap%Eig(inp%BMIN:inp%BMAX, :)
      olap_sec%Dij = &
          olap%Dij(inp%BMIN:inp%BMAX, inp%BMIN:inp%BMAX, :)

    else if (inp%SOCTYPE==2) then

      nb = inp%NBANDS
      nbasU = inp%BMAXU - inp%BMINU + 1
      nbasD = inp%BMAXD - inp%BMIND+ 1
      nbas = nbasU + nbasD
      allocate(basisU(nbasU))
      allocate(basisD(nbasD))
      allocate(basis(nbas))
      basisU = (/(i, i=inp%BMINU,inp%BMAXU)/)
      basisD = (/(i, i=inp%BMIND,inp%BMAXD)/)
      basis = (/basisU,basisD+nb/)

      olap_sec%Dij = cero
      olap_sec%Eig = olap%Eig(basis,:)
      olap_sec%Dij(1:nbasU,1:nbasU,:) = olap%Dij(basisU,basisU,:)
      olap_sec%Dij(nbasU+1:nbas,nbasU+1:nbas,:) = &
          olap%Dij(basisD,basisD+nb,:)
      ! olap_sec%Sij = olap%Sij(basis, basis, :)

      deallocate(basisU, basisD, basis)

    end if

  end subroutine

  subroutine writeNaEig(olap, inp)
    implicit none

    type(overlap), intent(in) :: olap
    type(namdInfo), intent(in) :: inp
    integer :: i, j, it, Nt, nbas

    open(unit=22, file='EIGTXT', status='unknown', action='write')
    open(unit=23, file='NATXT', status='unknown', action='write')

    Nt = inp%NSW - 1
    nbas = olap%NBANDS

    do it=1, Nt
      write(unit=22, fmt='(*(G26.17))') (olap%Eig(i,it), i=1,nbas)
      write(unit=23, fmt='(*(G26.17))') &
          ((olap%Dij(i,j, it), j=1,nbas), i=1,nbas)
    end do

    close(unit=22)
    close(unit=23)
  end subroutine

  subroutine readNaEig(olap, inp)
    implicit none

    type(overlap), intent(inout) :: olap
    type(namdInfo), intent(in) :: inp
    integer :: i, j, it, Nt, nbas, ierr

    open(unit=22, file='EIGTXT', status='unknown', action='read', iostat=ierr)
    if (ierr /= 0) then
      write(*,*) "EIGTXT does NOT exist!"
      stop
    end if
    open(unit=23, file='NATXT', status='unknown', action='read', iostat=ierr)
    if (ierr /= 0) then
      write(*,*) "NATXT does NOT exist!"
      stop
    end if

    Nt = inp%NSW - 1
    nbas = olap%NBANDS

    do it=1, Nt
      read(unit=22, fmt=*) (olap%Eig(i,it), i=1,nbas)
      read(unit=23, fmt='(*(G26.17))') &
          ((olap%Dij(i,j, it), j=1,nbas), i=1,nbas)
    end do

    close(unit=22)
    close(unit=23)
  end subroutine

end module couplings
