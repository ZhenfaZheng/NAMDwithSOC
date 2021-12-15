module calcsoc
contains

subroutine readSocCar(RUNPATH, NPROJ, SOC)
!===============================================================================
! Read SOC from SocCar
!   in arg(s): 
!     - RUNPATH: string of the path where SocCar lies in
!     - NPROJ: Number of projectors, can be got from NormalCAR
!
!  out arg(s):
!     - SOC(NPROJ, NPROJ, 4)
!
! Author: ionizing
! Date: Sept 15, 2019
!===============================================================================

  implicit none
  character(*), intent(in) :: RUNPATH
  complex(8), allocatable, intent(out) :: SOC(:,:,:)
  integer, intent(in) :: NPROJ
  integer, parameter :: IU = 18
  integer :: IERR

  integer :: i, j, k

  open(unit=IU, file=RUNPATH//'SocCar', action='read', form='formatted', status='old',&
    & iostat=IERR)
  if (0 /= IERR) stop 'Error opening file "SocCar"'

  allocate(SOC(NPROJ, NPROJ, 4))
  do k=1, 4
    do i=1, NPROJ
      do j=1, NPROJ
#ifdef ZZF_VERSION
        read(unit=IU, fmt='(D14.7,D14.7)', iostat=IERR, advance='no') SOC(i, j, k)
#else
        read(unit=IU, fmt='(F22.16,F22.16)', iostat=IERR, advance='no') SOC(i, j, k)
#endif
      enddo
      read(unit=IU, fmt='(A1)', iostat=IERR)
    enddo
    read(unit=IU, fmt='(A1)', iostat=IERR)
  enddo

  if (0 /= IERR) then
    stop 'Error reading SOC from SocCar'
  endif

  close(IU)
end subroutine readSocCar



subroutine readNormalCAR(RUNPATH, NBANDS, NKPTS, IKPT, CPROJ, NPROJ)
!===============================================================================
! Read NormalCAR and get CPRJ and NPROJ
!   in arg(s):
!     - RUNPATH: string of the path where NormalCAR lies in
!     - NBANDS: number of total bads
!     - NKPTS:  number of total k-points
!     - IKPT:   which k-point to deal with, start from 1
!
!  out arg(s):
!     - CPROJ(NPROJ, NBANDS, 2), *NOT ALLOCATED IN ADVANCE*, 2 = NSPIN
!     - NRPOJ: number of projectors
!
! Author: ionizing
! Date: Sept 16, 2019
!===============================================================================
#ifdef __INTEL_COMPILER
  use ifport
#endif

  implicit none
  character(*), intent(in) :: RUNPATH
  integer, intent(in)  :: NBANDS, NKPTS, IKPT
  integer, intent(out) :: NPROJ
  complex(8), allocatable, intent(out) :: CPROJ(:, :, :)

  integer, parameter :: IU = 16
  integer :: IERR, statb(13)
  integer :: rec_l, recl_

  integer :: LMDIM, NIONS, NRSPINORS
  integer :: NPROD, NPRO, NTYP
  real(8), allocatable :: CQIJ(:,:,:,:)
  complex(8), allocatable :: CPROJ_T(:,:,:,:)
  integer, allocatable :: LMMAX(:), NITYP(:)
  

  integer :: itype
  integer :: ispin, iikpt, iband

  open(unit=IU, file=RUNPATH//'NormalCAR', status='old', action='read', iostat=IERR,&
    & form='unformatted', access='stream')
  IERR = fstat(IU, statb)
  if (0 /= IERR) stop 'Error opening NormalCAR.'

  read(unit=IU, iostat=IERR) rec_l, LMDIM, NIONS, NRSPINORS, recl_
  if (rec_l /= recl_ .or. IERR /= 0) stop "88 wrong rec length"

  allocate(CQIJ(LMDIM, LMDIM, NIONS, NRSPINORS))
  read(unit=IU, iostat=IERR) rec_l, CQIJ(1:LMDIM, 1:LMDIM, 1:NIONS, 1:NRSPINORS), recl_
  if (rec_l /= recl_ .or. IERR /= 0) stop "92 wrong rec length"
  deallocate(CQIJ)

  read(unit=IU, iostat=IERR) rec_l, NPROD, NPRO, NTYP, recl_
  if (rec_l /= recl_ .or. IERR /= 0) stop "96 wrong rec length"

  allocate(LMMAX(NTYP))
  allocate(NITYP(NTYP))
  do itype=1, NTYP
    read(unit=IU, iostat=IERR) rec_l, LMMAX(itype), NITYP(itype), recl_
  enddo
  if (rec_l /= recl_) stop "103 wrong rec length"
  deallocate(NITYP)
  deallocate(LMMAX)

  read(unit=IU, iostat=IERR) rec_l
  NPROJ = rec_l / 16
  ! print *, 'NPROJ = ', NPROJ
#ifdef __INTEL_COMPILER
  IERR = fseek(lunit=IU, offset=-4, from=1)
#else
  call fseek(unit=IU, offset=-4, whence=1, status=IERR)
#endif
  allocate(CPROJ_T(NPROJ, NBANDS, NKPTS, 2))
  do ispin=1, 2
    do iikpt=1, NKPTS
      do iband=1, NBANDS
        read(unit=IU, iostat=IERR) rec_l, CPROJ_T(1:NPROJ, iband, iikpt, ispin), recl_
        if (rec_l /= recl_) stop '116 wrong rec length'
      enddo
    enddo
  enddo

  allocate(CPROJ(NPROJ, NBANDS, 2))
  CPROJ(:, :, :) = CPROJ_T(:, :, ikpt, :)
  deallocate(CPROJ_T)
  close(IU)
end subroutine readNormalCAR

subroutine calcHmm(RUNPATH, NBANDS, NKPTS, IKPT, HMM)
!==============================================================================
! calculate SOC matrix, needs readingSocCar and NormalCAR
!   in arg(s):
!     - RUNPATH: string of the folder where SocCar and NormalCAR be in
!     - NBANDS: number of bands in current system
!     - NKPTS:  number of k-points in current system
!     - IKPT:   which k-point to deal with, start from 1
!  out arg(s):
!     - HMM:    SOC matrix organized as the following shows
!
!                 uu | ud
!                 -------
!                 du | dd
!
! Author: ionizing
! Date: Sept 17, 2019
!==============================================================================

  implicit none
  character(*), intent(in)             :: RUNPATH
  integer, intent(in)                  :: NBANDS, NKPTS, IKPT
  complex(8), intent(inout) :: Hmm(:, :)

  integer :: NPROJ
  complex(8), allocatable :: CPROJ(:, :, :), SOC(:,:,:)
  ! CPROJ(NPROJ, NBANDS, 2)
  ! SOC(NPROJ, NPROJ, 4)
  complex(8), allocatable :: tmp(:, :)

  call readNormalCAR(RUNPATH, NBANDS, NKPTS, IKPT, CPROJ, NPROJ)
  call readSocCar(RUNPATH, NPROJ, SOC)

  ! print *, "NBANDS = ", NBANDS
  ! print *, "NPROJ = ", NPROJ

  if (rank(Hmm) /= 2 .or. any(shape(Hmm) /= (/2*NBANDS, 2*NBANDS/))) &
    &stop 'shape of pre-allocated Hmm not correct.'
  allocate(tmp(NBANDS, NBANDS))
  Hmm = (0, 0)

  tmp = (0, 0)
  tmp = matmul( conjg(transpose(CPROJ(:,:,1))), &
    &         matmul(SOC(:,:,1),CPROJ(:,:,1)) )
  Hmm(1:NBANDS, 1:NBANDS) = tmp

  tmp = (0, 0)
  tmp = matmul( conjg(transpose(CPROJ(:,:,2))), &
    &        matmul(SOC(:,:,4), CPROJ(:,:,2)))
  Hmm(NBANDS+1:2*NBANDS, NBANDS+1:2*NBANDS) = tmp

  tmp = (0, 0)
  tmp = matmul( conjg(transpose(CPROJ(:,:,1))), &
    &        matmul(SOC(:,:,2), CPROJ(:,:,2)))
  Hmm(1:NBANDS, NBANDS+1:2*NBANDS) = tmp

  tmp = (0, 0)
  tmp = matmul( conjg(transpose(CPROJ(:,:,2))), &
    &        matmul(SOC(:,:,3), CPROJ(:,:,1)))
  Hmm(NBANDS+1:2*NBANDS, 1:NBANDS) = tmp

  deallocate(tmp)
end subroutine calcHmm

end module calcsoc
