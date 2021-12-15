module fileio
  use prec

  implicit none

  type namdInfo
    integer :: BMIN
    integer :: BMAX
    integer :: BMINU
    integer :: BMAXU
    integer :: BMIND
    integer :: BMAXD
    integer :: NBASIS      ! No. of adiabatic states as basis
    integer :: NBANDS      ! No. of band of the system
    integer :: INIBAND     ! inititial adiabatic state of excited electron/hole
    integer :: NSW         ! No. of MD steps
    integer :: NAMDTINI    ! Initial time step of NAMD
    integer :: NAMDTIME    ! No. of steps of NAMD
    integer :: INISPIN     ! Spin component of initial state
    integer, allocatable, dimension(:) :: NAMDTINI_A    ! No. of steps of NAMD
    integer, allocatable, dimension(:) :: INIBAND_A     ! No. of steps of NAMD
    integer, allocatable, dimension(:) :: INISPIN_A     ! Initial spin of NAMD
    integer :: NTRAJ       ! No. of surface hopping trajectories
    integer :: NELM        ! No. of steps of electron wave propagation
    integer :: NSAMPLE     ! No. of steps of electron wave propagation
    real(kind=q) :: POTIM  ! Time step of MD run
    real(kind=q) :: TEMP   ! MD Temperature

    integer :: SOCTYPE     ! which spin representation to use, default 1.
                           ! 1: spin-adiabatic representation;
                           ! 2: spin-diabatic representation.

    ! hole or electron surface hopping
    logical :: LHOLE
    ! whether to perform surface hopping, right now the value is .TRUE.
    logical :: LSHP
    logical :: LCPTXT
    ! running directories
    character(len=256) :: RUNDIR
    character(len=256) :: TBINIT

  end type

  contains

    subroutine getUserInp(inp)
      implicit none

      type(namdInfo), intent(inout) :: inp

      ! local variables with the same name as those in "inp"
      integer :: bmin
      integer :: bmax
      integer :: bminU
      integer :: bmaxU
      integer :: bminD
      integer :: bmaxD
      integer :: nsw
      integer :: iniband
      integer :: nbands
      integer :: soctype
      integer :: namdtime
      ! integer :: namdtini
      integer :: ntraj
      integer :: nelm
      integer :: nsample
      real(kind=q) :: potim
      real(kind=q) :: temp

      ! hole or electron surface hopping
      logical :: lhole
      ! surface hopping?
      logical :: lshp
      logical :: lcpext
      ! running directories
      character(len=256) :: rundir
      character(len=256) :: tbinit


      namelist /NAMDPARA/ bmin, bmax, nsw,    &
                          bminU, bmaxU,       &
                          bminD, bmaxD,       &
                          soctype, nbands,    &
                          potim, ntraj, nelm, &
                          temp, rundir,       &
                          lhole, lshp, lcpext,&
                          namdtime,           &
                          nsample, tbinit

      integer :: ierr, i
      logical :: lext

      ! set default values for thos parameters
      rundir = 'run'
      tbinit = 'INICON'
      bmin = 0
      bmax = 0
      bminU = 0
      bmaxU = 0
      bminD = 0
      bmaxD = 0
      nbands = 0
      soctype = 1
      ! iniband = 0
      ntraj = 1000
      nelm = 1000
      lhole = .FALSE.
      lshp = .TRUE.
      ! namdtini = 1
      namdtime = 200
      potim = 1.0_q
      temp = 300.
      lcpext = .FALSE.

      open(file="inp", unit=8, status='unknown', action='read', iostat=ierr)
      if ( ierr /= 0 ) then
        write(*,*) "I/O error with input file: 'inp'"
      end if

      read(unit=8, nml=NAMDPARA)
      close(unit=8)

      if (soctype==2) allocate(inp%INISPIN_A(nsample))
      allocate(inp%INIBAND_A(nsample), inp%NAMDTINI_A(nsample))
      inquire(file=tbinit, exist=lext)
      if (.NOT. lext) then
        write(*,*) "File containing initial conditions does NOT exist!"
      else
        open(unit=9, file=tbinit, action='read')
        do i=1, nsample
          if (soctype==1) then
            read(unit=9,fmt=*) inp%NAMDTINI_A(i), inp%INIBAND_A(i)
          else if (soctype==2) then
            read(unit=9,fmt=*) &
             inp%NAMDTINI_A(i), inp%INIBAND_A(i), inp%INISPIN_A(i)
          end if
        end do
        close(9)
      end if

      ! do some checking...
      ! put the following checks in the future version
      ! if (bmin <= 0 .OR. bmax <= 0 .OR. bmin >= bmax) then
      !   write(*,*) "Please specify the correct BMIN/BMAX"
      !   stop
      ! end if

      ! if (iniband == 0 .OR. iniband < bmin .OR. iniband > bmax) then
      !   write(*,*) "Please specify the correct initial band!"
      !   stop
      ! end if

      ! if (nbands == 0) then
      !   write(*,*) "I need the No. of bands..."
      !   stop
      ! end if

      ! if (namdtini + namdtime - 1 > nsw) then
      !   write(*,*) "NAMDTIME too long..."
      !   stop
      ! end if
      ! here ends the currently simplified version of parameter checking

      ! assign the parameters
      inp%BMIN     = bmin
      inp%BMAX     = bmax
      inp%BMINU    = bminU
      inp%BMAXU    = bmaxU
      inp%BMIND    = bminD
      inp%BMAXD    = bmaxD
      inp%NBASIS   = bmax - bmin + 1
      inp%NSW      = nsw
      inp%NBANDS   = nbands
      inp%SOCTYPE  = soctype
      inp%NAMDTIME = namdtime
      inp%NTRAJ    = ntraj
      inp%NELM     = nelm
      inp%LHOLE    = lhole
      inp%LSHP     = lshp
      inp%RUNDIR   = trim(rundir)
      inp%TBINIT   = trim(tbinit)
      inp%NSAMPLE  = nsample
      inp%POTIM    = potim
      inp%LCPTXT   = lcpext
      inp%TEMP     = temp

      if (inp%SOCTYPE==2) then
        inp%NBASIS = bmaxU - bminU + bmaxD - bminD + 2
      end if

    end subroutine

    ! Need a subroutine to print out all the input parameters
    subroutine printUserInp(inp)
      implicit none
      type(namdInfo), intent(in) :: inp

      write(*,'(A)') "------------------------------------------------------------"

      if (inp%SOCTYPE==1) then
        write(*,'(A30,A3,I5)') 'BMIN',    ' = ', inp%BMIN
        write(*,'(A30,A3,I5)') 'BMAX',    ' = ', inp%BMAX
        write(*,'(A30,A3,I5)') 'NBANDS',  ' = ', inp%NBANDS
        write(*,'(A30,A3,I5)') 'INIBAND', ' = ', inp%INIBAND
        write(*,'(A30,A3,I5)') 'SOCTYPE', ' = ', inp%SOCTYPE
      else if (inp%SOCTYPE==2) then
        write(*,'(A30,A3,I5)') 'BMINU',   ' = ', inp%BMINU
        write(*,'(A30,A3,I5)') 'BMAXU',   ' = ', inp%BMAXU
        write(*,'(A30,A3,I5)') 'BMIND',   ' = ', inp%BMIND
        write(*,'(A30,A3,I5)') 'BMAXD',   ' = ', inp%BMAXD
        write(*,'(A30,A3,I5)') 'NBANDS',  ' = ', inp%NBANDS
        write(*,'(A30,A3,I5)') 'INIBAND', ' = ', inp%INIBAND
        write(*,'(A30,A3,I5)') 'INISPIN', ' = ', inp%INISPIN
        write(*,'(A30,A3,I5)') 'SOCTYPE', ' = ', inp%SOCTYPE
      end if

      write(*,'(A30,A3,I5)')   'NSW',     ' = ', inp%NSW
      write(*,'(A30,A3,F5.1)') 'POTIM',   ' = ', inp%POTIM
      write(*,'(A30,A3,F5.1)') 'TEMP',    ' = ', inp%TEMP

      write(*,'(A30,A3,I5)') 'NAMDTINI',  ' = ', inp%NAMDTINI
      write(*,'(A30,A3,I5)') 'NAMDTIME',  ' = ', inp%NAMDTIME
      write(*,'(A30,A3,I5)') 'NTRAJ',     ' = ', inp%NTRAJ
      write(*,'(A30,A3,I5)') 'NELM',      ' = ', inp%NELM

      write(*,'(A30,A3,L5)') 'LHOLE',     ' = ', inp%LHOLE
      write(*,'(A30,A3,L5)') 'LSHP',      ' = ', inp%LSHP
      write(*,'(A30,A3,L5)') 'LCPTXT',    ' = ', inp%LCPTXT
      write(*,'(A30,A3,A)')  'RUNDIR',    ' = ', TRIM(ADJUSTL(inp%rundir))

      write(*,'(A)') "------------------------------------------------------------"
    end subroutine

end module
