module shop
  use prec
  use fileio
  use hamil
  implicit none

  contains

  ! initialize the random seed from the system clock
  ! code from: http://fortranwiki.org/fortran/show/random_seed
  subroutine init_random_seed()
    implicit none
    integer :: i, n, clock
    integer, dimension(:), allocatable :: seed

    call random_seed(size = n)
    allocate(seed(n))

    call system_clock(count=clock)

    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    call random_seed(put = seed)

    deallocate(seed)
  end subroutine

  subroutine whichToHop(cstat, ks)
    implicit none

    integer, intent(inout) :: cstat
    type(TDKS), intent(in) :: ks

    integer :: i
    real(kind=q) :: lower, upper, r

    call random_number(r)

    do i=1, ks%ndim
      if (i == 1) then
        lower = 0
        upper = ks%sh_prop(cstat,i)
      else
        lower = upper
        upper = upper + ks%sh_prop(cstat,i)
      end if
      if (lower <= r .AND. r < upper) then
        cstat = i
        exit
      end if
    end do

  end subroutine

  subroutine calcprop(tion, cstat, ks, inp)
    implicit none

    type(TDKS), intent(inout) :: ks
    type(namdInfo), intent(in) :: inp
    integer, intent(in) :: tion
    integer, intent(in) :: cstat

    integer :: i, j
    real(kind=q) :: Akk
    real(kind=q) :: dE, kbT

    Akk = CONJG(ks%psi_a(cstat, tion)) * ks%psi_a(cstat, tion)

    if (inp%SOCTYPE==1) then
      ! P(k -> m) = 2 Re [ CONJG(C_k) * C_m * d_km ] / ( CONJG(C_k) * C_k )
      ks%Bkm = 2. * REAL(  CONJG(ks%psi_a(cstat, tion)) * ks%psi_a(:, tion) * &
                           ks%NAcoup(cstat, :, tion) )
    else if (inp%SOCTYPE==2) then
      ! P(k -> m) = 2 Re [ CONJG(C_k) * C_m * d_km ] / ( CONJG(C_k) * C_k ) &
      !    - 2 / hbar Im [ CONJG(C_k) * C_m * S_km ] / ( CONJG(C_k) * C_k )
      ks%Bkm = 2. * REAL(  CONJG(ks%psi_a(cstat, tion)) * ks%psi_a(:, tion) * &
                           ks%NAcoup(cstat, :, tion) ) &
             - 2. * AIMAG( CONJG(ks%psi_a(cstat, tion)) * ks%psi_a(:, tion) * &
                           ks%SOcoup(cstat, :, tion) / hbar )
    end if

    ks%sh_prop(cstat,:) = ks%Bkm / Akk * inp%POTIM

    kbT = inp%TEMP * BOLKEV

    if (inp%LHOLE) then
      do i=1, ks%ndim
        dE = ks%eigKs(cstat, tion) - ks%eigKs(i,tion)
        if (dE>0) then
          ks%sh_prop(cstat,i) = ks%sh_prop(cstat,i) * exp(-dE / kbT)
        end if
      end do
    else
      do i=1, ks%ndim
        dE = ks%eigKs(i,tion) - ks%eigKs(cstat, tion)
        if (dE>0) then
          ks%sh_prop(cstat,i) = ks%sh_prop(cstat,i) * exp(-dE / kbT)
        end if
      end do
    end if

    forall (i=1:ks%ndim, ks%sh_prop(cstat,i) < 0) ks%sh_prop(cstat,i) = 0

  end subroutine

  ! calculate surface hopping probabilities
  subroutine runSH(ks, inp)
    implicit none

    type(TDKS), intent(inout) :: ks
    type(namdInfo), intent(in) :: inp
    integer :: i, j, ibas, nbas, tion
    integer :: istat, cstat
    integer, allocatable :: cstat_all(:)

    ks%sh_pops = 0
    ks%sh_prop = 0
    nbas = ks%ndim
    allocate(cstat_all(inp%NTRAJ))

    if (inp%SOCTYPE==1) then
      istat = inp%INIBAND - inp%BMIN + 1
    else if (inp%SOCTYPE==2) then
      if (inp%INISPIN == 1) then
        istat = inp%INIBAND - inp%BMINU + 1
      else
        istat = inp%INIBAND - inp%BMIND + inp%BMAXU - inp%BMINU + 2
      end if
    end if

    ! initialize the random seed for ramdom number production
    call init_random_seed()

    cstat_all = istat

    do tion=1, inp%NAMDTIME

      do ibas=1, nbas
        call calcprop(tion, ibas, ks, inp)
      end do

      do i=1, inp%NTRAJ
        cstat = cstat_all(i)
        call whichToHop(cstat, ks)
        cstat_all(i) = cstat
        ks%sh_pops(cstat, tion) = ks%sh_pops(cstat, tion) + 1
      end do

    end do

    ks%sh_pops = ks%sh_pops / inp%NTRAJ

  end subroutine

  ! need a subroutine here to write the results we need
  subroutine printSH(ks, inp)
    implicit none
    type(TDKS), intent(in) :: ks
    type(namdInfo), intent(in) :: inp

    integer :: i, j, tion, ierr, io
    character(len=48) :: buf

    write(buf, *) inp%NAMDTINI
    open(unit=24, file='SHPROP.' // trim(adjustl(buf)), status='unknown', action='write', iostat=ierr)
    open(unit=25, file='PSICT.' // trim(adjustl(buf)), status='unknown', action='write', iostat=ierr)
    if (ierr /= 0) then
      write(*,*) "SHPROP file I/O error!"
      stop
    end if

    do io = 24, 25

      if (inp%SOCTYPE==1) then
        write(io,'(A,A12,A3,I5)') '#', 'BMIN',    ' = ', inp%BMIN
        write(io,'(A,A12,A3,I5)') '#', 'BMAX',    ' = ', inp%BMAX
        write(io,'(A,A12,A3,I5)') '#', 'NBANDS',  ' = ', inp%NBANDS
        write(io,'(A,A12,A3,I5)') '#', 'INIBAND', ' = ', inp%INIBAND
        write(io,'(A,A12,A3,I5)') '#', 'SOCTYPE', ' = ', inp%SOCTYPE
      else if (inp%SOCTYPE==2) then
        write(io,'(A,A12,A3,I5)') '#', 'BMINU',   ' = ', inp%BMINU
        write(io,'(A,A12,A3,I5)') '#', 'BMAXU',   ' = ', inp%BMAXU
        write(io,'(A,A12,A3,I5)') '#', 'BMIND',   ' = ', inp%BMIND
        write(io,'(A,A12,A3,I5)') '#', 'BMAXD',   ' = ', inp%BMAXD
        write(io,'(A,A12,A3,I5)') '#', 'NBANDS',  ' = ', inp%NBANDS
        write(io,'(A,A12,A3,I5)') '#', 'INIBAND', ' = ', inp%INIBAND
        write(io,'(A,A12,A3,I5)') '#', 'INISPIN', ' = ', inp%INISPIN
        write(io,'(A,A12,A3,I5)') '#', 'SOCTYPE', ' = ', inp%SOCTYPE
      end if

      write(io,'(A,A12,A3,I5)')   '#', 'NSW',     ' = ', inp%NSW
      write(io,'(A,A12,A3,F5.1)') '#', 'POTIM',   ' = ', inp%POTIM
      write(io,'(A,A12,A3,F5.1)') '#', 'TEMP',    ' = ', inp%TEMP

      write(io,'(A,A12,A3,I5)') '#', 'NAMDTINI',  ' = ', inp%NAMDTINI
      write(io,'(A,A12,A3,I5)') '#', 'NAMDTIME',  ' = ', inp%NAMDTIME
      write(io,'(A,A12,A3,I5)') '#', 'NTRAJ',     ' = ', inp%NTRAJ
      write(io,'(A,A12,A3,I5)') '#', 'NELM',      ' = ', inp%NELM

      write(io,'(A,A12,A3,L5)') '#', 'LHOLE',     ' = ', inp%LHOLE
      write(io,'(A,A12,A3,L5)') '#', 'LSHP',      ' = ', inp%LSHP
      write(io,'(A,A12,A3,L5)') '#', 'LCPTXT',    ' = ', inp%LCPTXT
      write(io,'(A,A12,A3,A)')  '#', 'RUNDIR',    ' = ', TRIM(ADJUSTL(inp%rundir))
    end do

    do tion=1, inp%NAMDTIME
      write(unit=24, fmt=*) tion * inp%POTIM, SUM(ks%eigKs(:,tion) * ks%sh_pops(:,tion)), &
                            (ks%sh_pops(i,tion), i=1, ks%ndim)
      write(unit=25, fmt=*) tion * inp%POTIM, SUM(ks%eigKs(:,tion) * ks%pop_a(:,tion)), &
                            (ks%psi_a(i,tion), i=1, ks%ndim)
                            ! (ks%pop_a(i,tion), i=1, ks%ndim)
    end do

    close(24)
    close(25)

  end subroutine

end module
