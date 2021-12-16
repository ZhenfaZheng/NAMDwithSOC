#include "symbol.inc"
!#define debug
      MODULE RELATIVISTIC
      USE prec
      
      REAL(q), ALLOCATABLE, PRIVATE, SAVE :: SPINORB_MATRIX_ELEMENTS(:,:,:)
      REAl(q), ALLOCATABLE, PRIVATE, SAVE :: SPINORB_ENERGY_PER_SITE(:)
      

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! ADDED BY QijingZheng
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! < f_i(r) | d V_KS(r) / dr * r^{-1} | f_j(r) >
      REAl(q), ALLOCATABLE, PRIVATE, SAVE :: SPINORB_VKS_RAD(:,:,:)
      ! < Y_i; \sigma_1 | S L | Y_j; \sigma_2 > < f_i(r) | d V_KS(r) / dr * r^{-1} | f_j(r) >
      COMPLEX(q), ALLOCATABLE, PRIVATE, SAVE :: SPINORB_VKS(:,:,:,:)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CONTAINS
      
!*******************************************************************

      SUBROUTINE SPINORB_STRENGTH(POT, RHOC, POTVAL, R, DLLMM, CHANNELS, L, W, Z, THETA, PHI)

!*******************************************************************
!
!  the potential is given by
!       V(r) =  \sum_lm pot_lm(r) * Y_lm(r)
!  where   pot_lm(r) is stored in POT(2l+1+m,..) l=0,..,LMAX, m=0,..,2*l
!  we use only the radial component pot_00(r)
!
!  the wavefunction psi(r) can be obtained from the stored
!  coefficients using:
!      psi(r) = \sum_lmn Y_lm(r) w_ln(r) r
!*******************************************************************
      USE prec
      USE constant
      USE radial
      IMPLICIT NONE

      REAL(q) POT(:)         ! spherical contribution to potential w.r.t reference potential
      REAL(q) RHOC(:)        ! electronic core charge
      REAL(q) POTVAL(:)      ! minus potential of atom
      TYPE(rgrid) :: R
      REAL(q) W(:,:)         ! wavefunctions phi(r,l)
      OVERLAP DLLMM(:,:,:)   ! contribution to H from so-coupling
      INTEGER CHANNELS, L(:) 
      REAL(q) THETA,PHI      ! Euler angle
      REAL(q) Z              ! charge of the nucleus
! local
      INTEGER I,J,LM,LMP,M,MP,CH1,CH2,LL,LLP
      REAL(q) APOT(R%NMAX)   ! average potential (up down)
      REAL(q) DPOT(R%NMAX)   ! radial derivative   of potential APOT
      REAL(q) RHOT(R%NMAX)   ! charge density
      REAL(q) ksi(R%NMAX)
      REAL(q), PARAMETER :: C = 137.037  ! speed of light in a.u.
      REAL(q), PARAMETER :: INVMC2=7.45596E-6
!                           invmc2=hbar^2/2(m_e c)^2 in A^2
      INTEGER, PARAMETER :: LMAX=3, MMAX=LMAX*2+1
      COMPLEX(q) DUMMY(MMAX,MMAX,3,LMAX)
      COMPLEX(q) LS(MMAX,MMAX,4,LMAX)
      REAL(q) SUM, SCALE

      LS=(0._q,0._q)
      CALL SETUP_LS(1,THETA,PHI,DUMMY(1:3,1:3,1:3,1),LS(1:3,1:3,1:4,1))
      CALL SETUP_LS(2,THETA,PHI,DUMMY(1:5,1:5,1:3,2),LS(1:5,1:5,1:4,2))
      CALL SETUP_LS(3,THETA,PHI,DUMMY(1:7,1:7,1:3,3),LS(1:7,1:7,1:4,3))
      
!     thats just Y_00
      SCALE=1/(2*SQRT(PI))
!     unfortunately the PAW method operates usually only with
!     difference potentials (compared to isolated atom)
!     we need to evaluate a couple of terms

!     lets first calculate the Hartree potential of the core electrons
      CALL RAD_POT_HAR(0, R, APOT, RHOC, SUM)
!     add the potential of the nucleus (delta like charge Z at origin)
      APOT=APOT*SCALE - FELECT/R%R*Z
!     subtract reference potential POTVAL (previously added to POT(:,:) (see RAD_POT)
!     this one contains essentially valence only contributions
      APOT=APOT-POTVAL
!     finally add the current potential (average spin up and down)
      APOT=APOT+POT(1:SIZE(APOT))* SCALE
!     gradient
      CALL GRAD(R,APOT,DPOT)
!     ksi(r)=  hbar^2/2(m_e c)^2 1/r d V(r)/d r
!     KSI(:)=INVMC2*DPOT(:)/ R%R 

!     ZORA: includes a radially variant relativistic mass correction 1/(1-V(r)/2/c^2)^2 in Hartree atomic units,
!     in accordance with for instance Eq. 42 of Lenthe et al., J. Chem. Phys. 99, 4597 (1993) and Eq. 6.16 in
!     "the ZORA equation", Thesis of E. van Lenthe.)
      DO I=1,R%NMAX
!        KSI(I)=INVMC2*(RYTOEV/(RYTOEV-0.5_q*APOT(I)/C/C))*DPOT(I)/R%R(I)
         KSI(I)=INVMC2*(1._q/(1._q-APOT(I)/(2._q*RYTOEV)/C/C/2._q)**2)*DPOT(I)/R%R(I)
      ENDDO

#ifdef debug
      DO I=1,R%NMAX 
         WRITE(78,'(4F20.7)') R%R(I),APOT(I)*R%R(I),POTVAL(I)*R%R(I),KSI(I)*R%R(I)
      ENDDO
      REWIND(78)
#endif

!     calculates the integral
!     D(ll,LM) =  D(ll,LM) 
!'        + \int dr  w_ln(r)  ksi(r)  w_ln'(r) 
!         * \int dOmega  Y_lm LS Y_lm
!     on the radial grid, inside the augmentation sphere only
 
      LM =1
      DO CH1=1,CHANNELS
      LMP=1
      DO CH2=1,CHANNELS
        DO I=1,R%NMAX
           RHOT(I)=W(I,CH1)*W(I,CH2)
        END DO
        LL = L(CH1)
        LLP= L(CH2)
!     calculation is made only for l=2 and l=1 orbitals
!     a spherical potential is assumed
        IF (LL == LLP .AND. LL>0 .AND. LL<=LMAX ) THEN
          SUM=0
          DO I=1,R%NMAX 
!      The integral is made only inside the augmentation sphere
!            IF(R%R(I) <= R%RMAX) THEN
              SUM= SUM+KSI(I)*RHOT(I)*R%SI(I)
!            ENDIF
          END DO
          SUM=SUM
!
! VASP uses a reverted notation (for efficiency reason)
!  D(lm,l'm',alpha+2*alpha') =  <alpha'| < y_l'm' | D | y_lm>  |alpha>
! therefore we need a little bit of reindexing (not too complicated)
          DO I=0,1
          DO J=0,1
          DO M =1,2*LL+1
          DO MP=1,2*LL+1
             DLLMM(LMP+MP-1,LM+M-1,J+2*I+1)=DLLMM(LMP+MP-1,LM+M-1,J+2*I+1)+ &
             SUM*LS(M,MP,I+2*J+1,LL)
          END DO
          END DO
          END DO
          END DO
        ENDIF

      LMP=LMP+(2*LLP+1)
      ENDDO
      LM= LM+ (2*LL+1)
      ENDDO

      END SUBROUTINE SPINORB_STRENGTH 

!*******************************************************************
!
!  the potential is given by
!       V(r) =  \sum_lm pot_lm(r) * Y_lm(r)
!  where   pot_lm(r) is stored in POT(2l+1+m,..) l=0,..,LMAX, m=0,..,2*l
!  we use only the radial component pot_00(r)
!
!  the wavefunction psi(r) can be obtained from the stored
!  coefficients using:
!      psi(r) = \sum_lmn Y_lm(r) w_ln(r) r
!*******************************************************************

      SUBROUTINE SPINORB_STRENGTH_KAKA(POT, RHOC, POTVAL, R, CHANNELS, &
                                       NI, NIONS, L, W, Z, THETA, PHI)
      USE prec
      USE constant
      USE radial
      USE wave
      IMPLICIT NONE

      INTEGER :: NI          ! the index of the atom
      INTEGER :: NIONS       ! number of atoms
      REAL(q) POT(:)         ! spherical contribution to potential w.r.t reference potential
      REAL(q) RHOC(:)        ! electronic core charge
      REAL(q) POTVAL(:)      ! minus potential of atom
      TYPE(rgrid) :: R
      REAL(q) W(:,:)         ! wavefunctions phi(r,l)
      INTEGER CHANNELS, L(:) 
      REAL(q) THETA,PHI      ! Euler angle
      REAL(q) Z              ! charge of the nucleus

      INTEGER I,J,LM,LMP,M,MP,CH1,CH2,LL,LLP
      REAL(q) APOT(R%NMAX)   ! average potential (up down)
      REAL(q) DPOT(R%NMAX)   ! radial derivative   of potential APOT
      REAL(q) RHOT(R%NMAX)   ! charge density
      REAL(q) ksi(R%NMAX)
      REAL(q), PARAMETER :: C = 137.037       ! speed of light in a.u.
      REAL(q), PARAMETER :: INVMC2=7.45596E-6 ! invmc2=hbar^2/2(m_e c)^2 in eV/A^2
      INTEGER, PARAMETER :: LMAX=3, MMAX=LMAX*2+1
      COMPLEX(q) DUMMY(MMAX,MMAX,3,LMAX)
      COMPLEX(q) LS(MMAX,MMAX,4,LMAX)
      REAL(q) SUM, SCALE

      LS=(0._q,0._q)
      CALL SETUP_LS(1,THETA,PHI,DUMMY(1:3,1:3,1:3,1),LS(1:3,1:3,1:4,1))
      CALL SETUP_LS(2,THETA,PHI,DUMMY(1:5,1:5,1:3,2),LS(1:5,1:5,1:4,2))
      CALL SETUP_LS(3,THETA,PHI,DUMMY(1:7,1:7,1:3,3),LS(1:7,1:7,1:4,3))

      ! allocate if necessary
      IF (.NOT. ALLOCATED(SPINORB_VKS_RAD)) THEN
         ALLOCATE(SPINORB_VKS_RAD(16,16,NIONS))
         SPINORB_VKS_RAD = 0._q
         ! write(*,*) "rott", "ALLOCATED"
      ENDIF
      IF (.NOT. ALLOCATED(SPINORB_VKS)) THEN
         ALLOCATE(SPINORB_VKS(64,64,4,NIONS))
         SPINORB_VKS = (0._q, 0._q)
         ! write(*,*) "kaka", "ALLOCATED"
      ENDIF
      SPINORB_VKS_RAD(:,:,NI) = 0._q
      SPINORB_VKS(:,:,:,NI)   = (0._q, 0._q)

      ! thats just Y_00
      SCALE=1/(2*SQRT(PI))
      ! unfortunately the PAW method operates usually only with
      ! difference potentials (compared to isolated atom)
      ! we need to evaluate a couple of terms

      ! lets first calculate the Hatree potential of the core electrons
      CALL RAD_POT_HAR(0, R, APOT, RHOC, SUM)
      ! add the potential of the nucleus (delta like charge Z at origin)
      APOT=APOT*SCALE - FELECT/R%R*Z
      ! subtract reference potential POTVAL (previously added to POT(:,:) (see RAD_POT)
      ! this one contains essentially valence only contributions
      APOT=APOT-POTVAL
      ! finally add the current potential (average spin up and down)
      APOT=APOT+POT(1:SIZE(APOT))* SCALE
      ! gradient
      CALL GRAD(R,APOT,DPOT)

!     ZORA: includes a radially variant relativistic mass correction 1/(1-V(r)/2/c^2)^2 in Hartree atomic units,
!     in accordance with for instance Eq. 42 of Lenthe et al., J. Chem. Phys. 99, 4597 (1993) and Eq. 6.16 in
!     "the ZORA equation", Thesis of E. van Lenthe.)
      DO I=1,R%NMAX
!        KSI(I)=INVMC2*(RYTOEV/(RYTOEV-0.5_q*APOT(I)/C/C))*DPOT(I)/R%R(I)
         KSI(I)=INVMC2*(1._q/(1._q-APOT(I)/(2._q*RYTOEV)/C/C/2._q)**2)*DPOT(I)/R%R(I)
      ENDDO

#ifdef debug
      DO I=1,R%NMAX 
         WRITE(78,'(4F20.7)') R%R(I),APOT(I)*R%R(I),POTVAL(I)*R%R(I),KSI(I)*R%R(I)
      ENDDO
      REWIND(78)
#endif

      !     calculates the integral
      !     D(ll,LM) =  D(ll,LM) 
      !         + \int dr  w_ln(r)  ksi(r)  w_ln'(r) 
      !         * \int dOmega  Y_lm LS Y_lm
      !     on the radial grid, inside the augmentation sphere only
 
      LM =1
      DO CH1=1,CHANNELS
        LMP=1
        LL = L(CH1)
        DO CH2=1,CHANNELS
          LLP= L(CH2)

          ! calculation is made only for l=2 and l=1 orbitals
          ! a spherical potential is assumed
          IF (LL == LLP .AND. LL>0 .AND. LL<=LMAX ) THEN

            DO I=1,R%NMAX
               RHOT(I)=W(I,CH1)*W(I,CH2)
            END DO

            SUM=0
            DO I=1,R%NMAX 
               ! The integral is made only inside the augmentation sphere
               ! IF(R%R(I) <= R%RMAX) &
                SUM= SUM+KSI(I)*RHOT(I)*R%SI(I)
            END DO
            ! write(*,*) "Qijing", NI, CH1, CH2, SUM
            SPINORB_VKS_RAD(CH1,CH2,NI) = SUM

            ! VASP uses a reverted notation (for efficiency reason)
            !  D(lm,l'm',alpha+2*alpha') =  <alpha'| < y_l'm' | D | y_lm>  |alpha>
            ! therefore we need a little bit of reindexing (not too complicated)

            DO I=0,1
              DO J=0,1
                DO M =1,2*LL+1
                  DO MP=1,2*LL+1
                     SPINORB_VKS(LMP+MP-1,LM+M-1,J+2*I+1,NI) = SPINORB_VKS(LMP+MP-1,LM+M-1,J+2*I+1,NI) + &
                                                               SUM * LS(M,MP,I+2*J+1,LL)
                     ! if (NI == 6) then
                     !   write(*,*) NI, LL, I, J
                     !   write(*,*) SUM * LS(:,:,I+2*J+1,LL)
                     ! end if
                  END DO
                END DO
              END DO
            END DO

          ENDIF

          LMP=LMP+(2*LLP+1)
        ENDDO
        LM= LM+ (2*LL+1)
      ENDDO

      END SUBROUTINE SPINORB_STRENGTH_KAKA 

      SUBROUTINE WRT_SPINORB_AEPWFC_KAKA(WDES, T_INFO, IO, PP)
        USE prec
        USE wave
        USE base
        USE poscar
        USE pseudo

        IMPLICIT NONE

        TYPE (wavedes)    :: WDES
        TYPE (type_info)  :: T_INFO
        TYPE (in_struct)  :: IO
        TYPE (potcar)     :: PP(T_INFO%NTYP)

        ! local variables
        INTEGER :: NI, I, J, K
        INTEGER :: NPROJ
        INTEGER :: ILOWER, IUPPER, ITP, LMMAX, LMAX
        COMPLEX(q), ALLOCATABLE :: DUMP(:,:,:)
        CHARACTER(len=1024)     :: fbuff

        ! allocate if necessary
        IF (.NOT.ALLOCATED(SPINORB_VKS_RAD)) THEN
           ALLOCATE(SPINORB_VKS_RAD(16,16,T_INFO%NIONS))
           SPINORB_VKS_RAD = 0._q
        ENDIF
        IF (.NOT.ALLOCATED(SPINORB_VKS)) THEN
           ALLOCATE(SPINORB_VKS(64,64,4,T_INFO%NIONS))
           SPINORB_VKS = (0._q, 0._q)
        ENDIF

        ! communicate
        ! write(*,*) WDES%COMM%NODE_ME, "======", IO%IU6
        CALLMPI( M_sum_z(WDES%COMM, SPINORB_VKS,     64*64*4*T_INFO%NIONS) )
        CALLMPI( M_sum_d(WDES%COMM, SPINORB_VKS_RAD, 16*16*T_INFO%NIONS) )

        IF (IO%IU6>0) THEN
          NPROJ = 0
          DO I=1, T_INFO%NTYP
            NPROJ = NPROJ + PP(I)%LMMAX * T_INFO%NITYP(I)
            ! write(*,*) "kaka ", I, PP(I)%LMMAX, T_INFO%TYPE(I), NPROJ
          ENDDO
          write(fbuff,*) "(", NPROJ * 2, "F22.16)"
          ALLOCATE(DUMP(NPROJ,NPROJ,4))
          
          DUMP = (0.0_q, 0.0_q)

          do K=1, 4
            ILOWER = 1
            do I=1, T_INFO%NIONS
              ITP    = T_INFO%ITYP(I)
              LMMAX  = PP(ITP)%LMMAX
              IUPPER = ILOWER + LMMAX - 1
              ! write(*,*) "kaka", ITP, LMMAX, ILOWER, IUPPER
              DUMP(ILOWER:IUPPER,ILOWER:IUPPER,K) = DUMP(ILOWER:IUPPER,ILOWER:IUPPER,K) + &
                                                    SPINORB_VKS(1:LMMAX,1:LMMAX,K,I)
              ILOWER = IUPPER + 1
            end do
          end do

          ! spin-orbit matrix of shape (nproj, nproj, 4)
          OPEN(FILE='SocCar', unit=2000, status='unknown')
          DO K=1, 4
            DO I=1, NPROJ
              write(2000, TRIM(fbuff)) (DUMP(J,I,K),J=1, NPROJ)
            END DO
            write(2000,*)
          END DO
          close(2000)

          ! spin-orbit matrix of shape (64, 64, 4, NIONS)
          OPEN(FILE='SocAllCar', unit=2100, status='unknown')
          DO NI=1, T_INFO%NIONS
            DO K=1,4
              DO I=1, 64
                write(2100, '(128F22.16)') (SPINORB_VKS(I,J,K,NI), J=1,64)
              END DO
              write(2100, *) 
            END DO
            write(2100, *) 
          END DO
          close(2100)

          ! radial spin-orbit matrix of shape (16, 16, 4, NIONS)
          OPEN(FILE='SocRadCar', unit=2200, status='unknown')
          
          write(2200, '(A)') "L_START"
          ! Number_of_Atom_Types, Number_of_Atoms
          write(2200, '(3I5)') T_INFO%NTYP, T_INFO%NIONS, NPROJ
          LMAX = 0
          do ITP=1, T_INFO%NTYP
            if (LMAX < PP(ITP)%LMAX) then
              LMAX = PP(ITP)%LMAX  ! Number of projectors for each elements
            end if
            ! write quantum number L of the projector function for each
            ! type of atom
            write(fbuff,*) "(", PP(ITP)%LMAX, "I3)"
            write(2200, TRIM(fbuff)) (PP(ITP)%LPS(J), J=1, PP(ITP)%LMAX)
          end do
          ! Chemical symbol for each atom
          write(fbuff,*) "(1X,", T_INFO%NIONS, "A3)"
          write(2200, TRIM(fbuff)) (T_INFO%TYPE(T_INFO%ITYP(J)), J=1, T_INFO%NIONS)
          ! Type index for each atom
          write(fbuff,*) "(", T_INFO%NIONS, "I3)"
          write(2200, TRIM(fbuff)) (T_INFO%ITYP(J), J=1, T_INFO%NIONS)
          write(2200, '(A)') "L_END"

          DO NI=1, T_INFO%NIONS
            DO I=1, LMAX
              ! write the radial soc matrix elements
              write(fbuff,*) "(", 16, "F22.16)"
              write(2200, fbuff) (SPINORB_VKS_RAD(I,J, NI), J=1, LMAX)
            END DO
            write(2200, *) 
          END DO
          close(2200)

        END IF

      END SUBROUTINE WRT_SPINORB_AEPWFC_KAKA

!**********************************************************************
!
!**********************************************************************
      SUBROUTINE CALC_SPINORB_MATRIX_ELEMENTS(WDES,PP,T_INFO,NI,CSO,COCC)
      USE prec
      USE wave
      USE pseudo
      USE poscar
      IMPLICIT NONE
      TYPE (wavedes)  WDES
      TYPE (potcar) PP
      TYPE (type_info) T_INFO
      INTEGER NI
      OVERLAP CSO(:,:,:),COCC(:,:,:)
      ! local variables
      INTEGER CH1,CH2,LM,LMP,LL,LLP,M,MP
      INTEGER ISPINOR1,ISPINOR2
           
      ! quick return
      IF (WDES%NCDIJ/=4.OR.(.NOT.WDES%LSORBIT)) RETURN
      
      ! allocate if necessary
      IF (.NOT.ALLOCATED(SPINORB_MATRIX_ELEMENTS)) THEN
         ALLOCATE(SPINORB_MATRIX_ELEMENTS(16,16,T_INFO%NIONS))
         SPINORB_MATRIX_ELEMENTS=0
      ENDIF
      IF (.NOT.ALLOCATED(SPINORB_ENERGY_PER_SITE)) THEN
         ALLOCATE(SPINORB_ENERGY_PER_SITE(T_INFO%NIONS))
         SPINORB_ENERGY_PER_SITE=0
      ENDIF
      
      SPINORB_MATRIX_ELEMENTS(:,:,NI)=0
      SPINORB_ENERGY_PER_SITE(NI)=0
      
      LM=1      
      DO CH1=1,PP%LMAX
      LMP=1
      DO CH2=1,PP%LMAX

      LL=PP%LPS(CH1); LLP=PP%LPS(CH2)
         
      IF (LL==LLP) THEN 
      
         DO ISPINOR1=0,1
         DO ISPINOR2=0,1
         DO M=1,2*LL+1
         DO MP=1,2*LLP+1
            SPINORB_MATRIX_ELEMENTS(LL*LL+M,LLP*LLP+MP,NI)= &
           &   SPINORB_MATRIX_ELEMENTS(LL*LL+M,LLP*LLP+MP,NI)+ &
           &   REAL(CSO(LMP+MP-1,LM+M-1,ISPINOR2+2*ISPINOR1+1)* &
#ifdef realmode
           &   COCC(LMP+MP-1,LM+M-1,ISPINOR2+2*ISPINOR1+1),q)
#else
           &   CONJG(COCC(LMP+MP-1,LM+M-1,ISPINOR2+2*ISPINOR1+1)),q)
#endif
            SPINORB_ENERGY_PER_SITE(NI)=SPINORB_ENERGY_PER_SITE(NI)+ &
           &   REAL(CSO(LMP+MP-1,LM+M-1,ISPINOR2+2*ISPINOR1+1)* &
#ifdef realmode
           &   COCC(LMP+MP-1,LM+M-1,ISPINOR2+2*ISPINOR1+1),q)
#else
           &   CONJG(COCC(LMP+MP-1,LM+M-1,ISPINOR2+2*ISPINOR1+1)),q)
#endif            
         ENDDO
         ENDDO
         ENDDO
         ENDDO
         
      ENDIF
         
      LMP=LMP+(2*LLP+1)
      ENDDO
      LM=LM+(2*LL+1)
      ENDDO
      
      RETURN
      END SUBROUTINE CALC_SPINORB_MATRIX_ELEMENTS


!**********************************************************************
!
!**********************************************************************
      SUBROUTINE WRITE_SPINORB_MATRIX_ELEMENTS(WDES,T_INFO,IO)
      USE prec
      USE wave
      USE base
      USE poscar
      IMPLICIT NONE
      TYPE (wavedes) WDES
      TYPE (type_info) T_INFO
      TYPE (in_struct) IO
      ! local variables
      INTEGER NI,L,M
      
      ! quick return
      IF (WDES%NCDIJ/=4.OR.(.NOT.WDES%LSORBIT)) RETURN
      IF (.NOT.ALLOCATED(SPINORB_MATRIX_ELEMENTS)) THEN
         ALLOCATE(SPINORB_MATRIX_ELEMENTS(16,16,T_INFO%NIONS))
         SPINORB_MATRIX_ELEMENTS=0
      ENDIF
      IF (.NOT.ALLOCATED(SPINORB_ENERGY_PER_SITE)) THEN
         ALLOCATE(SPINORB_ENERGY_PER_SITE(T_INFO%NIONS))
         SPINORB_ENERGY_PER_SITE=0
      ENDIF
      
      ! communicate
      CALLMPI( M_sum_d(WDES%COMM,SPINORB_MATRIX_ELEMENTS,16*16*T_INFO%NIONS) )
      CALLMPI( M_sum_d(WDES%COMM,SPINORB_ENERGY_PER_SITE,T_INFO%NIONS) )
      
      IF (IO%IU6>0) THEN
         WRITE(IO%IU6,'(/A)') ' Spin-Orbit-Coupling matrix elements'
         DO NI=1,T_INFO%NIONS
            WRITE(IO%IU6,'(/A,I4,A,F14.7)') ' Ion: ',NI,'  E_soc: ',SPINORB_ENERGY_PER_SITE(NI)
            DO L=1,3
               WRITE(IO%IU6,'(A,I4)') ' l=',L
               DO M=L*L+1,L*L+2*L+1
                  WRITE(IO%IU6,'(7F14.7)') SPINORB_MATRIX_ELEMENTS(L*L+1:L*L+2*L+1,M,NI)
               ENDDO
            ENDDO
         ENDDO
      ENDIF
      
      RETURN
      END SUBROUTINE WRITE_SPINORB_MATRIX_ELEMENTS


!**********************************************************************
!
! calculate the LS operator for arbitrary l-quantum number L 
! assuming a spin quantization axis rotated by \theta and \phi
! with respect to the z-axis (n.b. first LS is calculated assuming
! a quantization axis parallel to z, and then the matrix is
! rotated)
!
! LS(m1,m2,alpha1+alpha2*2+1)= <alpha1| <y_lm1| LS |y_lm2> |alpha2>
!
! with 
!  
! alpha1, alpha2 either 0 (=spinor up comp.) or 1 (=spinor down comp)
!
! N.B.: be aware that the storage layout with respect to m1, m2, 
! alpha1, and alpha2, is lateron changed to comply with the more 
! efficient reversed storage layout used in VASP.
!
! Presently also the L operator is passed on, for use in the
! orbital moment calculations in the module in orbmom.F
!
!**********************************************************************

      SUBROUTINE SETUP_LS(L,THETA,PHI,L_OP_R,LS)

      USE prec
      
      IMPLICIT NONE
      
      INTEGER L,M,M_,I,J,K
      
      REAL(q) C_UP,C_DW
      REAL(q) THETA,PHI
      
      COMPLEX(q) U_C2R(2*L+1,2*L+1),U_R2C(2*L+1,2*L+1),TMP(2*L+1,2*L+1)
      COMPLEX(q) L_OP_C(2*L+1,2*L+1,3),L_OP_R(2*L+1,2*L+1,3)
      COMPLEX(q) LS(2*L+1,2*L+1,4),LS_TMP(2*L+1,2*L+1,4)
      COMPLEX(q) ROTMAT(0:1,0:1)

! set up L operator (in units of h_bar) for complex spherical harmonics y_lm     
!
!   |y_lm1> L_k <y_lm2| = |y_lm1> L_OP_C(m1,m2,k) <y_lm2| , where k=x,y,z
!
      L_OP_C=(0._q,0._q)
      
      DO M=1,2*L+1
         M_=M-L-1   
         C_UP=SQRT(REAL((L-M_)*(L+M_+1),q))/2.0_q
         C_DW=SQRT(REAL((L+M_)*(L-M_+1),q))/2.0_q
         ! fill x-component
         IF ((M_+1)<= L) L_OP_C(M+1,M,1)=C_UP
         IF ((M_-1)>=-L) L_OP_C(M-1,M,1)=C_DW
         ! fill y-component
         IF ((M_+1)<= L) L_OP_C(M+1,M,2)=-CMPLX(0,C_UP,q)
         IF ((M_-1)>=-L) L_OP_C(M-1,M,2)= CMPLX(0,C_DW,q)
         ! fill z-component
         L_OP_C(M,M,3)=M_
      ENDDO
      
#ifdef debug
      WRITE(*,*)
      WRITE(*,*) 'L operator with respect to complex spherical harmonics'
      WRITE(*,*) 'L=',L
      DO I=1,3
      WRITE(*,*) 'component',I     
      WRITE(*,*) 'real part'
      DO M=1,2*L+1
         WRITE(*,'(10F8.4)') (REAL(L_OP_C(M,M_,I)),M_=1,2*L+1)
      ENDDO
      WRITE(*,*) 'imaginary part'
      DO M=1,2*L+1
         WRITE(*,'(10F8.4)') (AIMAG(L_OP_C(M,M_,I)),M_=1,2*L+1)
      ENDDO
      ENDDO
#endif
      
! set up transformation matrix real->complex spherical harmonics
!
!  <y_lm1|Y_lm2> = U_R2C(m1,m2) 
! 
! where y_lm and Y_lm are, respectively, the complex and real 
! spherical harmonics
!
      U_R2C=(0._q,0._q)
          
      DO M=1,2*L+1
         M_=M-L-1
         IF (M_>0) THEN
            U_R2C( M_+L+1,M)=(-1)**M_/SQRT(2._q)
            U_R2C(-M_+L+1,M)=1/SQRT(2._q)
         ENDIF
         IF (M_==0) THEN
            U_R2C(L+1,L+1)=1
         ENDIF
         IF (M_<0) THEN
            U_R2C( M_+L+1,M)= CMPLX(0,1/SQRT(2._q),q)
            U_R2C(-M_+L+1,M)=-CMPLX(0,(-1)**M_/SQRT(2._q),q)
         ENDIF
      ENDDO

! set up transformation matrix complex->real spherical harmonics
!
!  <Y_lm1|y_lm2> = U_C2R(m1,m2)
! 
! where y_lm and Y_lm are, respectively, the complex and real 
! spherical harmonics
!
      U_C2R=(0._q,0._q)
      
      DO M=1,2*L+1
         M_=M-L-1
         IF (M_>0) THEN
            U_C2R( M_+L+1,M)=(-1)**M_/SQRT(2._q)
            U_C2R(-M_+L+1,M)=CMPLX(0,(-1)**M_/SQRT(2._q),q)
         ENDIF
         IF (M_==0) THEN
            U_C2R(L+1,L+1)=1
         ENDIF
         IF (M_<0) THEN
            U_C2R( M_+L+1,M)=-CMPLX(0,1/SQRT(2._q),q)
            U_C2R(-M_+L+1,M)=1/SQRT(2._q)
         ENDIF
      ENDDO

#ifdef debug
      WRITE(*,*) 'U_CR: transformation matrix complex Y_lm to real Y_lm'
      WRITE(*,*) 'real part'
      DO M=1,2*L+1
         WRITE(*,'(10F8.4)') (REAL(U_C2R(M,M_)),M_=1,2*L+1)
      ENDDO
      WRITE(*,*) 'imaginary part'
      DO M=1,2*L+1
         WRITE(*,'(10F8.4)') (AIMAG(U_C2R(M,M_)),M_=1,2*L+1)
      ENDDO

      WRITE(*,*)
      WRITE(*,*) 'U_RC: transformation matrix real Y_lm to complex Y_lm'
      WRITE(*,*) 'real part'
      DO M=1,2*L+1
         WRITE(*,'(10F8.4)') (REAL(U_R2C(M,M_)),M_=1,2*L+1)
      ENDDO
      WRITE(*,*) 'imaginary part'
      DO M=1,2*L+1
         WRITE(*,'(10F8.4)') (AIMAG(U_R2C(M,M_)),M_=1,2*L+1)
      ENDDO

      TMP=(0._q,0._q)
      DO M=1,2*L+1
      DO M_=1,2*L+1
         DO I=1,2*L+1
            TMP(M,M_)=TMP(M,M_)+U_C2R(M,I)*U_R2C(I,M_)
         ENDDO
      ENDDO
      ENDDO

      WRITE(*,*)
      WRITE(*,*) 'U_C2R \dot U_R2C: should be an identity'
      WRITE(*,*) 'real part'
      DO M=1,2*L+1
         WRITE(*,'(10F8.4)') (REAL(TMP(M,M_)),M_=1,2*L+1)
      ENDDO
      WRITE(*,*) 'imaginary part'
      DO M=1,2*L+1
         WRITE(*,'(10F8.4)') (AIMAG(TMP(M,M_)),M_=1,2*L+1)
      ENDDO                  
#endif

! Calculate L operator (in units of h_bar) with respect to 
! the real spherical harmonics Y_lm
!
!    |Y_lm1> L_k <Y_lm2| = |Y_lm1> L_OP_R(m1,m2,k) <Y_lm2| , where k=x,y,z
!
! n.b. L_OP_R(m1,m2,k)= \sum_ij U_C2R(m1,i) L_OP_C(i,j,k) U_R2C(j,m2)
!
      L_OP_R=(0._q,0._q)

      DO M=1,2*L+1
      DO M_=1,2*L+1
         DO I=1,2*L+1
         DO J=1,2*L+1
            L_OP_R(M,M_,:)=L_OP_R(M,M_,:)+U_C2R(M,I)*L_OP_C(I,J,:)*U_R2C(J,M_)
         ENDDO
         ENDDO      
      ENDDO
      ENDDO
      
#ifdef debug
      WRITE(*,*)
      WRITE(*,*) 'L operator with respect to real spherical harmonics'
      WRITE(*,*) 'L=',L
      DO I=1,3
      WRITE(*,*) 'component',I     
      WRITE(*,*) 'real part'
      DO M=1,2*L+1
         WRITE(*,'(10F8.4)') (REAL(L_OP_R(M,M_,I)),M_=1,2*L+1)
      ENDDO
      WRITE(*,*) 'imaginary part'
      DO M=1,2*L+1
         WRITE(*,'(10F8.4)') (AIMAG(L_OP_R(M,M_,I)),M_=1,2*L+1)
      ENDDO
      ENDDO
#endif

! Calculate the SO (L \dot S) operator (in units of h_bar^2)
! <up| SO |up>
!     LS(:,:,1)= -L_OP_R(:,:,3)/2
! <up| SO |down>
!     LS(:,:,2)= -L_OP_R(:,:,1)/2 + (0._q,1._q)*L_OP_R(:,:,2)/2
! <down| SO |up>
!     LS(:,:,3)= -L_OP_R(:,:,1)/2 - (0._q,1._q)*L_OP_R(:,:,2)/2
! <down|SO|down>
!     LS(:,:,4)=  L_OP_R(:,:,3)/2

! Calculate the SO (L \dot S) operator (in units of h_bar^2)
! <up| SO |up>
      LS(:,:,1)=  L_OP_R(:,:,3)/2
! <up| SO |down>
      LS(:,:,2)=  L_OP_R(:,:,1)/2 + (0._q,1._q)*L_OP_R(:,:,2)/2
! <down| SO |up>
      LS(:,:,3)=  L_OP_R(:,:,1)/2 - (0._q,1._q)*L_OP_R(:,:,2)/2
! <down|SO|down>
      LS(:,:,4)= -L_OP_R(:,:,3)/2


#ifdef debug
      WRITE(*,*)
      WRITE(*,*) 'SO operator with respect to real spherical harmonics'
      WRITE(*,*) 'L=',L
      DO I=1,4
      WRITE(*,*) 'component',I     
      WRITE(*,*) 'real part'
      DO M=1,2*L+1
         WRITE(*,'(10F8.4)') (REAL(LS(M,M_,I)),M_=1,2*L+1)
      ENDDO
      WRITE(*,*) 'imaginary part'
      DO M=1,2*L+1
         WRITE(*,'(10F8.4)') (AIMAG(LS(M,M_,I)),M_=1,2*L+1)
      ENDDO
      ENDDO
#endif

! Rotate the LS operator by \theta and \phi

      ROTMAT(0,0)= COS(THETA/2)*EXP(-(0._q,1._q)*PHI/2)
      ROTMAT(0,1)=-SIN(THETA/2)*EXP(-(0._q,1._q)*PHI/2)
      ROTMAT(1,0)= SIN(THETA/2)*EXP( (0._q,1._q)*PHI/2)
      ROTMAT(1,1)= COS(THETA/2)*EXP( (0._q,1._q)*PHI/2)
! this rotation matrix is consistent with a rotation
! of a magnetic field by theta and phi according to
!
!                       cos \theta \cos phi    - sin \phi   cos \phi \sin \theta
! U(\theta, \phi) =   ( cos \theta \sin phi      cos \phi   sin \phi \sin \theta )
!                        - sin \theta               0             cos \theta
!      
! (first rotation by \theta and then by \phi)
! unfortunately this rotation matrix does not have
! the property U(\theta,\phi) = U^T(-\theta,-\phi)

      ! LS_TMP(m1,m2,J+2I+1) = \sum_K LS(m1,m2,J+2K+1)*ROTMAT(K,I)
      LS_TMP=(0._q,0._q)
      DO I=0,1
         DO J=0,1
            DO K=0,1
               LS_TMP(:,:,J+I*2+1)=LS_TMP(:,:,J+I*2+1)+LS(:,:,J+K*2+1)*ROTMAT(K,I)
            ENDDO
         ENDDO
      ENDDO

      ! LS(m1,m2,J+2I+1) = \sum_M LS_TMP(m1,m2,K+2I+1)*transpose(ROTMAT(J,K))
      LS=(0._q,0._q)
      DO I=0,1
         DO J=0,1
            DO K=0,1
               LS(:,:,J+I*2+1)=LS(:,:,J+I*2+1)+CONJG(ROTMAT(K,J))*LS_TMP(:,:,K+I*2+1)
            ENDDO
         ENDDO
      ENDDO

#ifdef debug
      WRITE(*,*) 'theta=',THETA,'phi=',PHI
      WRITE(*,*)
      WRITE(*,*) 'Rotated SO operator with respect to real spherical harmonics'
      WRITE(*,*) 'L=',L
      DO I=1,4
      WRITE(*,*) 'component',I     
      WRITE(*,*) 'real part'
      DO M=1,2*L+1
         WRITE(*,'(10F8.4)') (REAL(LS(M,M_,I)),M_=1,2*L+1)
      ENDDO
      WRITE(*,*) 'imaginary part'
      DO M=1,2*L+1
         WRITE(*,'(10F8.4)') (AIMAG(LS(M,M_,I)),M_=1,2*L+1)
      ENDDO
      ENDDO
#endif      

      END SUBROUTINE SETUP_LS

      END MODULE


!**********************************************************************
!
! write the Euler angles and the transformation matrix
!
!**********************************************************************


      SUBROUTINE WRITE_EULER( IU6, LSORBIT, SAXIS) 
        USE prec
        IMPLICIT NONE
        INTEGER IU6
        LOGICAL LSORBIT
        REAL(q) SAXIS(3)
      ! local
        REAL(q) ALPHA, BETA
        
        IF ( LSORBIT ) THEN
           CALL EULER(SAXIS, ALPHA, BETA)
           IF (IU6>=0) THEN
              WRITE(IU6,10) ALPHA, BETA,COS(BETA)*COS(ALPHA),-SIN(ALPHA),SIN(BETA)*COS(ALPHA), &
                   COS(BETA)*SIN(ALPHA),COS(ALPHA),SIN(BETA)*SIN(ALPHA), &
                   -SIN(BETA),0.0D0,COS(BETA), &
                   COS(BETA)*COS(ALPHA),COS(BETA)*SIN(ALPHA),-SIN(BETA), &
                   -SIN(ALPHA),COS(ALPHA),0.0D0, & 
                   SIN(BETA)*COS(ALPHA),SIN(BETA)*SIN(ALPHA),COS(BETA)
10            FORMAT(' Euler angles ALPHA=',F14.7,'  BETA=',F14.7// &
                   ' transformation matrix from SAXIS to cartesian coordinates',/ &
                   ' ---------------------------------------------------------',/ &
                   3(F14.7,' m_x',F14.7,' m_y',F14.7,' m_z' /),/ &
                   ' transformation matrix from cartesian coordinates to SAXIS',/ &
                   ' ---------------------------------------------------------',/ &
                   3(F14.7,' m_x',F14.7,' m_y',F14.7,' m_z' /)/)
           ENDIF
        ENDIF
        
      END SUBROUTINE WRITE_EULER



