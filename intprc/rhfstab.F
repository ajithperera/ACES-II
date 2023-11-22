
C STABILITY TEST FOR RHF WAVEFUNCTIONS

#include "flags.h"

      SUBROUTINE RHFSTAB(ICORE,MAXCOR,EIGVAL,NSIZ,IUHF)
      IMPLICIT INTEGER (A-Z)

      CHARACTER*34 TYPE
      CHARACTER*9 DIFSAME
      DOUBLE PRECISION ONE,ONEM,TWO,FACT,SNRM2,X,ZILCH
      DOUBLE PRECISION EIGVAL,EMIN,TOLER2,TOLER
      INTEGER SCR_SIZE, INFO
      logical bUHF, bStop, bWarn
C
      PARAMETER (TOLER=1.D-5, TOLER2=1.D-8)
      INTEGER IONE
      PARAMETER (IONE=1)
      LOGICAL FULLDIAG
C
      DIMENSION ICORE(MAXCOR),EIGVAL(2*NSIZ)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NF(2,2)
      COMMON /FLAGS/  IFLAGS(100)
      COMMON /FLAGS2/ IFLAGS2(500)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPY(255,2),DIRPRD(8,8)
      DATA TWO /2.0/
      DATA ONE /1.0/
      DATA ZILCH/0.0/
      DATA ONEM/-1.0/
      DATA FACT/1.0/
C
C The RPA solve default is set to full_diag. For large no-symmetry choose
C the  Davidson scheme. Ajith Perera, 01/2017. 
C
      FULLDIAG = .TRUE.
      IF (IFLAGS2(h_IFLAGS2_rpa_solve) .EQ. 0) Then
         FULLDIAG = .TRUE.
      ELSE IF  (IFLAGS2(h_IFLAGS2_rpa_solve) .EQ. 1) Then
         FULLDIAG = .FALSE.
      ENDIF

      WRITE(*,*) 'Performing stability analysis of RHF wavefunction.'
      WRITE(*,*)
      WRITE(*,*) 'Orbital rotation parameters will be evaluated for '
      WRITE(*,*) 'each symmetry block.'

C LOOP OVER IRREPS OF SPIN ORBITAL A MATRIX

C GET SINGLE EXCITATION DENOMINATORS FOR DIAGONAL ELEMENTS

      ISIZEA = NT(1)
      IDENA = MAXCOR+1-ISIZEA*IINTFP
      IF (MOD(IDENA,2).NE.1) IDENA=IDENA-1
      EMIN = 1.d0
      MXCOR = IDENA-1

      CALL UPDMOI(1,NT(1),3,90,0,0)

      bStop = .false.
      bWarn = .false.
      DO IRREP = 1, NIRREP
         NUMDIS=IRPDPD(IRREP,ISYTYP(1,25))
         DISSIZ=IRPDPD(IRREP,ISYTYP(2,25))
         SIZE=NUMDIS*DISSIZ
         I000=1
         I010=I000+4*SIZE*IINTFP
         I020=I010+MAX(DISSIZ,NUMDIS)*2*IINTFP
         I030=I020+MAX(DISSIZ,NUMDIS)*IINTFP
         IF (I030.GT.MAXCOR) CALL INSMEM('RHFSTAB',I030, MAXCOR)
         SCR_SIZE = MAXCOR + 1 - I030
         CALL MKARHF(ICORE(I000),ICORE(I010),IRREP)
         CALL SINDEN(ICORE(I010),ICORE(IDENA),IRREP,1)
#ifdef _DEBUG_LVLM1
         ROWDIM=2*DISSIZ
            CALL OUTPUT(ICORE(I000), 1, ROWDIM, 1, ROWDIM, ROWDIM,
     &                  ROWDIM, 1)
#endif
         DAISIZ=DISSIZ
         ROWDIM=2*DISSIZ
         IOFFAA=1
         IOFFBB=(SIZE*2+DISSIZ)*IINTFP+IOFFAA
         CALL XAXPY(DAISIZ,ONE,ICORE(I010),1,ICORE(IOFFAA),ROWDIM+1)
         CALL XAXPY(DAISIZ,ONE,ICORE(I010),1,ICORE(IOFFBB),ROWDIM+1)
         IF (FULLDIAG) THEN
#ifdef _DEBUG_LVLM1
            CALL OUTPUT(ICORE(I000), 1, ROWDIM, 1, ROWDIM, ROWDIM,
     &                  ROWDIM, 1)
#endif
            NROOT=ROWDIM
            CALL DSYEV('V','L',ROWDIM,ICORE(I000),ROWDIM,
     &                 EIGVAL,ICORE(I030),(SCR_SIZE/IINTFP),INFO)
#ifdef _DEBUG_LVL0
            Write(6,*) "The eigenvalues of RPA matrix"
            Write(6, "(4F10.5)") (EIGVAL(I), I=1, ROWDIM)
CSSS            Write(6, "(4F10.5)") (EIGVAL(I)*27.2113961, I=1, ROWDIM)
#endif
            IF (INFO.NE.0) THEN
               WRITE(*,*)
     & '@RHFSTAB: There was a problem diagonalizing the matrix.'
               CALL ERREX
            END IF
            ILOCEVC=I000
            ILOCSCR=I010
         ELSE
            NROOT=MAX(48/NIRREP,IFLAGS(75)+5)
            NROOT=MIN(ROWDIM,NROOT,15)
            MAXIT=IFLAGS2(105)
            I015=I010+MAXIT*MAXIT*IINTFP
            I020=I015+ROWDIM*NROOT*IINTFP
            I030=I020+ROWDIM*(MAXIT+NROOT)*IINTFP
            I040=I030+MAX(ROWDIM,NROOT)*IINTFP
            I050=I040+ROWDIM*IINTFP
            I060=I050+MAXIT*MAX(ROWDIM,MAXIT)*IINTFP
            I070=I060+ROWDIM*IINTFP
            CALL DAVID1(ICORE(I000), ICORE(I010),
     &                  EIGVAL,      ICORE(I015),
     &                  ICORE(I020), ICORE(I030),
     &                  ICORE(I040), ICORE(I050),
     &                  ICORE(I060), MAXIT,
     &                  NROOT,       IROOT,
     &                  ROWDIM,      DISSIZ,
     &                  TOLER,       IERR)
C
C Correctly Handle the cases where iterative diagonalization failed
c Ajith Perera, 04/2013.
C    
            If (IERR .GT. 0) Return
    
            ILOCEVC=I015
            ILOCSCR=I030
            WRITE(*,*) 'Davidson procedure found ',IROOT,' roots.'
            CALL SRTEVAL(IROOT,         ROWDIM,     EIGVAL,
     &                   ICORE(ILOCEVC),ICORE(I020),ICORE(I030))
            NROOT=IROOT
         END IF

C ANALYZE EIGENVALUES AND EIGENVECTORS REGARDING STABILITY TYPE
         CALL WHENFLT(NROOT,EIGVAL,1,ZILCH,ICORE(I020),NVAL)
         IF (FULLDIAG) THEN
            WRITE(*,*) 'There are ',NVAL,
     &                 ' instabilities within irrep ',IRREP
         ELSE
            WRITE(*,*) 'There are at least ',NVAL,
     &                 ' instabilities within irrep ',IRREP
         END IF

         DO I = 1, NVAL
            I001=I010+IINTFP*DISSIZ
            ILOC=ICORE(I020-1+I)
            IOFFA=ILOCEVC+(ILOC-1)*ROWDIM
            IOFFB=IOFFA+IINTFP*DISSIZ
            CALL SCOPY(2*DISSIZ,ICORE(IOFFA),1,ICORE(ILOCSCR),1)
            CALL SAXPY(DISSIZ,ONEM,ICORE(IOFFB),1,ICORE(ILOCSCR),1)
            X=SNRM2(DISSIZ,ICORE(ILOCSCR),1)
            IF (IRREP.EQ.1) THEN
               IF (X.LT.TOLER2) THEN
                  TYPE=' (RHF -> RHF of identical symmetry)'
                  bUHF = .false.
               ELSE
                  TYPE=' (RHF -> UHF of identical symmetry)'
                  bUHF = .true.
               END IF
            ELSE
               IF (X.LT.TOLER2) THEN
                  TYPE=' (RHF -> RHF with broken symmetry)'
                  bUHF = .false.
               ELSE
                  TYPE=' (RHF -> UHF with broken symmetry)'
                  bUHF = .true.
               END IF
            END IF
            WRITE(*,'(T4,A,I2,A,F15.10,A)')
     &              'Eigenvalue #',I,' = ',EIGVAL(ILOC),TYPE

            if (irrep.eq.1) then
               if (iFlags(h_IFLAGS_rot_evec).eq.0) then
c               o save the most significant symmetry-allowed instability
                  IF ((EIGVAL(ILOC).LT.EMIN).AND..NOT.bUHF) THEN
c                  o dump the RHF->RHF rotation vector
                     CALL PUTLST(ICORE(IOFFA),1,1,1,3,90)
                     EMIN=EIGVAL(ILOC)
                  END IF
               else
c               o follow the user-specified root
                  if (i.eq.iFlags(h_IFLAGS_rot_evec)) then
c                  o dump the rotation vector
                     CALL PUTLST(ICORE(IOFFA),1,1,1,3,90)
                     EMIN=EIGVAL(ILOC)
                     bStop=bUHF
                  end if
               end if
            else
               IF (EIGVAL(ILOC).LT.EMIN) bWarn = .true.
            end if

c        END DO I = 1, NVAL
         END DO

         if (irrep.eq.1) then
            i = iFlags(h_IFLAGS_rot_evec)
            if (i.lt.0.or.nval.lt.i) then
               write(*,*) '@RHFSTAB: Assertion failed.'
               write(*,*) '          ROT_EVEC = ',i
               write(*,*) '          # roots  = ',nval
               call errex
            end if
         end if

c     END DO IRREP=1,NIRREP
      END DO

c   o warn the user if there is a more significant instability outside irrep 1
      if (bWarn) then
         write(*,*)
         write(*,*) '@RHFSTAB: WARNING: There is a more significant',
     &              ' instability that breaks the'
         write(*,*) '          symmetry of the wavefunction. To follow',
     &              ' that instability, the'
         write(*,*) '          computational point group symmetry must',
     &              ' be lowered such that the'
         write(*,*) '          instability shows up in irrep 1.'
      end if

c   o stop if the user selected a UHF instability
      if (bStop) then
         write(*,*)
         write(*,*) '@RHFSTAB: ERROR: Only RHF->RHF instabilities are',
     &              ' allowed with an RHF'
         write(*,*) '          reference wavefunction. Either select a',
     &              ' different root with the'
         write(*,*) '          ROT_EVEC keyword or change the ',
     &              'reference to UHF with REF=UHF.'
         call errex
      end if

c   o tell the world that a rotation exists
      IF (EMIN.NE.1.d0) THEN
         WRITE(*,*) 'Eigenvector for EVAL = ',EMIN,
     &              ' has been written to MOINTS file.'
         CALL PUTREC(20,'JOBARC','SCFKICK ',IONE,1)
         if (iFlags(h_IFLAGS_hfstab).eq.2) CALL ROTT1(ICORE,MAXCOR,IUHF)
      END IF

      RETURN
      END

