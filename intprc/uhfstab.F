#include "flags.h"
C UHF->UHF STABILITY TEST FOR HF WAVEFUNCTIONS

      SUBROUTINE UHFSTAB(ICORE,MAXCOR,EIGVAL,NSIZ,IUHF)
      IMPLICIT INTEGER (A-Z)

      CHARACTER*34 TYPE
      CHARACTER*9 DIFSAME
      DOUBLE PRECISION ONE,ONEM,TWO,FACT,SNRM2,X,ZILCH
      DOUBLE PRECISION EIGVAL,EMIN,TOLER2,TOLER
      INTEGER SCR_SIZE, INFO
      logical bWarn

      PARAMETER (TOLER=1.D-5, TOLER2=1.0D-8 )
      INTEGER IONE
      PARAMETER (IONE=1)
      LOGICAL FULLDIAG

      DIMENSION ICORE(MAXCOR),EIGVAL(2*NSIZ)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NF(2,2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPY(255,2),DIRPRD(8,8)
      COMMON /FLAGS/  IFLAGS(100)
      COMMON /FLAGS2/ IFLAGS2(500)
      DATA TWO /2.0/
      DATA ONE /1.0/
      DATA ZILCH/0.0/
      DATA ONEM/-1.0/
      DATA FACT/1.0/
C
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

      WRITE(6,1000)
1000  FORMAT(T3,'@UHFSTAB-I, Performing stability analysis of ',
     &          'UHF wavefunction.')
      WRITE(6,1001)
1001  FORMAT(T3,'Orbital rotation parameters will be evaluated for ',
     &          'each symmetry block.')
C
C LOOP OVER IRREPS OF SPIN ORBITAL A MATRIX
C
C GET SINGLE EXCITATION DENOMINATORS FOR DIAGONAL ELEMENTS
C
      ISIZEA=NT(1)
      ISIZEB=NT(2)
      IDENA=MAXCOR-ISIZEA*IINTFP+1
      IF(MOD(IDENA,2).NE.1)IDENA=IDENA-1
      IDENB=IDENA-ISIZEB*IINTFP
      EMIN=1.d0
      CALL UPDMOI(1,NT(1),3,90,0,0)
      CALL UPDMOI(1,NT(2),4,90,0,0)
      bWarn = .false.
      DO 10 IRREP=1,NIRREP
       SIZEAA=IRPDPD(IRREP,ISYTYP(1,19))
       SIZEBB=IRPDPD(IRREP,ISYTYP(1,20))
       ROWDIM=SIZEAA+SIZEBB
       IOFFAAAA=1
       IOFFBBBB=IOFFAAAA+IINTFP*(SIZEAA*ROWDIM+SIZEAA)
       IOFFBBAA=IOFFAAAA+IINTFP*ISIZEA
       IOFFAABB=IOFFAAAA+IINTFP*(SIZEAA*ROWDIM)
       SIZTOT=ROWDIM*ROWDIM
       I000=1
       I010=I000+SIZTOT*IINTFP
       I020=I010+2*ROWDIM*IINTFP
       I030=I020+ROWDIM*IINTFP
C 07/96 bug fix: Ajith
       IF(I030.GT.MAXCOR)CALL INSMEM('UHFSTAB', I030, MAXCOR)
       SCR_SIZE=MAXCOR+1-I030
       CALL MKAUHF(ICORE(I000),ICORE(I010),IRREP)
       IOFFA=I010
       IOFFB=I010+SIZEAA*IINTFP
C
       CALL SINDEN(ICORE(IOFFA),ICORE(IDENA),IRREP,1)
       CALL SINDEN(ICORE(IOFFB),ICORE(IDENB),IRREP,2)
#ifdef _DEBUG_LVLM1
       CALL OUTPUT(ICORE(IOFFA), 1,SIZEAA, 1, 1, 1, SIZEAA, 1)
       CALL OUTPUT(ICORE(IOFFB), 1,SIZEBB, 1, 1, 1, SIZEAA, 1)
#endif

#ifdef _DEBUG_LVLM1
            CALL OUTPUT(ICORE(I000), 1, ROWDIM, 1, ROWDIM, ROWDIM,
     &                  ROWDIM, 1)
#endif
C
       CALL SAXPY(SIZEAA,ONE,ICORE(IOFFA),1,ICORE(IOFFAAAA),ROWDIM+1)
       CALL SAXPY(SIZEBB,ONE,ICORE(IOFFB),1,ICORE(IOFFBBBB),ROWDIM+1)
C
       IF(FULLDIAG)THEN
        NROOT=ROWDIM
#ifdef _DEBUG_LVLM
            CALL OUTPUT(ICORE(I000), 1, ROWDIM, 1, ROWDIM, ROWDIM,
     &                  ROWDIM, 1)
#endif
        CALL DSYEV('V','L',ROWDIM,ICORE(I000),ROWDIM,
     &             EIGVAL,ICORE(I030),(SCR_SIZE/IINTFP),INFO)
#ifdef _DEBUG_LVL0
            Write(6,*) "The eigenvalues of RPA matrix"
CSSS            Write(6, "(4F10.5)") (EIGVAL(I)*27.2113961, I=1, ROWDIM)
            Write(6, "(4F10.5)") (EIGVAL(I), I=1, ROWDIM)
#endif

        IF (INFO.NE.0) THEN
           WRITE(*,*)
     & '@UHFSTAB: There was a problem diagonalizing the matrix.'
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
        I040=I030+MAX(NROOT,ROWDIM)*IINTFP
        I050=I040+ROWDIM*IINTFP
        I060=I050+MAXIT*MAX(ROWDIM,MAXIT)*IINTFP
        I070=I060+ROWDIM*IINTFP
        CALL DAVID1(ICORE(I000),ICORE(I010),EIGVAL,ICORE(I015),
     &              ICORE(I020),ICORE(I030),ICORE(I040),ICORE(I050),
     &              ICORE(I060),MAXIT,NROOT,IROOT,ROWDIM,ROWDIM,
     &              TOLER,IERR)
C
C Correctly Handle the cases where iterative diagonalization failed
c Ajith Perera, 04/2013.
C
        If (IERR .GT. 0) Return

        ILOCEVC=I015
        ILOCSCR=I030
        WRITE(6,1010)IROOT
1010    FORMAT(T3,'Davidson procedure found ',I5,' roots.')
        CALL SRTEVAL(IROOT,ROWDIM,EIGVAL,ICORE(ILOCEVC),
     &               ICORE(I020),ICORE(I030))
        NROOT=IROOT
       ENDIF
C
C ANALYZE EIGENVALUES AND EIGENVECTORS REGARDING STABILITY TYPE
C
       CALL WHENFLT(NROOT,EIGVAL,1,ZILCH,ICORE(I020),NVAL)
C
       WRITE(6,100)NVAL,IRREP
100    FORMAT(T3,'There are ',I3,' instabilities within irrep ',
     &        i3,'.')
       DO 12 I=1,NVAL
        ILOC=ICORE(I020+I-1)
        IOFFA=ILOCEVC+(ILOC-1)*ROWDIM
        IF(IRREP.EQ.1)THEN
         WRITE(6,104)I,EIGVAL(ILOC)
        ELSE
         WRITE(6,105)I,EIGVAL(ILOC)
        ENDIF
104     FORMAT(T4,'Eigenvalue #',I2,' = ',F15.10,
     &  ' (UHF -> UHF of identical symmetry)')
105     FORMAT(T4,'Eigenvalue #',I2,' = ',F15.10,
     &  ' (UHF -> UHF of broken symmetry)')
        IF(IRREP.EQ.1)THEN
         IF(EIGVAL(ILOC).LT.EMIN.AND.IFLAGS(75).EQ.0.OR.
     &      I.EQ.IFLAGS(75))THEN
C
C FOR IRREP=1, CHECK TO SEE IF THE USER HAS SPECIFIED A SPECIFIC
C VECTOR TO FOLLOW AND DUMP THIS ONE TO THE MOINTS FILE (DEFAULTS TO
C MINIMUM EIGENVALUE).
C
C DUMP THE ORBITAL ROTATION VECTOR TO MOINTS FILE
C
          IOFFB=NT(1)*IINTFP+IOFFA
          CALL PUTLST(ICORE(IOFFA),1,1,1,3,90)
          CALL PUTLST(ICORE(IOFFB),1,1,1,4,90)
          EMIN=EIGVAL(ILOC)
         ENDIF
        ELSE
         IF(EIGVAL(ILOC).LT.EMIN) bWarn=.true.
        ENDIF
C   
12     CONTINUE
       if (irrep.eq.1) then
          i = iFlags(75)
          if (i.lt.0.or.nval.lt.i) then
             write(*,*) '@UHFSTAB: Assertion failed.'
             write(*,*) '          ROT_EVEC = ',i
             write(*,*) '          # roots  = ',nval
             call errex
          end if
       end if
10    CONTINUE
c   o warn the user if there is a more significant instability outside irrep 1
      if (bWarn) then
         write(*,*)
         write(*,*) '@UHFSTAB: WARNING: There is a more significant',
     &              ' instability that breaks the'
         write(*,*) '          symmetry of the wavefunction. To follow',
     &              ' that instability, the'
         write(*,*) '          computational point group symmetry must',
     &              ' be lowered such that the'
         write(*,*) '          instability shows up in irrep 1.'
      end if
C
      IF(EMIN.NE.1.d0)THEN
       WRITE(6,2000)EMIN
       CALL PUTREC(20,'JOBARC','SCFKICK ',IONE,1) 
       IF (IFLAGS(74).EQ.2) CALL ROTT1(ICORE,MAXCOR,IUHF)
      ENDIF
2000  FORMAT(T3,'Eigenvector for EVAL= ',F10.5,' has been ',
     &          'written to MOINTS file.')
      RETURN
      END
