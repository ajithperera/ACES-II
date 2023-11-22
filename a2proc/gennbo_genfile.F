
      SUBROUTINE GENNBO_GENFILE(NATOMS,NBAS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

#include "mxatms.par"
#include "maxbasfn.par"
#include "machsp.com"

      DIMENSION COORD(3,Mxatms),NUCCHG(Mxatms),OVLAP(NBAS*NBAS)
      DIMENSION D(NBAS*NBAS),C(NBAS*NBAS),CSC(NBAS*NBAS),F(NBAS*NBAS)
      CHARACTER*80 FNAME
      CHARACTER*4 ATMLBL(MAXBASFN),LABEL(20),ANGLBL(MAXBASFN)
      CHARACTER*1 ALPH(26),ALPH2(26)
      INTEGER CENTERBF(MAXBASFN)
      LOGICAL PURED,PUREF

      DATA LABEL/'S','X','Y','Z','XX','XY','XZ','YY','YZ','ZZ',
     &           'F300','F210','F201','F120','F111','F102',
     &           'F030','F021','F012','F003'/
      DATA ALPH/'A','B','C','D','E','F','G','H','I','J','K',
     &          'L','M','N','O','P','Q','R','S','T','U','V',
     &          'W','X','Y','Z'/
      DATA ALPH2/'a','b','c','d','e','f','g','h','i','j','k',
     &           'l','m','n','o','p','q','r','s','t','u','v',
     &           'w','x','y','z'/

      OPEN(UNIT=47,FILE='FILE47')

      WRITE(47,*)'$GENNBO'
      WRITE(47,'(T2,A,I3)')'NATOMS=',NATOMS
      WRITE(47,'(T2,A,I3)')'NBAS=',NBAS
      WRITE(47,*)'BODM'
      WRITE(47,*)'BOHR'
      WRITE(47,*)'$END'
      WRITE(47,*)
      WRITE(47,*)'$NBO'
c      WRITE(47,*)'3CBOND'
      WRITE(47,*)'PRINT=4'
      WRITE(47,*)'FNLMO=W49'
      WRITE(47,*)'NLMO'
      WRITE(47,*)'$END'
      WRITE(47,*)


      WRITE(47,*)'$COORD'
      WRITE(47,*)'NATURAL ORBITAL ANALYSIS USING ACES2 INFORMATION'
      CALL GETREC(-1,'JOBARC','COORD',3*NATOMS,COORD)
      CALL GETREC(-1,'JOBARC','ATOMCHRG',NATOMS,NUCCHG)

      DO I = 1, NATOMS
         WRITE(47,*)NUCCHG(I),NUCCHG(I),COORD(1,I),COORD(2,I),COORD(3,I)
      END DO
      WRITE(47,*)'$END'


c      CALL GETREC(-1,'JOBARC','CENTERBF',NBAS,CENTERBF)

      open(100,status='scratch')
      open(200,status='scratch')
      open(300,status='scratch')

      CALL GFNAME('IIII',FNAME,ILENGTH)
      OPEN(UNIT=10,FILE='IIII',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      CALL LOCATE(10,'LABBASIS')

      DO I = 1, NBAS
         READ(10)J,ATMLBL(I),ANGLBL(I),IJUNK
         write(100,*) ATMLBL(I)
         backspace(100)
         read(100,110) ATMLBL(I)
110     FORMAT (4X,A1)
         IF (ATMLBL(I).EQ.ALPH(1))THEN
            write(200,*)11
         ELSEIF (ATMLBL(I).EQ.ALPH(2))THEN
            write(200,*)12
         ELSEIF (ATMLBL(I).EQ.ALPH(3))THEN
            write(200,*)13
         ELSEIF (ATMLBL(I).EQ.ALPH(4))THEN
            write(200,*)14
         ELSEIF (ATMLBL(I).EQ.ALPH(5))THEN
            write(200,*)15
         ELSEIF (ATMLBL(I).EQ.ALPH(6))THEN
            write(200,*)16
         ELSEIF (ATMLBL(I).EQ.ALPH(7))THEN
            write(200,*)17
         ELSEIF (ATMLBL(I).EQ.ALPH(8))THEN
            write(200,*)18
         ELSEIF (ATMLBL(I).EQ.ALPH(9))THEN
            write(200,*)19
         ELSEIF (ATMLBL(I).EQ.ALPH(10))THEN
            write(200,*)20
         ELSEIF (ATMLBL(I).EQ.ALPH(11))THEN
            write(200,*)21
         ELSEIF (ATMLBL(I).EQ.ALPH(12))THEN
            write(200,*)22
         ELSEIF (ATMLBL(I).EQ.ALPH(13))THEN
            write(200,*)23
         ELSEIF (ATMLBL(I).EQ.ALPH(14))THEN
            write(200,*)24
         ELSEIF (ATMLBL(I).EQ.ALPH(15))THEN
            write(200,*)25
         ELSEIF (ATMLBL(I).EQ.ALPH(16))THEN
            write(200,*)26
         ELSEIF (ATMLBL(I).EQ.ALPH(17))THEN
            write(200,*)27
         ELSEIF (ATMLBL(I).EQ.ALPH(18))THEN
            write(200,*)28
         ELSEIF (ATMLBL(I).EQ.ALPH(19))THEN
            write(200,*)29
         ELSEIF (ATMLBL(I).EQ.ALPH(20))THEN
            write(200,*)30
         ELSEIF (ATMLBL(I).EQ.ALPH(21))THEN
            write(200,*)31
         ELSEIF (ATMLBL(I).EQ.ALPH(22))THEN
            write(200,*)32
         ELSEIF (ATMLBL(I).EQ.ALPH(23))THEN
            write(200,*)33
         ELSEIF (ATMLBL(I).EQ.ALPH(24))THEN
            write(200,*)34
         ELSEIF (ATMLBL(I).EQ.ALPH(25))THEN
            write(200,*)35
         ELSEIF (ATMLBL(I).EQ.ALPH(26))THEN
            write(200,*)36
         ELSEIF (ATMLBL(I).EQ.ALPH2(1))THEN
            write(200,*)37
         ELSEIF (ATMLBL(I).EQ.ALPH2(2))THEN
            write(200,*)38
         ELSEIF (ATMLBL(I).EQ.ALPH2(3))THEN
            write(200,*)39
         ELSEIF (ATMLBL(I).EQ.ALPH2(4))THEN
            write(200,*)40
         ELSEIF (ATMLBL(I).EQ.ALPH2(5))THEN
            write(200,*)41
         ELSEIF (ATMLBL(I).EQ.ALPH2(6))THEN
            write(200,*)42
         ELSEIF (ATMLBL(I).EQ.ALPH2(7))THEN
            write(200,*)43
         ELSEIF (ATMLBL(I).EQ.ALPH2(8))THEN
            write(200,*)44
         ELSEIF (ATMLBL(I).EQ.ALPH2(9))THEN
            write(200,*)45
         ELSEIF (ATMLBL(I).EQ.ALPH2(10))THEN
            write(200,*)46
         ELSEIF (ATMLBL(I).EQ.ALPH2(11))THEN
            write(200,*)47
         ELSEIF (ATMLBL(I).EQ.ALPH2(12))THEN
            write(200,*)48
         ELSEIF (ATMLBL(I).EQ.ALPH2(13))THEN
            write(200,*)49
         ELSEIF (ATMLBL(I).EQ.ALPH2(14))THEN
            write(200,*)50
         ELSEIF (ATMLBL(I).EQ.ALPH2(15))THEN
            write(200,*)51
         ELSEIF (ATMLBL(I).EQ.ALPH2(16))THEN
            write(200,*)52
         ELSEIF (ATMLBL(I).EQ.ALPH2(17))THEN
            write(200,*)53
         ELSEIF (ATMLBL(I).EQ.ALPH2(18))THEN
            write(200,*)54
         ELSEIF (ATMLBL(I).EQ.ALPH2(19))THEN
            write(200,*)55
         ELSEIF (ATMLBL(I).EQ.ALPH2(20))THEN
            write(200,*)56
         ELSEIF (ATMLBL(I).EQ.ALPH2(21))THEN
            write(200,*)57
         ELSEIF (ATMLBL(I).EQ.ALPH2(22))THEN
            write(200,*)58
         ELSEIF (ATMLBL(I).EQ.ALPH2(23))THEN
            write(200,*)59
         ELSEIF (ATMLBL(I).EQ.ALPH2(24))THEN
            write(200,*)60
         ELSEIF (ATMLBL(I).EQ.ALPH2(25))THEN
            write(200,*)61
         ELSEIF (ATMLBL(I).EQ.ALPH2(26))THEN
            write(200,*)62
         ELSE
            write(200,*)ATMLBL(I)
         END IF
         backspace(200)
         read(200,*) J
         IF (J==0) THEN
            write(300,*) 10
         ELSE
            write(300,*) J
         END IF
         backspace(300)
         read(300,*) ATMLBL(I)
         write(*,*) I,ATMLBL(I),ANGLBL(I)
      END DO

      close(100)
      close(200)
      close(300)

      WRITE(47,*)
      WRITE(47,*)'$BASIS'
      WRITE(47,*)'CENTER ='

      DO I = 1, NBAS
c         WRITE(47,*)CENTERBF(I)
          WRITE(47,*)ATMLBL(I)
      END DO

      WRITE(47,*)
      WRITE(47,*)'LABEL ='
      IDTYPE=0
      PURED=.TRUE.
      PUREF=.TRUE.

      DO I = 1, NBAS
         IF(ANGLBL(I).EQ.LABEL(10))THEN
            PURED=.FALSE.
         END IF
         IF(ANGLBL(I).EQ.LABEL(20))THEN
            PUREF=.FALSE.
         END IF
      END DO

! VMOL LABELING IS WRONG  !!!

      DO I = 1, NBAS
         IF(ANGLBL(I).EQ.LABEL(1))THEN
            WRITE(47,*)'1'
         ELSE IF(ANGLBL(I).EQ.LABEL(2))THEN
            WRITE(47,*)'101'
         ELSE IF(ANGLBL(I).EQ.LABEL(3))THEN
            WRITE(47,*)'102'
         ELSE IF(ANGLBL(I).EQ.LABEL(4))THEN
            WRITE(47,*)'103'
         ELSE IF(ANGLBL(I).EQ.LABEL(5))THEN
            IF (PURED) THEN
               WRITE(47,*)'255'
            ELSE
               WRITE(47,*)'201'
            END IF
         ELSE IF(ANGLBL(I).EQ.LABEL(6))THEN
            IF (PURED) THEN
               WRITE(47,*)'251'
            ELSE
               WRITE(47,*)'202'
            END IF
         ELSE IF(ANGLBL(I).EQ.LABEL(7))THEN
            IF (PURED) THEN
               WRITE(47,*)'253'
            ELSE
               WRITE(47,*)'205'
            END IF
         ELSE IF(ANGLBL(I).EQ.LABEL(8))THEN
            IF (PURED) THEN
               WRITE(47,*)'254'
            ELSE
               WRITE(47,*)'204'
            END IF
         ELSE IF(ANGLBL(I).EQ.LABEL(9))THEN
            IF (PURED) THEN
               WRITE(47,*)'252'
            ELSE
               WRITE(47,*)'203'
            END IF
         ELSE IF(ANGLBL(I).EQ.LABEL(10))THEN
               WRITE(47,*)'206'

         ELSE IF(ANGLBL(I).EQ.LABEL(11))THEN
            IF (PUREF) THEN
               WRITE(47,*) '351'
            ELSE
               WRITE(47,*) '301'
            END IF
         ELSE IF(ANGLBL(I).EQ.LABEL(12))THEN
            IF (PUREF) THEN
               WRITE(47,*) '352'
            ELSE
               WRITE(47,*) '303'
            END IF
         ELSE IF(ANGLBL(I).EQ.LABEL(13))THEN
            IF (PUREF) THEN
               WRITE(47,*) '353'
            ELSE
               WRITE(47,*) '304'
            END IF
         ELSE IF(ANGLBL(I).EQ.LABEL(14))THEN
            IF (PUREF) THEN
               WRITE(47,*) '354'
            ELSE
               WRITE(47,*) '302'
            END IF
         ELSE IF(ANGLBL(I).EQ.LABEL(15))THEN
            IF (PUREF) THEN
               WRITE(47,*) '355'
            ELSE
               WRITE(47,*) '305'
            END IF
         ELSE IF(ANGLBL(I).EQ.LABEL(16))THEN
            WRITE(47,*) '308'
         ELSE IF(ANGLBL(I).EQ.LABEL(17))THEN
            IF (PUREF) THEN
               WRITE(47,*) '356'
            ELSE
               WRITE(47,*) '307'
            END IF
         ELSE IF(ANGLBL(I).EQ.LABEL(18))THEN
            IF (PUREF) THEN
               WRITE(47,*) '357'
            ELSE
               WRITE(47,*) '309'
            END IF
         ELSE IF(ANGLBL(I).EQ.LABEL(19))THEN
            WRITE(47,*) '306'
         ELSE IF(ANGLBL(I).EQ.LABEL(20))THEN
            WRITE(47,*) '310'
         END IF
      END DO
      WRITE(47,*)'$END'
      WRITE(47,*)


c     The contract datalist is not necessary since it is used only for plotting data.



      CALL GETREC(-1,'JOBARC','AOOVRLAP',NBAS*NBAS,OVLAP)
      CALL GETREC(-1,'JOBARC','SCFDENSA',NBAS*NBAS,D)
      CALL GETREC(-1,'JOBARC','SCFEVCA0',NBAS*NBAS,C)
      CALL GETREC(-1,'JOBARC','FOCKA',NBAS*NBAS,F)
CCC BAZANTE
CCC for testing purposes
c      CALL NORMALIZE_OVERLAP_DENSITY_COEFFS (NBAS,C,D,F,OVLAP)
c      CALL PUTREC(-1,'JOBARC','AOOVRLAP',NBAS*NBAS,OVLAP)
c      CALL PUTREC(-1,'JOBARC','SCFDENSA',NBAS*NBAS,D)
c      CALL PUTREC(-1,'JOBARC','SCFEVCA0',NBAS*NBAS,C)
c      CALL PUTREC(-1,'JOBARC','FOCKA',NBAS*NBAS,F)
CCC

      WRITE(47,*)'$OVERLAP'
      WRITE(47,*)OVLAP
      WRITE(47,*)'$END'

      WRITE(47,*)
      WRITE(47,*)'$DENSITY'
      WRITE(47,*)D
      WRITE(47,*)'$END'

      WRITE(47,*)
      WRITE(47,*)'$FOCK'
      WRITE(47,*)F
      WRITE(47,*)'$END'

      WRITE(47,*)
      WRITE(47,*)'$LCAOMO'
      WRITE(47,*)C
      WRITE(47,*)'$END'



      CLOSE(47)

      RETURN
      END

