#include "flags.h"

      SUBROUTINE QRHFSET
C
C  THIS ROUTINE SETS UP THE QRHFINF COMMON BLOCK, WHICH IS
C  VERY USEFUL FOR THE EVALUATION OF THE QRHF Z-VECTOR EQUATIONS
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INTEGER POP,VRT,POPRHF,VRTRHF,POPDOC,VRTDOC,DIRPRD
      LOGICAL RESPONSE 
C
c&line mod
#include "maxbasfn.par"
      DIMENSION ISET(255),ISET2(255),SCR(2*MAXBASFN)
C
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NF1(2),NF2(2) 
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/QRHFINF/POPRHF(8),VRTRHF(8),NOSH1(8),NOSH2(8),
     &               POPDOC(8),VRTDOC(8),NAI,N1I,N2A,
c&line mod
     &               NUMISCF,NUMASCF,ISPINP,ISPINM,IQRHF
      COMMON/SYMINF/NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      COMMON /FLAGS/ IFLAGS(100)
C
      IONE=1

      RESPONSE = .FALSE.
      RESPONSE = (IFLAGS(h_IFLAGS_density) .EQ.1)
      IF (RESPONSE) THEN
         RETURN
      ENDIF 
C
      CALL GETREC(20,'JOBARC','QRHFTOT ',IONE,NMOD)
      if (nmod.gt.255) then
         print *, '@QRHFSET: Assertion failed.'
         print *, '          iset dimension = 255'
         print *, '          nmod = ',nmod
         call errex
      end if
      CALL GETREC(20,'JOBARC','QRHFIRR ',NMOD,ISET)
c&new line
      CALL GETREC(-1,'JOBARC','QRHFSPN ',NMOD,ISET2)
c
c this is a temporary kludge until qrhfx, etc routines are
c available
c
c      if(iset(1).lt.0)iflags(32)=-iset(1)
c      if(iset(1).gt.0)iflags(33)=iset(1)
c
c I have no idea who left this uncommented. Unless this is commented
c first order properties with QRHF reference will not work. The
c first order properties with QRHF was tested numericaly just to
c make sure that they are correct, and they are correct. 
c 04/2006, Ajith Perera.
c      iflags(38)=0
c&new lines
C
C     DETERMINE SPIN
C
      ISPINP=0
      ISPINM=0
      DO 11 I=1,NMOD
         IF(ISET(I).LT.0) THEN
            IF(ISET2(I).EQ.1) THEN
               IF(ISPINP.EQ.2) THEN
                  CALL QRHFERR('QRHFSET-P')
                  CALL ERREX
               ELSE
                  ISPINP=1
               ENDIF
            ELSE
               IF(ISPINP.EQ.1) THEN
                  CALL QRHFERR('QRHFSET-P')
                  CALL ERREX
               ELSE
                  ISPINP=2
               ENDIF
            ENDIF
         ELSE
            IF(ISET2(I).EQ.2) THEN
               IF(ISPINM.EQ.1) THEN
                  CALL QRHFERR('QRHFSET-M')
                  CALL ERREX
               ELSE
                  ISPINM=2
               ENDIF
            ELSE
               IF(ISPINM.EQ.2) THEN
                  CALL QRHFERR('QRHFSET-M')
                  CALL ERREX
               ELSE
                  ISPINM=1
               ENDIF
            ENDIF
         ENDIF
11    CONTINUE
C
C     NOW DETERMINE THE IRREP WHICH OCCUPITION MATCHES THE RHF'S ONE
C
      IQRHFP=0
      IQRHFM=0
      IF(ISPINP.NE.0) IQRHFP=3-ISPINP
      IF(ISPINM.NE.0) IQRHFM=3-ISPINM
      IF(IQRHFP.NE.0)THEN
         IF(IQRHFM.NE.0.AND.IQRHFM.NE.IQRHFP) THEN
            CALL QRHFERR('QRHFSET-3')
            CALL ERREX
         ELSE
            IQRHF=IQRHFP
         ENDIF
      ELSEIF(IQRHFM.NE.0) THEN
         IQRHF=IQRHFM
      ELSE
         CALL QRHFERR('QRHFSET-4')
         CALL ERREX
      ENDIF
c&end new lines
C
c&two lines mod
      CALL ICOPY(NIRREP,POP(1,IQRHF),1,POPRHF,1)
      CALL ICOPY(NIRREP,VRT(1,IQRHF),1,VRTRHF,1)
c&new lines
      CALL IZERO(NOSH1,NIRREP)
      CALL IZERO(NOSH2,NIRREP)
c&new lines end
C 
      DO 10 I=1,NMOD
       INDEX=ISET(I)
       IF(INDEX.GT.0)THEN
c&del        POPRHF(INDEX)=POPRHF(INDEX)-1
c&new line
        NOSH2(INDEX)=NOSH2(INDEX)+1
       ELSE
c&del        VRTRHF(-INDEX)=VRTRHF(-INDEX)-1
c&new line
        NOSH1(-INDEX)=NOSH1(-INDEX)+1
       ENDIF
10    CONTINUE
C
      NAI=0
      N1I=0
      N2A=0
      NUMISCF=0
      NUMASCF=0
      DO 20 IRREP=1,NIRREP
c&lines modified
       POPDOC(IRREP)=POPRHF(IRREP)-NOSH1(IRREP)
       VRTDOC(IRREP)=VRTRHF(IRREP)-NOSH2(IRREP)
c&end modif
       NAI=NAI+POPRHF(IRREP)*VRTRHF(IRREP)
       N1I=N1I+NOSH1(IRREP)*POPDOC(IRREP)
       N2A=N2A+NOSH2(IRREP)*VRTDOC(IRREP)
       NUMISCF=NUMISCF+POPRHF(IRREP)
       NUMASCF=NUMASCF+VRTRHF(IRREP)
20    CONTINUE
      NBAS=NUMISCF+NUMASCF
C
      IF(IFLAGS(1).GE.10)THEN
       WRITE(6,1000)
1000   FORMAT(T3,'@QRHFSET-I, Contents of QRHFINF common block:')
       WRITE(6,1001)(POPRHF(I),I=1,NIRREP)
       WRITE(6,1001)(VRTRHF(I),I=1,NIRREP)
       WRITE(6,1001)(NOSH1(I),I=1,NIRREP)
       WRITE(6,1001)(NOSH2(I),I=1,NIRREP)
       WRITE(6,1001)(POPDOC(I),I=1,NIRREP)
       WRITE(6,1001)(VRTDOC(I),I=1,NIRREP)
       WRITE(6,1002)NAI
       WRITE(6,1003)N1I
       WRITE(6,1004)N2A
       WRITE(6,1005)NUMISCF
       WRITE(6,1006)NUMASCF
c&new lines
       WRITE(6,1007)ISPINP
       WRITE(6,1008)ISPINM
       WRITE(6,1009)IQRHF
c&end new lines
1001   FORMAT(8I5)
1002   FORMAT(T3,'Length of AI vector    : ',I5)
1003   FORMAT(T3,'Length of 1I vector    : ',I5)
1004   FORMAT(T3,'Length of 2A vector    : ',I5)
1005   FORMAT(T3,'Occupied RHF orbitals  : ',I5)
1006   FORMAT(T3,'Virtual  RHF orbitals  : ',I5)
c&new lines
 1007  FORMAT(T3,'Spin of QRHF_P orbitals: ',I5)
 1008  FORMAT(T3,'Spin of QRHF_M orbitals: ',I5)
 1009  FORMAT(T3,'Spin of RHF occupation : ',I5)
c&end new lines
      ENDIF
C 
C REORDER SCF EIGENVALUES BASED ON *RHF OCCUPATION*
C
      CALL GETREC(20,'JOBARC','QRHFEVAL',NBAS*IINTFP,SCR)
      IOFFIRR =1
      IOFFOCC2=NBAS+1
      IOFFVRT2=NBAS+NUMISCF+1
      DO 100 IRREP=1,NIRREP
       NOCC=POPRHF(IRREP)
       NVRT=VRTRHF(IRREP)
       IOFFOCC=IOFFIRR
       IOFFVRT=IOFFIRR+NOCC
       CALL SCOPY(NOCC,SCR(IOFFOCC),1,SCR(IOFFOCC2),1)
       CALL SCOPY(NVRT,SCR(IOFFVRT),1,SCR(IOFFVRT2),1)
       IOFFIRR=IOFFIRR+NOCC+NVRT
       IOFFOCC2=IOFFOCC2+NOCC
       IOFFVRT2=IOFFVRT2+NVRT
100   CONTINUE
      CALL PUTREC(20,'JOBARC','RHFEVAL ',NBAS*IINTFP,SCR(NBAS+1))
C
      RETURN
      END
