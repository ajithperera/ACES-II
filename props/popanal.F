      SUBROUTINE POPANAL(NAO,NMO,DMO,DAO,SCR,MAXCOR,IUHF,ISPIN,
     &                   LABEL,ILENGTH,ITYPE)
C
C THIS ROUTINE CARRIES OUT A POPULATION ANALYSIS ON THE INPUT
C ONE-PARTICLE DENSITY MATRIX DMO.  IT IS ASSUMED THAT DMO IS
C EXPRESSED IN THE MOLECULAR ORBITAL BASIS ON INPUT
C
C A bug fix for large basis sets. The orginal ISCR was
C declared as ISCR(100, 3) and was not enough to handle large
C basis set calcs. Now it is to what it supposed to be.
C Ajith Perera, 03/2007. 
C 
C
CEND
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include "maxbasfn.par"
C
      PARAMETER (MXATMS=150)
      CHARACTER*8 SCFVC1(2)
      CHARACTER*8 SCFVC2(2)
      CHARACTER*5 SPTYPE(2)
      CHARACTER*1 ORBTYP(7)
      CHARACTER*80 LABEL
      DIMENSION DMO(NMO,NMO),DAO(NAO,NAO),SCR(MAXCOR),
     &          ISCR(MAXBASFN,3),IJUNK(MXATMS)
C
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
C
      DATA ONE,ONEM,ZILCH/1.0D0,-1.0D0,0.0D0/
      DATA SCFVC1/'SCFEVECA','SCFEVECB'/
      DATA SCFVC2/'SCFEVCA0','SCFEVCB0'/
      DATA ORBTYP/'S','P','D','F','G','H','I'/
      DATA SPTYPE/'Alpha','Beta '/
C
      write(6,*)
      IF(ILENGTH.EQ.3) THEN
       write(6,9903) LABEL
9903   FORMAT(T3,'Mulliken population analysis of ',A3,
     &          ' density.')
      ELSE IF(ILENGTH.EQ.4) THEN
       write(6,9904) LABEL
9904   FORMAT(T3,'Mulliken population analysis of ',A4,
     &          ' density.')
      ELSE IF(ILENGTH.EQ.5) THEN
       write(6,9903) LABEL
9905   FORMAT(T3,'Mulliken population analysis of ',A5,
     &          ' density.')
      ELSE IF(ILENGTH.EQ.6) THEN
       write(6,9906) LABEL
9906   FORMAT(T3,'Mulliken population analysis of ',A6,
     &          ' density.')
      ELSE IF(ILENGTH.EQ.7) THEN
       write(6,9907) LABEL
9907   FORMAT(T3,'Mulliken population analysis of ',A7,
     &          ' density.')
      ELSE IF(ILENGTH.EQ.8) THEN
       write(6,9908) LABEL
9908   FORMAT(T3,'Mulliken population analysis of ',A8,
     &          ' density.')
      ELSE IF(ILENGTH.EQ.9) THEN
       write(6,9909) LABEL
9909   FORMAT(T3,'Mulliken population analysis of ',A9,
     &          ' density.')
      ELSE IF(ILENGTH.EQ.10) THEN
       write(6,9910) LABEL
9910   FORMAT(T3,'Mulliken population analysis of ',A10,
     &          ' density.')
      ELSE IF(ILENGTH.EQ.11) THEN
       write(6,9911) LABEL
9911   FORMAT(T3,'Mulliken population analysis of ',A11,
     &          ' density.')
      ELSE IF(ILENGTH.EQ.12) THEN
       write(6,9912) LABEL
9912   FORMAT(T3,'Mulliken population analysis of ',A12,
     &          ' density.')
      ELSE IF(ILENGTH.EQ.13) THEN
       write(6,9913) LABEL
9913   FORMAT(T3,'Mulliken population analysis of ',A13,
     &          ' density.')
      ELSE IF(ILENGTH.EQ.14) THEN
       write(6,9914) LABEL
9914   FORMAT(T3,'Mulliken population analysis of ',A14,
     &          ' density.')
      ELSE IF(ILENGTH.EQ.15) THEN
       write(6,9915) LABEL
9915   FORMAT(T3,'Mulliken population analysis of ',A15,
     &          ' density.')
      ENDIF
      IF(IUHF.EQ.1) THEN
       write(6,9999) SPTYPE(ISPIN)
9999   FORMAT(T3,A4,' density is analyzed.')
      ELSE
       write(6,9998)
9998   FORMAT(T3,' Total density is analyzed.')
      ENDIF
      write(6,*)  
      IONE=1
      CALL GETREC(20,'JOBARC','NAOBASFN',IONE,NAOTOT)
      CALL GETREC(20,'JOBARC','NATOMS  ',IONE,NATOMS)
      CALL GETREC(20,'JOBARC','MAP2ZMAT',NATOMS,IJUNK)
      CALL GETREC(20,'JOBARC','CENTERBF',NAOTOT,ISCR(1,1))
      CALL GETREC(20,'JOBARC','ANGMOMBF',NAOTOT,ISCR(1,2))
C
C TRANSFORM DENSITY MATRIX TO THE SYMMETRY ADAPTED AO BASIS
C
      I000=1
      I010=I000+NAOTOT*NAOTOT
      I020=I010+NAOTOT*NAOTOT
      CALL MO2AO3(DMO,DAO,SCR(I000),SCR(I010),NAO,NMO,
     &            ISPIN,ITYPE)

C
C EVALUATE D*S
C
      CALL GETREC(20,'JOBARC','AOOVRLAP',IINTFP*NAO*NAO,SCR)
      CALL XGEMM ('N','N',NAO,NAO,NAO,ONE,DAO,NAO,SCR,NAO,
     &            ZILCH,SCR(I010),NAO)
C
C TRANSFORM THIS TO THE FULL BASIS
C
      CALL GETREC(20,'JOBARC','AO2SOINV',IINTFP*NAO*NAOTOT,SCR(I000))
      CALL XGEMM ('N','N',NAO,NAOTOT,NAO,ONE,SCR(I010),NAO,SCR(I000),
     &            NAO,ZILCH,SCR(I020),NAO)
      CALL GETREC(20,'JOBARC','AO2SO   ',IINTFP*NAO*NAOTOT,SCR(I000))
      CALL XGEMM ('N','N',NAOTOT,NAOTOT,NAO,ONE,SCR(I000),NAOTOT,
     &            SCR(I020),NAO,ZILCH,SCR(I010),NAOTOT)

C
C COPY DIAGONAL ELEMENTS TO SCRATCH VECTOR AND THEN TRANSFORM IT
C TO THE NON-SYMMETRY ADAPTED AO BASIS
C
      CALL SCOPY(NAOTOT,SCR(I010),NAOTOT+1,SCR,1)
C
      WRITE(6,1005)
      WRITE(6,1002)
      WRITE(6,1000)
      WRITE(6,1002)
C
      DO 101 I=1,NAOTOT
       WRITE(6,1001)IJUNK(ISCR(I,1)),ORBTYP(1+ISCR(I,2)),SCR(I)        
101   CONTINUE
      WRITE(6,1002)
C
C CALCULATE POPULATIONS PER ATOM
C
      WRITE(6,1004)
      WRITE(6,1002)
      WRITE(6,1000)
      WRITE(6,1002)
      DO 102 IATOM=1,NATOMS
       CALL WHENEQ(NAOTOT,ISCR(1,1),1,IATOM,ISCR(1,3),NMATCH)
       IF(NMATCH.NE.0)THEN
        CALL GATHER(NMATCH,SCR(I010),SCR(I000),ISCR(1,3))
        X=SSUM(NMATCH,SCR(I010),1)
        WRITE(6,1003)IJUNK(IATOM),X
       ENDIF
102   CONTINUE
      WRITE(6,1002)
C
      RETURN
1000  FORMAT(T6,'Z-matrix',/,T7,'Center',T23,'Function',
     &       T41,'Population')
1001  FORMAT(T8,I3,T27,A,T40,F12.8)
1002  FORMAT(' ',55('-'))
1003  FORMAT(T8,I3,T40,F12.8)
1004  FORMAT(T3,'Population analysis by atoms (atomic charges).')
1005  FORMAT(T3,'Population analysis by orbitals.')
      END
