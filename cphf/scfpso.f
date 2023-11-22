      SUBROUTINE SCFPSO(IRREP,NPERTI,UAIA,UAIB,
     &                  BAIA,BAIB,JPSO,JPSO2,SCR1,
     &                  SCR2,NCOORD,IPERT)
C
C THIS ROUTINE CALCULATES THE PARAMAGNETIC SPIN-ORBIT  
C CONTRIBUTION TO THE INDIRECT SPIN-SPIN COUPLING CONSTANTS
C
CEND
C
C JG 4/93
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL SCF,NONHF,LAST
      INTEGER DIRPRD
      DOUBLE PRECISION JPSO,JPSO2,MPROTON
C
      DIMENSION UAIA(1),UAIB(1),BAIA(1),BAIB(1),
     &          JPSO(NCOORD,NCOORD),JPSO2(NCOORD,NCOORD),
     &          SCR1(1),SCR2(1)
C
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),NJUNK(18)
      COMMON/METHOD/IUHF,SCF,NONHF
      COMMON/PERT/NTPERT,NNPERT(8),IIPERT(8),IXPERT,IYPERT,IZPERT,
     &            IYZPERT,IXZPERT,IXYPERT,ITRANSX,ITRANSY,ITRANSZ,
     &            NUCIND
      COMMON/FLAGS/IFLAGS(100)
      COMMON/SYMINF/NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
C
      DATA IONE/1/
      DATA HALFM,ONE /-0.5D0,1.0D0/
      DATA TWO,TWOM,FOUR /2.D0,-2.D0,4.D0/
C
C C: VELOCITY OF LIGHT IN A.U.
C
      DATA C /137.035987D0/
C
C MPROTON: MASS OF PROTON IN A.U.
C
      DATA MPROTON /1836.152736D0/
C
C CONVERT: CONVERSION FROM A.U. TO SEC  (FOR TIME UNITS)
C
      DATA  CONVERT /2.4188843D-17/
C
      DATA PI /3.141592654D0/
C
C CONVERSION OF CALCULATED COUPLING TENSOR TO THE USUAL UNITS (HZ)
C
C   A FACTOR OF TWO ACCOUNTS FOR THE FACT THAT THE SUM RUNS ONLY OVER 
C   N > N'. THE DIAGONAL ELEMENTS ARE THEREFORE NOT CORRECT BUT
C   ALSO OF NO INTEREST (AT LEAST AT THE MOMENT)
C
      FACT=TWOM/(FOUR*(MPROTON**2)*(C**4)*CONVERT*TWO*PI)
C
C DETERMINE LENGTHS OF UAI AND BAI VECTORS
C
      NAA=IRPDPD(IRREP,9)
      NAB=IRPDPD(IRREP,10)
C 
C  LOOP OVER ALL NUCLEAR SPINS WITHIN THIS IRREP
C
      DO 1000 IPERT1=1,NPERTI
C
       IND1=IPERT+IPERT1
C
       IOFF1=1+(IPERT1-1)*NAA
       IOFF1B=1+(IPERT1-1)*NAB
C
C  LOOP OVER ALL PERTURBATIONS IN THIS IRREP
C
      DO 1000 IPERT2=1,NPERTI
C
       IND2=IPERT+IPERT2
C
       IOFF4=1+(IPERT2-1)*NAA
       IOFF4B=1+(IPERT2-1)*NAB
C
       JPSO(IND1,IND2)=
     &        -FACT*SDOT(NAA,UAIA(IOFF4),1,BAIA(IOFF1),1)
C
       IF(IUHF.EQ.0) THEN
        JPSO(IND1,IND2)=JPSO(IND1,IND2)*TWO
       ELSE
        JPSO(IND1,IND2)=JPSO(IND1,IND2)
     &        -FACT*SDOT(NAB,UAIB(IOFF4B),1,BAIB(IOFF1B),1)
C
       ENDIF
1000  CONTINUE
C
      IF(IND1.EQ.NCOORD) THEN
C
C FOR CASES WITH SYMMETRY, TRANSFORM THE COUPLING TENSOR BACK TO
C THE NON-SYMMETRIC REPRESENTATION
C
       IF(NIRREP.NE.1) THEN
C
        CALL ZERO(JPSO2,NCOORD*NCOORD)
        CALL TRAHES(JPSO,JPSO2,SCR1,SCR2,NCOORD)
c YAU : old
c       CALL ICOPY(NCOORD*NCOORD*IINTFP,JPSO2,1,JPSO,1)
c YAU : new
        CALL DCOPY(NCOORD*NCOORD,JPSO2,1,JPSO,1)
c YAU : end
C
       ENDIF
C
C PRINT PARAMAGNETIC SO CONTRIBUTION TO J
C
       CALL GETREC(20,'JOBARC','NINDATOM',IONE,NUCIND)
       CALL GETREC(20,'JOBARC','MULTATOM',NUCIND,SCR1)
       CALL GETREC(20,'JOBARC','NAMCOORD',IINTFP*3*NUCIND,
     &             SCR1(1+NUCIND))
C
C 10/13 Ajith change the print level
C
       IF (IFLAGS(2) .GE. 20) THEN
          IF(SCF) THEN
             CALL HEADER
     &  (' Paramagnetic part of the SO contribution to J (in Hz)',-1)
      ELSE
         CALL HEADER(
     & ' Paramagnetic SCF part of the SO contribution to J (in Hz)',-1)
         ENDIF
C
C JPRI ALSO CONVERTS THE COUPLINGS TO THE REAL J (GIVEN IN HZ)
C BY MULTIPLYING WITH THE APPROPRIATE NUCLEAR G-FACTORS
C

         CALL JPRI(JPSO,NCOORD,NUCIND,SCR1,SCR1(1+NUCIND),
     &             SCR2,.TRUE.) 
      ENDIF
C
C SAVE JSO ON THE FILE `JSO' FOR LATER USE 
C
      OPEN(UNIT=82,FILE='JSO',STATUS='UNKNOWN',FORM='FORMATTED')
C
C for later use, add pso contribution to dso contribution
C
c      OPEN(UNIT=82,FILE='JSO',STATUS='OLD',FORM='FORMATTED')
c      DO 2000 ICOORD=1,NCOORD
c       READ(82,3001)(HESS2(ICOORD,J),J=1,NCOORD)
2000  CONTINUE
3001  FORMAT(3F20.10)
c      CALL SAXPY(NCOORD*NCOORD,ONE,HESS2,1,HESS,1)
      REWIND(82)
      DO 2100 ICOORD=1,NCOORD
       WRITE(82,3001) (JPSO(ICOORD,J),J=1,NCOORD) 
2100  CONTINUE
      CLOSE(UNIT=82,STATUS='KEEP') 
C
      IF(.NOT.SCF) THEN
C
C SAVE JSO ON THE FILE `JSO' FOR LATER USE 
C
       OPEN(UNIT=82,FILE='JSOSCF',STATUS='UNKNOWN',FORM='FORMATTED')
c       DO 2001 ICOORD=1,NCOORD
c        READ(82,3001)(HESS2(ICOORD,J),J=1,NCOORD)
2001   CONTINUE
c       CALL SAXPY(NCOORD*NCOORD,ONE,HESS2,1,HESS,1)
       REWIND(82)
       DO 2101 ICOORD=1,NCOORD
        WRITE(82,3001) (JPSO(ICOORD,J),J=1,NCOORD) 
2101   CONTINUE
       CLOSE(UNIT=82,STATUS='KEEP') 
      ENDIF
      ENDIF
      RETURN
      END