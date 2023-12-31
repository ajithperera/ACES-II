      
      SUBROUTINE FIXSYM(ISYM,NBAS,SCR,ISCR)
C
C THIS ROUTINE IS A CRUDE FIRST EFFORT TO ENSURE THAT CALCULATION
C TYPES INVOLVING CHANGES IN THE POINT GROUP SYMMETRY (PARTICULARLY
C FINITE DIFFERENCE FREQUENCY CALCULATIONS) CAN BE CARRIED OUT IN
C A BLACK-BOX FASHION.  THE ROUTINE TRIES TO PICK THE APPROPRIATE
C IRREDUCIBLE REPRESENTATION OF THE POINT GROUP WHICH CORRESPONDS
C TO THE SYMMETRY OF THE DOMINANT
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIRPRD,POP,VRT
      LOGICAL PRINT
      CHARACTER*8 LABEVL(2),LABIRR(2)
      DIMENSION SCR(NBAS),ISCR(NBAS),IOFFO(8,2),IOFFV(8,2)
C
      DATA LABEVL /'SCFEVALA','SCFEVALB'/
      DATA LABIRR /'IRREPALP','IRREPBET'/
C
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYMINF/NSTART,NIRREP,IRREPY(255,2),DIRPRD(8,8)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/FLAGS/IFLAGS(100)
      COMMON/INFO/NOCCO(2),NVRTO(2)
      COMMON /SYMLOC/ ISYMOFF(8,8,25)
C
      INDX(I,J,N)=I+(J-1)*N
C
      PRINT=IFLAGS(1).GE.0
C
C CALCULATE OFFSETS
C
      DO 1 ISPIN=1,2
       IOFFO(1,ISPIN)=0
       IOFFV(1,ISPIN)=0
       DO 2 IRREP=1,NIRREP-1
        IOFFO(IRREP+1,ISPIN)=IOFFO(IRREP,ISPIN)+POP(IRREP,ISPIN)
        IOFFV(IRREP+1,ISPIN)=IOFFV(IRREP,ISPIN)+VRT(IRREP,ISPIN)
2      CONTINUE
1     CONTINUE
C
C READ EIGENVALUES OF PRINCIPAL DEPOPULATED AND POPULATED ORBITALS
C FROM JOBARC -- THESE CORRESPOND TO EIGENVALUES AT THE REFERENCE
C GEOMETRY
C
      IONE=1
      CALL GETREC(20,'JOBARC','PRINFROM',IINTFP,EVALI)
      CALL GETREC(20,'JOBARC','PRININTO',IINTFP,EVALA)
      CALL GETREC(20,'JOBARC','PRINSPIN',IONE,ISPIN)
C
C FIND CLOSEST CORRESPONDENCE BETWEEN THESE EIGENVALUES AND THOSE
C OF THE PRESENT STRUCTURE
C
      CALL GETREC(20,'JOBARC',LABEVL(ISPIN),NBAS*IINTFP,SCR)
      DIFMINI=1.D+30
      DIFMINA=1.D+30
      DO 10 IORB=1,NBAS
       XI=ABS(SCR(IORB)-EVALI)
       XA=ABS(SCR(IORB)-EVALA)
       IF(XI.LT.DIFMINI)THEN
        ILOCI=IORB
        DIFMINI=XI
       ENDIF
       IF(XA.LT.DIFMINA)THEN
        ILOCA=IORB
        DIFMINA=XA
       ENDIF
10    CONTINUE
C
C GET SYMMETRIES OF THESE ORBITALS - AND OFFSET OF VECTOR 
C IN SYMMETRY PACKED LIST
C
      CALL GETREC(20,'JOBARC',LABIRR(ISPIN),NBAS,ISCR)
      IRREPI=ISCR(ILOCI)
      IRREPA=ISCR(ILOCA)
      ISYM=DIRPRD(IRREPA,IRREPI)
C
      INDEXA=ILOCA-NOCCO(ISPIN)-IOFFV(IRREPA,ISPIN)
      INDEXI=ILOCI-IOFFO(IRREPI,ISPIN)
      IPOSPCK=INDX(INDEXA,INDEXI,VRT(IRREPA,ISPIN))+
     &             ISYMOFF(IRREPI,ISYM,8+ISPIN)-1
C
      CALL PUTREC(20,'JOBARC','TDAGUESS',IONE,IPOSPCK)
C
      IF(PRINT)THEN
       WRITE(6,1000)ISYM
1000   FORMAT(T3,'@FIXSYM-I, TRANSITION APPEARS TO ',
     &           'BELONG TO SYMMETRY BLOCK ',I3,'.')
       WRITE(6,1001)ILOCI,ILOCA
1001   FORMAT(T3,'DOMINANT CONTRIBUTION ',I3,'->',I3)
      ENDIF
      RETURN
      END
