      SUBROUTINE SCFMG2(IRREP,NPERTB,NPERTI,UAIA,UAIB,SIJA,SIJB,
     &                  BAIA,BAIB,BIJA,BIJB,MSZ,MSZ2,
     &                  CSH,CSH2,SCR1,SCR2,ISCR,NCOORD,IPERT,
     &                  IXYZSYM,LAST)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL SCF,NONHF,LAST
      INTEGER DIRPRD
      DOUBLE PRECISION MSZ,MSZ2
C
      DIMENSION UAIA(1),UAIB(1),BAIA(1),BAIB(1),
     &          SIJA(1),SIJB(1),BIJA(1),BIJB(1),
     &          MSZ(3,3),MSZ2(3,3),
     &          CSH(3,NCOORD),CSH2(3,NCOORD),
     &          SCR1(1),SCR2(1),IXYZSYM(3)
C
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),NJUNK(18)
      COMMON/METHOD/IUHF,SCF,NONHF
      COMMON/PERT/NTPERT,NNPERT(8),IIPERT(8),IXPERT,IYPERT,IZPERT,
     &            IYZPERT,IXZPERT,IXYPERT,ITRANSX,ITRANSY,ITRANSZ,
     &            NUCIND
      COMMON/SYMINF/NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
C
      DATA IONE/1/
      DATA HALFM,ONE /-0.5D0,1.0D0/
      DATA TWO,FOUR /2.D0,4.D0/
      DATA CONVERT /26.62566914D+00/
C
C DETERMINE LENGTHS OF UAI AND BAI VECTORS
C
      NAA=IRPDPD(IRREP,9)
      NAB=IRPDPD(IRREP,10)
      NIA=IRPDPD(IRREP,21)
      NIB=IRPDPD(IRREP,22)
C 
      IND=0
      ISTART1=1
C
C  LOOP OVER ALL MAGNETIC FIELD COMPONENTS WITHIN
C  THIS IRREP
C
      DO 1000 IFIELD=1,NPERTB
C
       IOFF1=1+(IFIELD-1)*NAA
       IOFF2=1+(IFIELD-1)*NIA
       IOFF1B=1+(IFIELD-1)*NAB
       IOFF2B=1+(IFIELD-1)*NIB
C
C DETERMINE POSITION IN STANDARD ORDERING X,Y,Z
C
       DO 10 IXYZ1=ISTART1,3
        IF(IXYZSYM(IXYZ1).EQ.IRREP-1) THEN
         IND1=IXYZ1
         GO TO 15
        ENDIF
10     CONTINUE
       CALL ERREX
15     CONTINUE      
       ISTART1=IND1+1
C
C  LOOP OVER ALL PERTURBATIONS IN THIS IRREP
C
      DO 1000 IPERT1=1,NPERTI
C
       IND2=IPERT+IPERT1
C
       IOFF4=1+(IPERT1-1)*NAA
       IOFF5=1+(IPERT1-1)*NIA
       IOFF4B=1+(IPERT1-1)*NAB
       IOFF5B=1+(IPERT1-1)*NIB
C
       CSH(IND1,IND2)=
     &        -SDOT(NAA,UAIA(IOFF4),1,BAIA(IOFF1),1)
     &        -HALFM*SDOT(NIA,BIJA(IOFF5),1,SIJA(IOFF2),1)
C
       IF(IUHF.EQ.0) THEN
        CSH(IND1,IND2)=CSH(IND1,IND2)*TWO
       ELSE
        CSH(IND1,IND2)=CSH(IND1,IND2)
     &        -SDOT(NAB,UAIB(IOFF4B),1,BAIB(IOFF1B),1)
     &        -HALFM*SDOT(NIB,BIJB(IOFF5B),1,SIJB(IOFF2B),1)
C
       ENDIF
       CSH(IND1,IND2)=TWO*CSH(IND1,IND2)
1000  CONTINUE
C
C UPDATE THE CSHIFT FILE (ONLY IN THE LAST STEP !)
C
      IF(LAST) THEN
C
C CONVERT FROM A.U. TO PPM
C
       CALL SSCAL(3*NCOORD,CONVERT,CSH,1)
C
C FOR CASES WITH SYMMETRY, TRANSFORM THE CHEMICAL SHIFTS BACK TO
C THE NON-SYMMETRIC REPRESENTATION
C
       IF(NIRREP.NE.1) THEN
C
        CALL ZERO(CSH2,NCOORD*3)
        CALL TRACSH(CSH,CSH2,SCR1,NCOORD)
c YAU : old
c       CALL ICOPY(NCOORD*3*IINTFP,CSH2,1,CSH,1)
c YAU : new
        CALL DCOPY(NCOORD*3,CSH2,1,CSH,1)
c YAU : end
C
       ENDIF
C
C  PRINT PARAMAGNETIC PART OF SHIELDING TENSOR
C 
       CALL GETREC(20,'JOBARC','NINDATOM',IONE,NUCIND)
       CALL GETREC(20,'JOBARC','MULTATOM',NUCIND,SCR1)
       CALL GETREC(20,'JOBARC','NAMCOORD',IINTFP*3*NUCIND,    
     &             SCR1(1+NUCIND))
       IF(SCF) THEN
        CALL HEADER(' Paramagnetic part of shielding tensor',-1)
       ELSE
        CALL HEADER(
     & 'SCF contribution to paramagnetic part of shielding tensor',-1)
       ENDIF
       CALL CSHPRI(CSH,NCOORD,NUCIND,SCR1,SCR1(1+NUCIND),0)
C
       OPEN(UNIT=82,FILE='CSHIFT',STATUS='OLD',FORM='FORMATTED')
       DO 2000 IATOM=1,NCOORD/3
        DO 2000 IFIELD=1,3
         READ(82,3001) CSH2(IFIELD,3*IATOM-2),
     &                 CSH2(IFIELD,3*IATOM-1),
     &                 CSH2(IFIELD,3*IATOM)
2000   CONTINUE
C
3001   FORMAT(4F20.10)
C
       CALL SAXPY(NCOORD*3,ONE,CSH,1,CSH2,1)
C
       REWIND(82)
C
       DO 2100 IATOM=1,NCOORD/3
        DO 2100 IFIELD=1,3
         WRITE(82,3001) CSH2(IFIELD,3*IATOM-2),
     &                  CSH2(IFIELD,3*IATOM-1),
     &                  CSH2(IFIELD,3*IATOM)
2100   CONTINUE
C
       CLOSE(UNIT=82,STATUS='KEEP')
C
       IF(.NOT.SCF) THEN
        OPEN(UNIT=82,FILE='CSHIFTSCF',STATUS='OLD',
     &       FORM='FORMATTED')
        DO 3000 IATOM=1,NCOORD/3
         DO 3000 IFIELD=1,3
          READ(82,3001) CSH2(IFIELD,3*IATOM-2),
     &                  CSH2(IFIELD,3*IATOM-1),
     &                  CSH2(IFIELD,3*IATOM)
3000    CONTINUE
C
C
        CALL SAXPY(NCOORD*3,ONE,CSH,1,CSH2,1)
C
        REWIND(82)
C
        DO 3100 IATOM=1,NCOORD/3
         DO 3100 IFIELD=1,3
          WRITE(82,3001) CSH2(IFIELD,3*IATOM-2),
     &                   CSH2(IFIELD,3*IATOM-1),
     &                   CSH2(IFIELD,3*IATOM)
3100    CONTINUE
C
        CLOSE(UNIT=82,STATUS='KEEP')
       ENDIF
C
C PRINT RESULTS
C
      IF(SCF) THEN
       CALL HEADER(' Total shielding tensor',-1)
       CALL CSHPRI(CSH2,NCOORD,NUCIND,SCR1,SCR1(1+NUCIND),0)
      ENDIF 
C
      ENDIF 
C
      RETURN
      END
