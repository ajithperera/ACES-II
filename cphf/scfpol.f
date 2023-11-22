        
      SUBROUTINE SCFPOL(IRREP,NPERT,UAIA,UAIB,BAIA,BAIB,
     *                  IXYZSYM)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL SCF,NONHF
C
      DIMENSION IXYZSYM(3)
      DIMENSION UAIA(1),UAIB(1),BAIA(1),BAIB(1)
C
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),NJUNK(18)
      COMMON/METHOD/IUHF,SCF,NONHF
C
      COMMON/POLAR/POLARS(3,3), POLFLT(3,3)
      DATA TWO,FOUR /2.D0,4.D0/
C
C DETERMINE LENGTHS OF UAI AND BAI VECTORS
C
      NAA=IRPDPD(IRREP,ISYTYP(1,18))
      NBB=IRPDPD(IRREP,ISYTYP(1,17))
C 
      IND=0
      ISTART1=1
C
C  LOOP OVER ALL PERTURBATIONS IN THIS IRREP
C
      DO 1000 IPERT1=1,NPERT
C
       DO  10 IXYZ1=ISTART1,3
        IF(IXYZSYM(IXYZ1).EQ.IRREP-1) THEN
         IND1=IXYZ1
         GO TO 15
        ENDIF
10     CONTINUE
       CALL ERREX
15     CONTINUE
       ISTART1=IND1+1
       ISTART2=1
       IOFFDAA=(IPERT1-1)*NAA+1
       IOFFDBB=(IPERT1-1)*NBB+1
C
       DO 1000 IPERT2=1,IPERT1
C 
        DO 20 IXYZ2=ISTART2,3
         IF(IXYZSYM(IXYZ2).EQ.(IRREP-1)) THEN
          IND2=IXYZ2
          GO TO 25
         ENDIF
20      CONTINUE
        CALL ERREX 
25      CONTINUE
        ISTART2=IND2+1
        IOFFBAA=(IPERT2-1)*NAA+1
        IOFFBBB=(IPERT2-1)*NBB+1
C
        IND=IND+1
C
        POLARS(IND1,IND2)=SDOT(NAA,BAIA(IOFFBAA),1,UAIA(IOFFDAA),1)
C
        IF(IUHF.EQ.0) THEN
         POLARS(IND1,IND2)=-FOUR*POLARS(IND1,IND2)
        ELSE
         POLARS(IND1,IND2)=POLARS(IND1,IND2)+SDOT(NBB,BAIB(IOFFBBB),1,
     &                              UAIB(IOFFDBB),1)
C
         POLARS(IND1,IND2)=-TWO*POLARS(IND1,IND2)
C
        ENDIF
        POLARS(IND2,IND1)=POLARS(IND1,IND2)
1000  CONTINUE
C
      RETURN
      END
