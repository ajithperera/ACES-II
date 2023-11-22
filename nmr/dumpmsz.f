      SUBROUTINE DUMPMSZ(MSZ,SCR)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER DIRPRD 
      DOUBLE PRECISION MSZ
C
      DIMENSION MSZ(3,3),SCR(1)
C
      COMMON/SYMINF/NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/PERT3/IXYZSYM(3)
C
      DATA IONE/1/
      DATA HALF,FORUTH/0.5D0,0.25D0/
C
      OPEN(UNIT=82,FILE='CHI',STATUS='OLD',FORM='FORMATTED')
      REWIND(82)
C
C FOR CASES WITH SYMMETRY, TRANSFORM THE CHEMICAL SHIFT TENSORS
C BACK TO THE NON-SYMMETRIC REPRESENTATION
C
      CALL SSCAL(9,HALF,MSZ,1)
C
      IF(NIRREP.NE.1) THEN
C
C REORDER ALSO ALL FIELD DIRECTIONS
C
       DO 10 IFIELD1=1,3
        INEW1=IXYZSYM(IFIELD1)
        DO 10 IFIELD2=1,3
         INEW2=IXYZSYM(IFIELD2)
         SCR(INEW1+(INEW2-1)*3)=MSZ(IFIELD1,IFIELD2)
10     CONTINUE
C
       CALL SCOPY(9,SCR,1,MSZ,1)
C
      ENDIF
C
      CALL GETREC(20,'JOBARC','NINDATOM',IONE,NUCIND)
      CALL GETREC(20,'JOBARC','MULTATOM',NUCIND,SCR)
      CALL GETREC(20,'JOBARC','NAMCOORD',IINTFP*3*NUCIND,SCR(1+NUCIND))
C
      write(6,6001)
6001  FORMAT(/,' MBPT(2) contribution to the paramagnetic part',
     &       ' of the magnetazibility tensor:',/)
      CALL MSZPRI(MSZ,SCR)
C
      RETURN
      END
