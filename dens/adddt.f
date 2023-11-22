      SUBROUTINE ADDDT(DENS,SCR,NSIZE,NIRREP,NUM,LIST,ISPIN,CANON)
C
C  THIS ROUTINE ADDS THE TRIPLE CONTRIBUTION
C  TO THE OCC-OCC (VIRT-VIRT) PART OF THE
C  DENSITY.
C
CEND
C
C CODED 6/91 JG
C Standard perturbed orbitals added 7/5/97 JDW.
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL CANON
      COMMON /FLAGS/ IFLAGS(100)
      DIMENSION DENS(1),SCR(1),NUM(8)
C
      CALL GETLST(SCR,1,1,1,ISPIN,LIST)
C
ccc      write(6,*) ' in ADDDT. nsize,nirrep,num,list,ispin ',
ccc     &                       nsize,nirrep,num,list,ispin
ccc      write(6,*) ' scr,dens on entry '
ccc      do 10 i=1,nsize
ccc       write(6,*) scr(i),dens(i)
ccc   10 continue
C
CJDW 7/5/97. Attempting standard orbitals.
ccc      IF(IFLAGS(64).EQ.0)THEN
ccc       write(6,*) ' *** adddt-i, attempting standard orbitals '
ccc       CALL VADD(DENS,DENS,SCR,NSIZE,1.0D+00)
ccc       write(6,*) ' dens after adding scr '
ccc       do 20 i=1,nsize
ccc       write(6,*) dens(i)
ccc   20  continue
ccc       RETURN
ccc      ENDIF
C
      IF(CANON)THEN
C
C  LOOP OVER ALL IRREPS
C
        IOFF1=0
        IOFF2=0
        DO 100 IRREP=1,NIRREP
C
        N=NUM(IRREP)
C
        DO 150 I=1,N
C
        DENS(IOFF1+I+(I-1)*N)=DENS(IOFF1+I+(I-1)*N)+SCR(IOFF2+I)
C
  150   CONTINUE 
C
        IOFF1=IOFF1+N*N
        IOFF2=IOFF2+N
C
  100   CONTINUE
      ELSE
        CALL VADD(DENS,DENS,SCR,NSIZE,1.0D+00)
      ENDIF
C
      RETURN
      END
