      SUBROUTINE ADDFX(XIA,XIANEW,FAB,FIJ,FAJ,LENAB,LENIJ,LENAJ,
     &                 NA,NB)
C
C  This routine adds the fock matrix contributions into the product
C
C      XIANEW = F*XIA
C
C  which arise from the additional Fock matrix terms which need to
C  be included in the A matrix when doing an ROHF calculation.
C
C  The total A matrix formulas are:
C
C   A(AI,BJ) = <AB//IJ> + <AJ//IB> + DEL(IJ)*F(AB) - DEL(AB)*F(IJ)  AA-BLOCK
C
C   A(ai,bj) = <ab//ij> + <aj//ib> + DEL(ij)*f(ab) - DEL(ab)*f(ij)  BB-BLOCK
C
C   A(AI,bj) = 2  <Ab//Ij> + DEL(Ib)*f(Aj)                          AB-BLOCK
C
C   A(ai,BJ) = 2  <aB//iJ> + DEL(Aj)*f(Ib)                          BA-BLOCK
C                                                      (TRANPOSE OF AB-BLOCK)
C  
C
CEND
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER DIRPRD,POP,VRT
      DIMENSION XIA(1),XIANEW(1)
      DIMENSION FAB(LENAB),FIJ(LENIJ),FAJ(LENAJ)
      COMMON/INFO/NOCCO(2),NVRTO(2)
      COMMON/SYMINF/NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
CJDW KKB stuff
C     COMMON/SYM/POP(8,2),VRT(8,2),NJUNK(6)
      COMMON/SYM2/POP(8,2),VRT(8,2),NJUNK(6)
CJDW END
C
      DATA ONE,ONEM,TWO /1.0D+0,-1.0D+0,2.0D+0/
C
C  First add the contributions arising from the Fock matrix
C  contributions into the AAAA and BBBB  A matrices.
C
C  LOOP OVER SPINS
C
      IOFFX1=1
      IOFFX2=1
      IOFFAB=1
      IOFFIJ=1
      DO 1000 ISPIN=1,2
C
C  Now we add the DEL(ij)*F(ab)  and -DEL(ab)*F(ij) contributions
C
        DO 300 IRREP=1,NIRREP
C
          NOCC=POP(IRREP,ISPIN)
          NVRT=VRT(IRREP,ISPIN)
C
C  add fock matrix elements F(a,b) contracted with XIA
C
          IF(MIN(NOCC,NVRT).NE.0) THEN
            CALL XGEMM('N','N',NVRT,NOCC,NVRT,ONE,FAB(IOFFAB),
     &                 NVRT,XIA(IOFFX1),NVRT,ONE,XIANEW(IOFFX1),NVRT)
          ENDIF
          IOFFX1=IOFFX1+NVRT*NOCC
          IOFFAB=IOFFAB+NVRT*NVRT
C
C  Add fock matrix elements F(i,j) contracted with XIA
C
          IF(MIN(NOCC,NVRT).NE.0) THEN
            CALL XGEMM('N','T',NVRT,NOCC,NOCC,ONEM,XIA(IOFFX2),NVRT,
     &                 FIJ(IOFFIJ),NOCC,ONE,XIANEW(IOFFX2),NVRT)
          ENDIF
          IOFFX2=IOFFX2+NVRT*NOCC
          IOFFIJ=IOFFIJ+NOCC*NOCC
C
  300   CONTINUE
C
 1000 CONTINUE
C
      IOFFRPA=0
      IOFFRPB=0
      DO 400 IRREP=1,NIRREP
        NOCCA=POP(IRREP,1)
        NOCCB=POP(IRREP,2)
        NVRTA=VRT(IRREP,1)
        NVRTB=VRT(IRREP,2)
        NOPEN=NOCCA-NOCCB
C
        IF(NOPEN.LE.0) GOTO 410
C
        IOFFX1=NA+IOFFRPB+1
        IOFFX2=IOFFRPA+NVRTA*NOCCB+1
        IOFFAJ=IOFFRPB+NOPEN+1
C
C  The contribution we need to add is DEL(Ib)*f(Aj).
C
        IF(MIN(NVRTA,NOCCB,NVRTB).NE.0) THEN
          CALL XGEMM('N','T',NVRTA,NOPEN,NOCCB,ONE,FAJ(IOFFAJ),NVRTB,
     &               XIA(IOFFX1),NVRTB,ONE,XIANEW(IOFFX2),NVRTA)
        ENDIF
C
C  The contribution we add now is DEL(aJ)*f(iB)
C
        IF(MIN(NOCCB,NVRTA,NVRTB).NE.0) THEN
          CALL XGEMM('T','N',NOPEN,NOCCB,NVRTA,ONE,XIA(IOFFX2),NVRTA,
     &               FAJ(IOFFAJ),NVRTB,ONE,XIANEW(IOFFX1),NVRTB)
        ENDIF
C
  410   CONTINUE
C
        IOFFRPB=IOFFRPB+NOCCB*NVRTB
        IOFFRPA=IOFFRPA+NOCCA*NVRTA
C
  400 CONTINUE
C
      RETURN
      END
