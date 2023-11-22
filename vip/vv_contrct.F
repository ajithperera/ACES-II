      SUBROUTINE VV_CONTRCT(S2IN,HBAR,S2EX,NROW,NCOL,NSUM)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      INTEGER POP,VRT,SIRREP
      DIMENSION S2IN(NROW*NSUM),HBAR(NCOL*NSUM) 
      DIMENSION S2EX(NROW*NCOL)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD 
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SINFO/ NS(8), SIRREP
      COMMON/SLISTS/LS1IN, LS1OUT, LS2IN(2,2), LS2OUT(2,2)
C
      DATA ONE,ONEM  /1.0D0,-1.0D0/
      DATA ZILCH/0.0/

#ifdef _DEBUG_LVLM
      write(6,"(a)") "HBAR"
      call output(HBAR,1,nsum,1,ncol,nsum,ncol,1)
      write(6,"(a)") "S2IN"
      call output(s2in,1,nrow,1,nsum,nrow,nsum,1)
#endif
      IF (NSUM*NROW*NCOL .NE. 0) THEN
         CALL XGEMM('N','N',NROW,NCOL,NSUM,ONE,S2IN,
     +               NROW,HBAR,NSUM,ONE,S2EX,NROW)

#ifdef _DEBUG_LVLM
      write(6,"(a)") "S2EX"
      call output(s2ex,1,nrow,1,ncol,nrow,ncol,1)
#endif
      ENDIF 
      RETURN
      END 
