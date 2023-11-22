      SUBROUTINE E4LNHF(ICORE,MAXCOR,IUHF,E4L,ET1SQ)
C
C  DRIVER FOR CALCULATION OF THE "LINEAR" PART OF THE MBPT(4)
C  ENERGY FOR NON-HF REFERENCE CASES 
C
C
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ONEM,ZILCH,Z12,Z21,ZLS,ZLD,E4L,ET1SQ
      LOGICAL LAMBDA
      DIMENSION ICORE(MAXCOR)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      DATA ONE  /1.0/   
      DATA ONEM /-1.0/
      DATA ZILCH /0.0/
C
C FIRST GET THE F*T1[2]*T2[1] CONTRIBUTION
C
      CALL T1T2F(ICORE,MAXCOR,IUHF,2,Z21)
C
C NOW WE NEED TO PRODUCE RESORTED LISTS OF T2[2] TO DO THE NEXT PIECE.
C  DO THAT NOW.
C
      CALL RNABIJ(ICORE,MAXCOR,IUHF,'R')
C
C NOW DO THE F*T1[1]*T2[2] CONTRIBUTION
C
      CALL T1T2F(ICORE,MAXCOR,IUHF,0,Z12)
C
C NOW GET THE T1[2]*D*T1[2] CONTRIBUTION
C
      CALL GETS4(ICORE,MAXCOR,ZLS,IUHF)
C
C NOW GET THE T2[2]*D*T2[2] CONTRIBUTION
C
      CALL GETD4(ICORE,MAXCOR,ZLD,IUHF)
C
C COMPUTE THE TOTAL LINEAR CONTRIBUTION TO THE FOURTH-ORDER
C  ENERGY
C
      E4L=Z12-Z21+ZLS+ZLD
C
C NOW COMPUTE THE T1*T1*W PART OF THE ENERGY
C
      CALL ROHFT1SQ(ICORE,MAXCOR,2,IUHF,ET1SQ)
C
      WRITE(6,200)E4L,ET1SQ
200   FORMAT(T3,' Linear fourth-order contribution ',F15.10,'.',/,
     &       T3,' T1[1]*T1[2]*W contribution       ',F15.10,'.')
C
C NOW PUT THE T2[1] RESORTED LISTS BACK ON DISK
C
      CALL RNABIJ(ICORE,MAXCOR,IUHF,'T')
      RETURN
      END
