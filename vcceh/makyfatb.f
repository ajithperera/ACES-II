C
      SUBROUTINE MAKYFATB(ZFN, YFA, T1, ICORE, NFVOZ, NFVVY, NFVOT,
     &                    MAXCOR, IUHF, IRREPX, ISPIN)  
C
C Take the product Y(F,A) = Z(F,N)*T(N,A) in the calculation
C of three body contribution to the quadratic term.
C
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ONEM,ZILCH,ZFN,YFA,T1
      DIMENSION ICORE(MAXCOR),ZFN(NFVOZ),YFA(NFVVY),T1(NFVOT),
     &          IZOFF(8)
C      
      COMMON /INFO/ NOCA,NOCB,NVRTA,NVRTB
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON /SYM/ POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF1BB,
     &             NF2AA,NF2BB
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
C   
C Common blocks used in the quadratic term
C 
      COMMON /QADINTI/ INTIMI, INTIAE, INTIME
      COMMON /QADINTG/ INGMNAA, INGMNBB, INGMNAB, INGABAA,
     &                 INGABBB, INGABAB, INGMCAA, INGMCBB,
     &                 INGMCABAB, INGMCBABA, INGMCBAAB,
     &                 INGMCABBA
      COMMON /APERTT1/ IAPRT1AA             
      COMMON /BPERTT1/ IBPRT1AA
      COMMON /APERTT2/ IAPRT2AA1, IAPRT2BB1, IAPRT2AB1, IAPRT2AB2, 
     &                 IAPRT2AB3, IAPRT2AB4, IAPRT2AA2, IAPRT2BB2,
     &                 IAPRT2AB5 
C
      DATA ONE /1.0D+00/ 
      DATA ONEM /-1.0D+00/
      DATA ZILCH /0.0D+00/
C
      IOFFT = 1
      IOFFZ = 1
      IZOFF(1) = 1
C
      DO 5000 IRREPAO = 2, NIRREP
         IRREP = IRREPAO - 1
         IRREPF = IRREP
C
         IZOFF(IRREPAO) = IZOFF(IRREPAO - 1) + VRT(IRREP, ISPIN)*
     &                    VRT(IRREPF, ISPIN)
C
 5000 CONTINUE
C
C Carry out the multiplication
C
      CALL ZERO (YFA, NFVVY)
C
      DO 10 IRREPN = 1, NIRREP
C
         IRREPZN = IRREPN
         IRREPZF = DIRPRD(IRREPZN, IRREPX)
         IRREPTN = IRREPZN
         IRREPTA = DIRPRD(IRREPTN, IRREPX)
C
         NVRTZF = VRT(IRREPZF, ISPIN)
         NOCCZN = POP(IRREPZN, ISPIN)
         NOCCTN = POP(IRREPTN, ISPIN)
         NVRTTA = VRT(IRREPTA, ISPIN)
C
         IOFFY = IZOFF(IRREPTA)
C
         IF (MIN(NVRTZF, NOCCZN, NOCCTN, NVRTTA) .GT. 0) THEN
C     
            CALL XGEMM('N', 'T', NVRTZF, NVRTTA, NOCCZN, ONE,
     &                  ZFN(IOFFZ), NVRTZF, T1(IOFFT), NVRTTA, ZILCH,
     &                  YFA(IOFFY), NVRTZF)
         ENDIF
C
C Update the offsets
C
         IOFFT = IOFFT + NVRTTA*NOCCTN
         IOFFZ = IOFFZ + NVRTZF*NOCCZN
C
 10   CONTINUE
C
C Return to the calling program
C
      RETURN
      END
