C
      SUBROUTINE MAKYNITB(ZFN, YNI, T1, ICORE, NFVOZ, NFOOY, NFVOT,
     &                    MAXCOR, IUHF, IRREPX, ISPIN)  
C
C Take the product Y(N,I) = Z(F,N)*T(F,I) in the calculation
C of three body contribution to the quadratic term.
C
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ONEM,ZILCH,ZFN,YNI,T1
      DIMENSION ICORE(MAXCOR),ZFN(NFVOZ),YNI(NFOOY),T1(NFVOT),
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
      IOFFY = 1
      IZOFF(1) = 1
C
      DO 5000 IRREPN = 2, NIRREP
         IRREP = IRREPN - 1
         IRREPF = DIRPRD(IRREP, IRREPX)
C
         IZOFF(IRREPN) = IZOFF(IRREPN - 1) + POP(IRREP, ISPIN)*
     &                   VRT(IRREPF, ISPIN)
 5000 CONTINUE
C
C Carry out the multiplication
C
      CALL ZERO (YNI, NFOOY)
C
      DO 10 IRREPI = 1, NIRREP
C
         IRREPTI = IRREPI
         IRREPTF = DIRPRD(IRREPTI, IRREPX)
         IRREPZF = IRREPTF
         IRREPZN = DIRPRD(IRREPZF, IRREPX)

         NVRTZF = VRT(IRREPZF, ISPIN)
         NOCCZN = POP(IRREPZN, ISPIN)
         NOCCTI = POP(IRREPTI, ISPIN)
         NVRTTF = VRT(IRREPTF, ISPIN)
C
         IOFFZ = IZOFF(IRREPZN)
C
         IF (MIN(NVRTZF, NOCCZN, NOCCTI, NVRTTF) .GT. 0) THEN
C
            CALL XGEMM('T', 'N', NOCCZN, NOCCTI, NVRTZF, ONE,
     &                  ZFN(IOFFZ), NVRTZF, T1(IOFFT), NVRTTF, ZILCH,
     &                  YNI(IOFFY), NOCCZN)
         ENDIF
C
C Update the offsets
C
         IOFFT = IOFFT + NVRTTF*NOCCTI
         IOFFY = IOFFY + NOCCTI*NOCCZN
C
 10   CONTINUE
C
C Return to the calling program
C
      RETURN
      END
