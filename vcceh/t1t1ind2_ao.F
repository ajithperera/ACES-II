C
      SUBROUTINE T1T1IND2_AO(ICORE, MAXCOR, IRREPX, IOFFSET, IUHF)
C
C This subroutine calculates 
C               E        F
C + P_(AB) SUM T(ALP) * T(BET) * Hbar(AB, EF) to the
C          E<F  I        J
C doubles equation when Hbar(AB,EF) is not formed on the disk.
C Spin integrated formulas are as follows. Here alpha and beta is implicit.
C
C RHF
C      
C  Z(Ij,Ab) = T(I,E)*T(j,f)*Hbar(Ab,Ef) [ABAB]
C
C UHF
C
C  Z(IJ,AB) = [T(I,E)*T(J,F) - T(J,E)*T(I,F)]*Hbar(AB,EF) [AAAA]
C
C  Z(IJ,AB) = [T(i,e)*T(j,f) - T(j,e)*T(i,f)]*Hbar(ab,ef) [BBBB]
C
C  Z(Ij,Ab) = T(I,E)*T(j,f)*Hbar(Ab,Ef) [ABAB]
C 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      INTEGER DIRPRD, POP, VRT, DISSYH, DISSYQ, DISTMP
      INTEGER DISSYQAB
C
      DIMENSION ICORE(MAXCOR), IOFFT1A(8,2), IOFFT1B(8,2)
C
      COMMON /SYM/ POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF1BB,
     &             NF2AA,NF2BB
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
      COMMON /FLAGS/ IFLAGS(100)
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
      COMMON /AOPARTT1/ IAOTT1AA1, IAOTT1BB1, IAOTT1AB1, IAOTT1AA2,
     &                  IAOTT1BB2, IAOTT1AB2, IAOTT1AA3, IAOTT1BB3,
     &                  IAOTT1AB3

      DATA ONE, HALF, HALFM, ZILCH, ONEM / 1.00D+00, 0.5D+00, -0.5D+00,
     &                               0.00D+00, -1.00D00/
C 
C Allocate  memory for T1(ALPHA) and T1(BETA) amplitudes and get
C them from the disk.
C     
      LISTA = IAPRT1AA
      LISTB = IBPRT1AA
C     
      CALL GETPERT1(ICORE, MAXCOR, INTCOR, IUHF, IRREPX, LISTA,
     &              IOFFSET, IOFFT1A)
C     
      CALL GETPERT1(ICORE, INTCOR, MXCOR, IUHF, IRREPX, LISTB,
     &              IOFFSET, IOFFT1B)
C
      IRREPT1T1 = 1

C AAAA and BBBB spin cases (UHF only)
C     
      IF (IUHF .EQ. 1) THEN
C
         DO 10 ISPIN = 1, 2
C             
            LISTQ = (IAOTT1AA1 - 1) + ISPIN
C
C Loop over irreps
C
            DO 20 IRREPQIJ = 1, NIRREP

               IRREPQAB = DIRPRD(IRREPQIJ, IRREPT1T1)
C     
               DISSYQ   = IRPDPD(IRREPQAB, ISYTYP(1, LISTQ))
               NUMSYQ   = IRPDPD(IRREPQIJ, ISYTYP(2, LISTQ))
C 
               I000 = 1
               I010 = I000 + DISSYQ*NUMSYQ*IINTFP

               IF (MIN(DISSYQ, NUMSYQ) .NE. 0) THEN
C
                  IF (I010 .LT. MXCOR) THEN
C     
                     CALL MKT1TAU(ICORE(I000),
     &                            ICORE(IOFFT1A(IRREPX, ISPIN)),
     &                            ICORE(IOFFT1B(IRREPX, ISPIN)),
     &                            DISSYQ, NUMSYQ,
     &                            POP(1,ISPIN), POP(1, ISPIN),
     &                            VRT(1, ISPIN), VRT(1, ISPIN),
     &                            IRREPQIJ, IRREPX, IRREPX, 
     &                            ISPIN, ONE)
                  ELSE
                     CALL INSMEM('T1T1IND2_AO', I010, MXCOR)
                  ENDIF
#ifdef _DEBUG_LVLM
            call checksum("t1t1ind2_ao:AA",icore(i000),numsyq*dissyq)
#endif 
               CALL PUTLST(ICORE(I000), 1, NUMSYQ, 1, IRREPQIJ, LISTQ) 
               ENDIF

 20         CONTINUE
          
 10      CONTINUE
C
      ENDIF
C
C ABAB spin case. Have contributions in both UHF and RHF cases.
C 
      LISTQ = IAOTT1AB1
      ISPIN = 3
C
      DO 30 IRREPQIJ = 1, NIRREP
           
         IRREPQAB = DIRPRD(IRREPQIJ, IRREPT1T1)
C
         DISSYQ = IRPDPD(IRREPQAB, ISYTYP(1, LISTQ))
         NUMSYQ = IRPDPD(IRREPQIJ, ISYTYP(2, LISTQ))
C
         I000 = 1
         I010 = I000 + DISSYQ*NUMSYQ*IINTFP
C         
         IF (MIN(DISSYQ, NUMSYQ) .NE. 0) THEN
C
            IF (I010 .LT. MXCOR) THEN
C       
               CALL MKT1TAU(ICORE(I000), 
     &                      ICORE(IOFFT1A(IRREPX, 1)),
     &                      ICORE(IOFFT1B(IRREPX, 2)), 
     &                      DISSYQ, NUMSYQ, POP(1, 1),
     &                      POP(1, 2), VRT(1, 1), VRT(1, 2),
     &                      IRREPQIJ, IRREPX, IRREPX, 
     &                      ISPIN, ONE) 
            ELSE
               CALL INSMEM('T1T1IND2_AO', I010, MXCOR)
            ENDIF
#ifdef _DEBUG_LVLM
            call checksum("t1t1ind2_ao:AB",icore(i000),numsyq*dissyq)
#endif 
            CALL PUTLST(ICORE(I000), 1, NUMSYQ, 1, IRREPQIJ, LISTQ) 
         ENDIF
C     
 30   CONTINUE
C
      LISTT2  = IAOTT1AA1 - 1 
      LISTQ   = 60 
      LISTAOS = IAOTT1AA2 - 1
      LISTAOT = IAOTT1AA3 - 1

      CALL QDRAOLAD(ICORE, MXCOR, IUHF, .FALSE., .FALSE., IRREPT1T1,
     &              0, LISTT2, LISTQ, LISTAOS, LISTAOT) 

      CALL QW1LAD(ICORE, MXCOR, IUHF, LISTQ+1, LISTT2+1, 1, IRREPT1T1) 

      CALL MODAIBC(ICORE,MAXCOR,IUHF,ONEM)

      CALL T1W1(ICORE, MXCOR, IUHF, IRREPT1T1, LISTT2 , 90, 26, 
     &          LISTQ) 
      
      CALL MODAIBC(ICORE,MAXCOR,IUHF,ONE)

      RETURN
      END
