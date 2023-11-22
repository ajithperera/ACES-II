C
      SUBROUTINE GENQINT(ICORE, MAXCOR, IRREPX, IOFFSET, IUHF)
C
C  Driver for I and G intermediates. Individual expressions
C  are given in respective routines and they are slightly different 
C  from CC W and F intermidiates.
C
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION TCPU,TSYS,FACT
      DIMENSION ICORE(MAXCOR)
C
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT1AA,NT1BB,NF1AA,NF1BB,
     &             NF2AA,NF2BB
      COMMON /MACHSP/ IINTLN, IFLTLN, IINTFP, IALONE, IBITWD      
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
      COMMON /TIMEINFO/ TIMEIN, TIMENOW, TIMETOT, TIMENEW
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
      FACT = 1.00D+00
      INDX = 1
C
C Zero out the lists to which the contributions are being written.
C
      DO 10 ISPIN = 3, 3 - 2*IUHF, -1
C
         CALL ZEROLIST(ICORE, MAXCOR, (INGMNAA - 1) + ISPIN)
         CALL ZEROLIST(ICORE, MAXCOR, (INGABAA - 1) + ISPIN)
         CALL ZEROLIST(ICORE, MAXCOR, (INGMCAA - 1) + ISPIN)
C
 10   CONTINUE
C     
      IF (IUHF .EQ. 0) THEN
         CALL ZEROLIST(ICORE, MAXCOR, INGMCBAAB)
      ELSE
         CALL ZEROLIST(ICORE, MAXCOR, INGMCBABA)
         CALL ZEROLIST(ICORE, MAXCOR, INGMCABBA)
         CALL ZEROLIST(ICORE, MAXCOR, INGMCBAAB)
      ENDIF
C
      DO 20 ISPIN = 1, IUHF + 1
C
         NFVV = IRPDPD(IRREPX,  18 + ISPIN)
         NFOO = IRPDPD(IRREPX,  20 + ISPIN)
         NFVO = IRPDPD(IRREPX,   8 + ISPIN)
         NSIZE = MAX(NFVV, NFOO, NFVO)*IINTFP
C
         CALL IZERO(ICORE, NSIZE)
C
         CALL PUTLST(ICORE, 1, 1, 1, ISPIN, INTIMI)
         CALL PUTLST(ICORE, 1, 1, 1, ISPIN, INTIAE)
         CALL PUTLST(ICORE, 1, 1, 1, ISPIN, INTIME)
C
 20   CONTINUE
C
C Generate G intermediates.
C
      IF (IFLAGS(1) .GE. 10) THEN
         CALL HEADER('Calculating G and I intermediates', 0, LUOUT)
         CALL TIMER(1)
      ENDIF
C
C Hole-Hole ladder intermediate.
C
      CALL T2QGMNIJ(ICORE, MAXCOR, IRREPX, IOFFSET, IUHF, FACT)
C
C Ring type G intermediates. Form the T1 contributions to the 
C G(MB, EJ) intermediate.
C 
      CALL QT1RING(ICORE, MAXCOR, IRREPX, IOFFSET, IUHF)
C
      CALL T2QGMBEJ(ICORE, MAXCOR, 'ABAB', IUHF, IRREPX, IOFFSET)
      CALL T2QGMBEJ(ICORE, MAXCOR, 'ABBA', IUHF, IRREPX, IOFFSET)
C
      IF(IUHF .EQ. 1) THEN
         CALL T2QGMBEJ(ICORE, MAXCOR, 'AAAA', IUHF, IRREPX, IOFFSET)
         CALL T2QGMBEJ(ICORE, MAXCOR, 'BBBB', IUHF, IRREPX, IOFFSET)
         CALL T2QGMBEJ(ICORE, MAXCOR, 'BABA', IUHF, IRREPX, IOFFSET)
         CALL T2QGMBEJ(ICORE, MAXCOR, 'BAAB', IUHF, IRREPX, IOFFSET)
      ENDIF
C     
C For quadratic contributions G(MN,IJ) intermediate have to be modified
C by some terms including single excitations.
C
      CALL T1QGMNIJ(ICORE, MAXCOR, IRREPX, IOFFSET, IUHF)
C
      IF (IFLAGS(1) .GE. 10) THEN
         CALL TIMER (1)
         WRITE(LUOUT, 101) TIMENEW
 101     FORMAT(T3, '@GENQINT-I, G intermediates required ',F9.3,
     &           ' seconds.')
      ENDIF
C
C I(AE), I(MI) and I(ME) intermediates.
C
      IF(IFLAGS(1) .GE. 10) CALL TIMER(1)
C     
      CALL T1QIME(ICORE, MAXCOR, IRREPX, IOFFSET, IUHF, INDX)
      CALL T1QIAE(ICORE, MAXCOR, IRREPX, IOFFSET, IUHF, INDX)
      CALL T1QIMI(ICORE, MAXCOR, IRREPX, IOFFSET, IUHF, INDX)
C
      IF(IUHF .NE. 0) THEN
         INDX = 2
         CALL T1QIME(ICORE, MAXCOR, IRREPX, IOFFSET, IUHF, INDX)
         CALL T1QIAE(ICORE, MAXCOR, IRREPX, IOFFSET, IUHF, INDX)
         CALL T1QIMI(ICORE, MAXCOR, IRREPX, IOFFSET, IUHF, INDX)
      ENDIF
C
      CALL T2QIMI(ICORE, MAXCOR, IRREPX, IOFFSET, IUHF, FACT)
      CALL T2QIAE(ICORE, MAXCOR, IRREPX, IOFFSET, IUHF, FACT)
C
C Add - 1/2 SUM M T(M,A)*Hbar(M,E) to I(AE) and + 1/2 SUM E T(I,E)*Hbar(M,E)
C In RHF CC case these contributions are zero and not added.
C
      CALL ADDFME(ICORE, MAXCOR, IRREPX, IOFFSET, IUHF)
C
      IF (IFLAGS(1) .GE. 10) THEN
         CALL TIMER(1)
         WRITE(LUOUT, *) 
         WRITE(LUOUT, 102) TIMENEW
 102     FORMAT(T3, '@GENQINT-I, I intermediates required ',F9.3,
     &           ' seconds.')
      ENDIF
C
      RETURN
      END
