      SUBROUTINE DINTRN10(EVEC,AODINT,BUF,RHF,NSTART,NEND,
     &                    IRREPX,ISPIN)
C
C THIS ROUTINE TRANSFORMS DERIVATIVE INTEGRALS FROM THE AO TO THE MO
C BASIS. THIS IS THE IJIK AND IJKI VERSION.
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIRPRD,AOPOP,AODSZ,AODIS,POP,VRT
      LOGICAL RHF
C
      DIMENSION AODINT(1),EVEC(1),BUF(1)
      DIMENSION NSTART(8),NEND(8) 
      DIMENSION IREORD(8)
C
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYM/POP(8,2),VRT(8,2),NDUM(6)
      COMMON/SYMINF/NDUMMY,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON/AOPOPS/AOPOP(8)
      COMMON/AOOFST/INDOCC(8,2)
      COMMON/INFO/NOCCO(2),NVRTO(2)
      COMMON/INFOI/TIMEI,NWRITI
      COMMON /TIMEINFO/ TIMEIN, TIMENOW, TIMETOT, TIMENEW 
C
      DATA ONE/1.0D0/
      DATA ZILCH/0.0D0/
C
      CALL TIMER(1)
C
      IOFF1=1
      IOFF2=1
C
      DO 1000 IRREP1=1,NIRREP
C
       DO 1010 IRREP2=1,NIRREP
C
        IOFF4I=IOFF1
C
        IRREP3=DIRPRD(IRREP2,IRREPX)
        IF(IRREP3.EQ.IRREP1.OR.IRREP2.EQ.IRREP1) GO TO 1010
C
        NAO1=AOPOP(IRREP1)
        NMO1=AOPOP(IRREP1)
        NSTART1=NSTART(IRREP1)
        NEND1=NEND(IRREP1)
        IOFFC2=INDOCC(IRREP2,1)
        NMO2=AOPOP(IRREP2)
        NAO2=AOPOP(IRREP2)
        NMO3=AOPOP(IRREP3)
        NAO3=AOPOP(IRREP3)
C
        AODSZ=2*AOPOP(IRREP3)*AOPOP(IRREP1)
        MODSZ=2*AOPOP(IRREP3)*AOPOP(IRREP1)
C
C TRANSFORM THE LHS TO THE (ISPIN,ISPIN) SPIN CASE
C
C LOOP OVER (XX) BLOCKS ON LHS AND TRANSFORM THEM TO (PI)
C
        DO 20 IOCC=NSTART1,NEND1
         CALL XGEMM('T','N',NMO2,AODSZ,NAO2,ONE,EVEC(IOFFC2),
     &              NAO2,AODINT(IOFF1),NAO2,ZILCH,BUF,NMO2)
         CALL TRANSP(BUF,AODINT(IOFF2),AODSZ,NMO2)
         IOFF1=IOFF1+NAO2*AODSZ
         IOFF2=IOFF2+NMO2*AODSZ
20      CONTINUE
C
C NOW WE HAVE (PI,XX) INTEGRALS.  TRANSPOSE TO (XX,PI)
C
C NOW TRANSFORM NEW LHS (XX) TO (PQ) [SPIN CASE (ISPIN,ISPIN)]
C
        IOFF5=IOFF4I  
        IOFF6=IOFF4I
        DO 21 IOCC=NSTART1,NEND1
         DO 22 JMO=1,NMO2
          IOFFC1=INDOCC(IRREP1,1)
          IOFFC3=INDOCC(IRREP3,1)
          IX=NAO3*NAO1+1
c YAU : old
c         CALL ICOPY(NAO3*NAO1*IINTFP,AODINT(IOFF5),1,BUF,1)
c YAU : new
          CALL DCOPY(NAO3*NAO1,AODINT(IOFF5),1,BUF,1)
c YAU : end
          CALL GHTRAN(BUF,BUF,EVEC,IOFFC3,IOFFC1,1,1,
     &                BUF(IX),0,1,NAO3,NAO1,NMO3,NMO1)
c YAU : old
c         CALL ICOPY(NMO3*NMO1*IINTFP,BUF,1,AODINT(IOFF6),1)
c YAU : new
          CALL DCOPY(NMO3*NMO1,BUF,1,AODINT(IOFF6),1)
c YAU : end
          IOFF5=IOFF5+NAO3*NAO1
          IOFF6=IOFF6+NMO3*NMO1
c YAU : old
c         CALL ICOPY(NAO3*NAO1*IINTFP,AODINT(IOFF5),1,BUF,1)
c YAU : new
          CALL DCOPY(NAO3*NAO1,AODINT(IOFF5),1,BUF,1)
c YAU : end
          CALL GHTRAN(BUF,BUF,EVEC,IOFFC1,IOFFC3,1,1,
     &                BUF(IX),0,1,NAO1,NAO3,NMO1,NMO3)
c YAU : old
c         CALL ICOPY(NMO3*NMO1*IINTFP,BUF,1,AODINT(IOFF6),1)
c YAU : new
          CALL DCOPY(NMO3*NMO1,BUF,1,AODINT(IOFF6),1)
c YAU : end
          IOFF5=IOFF5+NAO3*NAO1
          IOFF6=IOFF6+NMO3*NMO1
22       CONTINUE
21      CONTINUE
C
C END LOOP OVER IRREPS
C
1010   CONTINUE
1000  CONTINUE
C
C ALLOCATE MEMORY
C
      IOFFBUF1=1
      DO 1011 IRREP1=1,NIRREP
       DO 1011 IRREP2=1,NIRREP
       IRREP3=DIRPRD(IRREP2,IRREPX)
       IF(IRREP3.NE.IRREP1.AND.IRREP2.NE.IRREP1) THEN
        IOFFBUF1=MAX(IOFFBUF1,1+AOPOP(IRREP1)*AOPOP(IRREP2)*
     &               AOPOP(IRREP3))
       ENDIF
1011  CONTINUE
C
      IOFFBUF2=IOFFBUF1
      DO 1020 IRREP1=1,NIRREP
       DO 1021 IRREP2=1,NIRREP
        LENGTH=MAX(POP(IRREP2,1)*POP(IRREP1,1),
     &             VRT(IRREP1,1)*VRT(IRREP2,1),
     &             VRT(IRREP2,1)*POP(IRREP1,1))
        IOFFBUF2=MAX(IOFFBUF2,IOFFBUF1+2*LENGTH)
1021   CONTINUE
1020  CONTINUE
C
      CALL DUMP10(AODINT,BUF,BUF(IOFFBUF1),BUF(IOFFBUF2),
     &           NSTART,NEND,IRREPX)
C
      CALL TIMER(1)
      TIMEI=TIMEI+TIMENEW
C
      RETURN
      END
