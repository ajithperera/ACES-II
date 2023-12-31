      SUBROUTINE DINTRN8(EVEC,AODINT,BUF,RHF,NSTART,NEND,
     &                   IRREPX,ISPIN)
C
C THIS ROUTINE TRANSFORMS DERIVATIVE INTEGRALS FROM THE AO TO THE MO
C BASIS. THIS IS THE IIJK VERSION.
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
       IOFF4I=IOFF2
C
       NAO1=AOPOP(IRREP1)
       NMO1=AOPOP(IRREP1)
       NSTART1=NSTART(IRREP1)
       NEND1=NEND(IRREP1)
       IOFFC1=INDOCC(IRREP1,ISPIN)
C
       AODSZ=0
       MODSZ=0
       DO 1001 IRREP3=1,NIRREP
        IRREP4=DIRPRD(IRREPX,IRREP3)
        IF(IRREP3.NE.IRREP1.AND.IRREP4.NE.IRREP1) THEN
         AODSZ=AODSZ+AOPOP(IRREP3)*AOPOP(IRREP4)
         MODSZ=MODSZ+AOPOP(IRREP3)*AOPOP(IRREP4)
        ENDIF
1001   CONTINUE
C
C TRANSFORM THE LHS TO THE (ISPIN,ISPIN) SPIN CASE
C
C
C LOOP OVER (XX) BLOCKS ON LHS AND TRANSFORM THEM TO (PI)
C
       DO 20 IOCC=NSTART1,NEND1
        CALL XGEMM('T','N',NMO1,AODSZ,NAO1,ONE,EVEC(IOFFC1),
     &             NAO1,AODINT(IOFF1),NAO1,ZILCH,BUF,NMO1)
        CALL TRANSP(BUF,AODINT(IOFF2),AODSZ,NMO1)
        IOFF1=IOFF1+NAO1*AODSZ
        IOFF2=IOFF2+NMO1*AODSZ
20     CONTINUE
C
C NOW WE HAVE (PI,XX) INTEGRALS.  TRANSPOSE TO (XX,PI)
C
C TRANSFORM NEW LHS (XX) TO (PQ) [SPIN CASE (ISPIN,ISPIN)]
C
       IOFF5=IOFF4I
       IOFF6=IOFF4I
       DO 21 IOCC=NSTART1,NEND1
        DO 22 JMO=1,NMO1
         DO 23 IRREP3=1,NIRREP
          IRREP4=DIRPRD(IRREPX,IRREP3)
          IF(IRREP3.NE.IRREP1.AND.IRREP4.NE.IRREP1) THEN
           NAO3=AOPOP(IRREP3)
           NAO4=AOPOP(IRREP4)
           NMO3=AOPOP(IRREP3)
           NMO4=AOPOP(IRREP4)
           IOFFC3=INDOCC(IRREP3,1)
           IOFFC4=INDOCC(IRREP4,1)
           IX=NAO3*NAO4+1
c YAU : old
c          CALL ICOPY(NAO3*NAO4*IINTFP,AODINT(IOFF5),1,BUF,1)
c YAU : new
           CALL DCOPY(NAO3*NAO4,AODINT(IOFF5),1,BUF,1)
c YAU : end
           CALL GHTRAN(BUF,BUF,EVEC,IOFFC4,IOFFC3,1,1,
     &                 BUF(IX),0,1,NAO4,NAO3,NMO4,NMO3)
c YAU : old
c          CALL ICOPY(NMO4*NMO3*IINTFP,BUF,1,AODINT(IOFF6),1)
c YAU : new
           CALL DCOPY(NMO4*NMO3,BUF,1,AODINT(IOFF6),1)
c YAU : end
           IOFF5=IOFF5+NAO4*NAO3
           IOFF6=IOFF6+NMO4*NMO3
          ENDIF
23       CONTINUE
22      CONTINUE
21     CONTINUE
C
C END LOOP OVER IRREPS
C
1000  CONTINUE
C
C ALLOCATE MEMORY
C
      IOFFBUF1=1
      DO 1010 IRREP1=1,NIRREP
       DO 1010 IRREP3=1,NIRREP
        IRREP4=DIRPRD(IRREPX,IRREP3)
        IF(IRREP3.NE.IRREP1.AND.IRREP4.NE.IRREP1) THEN
         IOFFBUF1=MAX(IOFFBUF1,1+AOPOP(IRREP1)*AOPOP(IRREP3)*
     &                AOPOP(IRREP4))
        ENDIF
1010  CONTINUE
C
      IOFFBUF2=IOFFBUF1
      DO 1020 IRREP1=1,NIRREP
       DO 1020 IRREP2=1,NIRREP
        LENGTH=MAX(POP(IRREP1,1)*POP(IRREP2,1),
     &             VRT(IRREP1,1)*VRT(IRREP2,1),
     &             VRT(IRREP1,1)*POP(IRREP2,1))
        IOFFBUF2=MAX(IOFFBUF2,IOFFBUF1+2*LENGTH)
1020  CONTINUE
C
      CALL DUMP8(AODINT,BUF,BUF(IOFFBUF1),BUF(IOFFBUF2),
     &           NSTART,NEND,IRREPX)
    
      CALL TIMER(1)

      TIMEI=TIMEI+TIMENEW
C
      RETURN
      END
