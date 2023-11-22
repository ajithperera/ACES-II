      SUBROUTINE Z2inIVV(AIVV,ICORE,MAXCOR,IUHF,TAU,IT1OFF,
     &                   IT1END, ITIME)
C
C
C THIS IS A MODIFIED INTV1 ROUTINE WHICH ADDS A Z2*TAU
C CONTRIBUTION TO I(ab). Z2 IS OBTAINED FROM T2 AND
C <AB||CD> IN ATOMIC BASIS.
C  
C CODED   JULY/90    JG
C MODIFIED 12/03  PR; complete massacre. Finaly corrected 
C by A. Perera. 05/05. 
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INTEGER DIRPRD,DISSYT,POP,VRT
      LOGICAL TAU, CC
      LOGICAL MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC,MBPT4
      LOGICAL DENS,GRAD,QRHF,NONHF,ROHF,SEMI,CANON,TRIP1,TRIP2,GABCD
C
      DIMENSION ICORE(MAXCOR),AIVV(1),IT1OFF(2)
      COMMON/SYM2/POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF1BB,NF2AA,NF2BB
      COMMON/INFO/NOCCO(2),NVRTO(2)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/METH/MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC
      COMMON/SYMINF/NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),NTOT(18) 
C
      DATA ONE,ONEM,HALFM /1.0D+0,-1.0D+0, -0.50D0/
C
      MBPT4 = M4DQ .OR. M4SDQ .OR. M4SDTQ
      CC    = CCD .OR. CCSD .OR. QCISD

      IF (MBPT4) THEN
      Print*, " The AO ladder gradients for MBPT(4) is not available"
        CALL ERREX
      ENDIF
C
      MXCOR=MAXCOR
      IF(TAU) THEN
        CALL GETLST(ICORE(IT1OFF(1)),1,1,1,1,90)
        IF(IUHF.EQ.1) CALL GETLST(ICORE(IT1OFF(2)),1,1,1,2,90)
      ENDIF

      DO ISPIN = 3 -2*IUHF, 3
C
C  LISTG : Z2 amplitude
C  LISTTL: Lambda (143+) or T2 (43+) amplitudes depending on ITIME=1
C          or 2 respectively.
C  FACT  : prefactor
C
         IF (ITIME .EQ. 1) THEN
             IF (MBPT3 .OR. UCC) THEN
                 LISTTL = 43 + ISPIN 
                 FACT  = ONEM
             ELSE IF (MBPT4) THEN
                 LISTTL = 143 + ISPIN
                 FACT=ONEM
             ELSE IF (CC) THEN
                 LISTTL = 143 + ISPIN
                 FACT   =  HALFM 
             ENDIF
         ELSE 
             LISTTL = 43 + ISPIN
             FACT  =  HALFM
         ENDIF
         LISTG = 63 + ISPIN
C
         DO IRREP=1,NIRREP
C
            DISSYT = IRPDPD(IRREP,ISYTYP(1,LISTTL))
            NUMSYT = IRPDPD(IRREP,ISYTYP(2,LISTTL))
            IF (ISPIN.EQ.3) THEN
               NVRTSQ = IRPDPD(IRREP,13)
            ELSE
               NVRTSQ = IRPDPD(IRREP,18+ISPIN)
            ENDIF
            I001 = IT1END
            I002 = I001 + IINTFP*NUMSYT*NVRTSQ
            I003 = I002 + IINTFP*NUMSYT*NVRTSQ
            I004 = I003 + IINTFP*NUMSYT*DISSYT
            I005 = I004 + IINTFP*MAX(NUMSYT,DISSYT)
            I006 = I005 + IINTFP*MAX(NUMSYT,DISSYT)
            I007 = I006 + IINTFP*MAX(NUMSYT,DISSYT)
            IF (I007 .GT. MXCOR)CALL INSMEM('Z2INIVV',I007,MXCOR)
C
            CALL GETLST(ICORE(I001),1,NUMSYT,1,IRREP,LISTTL) 
C
            IF (TAU .AND. ITIME .EQ.2) THEN
               IF (ISPIN .LT. 3) THEN
C
                   CALL FTAU(ICORE(I001),ICORE(IT1OFF(ISPIN)),
     &                       ICORE(IT1OFF(ISPIN)),DISSYT,NUMSYT,
     &                       POP(1,ISPIN),POP(1,ISPIN),
     &                       VRT(1,ISPIN),VRT(1,ISPIN),IRREP,
     &                       ISPIN,ONE)
               ELSE
C
                   CALL FTAU(ICORE(I001),ICORE(IT1OFF(1)),
     &                       ICORE(IT1OFF(2)),DISSYT,NUMSYT,
     &                       POP(1,1),POP(1,2),VRT(1,1),
     &                       VRT(1,2),IRREP,ISPIN,ONE)
C
              ENDIF
            ENDIF 
C
            CALL GETTRN(ICORE(I002), ICORE(I004), DISSYT, NUMSYT,
     &                  1, IRREP, LISTG)
C
            IF (ISPIN .LT. 3) CALL SYMEXP(IRREP, VRT(1,ISPIN), NUMSYT,
     &                                    ICORE(I002))
C
            IF (IUHF .EQ. 0) CALL SPINAD1(IRREP, POP(1,1), DISSYT, 
     &                                    ICORE(I001), ICORE(I004),
     &                                    ICORE(I005))
            CALL TRANSP(ICORE(I001), ICORE(I004), NUMSYT, DISSYT)
            CALL DCOPY(NUMSYT*DISSYT, ICORE(I004), 1, ICORE(I001), 1)
C           
            IF (ISPIN .LT. 3) THEN
               CALL SYMEXP(IRREP, VRT(1, ISPIN), NUMSYT, ICORE(I001))
               IORB1 = ISPIN
               IORB2 = ISPIN
               IOFFI = 1 + (ISPIN - 1)*NF2AA
            ELSE
               IORB1 = 1
               IORB2 = 2
               IOFFI = 1 + IUHF*NF2AA
            ENDIF
C
            IOFFV = 0
            DO IRREPA = 1, NIRREP
               IRREPE = DIRPRD(IRREPA, IRREP)      
               IRREPB = IRREPA
C                
               NUME = VRT(IRREPE, IORB1)
               NUMA = VRT(IRREPA, IORB2)
               NUMB = VRT(IRREPB, IORB2)
C
               CALL XGEMM("T", "N", NUMA, NUMB, NUMSYT*NUME, FACT,
     &                    ICORE(I002+IOFFV), NUMSYT*NUME, 
     &                    ICORE(I001+IOFFV), NUMSYT*NUME, ONE,
     &                    AIVV(IOFFI), NUMA)
               IOFFV = IINTFP*NUMSYT*NUME*NUMA + IOFFV
               IOFFI = NUMA*NUMB + IOFFI
            ENDDO
C
           IF (IUHF .EQ. 1 .AND. ISPIN .EQ. 3) THEN
C
               CALL SYMTR1(IRREP,VRT(1,1),VRT(1,2),NUMSYT,
     &                     ICORE(I001),ICORE(I004),ICORE(I005),
     &                     ICORE(I006))
               CALL SYMTR1(IRREP,VRT(1,1),VRT(1,2),NUMSYT,
     &                     ICORE(I002),ICORE(I004),ICORE(I005),
     &                     ICORE(I006))
               IORB1 = 2
               IORB2 = 1
               IOFFI = 1
               IOFFV = 0
           
               DO IRREPA = 1, NIRREP
                  IRREPE = DIRPRD(IRREPA, IRREP)
                  IRREPB = IRREPA
C
                  NUME = VRT(IRREPE, IORB1)
                  NUMA = VRT(IRREPA, IORB2)
                  NUMB = VRT(IRREPB, IORB2)
C
                  CALL XGEMM("T", "N", NUMA, NUMB, NUMSYT*NUME, FACT,
     &                        ICORE(I002+IOFFV), NUMSYT*NUME,
     &                        ICORE(I001+IOFFV), NUMSYT*NUME, ONE,
     &                        AIVV(IOFFI), NUMA)
                  IOFFV = IINTFP*NUMSYT*NUME*NUMA + IOFFV
                  IOFFI = NUMA*NUMB + IOFFI
               ENDDO
           ENDIF
C
        ENDDO
      ENDDO
C
      RETURN
      END
