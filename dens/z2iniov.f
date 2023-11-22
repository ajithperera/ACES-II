      SUBROUTINE Z2inIOV(AIOV,ICORE,MAXCOR,IUHF,TAU,IT1OFF,
     &                   IT1END,NTIMES)
C
C Devleopd by Ajith Perera, 04/2005 liberaly using rouitines developed 
C by Jurgen Gauss. 
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INTEGER DIRPRD,DISSYT,POP,VRT,DISSYW
      LOGICAL TAU, CC
      LOGICAL MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC,MBPT4
      LOGICAL DENS,GRAD,QRHF,NONHF,ROHF,SEMI,CANON,TRIP1,TRIP2,GABCD
C
      DIMENSION ICORE(MAXCOR),AIOV(1),IT1OFF(2)
      COMMON/SYM2/POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF1BB,NF2AA,NF2BB
      COMMON/INFO/NOCCO(2),NVRTO(2)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/METH/MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC
      COMMON/SYMINF/NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),NTOT(18) 
C
      DATA ONE,ONEM,HALFM /1.0D+0,-1.0D+0, -0.50D0/
C
      MBPT4 = M4DQ.OR.M4SDQ.OR.M4SDTQ 
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

      DO ISPIN = 4 - 3*IUHF, 4 
C
C  LISTG : Z2 amplitude
C  LISTTL and LISTTT: Lambda or T2 amplitudes respectively.
C  FACT  : prefactor
C
         ISPINL = MIN(3, ISPIN)
         IF (MBPT3 .OR. UCC) THEN
             LISTTL = 43 + ISPINL
             LISTT2 = 43 + ISPINL 
             FACT   = ONEM
         ELSE IF (MBPT4) THEN
             LISTTL = 43 + ISPINL
             LISTT2 = 143 + ISPINL 
             FACT   = HALFM 
         ELSE IF (CC) THEN
             LISTTL = 43 + ISPINL
             LISTT2 = 143 + ISPINL
             FACT   =  HALFM
         ENDIF
         LISTW = 26 + ISPIN
C
         DO IRREP = 1, NIRREP
C
            DISSYT = IRPDPD(IRREP,ISYTYP(1,LISTTL))
            NUMSYT = IRPDPD(IRREP,ISYTYP(2,LISTTL))
            DISSYW = IRPDPD(IRREP,ISYTYP(1,LISTW)) 
            NUMSYW = IRPDPD(IRREP,ISYTYP(2,LISTW)) 
            IF (ISPIN .GE. 3) THEN
               NVRTSQ = IRPDPD(IRREP,13)
            ELSE
               NVRTSQ = IRPDPD(IRREP,18+ISPIN)
            ENDIF
            I001 = IT1END
            I002 = I001 + (NTIMES-1)*IINTFP*NUMSYT*NVRTSQ
            I003 = I002 + (NTIMES-1)*IINTFP*NUMSYW*NUMSYT 
            I004 = I003 + IINTFP*NUMSYT*NVRTSQ 
            I005 = I004 + IINTFP*NUMSYW*NUMSYT
            IEND1= I005 + IINTFP*NUMSYT*DISSYT
            I006 = I005 
            I007 = I006 + IINTFP*MAX(NUMSYT,DISSYT)
            I008 = I007 + IINTFP*MAX(NUMSYT,DISSYT)
            IEND2= I008 + IINTFP*MAX(NUMSYT,DISSYT)
            IEND = MAX(IEND1, IEND2)
C
            MEMLEFT = MXCOR - I004
            IF (DISSYW .GT. 0) THEN
                MAXDIS=MEMLEFT/(DISSYW*IINTFP)
            ELSE
                MAXDIS=MAX(1,NUMSYW)
            ENDIF
            IF (MAXDIS .LT. 1 .OR. IEND .GT. MXCOR)
     &         CALL INSMEM("@Z2INAOOV", I009, MXCOR)
            ILEFTODO = NUMSYW
C
            CALL GETLST(ICORE(I001), 1, NUMSYT, 1, IRREP, LISTTL) 
            IF (NTIMES .EQ. 2) CALL GETLST(ICORE(I003),1,NUMSYT,1,
     &                                    IRREP,LISTT2)
C
            IF (TAU)  THEN
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
     &                       VRT(1,2),IRREP,3,ONE)
C
              ENDIF
            ENDIF 
C
            CALL ZERO(ICORE(I002), NUMSYW*NUMSYT)
            IF (NTIMES .EQ. 2) CALl ZERO(ICORE(I004), NUMSYW*NUMSYT)
C
            IOFFW = 1 
            IOFFQ = 0
 10         CONTINUE
            ICANREAD = MIN(ILEFTODO, MAXDIS)
            CALL GETLST(ICORE(I005), IOFFW, ICANREAD, 1, IRREP, LISTW)
            IOFFW    = IOFFW + ICANREAD 
            ILEFTODO = ILEFTODO - ICANREAD  
C
            CALL XGEMM("T", "N", NUMSYT, ICANREAD, DISSYW, ONE,
     &                  ICORE(I001), DISSYT, ICORE(I005), DISSYW,
     &                  ONE, ICORE(I002 + IOFFQ), NUMSYT)
C
            IF (NTIMES .EQ.2) CALL XGEMM("T", "N", NUMSYT, 
     &                                   ICANREAD, DISSYW, ONE,
     &                                   ICORE(I003), DISSYT, 
     &                                   ICORE(I005), DISSYW,
     &                                   ONE, ICORE(I004 + IOFFQ),
     &                                   NUMSYT)
            IOFFQ = IOFFQ + ICANREAD*NUMSYT*IINTFP
            IF (ILEFTODO .GT. 0) GO TO 10
C
            IF (IUHF .EQ. 0) THEN
               CALL SPINAD1(IRREP, POP(1,1), DISSYT,
     &                      ICORE(I001), ICORE(I006),
     &                      ICORE(I007))
C
               IF (NTIMES .EQ. 2) CALL SPINAD1(IRREP, POP(1,1),
     &                                        DISSYT, ICORE(I003),
     &                                        ICORE(I006),
     &                                        ICORE(I007))
            ENDIF
C
            CALL TRANSP(ICORE(I001), ICORE(I005), NUMSYT, DISSYT)
            CALL DCOPY(NUMSYT*DISSYT, ICORE(I005), 1, ICORE(I001), 1)
           IF (NTIMES .EQ. 2) THEN
               CALL TRANSP(ICORE(I003), ICORE(I005), NUMSYT, DISSYT)
               CALL DCOPY(NUMSYT*DISSYT, ICORE(I005), 1, ICORE(I003), 
     &                    1)
            ENDIF
C        
            IF (ISPIN .LT. 3) THEN
               CALL SYMEXP(IRREP, VRT(1, ISPIN), NUMSYT, ICORE(I001))
               IF (NTIMES .EQ. 2) CALL SYMEXP(IRREP, VRT(1, ISPIN), 
     &                                        NUMSYT, ICORE(I003))
               IORB1 = ISPIN
               IORB2 = ISPIN
               IOFFI = 1 + (ISPIN - 1)*NTAA
            ELSE IF (ISPIN .EQ. 3) THEN
               IORB1 = 2
               IORB2 = 1
               IOFFI = 1 
               CALL SYMTR1(IRREP,POP(1,1),VRT(1,2),NUMSYT,ICORE(I002),
     &                     ICORE(I006),ICORE(I007),ICORE(I008))
               CALL SYMTR1(IRREP,VRT(1,1),VRT(1,2),NUMSYT,ICORE(I001),
     &                     ICORE(I006),ICORE(I007),ICORE(I008))
               IF (NTIMES .EQ. 2) THEN 
                   CALL SYMTR1(IRREP,POP(1,1),VRT(1,2),NUMSYT,
     &                         ICORE(I004),ICORE(I006),ICORE(I007),
     &                         ICORE(I008))
                   CALL SYMTR1(IRREP,VRT(1,1),VRT(1,2),NUMSYT,
     &                         ICORE(I003), ICORE(I006),ICORE(I007),
     &                         ICORE(I008))
               ENDIF
            ELSE IF (ISPIN .EQ. 4) THEN
               IORB1 = 1
               IORB2 = 2 
               IOFFI = 1 + IUHF*NTAA 
            ENDIF
C
            IOFFT = 0
            IOFFQ = 0
            DO IRREPA = 1, NIRREP
               IRREPE = DIRPRD(IRREPA, IRREP)      
               IRREPI = IRREPA
C                
               NUME = VRT(IRREPE, IORB1)
               NUMA = VRT(IRREPA, IORB2)
               NUMI = POP(IRREPI, IORB2)
C
               CALL XGEMM("T", "N", NUMA, NUMI, NUMSYT*NUME, FACT,
     &                    ICORE(I003+IOFFT), NUMSYT*NUME, 
     &                    ICORE(I002+IOFFQ), NUMSYT*NUME, ONE,
     &                    AIOV(IOFFI), NUMA)
               IF (NTIMES .EQ. 2)  CALL XGEMM("T", "N", NUMA,
     &                                  NUMI, NUMSYT*NUME, FACT,
     &                                  ICORE(I001+IOFFT), 
     &                                  NUMSYT*NUME,
     &                                  ICORE(I004+IOFFQ), 
     &                                  NUMSYT*NUME, ONE,
     &                                  AIOV(IOFFI), NUMA)
C
               IOFFT = IINTFP*NUMSYT*NUME*NUMA + IOFFT
               IOFFQ = IINTFP*NUMSYT*NUME*NUMI + IOFFQ
               IOFFI = NUMA*NUMI + IOFFI
            ENDDO
C
         ENDDO
      ENDDO
C
      RETURN
      END
