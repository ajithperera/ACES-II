C
      SUBROUTINE MKT1GMBEJC2(ICORE, MAXCOR, IUHF, IRREPX, IOFFSET)
C
C This subroutine computes two T1 contributions to the 
C W(MB,EJ) intermediate.  
C
C     G(MbEj) = SUM T(j,f)*Hbar(Mb,Ef) - T(n,b)*Hbar(Mn,Ej)  [ABAB]
C     G(mBeJ) = SUM T(J,F)*Hbar(mB,eF) - T(N,B)*Hbar(mN,eJ)  [BABA]
C
C Contrast to the AAAA, BBBB, ABBA, BAAB contributions ABAB and BABA
C are calculated as positive contributions, so the spin integrated 
C expressions given above agree with the formulas given in JCP, 94,
C 4334, 1990.
C
C This routine uses out of core algorithm.
C
      IMPLICIT INTEGER (A-Z)
      LOGICAL CHANGE,INCORE,RHF
      DOUBLE PRECISION ONE,ONEM,ZILCH,ALPHA,BETA,FACTOR,SDOT
      CHARACTER*4 SSTSPN
      CHARACTER*4 SPCASE(2)
C
      DIMENSION ICORE(MAXCOR),IOFFT1(8,2),IOFFZLEJ(8,4),
     &          IOFFZLBM(8,2)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),
     &                DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYM/ POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF2AA,
     &             NF1BB,NF2BB
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /FILES/LUOUT,MOINTS
      COMMON /FLAGS/IFLAGS(100)
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
      DATA ONE /1.0/
      DATA ZILCH /0.0/
      DATA ONEM /-1.0/
      DATA SPCASE /'ABAB','BABA'/
C
      CHANGE = .TRUE.
      RHF = .FALSE.
      IF (IUHF .EQ. 0) RHF = .TRUE.
C
C First pick up T1 vector.
C
      TLIST = IAPRT1AA
C
      CALL GETPERT1(ICORE, MAXCOR, MXCOR, IUHF, IRREPX, TLIST, 
     &              IOFFSET, IOFFT1)
C
C Spin cases ABAB and BABA respectively.
C     
      DO 10 ISPIN = 1, 1 + IUHF
         LSTTAR = (INGMCABAB - 1) + ISPIN
         LSTSCR = 459 + ISPIN
         FACTOR = ONE
C
C Loop over irreps - In the first block of code this corresponds mB,
C while it is jE in the second block.
C
C   G(MbEj) = SUM T(j,f)*Hbar(Mb,Ef) - T(n,b)*Hbar(Mn,Ej)
C   G(mBeJ) = SUM T(J,F)*Hbar(mB,eF) - T(N,B)*Hbar(mN,eJ)
C
         DO 20 IRREPBM = 1 , NIRREP
C
            IRREPWBM = IRREPBM
            IRREPWEF = IRREPWBM
            IRREPGBM = IRREPWBM
            IRREPGEJ = DIRPRD(IRREPGBM, IRREPX)
            IRREPMBI = IRREPGBM
            IRREPMAJ = IRREPGEJ
C
C Computes offfset 
C
            IOFF1 = 0
            IOFF2 = 0
            IOFF3 = 0
            IOFF4 = 0
C
            DO 1001 IRREPA = 1, NIRREP
C
               IRREPE = IRREPA
               IRREPM = IRREPA
               IRREPJ = DIRPRD(IRREPE, IRREPGEJ)
               IRREPB = DIRPRD(IRREPM, IRREPGBM)
C
               IOFFZLEJ(IRREPA, 1) = IOFF1
               IOFFZLEJ(IRREPA, 2) = IOFF2
               IOFFZLBM(IRREPM, 1) = IOFF3
               IOFFZLBM(IRREPM, 2) = IOFF4
C
               IOFF1 = IOFF1 + POP(IRREPJ, 2)*VRT(IRREPE, 1)
               IOFF2 = IOFF2 + POP(IRREPJ, 1)*VRT(IRREPE, 2)
               IOFF3 = IOFF3 + VRT(IRREPB, 2)*POP(IRREPM, 1)
               IOFF4 = IOFF4 + VRT(IRREPB, 1)*POP(IRREPM, 2)
C
 1001       CONTINUE
C     
C Compute dimensions of target matrix.
C     
            LSTTMP = 40 - ISPIN
C
            DSZTAR = IRPDPD(IRREPGEJ, ISYTYP(1, LSTTAR))
            DISTAR = IRPDPD(IRREPGBM, ISYTYP(2, LSTTAR))
C
            DSZTMP = IRPDPD(IRREPMAJ, ISYTYP(1, LSTTMP))
            DISTMP = IRPDPD(IRREPMBI, ISYTYP(2, LSTTMP))
C
C First do    
C  W(MbEj) = SUM T(j,f)*Hbar(Mb,Ef) (ISPIN = 1)
C  W(mBeJ) = SUM T(J,F)*Hbar(mB,eF) (ISPIN = 2 OR RHF).
C     
C This product is initially packed as (jE,Mb) [(Je,Bm)].
C     
            LISTW1 = 28 + ISPIN
            IF (RHF) LISTW1 = 30
C
            DSZW = IRPDPD(IRREPWEF, ISYTYP(1, LISTW1))
            DISW = IRPDPD(IRREPWBM, ISYTYP(2, LISTW1))
C     
C I000 Holds the G(jE,Mb) [W(Je,mB)] TARGET.
C I010 Holds area eventually used as scratch SYMTR.
C I020 Holds the the (Ef,Mb) [(Fe,Bm)] integrals.
C     
            I000 = 1
            I010 = I000 + IINTFP*DISTMP*DSZTMP
C
            CALL IZERO(ICORE, IINTFP*DSZTMP*DISTMP)
C
            I011 = I010 + IINTFP*MAX(DISTMP, DSZTMP)
            I012 = I011 + IINTFP*MAX(DSZTMP, DISTMP)
            I020 = I012 + IINTFP*MAX(DSZTMP, DISTMP)
            I030 = I020 + IINTFP*DSZW*DISW
C
            IF (I030 .LE. MXCOR) THEN
               INCORE = .TRUE.
               CALL GETLST(ICORE(I020), 1, DISW, 2, IRREPWBM, LISTW1)
            ELSE
               INCORE = .FALSE.
               I030 = I020 + IINTFP*DSZW
            ENDIF
C
            DO 30 INUMBM = 1, DISW
               IF (INCORE) THEN
                  IOFFW1R = (INUMBM - 1)*DSZW*IINTFP + I020
               ELSE
                  CALL GETLST(ICORE(I020), INUMBM, 1, 2, IRREPWBM, 
     &                        LISTW1)
                  IOFFW1R = I020
               ENDIF
C
               IOFFW1L = 0
               IOFFZR  = (INUMBM - 1)*DSZTMP*IINTFP + I000
C     
               IF(ISPIN .EQ. 1 .AND. .NOT. RHF) THEN
C
                  DO 40 IRREPF = 1, NIRREP
C
                     IRREPE = DIRPRD(IRREPF, IRREPWEF) 
                     IRREPJ = DIRPRD(IRREPF, IRREPX)
C                     
                     IOFFT  = IOFFT1(IRREPF, 2)
                     IOFFW1 = IOFFW1R + IOFFW1L
                     IOFFZ  = IOFFZR + IOFFZLEJ(IRREPE, 1)*IINTFP
C
                     NROW = POP(IRREPJ, 2)
                     NCOL = VRT(IRREPE, 1)
                     NSUM = VRT(IRREPF, 2)
C
                     ALPHA = ONE*FACTOR
                     BETA  = ZILCH
C
                     IF (MIN (NROW, NCOL, NSUM) .GT. 0) THEN 
                        CALL XGEMM('T', 'T', NROW, NCOL, NSUM, ALPHA, 
     &                              ICORE(IOFFT), NSUM, ICORE(IOFFW1),
     &                              NCOL, BETA, ICORE(IOFFZ), NROW)
                     ENDIF
C
                     IOFFW1L = IOFFW1L + NCOL*NSUM*IINTFP
 40               CONTINUE
C
               ELSE IF (ISPIN .EQ. 2 .OR. RHF) THEN
C
                  DO 41 IRREPE = 1, NIRREP
C
                     IRREPF = DIRPRD(IRREPE, IRREPWEF)
                     IRREPJ = DIRPRD(IRREPF, IRREPX)
C
                     IOFFT  = IOFFT1(IRREPF, 1)
                     IOFFW1 = IOFFW1R + IOFFW1L
                     IOFFZ  = IOFFZR + IOFFZLEJ(IRREPE, 2)*IINTFP
C
                     NROW = POP(IRREPJ, 1)
                     NCOL = VRT(IRREPE, 2)
                     NSUM = VRT(IRREPF, 1) 
C
                     ALPHA = ONE*FACTOR
                     BETA  = ZILCH
C
                     IF (MIN (NROW, NCOL, NSUM) .GT. 0) THEN
                        CALL XGEMM('T', 'N', NROW, NCOL, NSUM, ALPHA,
     &                              ICORE(IOFFT), NSUM, ICORE(IOFFW1),
     &                              NSUM, BETA, ICORE(IOFFZ), NROW)
                     ENDIF
                     IOFFW1L = IOFFW1L + NCOL*NSUM*IINTFP
 41               CONTINUE
               ENDIF
 30         CONTINUE
C     
C Now we have a (jE,Mb) (ISPIN = 1) or (Je,Bm) (ISPIN = 2 or RHF)
C ordered quantity. The next piece will be ordered (bM,Ej) (ISPIN = 1) or
C (Bm,Je) (ISPIN = 2 or RHF), so we need to reorder what we have to match
C this, therby allowing accumulation in matrix multiply opertaions.
C     
            IF(ISPIN .EQ. 1 .AND. .NOT. RHF) THEN
C
               CALL SYMTR1(IRREPGBM, POP(1,1), VRT(1,2), DSZTMP,
     &                     ICORE(I000), ICORE(I010), ICORE(I011),
     &                     ICORE(I012))
               CALL SYMTR3(IRREPGEJ, POP(1,2), VRT(1,1), DSZTMP,
     &                     DISTMP, ICORE(I000), ICORE(I010), 
     &                     ICORE(I011), ICORE(I012))
            ENDIF
C
            I020 = I010 + DISTMP*DSZTMP*IINTFP
C
            IF (I020 .GE. MXCOR) CALL INSMEM('MKT1GMBEJC2', I020,
     &                           MXCOR)
C
            CALL TRANSP(ICORE(I000), ICORE(I010), DISTMP, DSZTMP)
c YAU : old
c           CALL ICOPY(DISTMP*DSZTMP*IINTFP,ICORE(I010),1,ICORE(I000),1)
c YAU : new
            CALL DCOPY(DISTMP*DSZTMP,ICORE(I010),1,ICORE(I000),1)
c YAU : end
C     
C Now do W(MbEj) = - T(n,b)*Hbar(Mn,Ej) (ISPIN = 1)
C        W(mBeJ) = - T(N,B)*Hbar(mN,eJ) (ISPIN = 2 or RHF)
C     
            LISTW2 = 8 + ISPIN
            IF (RHF) LISTW2 = 10
            LSTTMP = 37 + ISPIN + (1 - IUHF)
C
            IRREPWEJ = IRREPGEJ
            IRREPWMN = IRREPWEJ
C     
            DSZW = IRPDPD(IRREPWMN, ISYTYP(1, LISTW2))
            DISW = IRPDPD(IRREPWEJ, ISYTYP(2, LISTW2))
C            
            DSZTMP = IRPDPD(IRREPMBI, ISYTYP(1, LSTTMP))
            DISTMP = IRPDPD(IRREPMAJ, ISYTYP(2, LSTTMP))
C     
C I020 now holds the now holds the (Mn,Ej)and [(Nm,JE)] integrals.
C     
            I020 = I010 + IINTFP*DSZW
            I030 = I020 + IINTFP*DSZW*DISW
C
            IF (I030 .LE. MXCOR) THEN
               INCORE = .TRUE.
               CALL GETLST(ICORE(I020), 1, DISW, 2, IRREPWEJ, LISTW2)
            ELSE
               INCORE = .FALSE.
               I030 = I020 + IINTFP*DSZW
            ENDIF
C
            DO 50 INUMEJ = 1, DISW
               IF (INCORE) THEN
                  IOFFW2R = (INUMEJ - 1)*DSZW*IINTFP + I020
               ELSE
                  CALL GETLST(ICORE(I020), INUMEJ, 1, 2, IRREPWEJ,
     &                        LISTW2)
                  IOFFW2R = I020
               ENDIF
C
               IOFFW2L = 0
               IOFFZR = (INUMEJ - 1)*DSZTMP*IINTFP + I000
C
               IF(ISPIN .EQ. 1 .AND. .NOT. RHF) THEN
C
                  DO 60 IRREPN = 1, NIRREP
C
                     IRREPM = DIRPRD(IRREPN, IRREPWMN)
                     IRREPB = DIRPRD(IRREPN, IRREPX)
C
                     IOFFT  = IOFFT1(IRREPB, 2)
                     IOFFW2 = IOFFW2R + IOFFW2L
                     IOFFZ  = IOFFZR + IOFFZLBM(IRREPM, 1)*IINTFP
C
                     NROW = VRT(IRREPB, 2)
                     NCOL = POP(IRREPM, 1)
                     NSUM = POP(IRREPN, 2)
C     
                     ALPHA = ONEM*FACTOR
                     BETA  = ONE
C
                     IF (MIN (NROW,NCOL,NSUM) .GT. 0) THEN
                        CALL XGEMM('N', 'T', NROW, NCOL, NSUM, ALPHA,
     &                              ICORE(IOFFT), NROW, ICORE(IOFFW2),
     &                              NCOL, BETA, ICORE(IOFFZ), NROW)
                     ENDIF
C     
                     IOFFW2L= IOFFW2L + NCOL*NSUM*IINTFP
 60               CONTINUE
C
               ELSE IF (ISPIN .EQ. 2 .OR. RHF) THEN
C
                  DO 61 IRREPM = 1, NIRREP
C
                     IRREPN = DIRPRD(IRREPM, IRREPWMN)
                     IRREPB = DIRPRD(IRREPN, IRREPX)
C
                     IOFFT  = IOFFT1(IRREPB, 1)
                     IOFFW2 = IOFFW2R + IOFFW2L
                     IOFFZ  = IOFFZR + IOFFZLBM(IRREPM, 2)*IINTFP
C
                     NROW = VRT(IRREPB, 1)
                     NCOL = POP(IRREPM, 2)
                     NSUM = POP(IRREPN, 1)
C
                     ALPHA = ONEM*FACTOR
                     BETA  = ONE
C
                     IF (MIN(NROW, NCOL, NSUM) .GT. 0) THEN
                        CALL XGEMM('N', 'N', NROW, NCOL, NSUM, ALPHA,
     &                              ICORE(IOFFT), NROW, ICORE(IOFFW2),
     &                              NSUM, BETA, ICORE(IOFFZ), NROW)
                     ENDIF
                     IOFFW2L= IOFFW2L + NCOL*NSUM*IINTFP
 61               CONTINUE
               ENDIF
 50         CONTINUE
C     
C Reorder to 
C   (bM,Ej) ->  (bM,Ej) (ISPIN = 1)
C   (Bm,Je) ->  (Bm,eJ) (ISPIN = 2 or RHF)
C     
            IF (ISPIN .EQ. 1 .AND. .NOT. RHF) THEN
               CONTINUE
               SSTSPN = 'ABBA'
            ELSE
               CALL SYMTR1(IRREPWEJ, POP(1,1), VRT(1,2), DSZTMP,
     &                     ICORE(I000), ICORE(I010), ICORE(I011),
     &                     ICORE(I012))
               SSTSPN = 'BAAB'
            ENDIF
C
            I020 = I010 + DISTMP*DSZTMP*IINTFP
            CALL TRANSP(ICORE(I000), ICORE(I010), DISTMP, DSZTMP)
C     
C Now write these to disk for each irrep.
C (Ej,bM) (ISPIN = 1)
C (eJ,Bm) (ISPIN = 2 or RHF)
C     
            CALL PUTLST(ICORE(I010), 1, DSZTMP, 1, IRREPGBM, LSTSCR)
C
 20      CONTINUE
C     
C Now switch ordering
C   (Ej,bM) -> (EM,bj) (ISPIN = 1)
C   (eJ,Bm) -> (em,BJ) (ISPIN = 2 or RHF)
C     
         ISCSIZ = (NVRTO(1) + NVRTO(2))*(NOCCO(1) + NOCCO(2))
         INSIZ  = IDSYMSZ(IRREPX, ISYTYP(1, LSTSCR), ISYTYP(2, LSTSCR))
         TARSIZ = IDSYMSZ(IRREPX, ISYTYP(1, LSTTAR), ISYTYP(2, LSTTAR))
C         
         I000 = 1
         I010 = I000 + TARSIZ*IINTFP
         I020 = I010 + INSIZ*IINTFP
         I030 = I020 + ISCSIZ
C     
         IF(I030 .GT. MXCOR) CALL INSMEM('MKT1GMBEJC2', I030, MXCOR)
C
         CALL GETALL(ICORE(I010), INSIZ, IRREPX, LSTSCR)
         CALL ALTSYMPCK0(ICORE(I010), ICORE(I000), INSIZ, TARSIZ,
     &                   ICORE(I020), IRREPX, SSTSPN)
C
         CALL PUTALL(ICORE(I000), TARSIZ, IRREPX, LSTTAR)
C     
C
         IF (IFLAGS(1) .GE. 20) THEN
C
            NSIZE = TARSIZ
C     
            CALL HEADER('Checksum @-QT1RING3', 0, LUOUT)
C     
            WRITE(LUOUT, *) SPCASE(ISPIN), ' = ', SDOT(NSIZE,
     &                      ICORE(I000), 1, ICORE(I000), 1)
         ENDIF
C
 10   CONTINUE
C 
      RETURN
      END
