C
      SUBROUTINE MKT1GMBEJC1(ICORE, MAXCOR, IUHF, IRREPX, IOFFSET)
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
C This routine uses in core algorithm.
C
      IMPLICIT INTEGER (A-Z)
      LOGICAL CHANGE, RHF
      DOUBLE PRECISION ONE,ONEM,ZILCH,ALPHA,BETA,FACTOR, SDOT
      DIMENSION ICORE(MAXCOR),IOFFT1(8,2), IZOFF1(8), IZOFF2(8)
C
      CHARACTER*4 SSTSPN
      CHARACTER*4 SPCASE(2)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),
     &                DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /FILES/ LUOUT, MOINTS
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /SYM/ POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF2AA,
     &             NF1BB,NF2BB
      COMMON /INFO/ NOCCO(2),NVRTO(2)
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
      IF(IUHF .EQ. 0) RHF = .TRUE.
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
         DO  20 IRREPBM = 1, NIRREP
C
            IRREPWBM = IRREPBM
            IRREPWEF = IRREPWBM
            IRREPGBM = IRREPWBM
            IRREPGEJ = DIRPRD(IRREPGBM, IRREPX)
            IRREPMBI = IRREPGBM
            IRREPMAJ = IRREPGEJ
C
            LSTTMP = 37 + ISPIN + (1 - IUHF)
C
            DSZTAR = IRPDPD(IRREPGEJ, ISYTYP(1, LSTTAR))
            DISTAR = IRPDPD(IRREPGBM, ISYTYP(2, LSTTAR))
C
            DSZTMP = IRPDPD(IRREPMBI, ISYTYP(1, LSTTMP))            
            DISTMP = IRPDPD(IRREPMAJ, ISYTYP(2, LSTTMP))
C
C First do    
C  W(MbEj) = SUM T(j,f)*Hbar(Mb,Ef) (ISPIN = 1)
C  W(mBeJ) = SUM T(J,F)*Hbar(mB,eF) (ISPIN = 2 OR RHF).
C     
C This product is initially packed Mb,Ej [Bm,eJ].
C     
            CALL IZERO(ICORE, IINTFP*DSZTMP*DISTMP)
C
            ALPHA = ONE*FACTOR
            BETA  = ONE
C
            LISTW1 = 28 + ISPIN
            IF (RHF) LISTW1 = 30
C
            DSZW = IRPDPD(IRREPWEF, ISYTYP(1, LISTW1))
            DISW = IRPDPD(IRREPWBM, ISYTYP(2, LISTW1))
C
            I000 = 1
            I010 = I000 + IINTFP*DISTMP*DSZTMP
            I020 = I010 + IINTFP*DSZW*DISW
            I030 = I020 + IINTFP*DSZW
C     
C Read integrals I(Mb,Ef) [I(Bm,Fe)] matrix.
C     
            CALL GETTRN(ICORE(I010), ICORE(I020), DSZW, DISW, 2,
     &                  IRREPWBM, LISTW1)
            ITMP = DISW
            DISW = DSZW
            DSZW = ITMP
C     
C Transpose ket indices if ISPIN = 2
C     
            IF (ISPIN .EQ. 2 .OR. RHF) THEN
C
               I030 = I020 + IINTFP*MAX(DSZW, DISW)
               I040 = I030 + IINTFP*MAX(DSZW, DISW)
               I050 = I040 + IINTFP*MAX(DSZW, DISW)
C
               CALL SYMTR1(IRREPWEF, VRT(1,1), VRT(1,2), DSZW, 
     &                     ICORE(I010), ICORE(I020), ICORE(I030),
     &                     ICORE(I040))
            ENDIF
C     
C We have now
C  I(Mb,Ef)  [ISPIN = 1]
C  I(Bm,eF)  [ISPIN = 2 or RHF]
C     
C Setup and perform the matrix multiply.
C  Z(Mb,Ej) = I(MbE,f)*T(f,j) [ISPIN = 1]
C  Z(Bm,eJ) = I(Bme,F)*T(F,J) [ISPIN = 2 or RHF]
C     
            IOFFI = I010
            IZOFF1(1) = I000
C
C Set up an offset array for the target index. 
C     
            DO 5000 IRREPZJ = 2, NIRREP
               IRREP = IRREPZJ - 1
               IRREPZE = DIRPRD(IRREP, IRREPGEJ)
               IZOFF1(IRREPZJ) = IZOFF1(IRREPZJ - 1) +
     &                           DSZW*VRT(IRREPZE, ISPIN)*
     &                           POP(IRREP, 3 - ISPIN)*IINTFP
 5000       CONTINUE
C     
            DO 100 IRREPF = 1, NIRREP
C
               IRREPE = DIRPRD(IRREPF, IRREPWEF)
               IRREPJ = DIRPRD(IRREPF, IRREPX)
C
               NROWI = DSZW*VRT(IRREPE, ISPIN)
               NCOLI = VRT(IRREPF, 3 - ISPIN)
               NROWT = VRT(IRREPF, 3 - ISPIN)
               NCOLT = POP(IRREPJ, 3 - ISPIN)
C
               NROWZ = NROWI
               NCOLZ = NCOLT
C
               IOFFT = IOFFT1(IRREPF, 3 - ISPIN)
               IOFFZ= IZOFF1(IRREPJ)
C
               IF (MIN (NROWZ, NCOLZ, NCOLI) .GT. 0) THEN
                  CALL XGEMM('N', 'N', NROWZ, NCOLZ, NCOLI, ALPHA,
     &                        ICORE(IOFFI), NROWI, ICORE(IOFFT),
     &                        NROWT, BETA, ICORE(IOFFZ), NROWZ)
               ENDIF
C
               IOFFI = IOFFI + IINTFP*NROWI*NCOLI
C
 100        CONTINUE
C     
C Now we have a (Mb,Ej) (ISPIN = 1) or (Bm,eJ) (ISPIN = 2 or RHF) ordered
C quantity. The next piece will be ordered (Ej,Mb) (ISPIN = 1) or
C (Je,mB) (ISPIN = 2 or RHF), so we need to reorder what we have to match
C this, therby allowing accumulation in matrix multiply operations.
C     
            IF (ISPIN .EQ. 2 .OR. RHF) THEN
C
               I020 = I010 + IINTFP*MAX(DSZTMP, DISTMP)
               I030 = I020 + IINTFP*MAX(DSZTMP, DISTMP)
               I040 = I030 + IINTFP*MAX(DSZTMP, DISTMP)
C
               CALL SYMTR1(IRREPGEJ, VRT(1,2), POP(1,1), DSZTMP, 
     &                     ICORE(I000), ICORE(I010), ICORE(I020),
     &                     ICORE(I030))
               CALL SYMTR3(IRREPGBM, VRT(1,1), POP(1,2), DSZTMP,
     &                     DISTMP, ICORE(I000), ICORE(I010),
     &                     ICORE(I020),ICORE(I030))
            ENDIF
C
            I020 = I010 + IINTFP*DSZTMP*DISTMP
            CALL TRANSP(ICORE(I000), ICORE(I010), DISTMP, DSZTMP)
c YAU : old
c           CALL ICOPY(IINTFP*DISTMP*DSZTMP,ICORE(I010),1,ICORE(I000),1)
c YAU : new
            CALL DCOPY(DISTMP*DSZTMP,ICORE(I010),1,ICORE(I000),1)
c YAU : end
C     
C Now do W(MbEj) = - T(n,b)*Hbar(Mn,Ej) (ISPIN = 1)
C        W(mBeJ) = - T(N,B)*Hbar(mN,eJ) (ISPIN = 2 or RHF)
C
            ALPHA = ONEM*FACTOR
            BETA  = ONE
C
            LISTW2 = 8 + ISPIN
            IF (RHF) LISTW2 = 10
            LSTTMP = 40 - ISPIN
C
            IRREPWEJ = IRREPGEJ
            IRREPWMN = IRREPWEJ
C
            DSZTMP = IRPDPD(IRREPMAJ, ISYTYP(1, LSTTMP))
            DISTMP = IRPDPD(IRREPMBI, ISYTYP(2, LSTTMP))
C
            DSZW = IRPDPD(IRREPWMN, ISYTYP(1, LISTW2))
            DISW = IRPDPD(IRREPWEJ, ISYTYP(2, LISTW2))
C     
            I020 = I010 + IINTFP*DSZW*DISW
            I030 = I020 + IINTFP*DSZW
C
            CALL GETTRN(ICORE(I010), ICORE(I020), DSZW, DISW, 2,
     &                   IRREPWEJ, LISTW2)
            ITMP = DISW
            DISW = DSZW
            DSZW = ITMP
C     
C We now have
C  I(Ej,Mn) [ISPIN = 1]
C  I(Je,Nm) [ISPIN = 2 or RHF]
C     
            IF(ISPIN .EQ. 2 .OR. RHF) THEN
C
               I030 = I020 + IINTFP*MAX(DSZW, DISW)
               I040 = I030 + IINTFP*MAX(DSZW, DISW)
               I050 = I040 + IINTFP*MAX(DSZW, DISW)
C
               CALL SYMTR1(IRREPWMN, POP(1,1), POP(1,2), DSZW,
     &                     ICORE(I010), ICORE(I020), ICORE(I030),
     &                     ICORE(I040))
            ENDIF
C     
C We now have
C  I(Ej,Mn) [ISPIN = 1]
C  I(Je,mN) [ISPIN = 2 or RHF]
C
C And we can do the matrix multiply
C  Z(Ej,Mb) = I(Ej,Mn)*T(b,n)  [ISPIN = 1]
C  Z(Je,mB) = I(Je,mN)*T(B,N)  [ISPIN = 2 or RHF]
C     
            IOFFI = I010
C
C Set up offsets for the target.
C            
            IZOFF2(1) = I000
C
            DO 5100 IRREPZB = 2, NIRREP
               IRREP = IRREPZB - 1
               IRREPZM = DIRPRD(IRREP, IRREPGBM)
               IZOFF2(IRREPZB) = IZOFF2(IRREPZB - 1) + 
     &                           DSZW*POP(IRREPZM, ISPIN)*
     &                           VRT(IRREP, 3 - ISPIN)*IINTFP
 5100       CONTINUE
C
            DO 200 IRREPN = 1, NIRREP
C
               IRREPM = DIRPRD(IRREPN, IRREPWMN)
               IRREPB = DIRPRD(IRREPN, IRREPX)
C     
               NROWI = DSZW*POP(IRREPM, ISPIN)
               NCOLI = POP(IRREPN, 3 - ISPIN)
               NROWT = VRT(IRREPB, 3 - ISPIN)
               NCOLT = POP(IRREPN, 3 - ISPIN)
C
               NROWZ = NROWI
               NCOLZ = NROWT
C
               IOFFT = IOFFT1(IRREPB, 3 - ISPIN)
               IOFFZ = IZOFF2(IRREPB)
C
               IF (MIN (NROWZ, NCOLZ, NCOLI) .GT. 0) THEN
C
                  CALL XGEMM('N', 'T', NROWZ, NCOLZ, NCOLI, ALPHA,
     &                        ICORE(IOFFI), NROWI, ICORE(IOFFT),
     &                        NROWT, BETA, ICORE(IOFFZ), NROWZ)
               ENDIF
C
               IOFFI = IOFFI + IINTFP*NROWI*NCOLI
C
 200        CONTINUE
C     
C Reorder to 
C  (Ej,Mb) ->  (Ej,bM) [ISPIN = 1]
C  (Je,mB) ->  (eJ,Bm) [ISPIN = 2 or RHF]
C     
            I020 = I010 + IINTFP*MAX(DSZTMP, DISTMP)
            I030 = I020 + IINTFP*MAX(DSZTMP, DISTMP)
            I040 = I030 + IINTFP*MAX(DSZTMP, DISTMP)
C
            CALL SYMTR1(IRREPGBM, POP(1,ISPIN), VRT(1,3-ISPIN),
     &                  DSZTMP, ICORE(I000), ICORE(I010), ICORE(I020),
     &                  ICORE(I030))
C
            IF (ISPIN .EQ. 2 .OR. RHF) THEN
C
               I020 = I010 + IINTFP*MAX(DSZTMP, DISTMP)
               I030 = I020 + IINTFP*MAX(DSZTMP, DISTMP)
               I040 = I030 + IINTFP*MAX(DSZTMP, DISTMP)
C
               CALL SYMTR3(IRREPGEJ, POP(1,1), VRT(1,2), DSZTMP,
     &                     DISTMP, ICORE(I000), ICORE(I010),
     &                     ICORE(I020), ICORE(I030))
            ENDIF
C
            IF (ISPIN .EQ. 1 .AND. .NOT. RHF) THEN
               SSTSPN = 'ABBA'
            ELSE
               SSTSPN = 'BAAB'
            ENDIF 
C
C Now write these to disk for each irrep.     
C
            CALL PUTLST(ICORE(I000), 1, DISTMP, 1, IRREPGBM, LSTSCR)
C
 20      CONTINUE
C     
C Now switch ordering.
C  (Ej,bM) ->  (EM,bj) [ISPIN = 1]
C  (eJ,Bm) ->  (em,BJ) [ISPIN = 2 or RHF]
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
         IF (I030 .GT. MXCOR) CALL INSMEM('MKT1GMBEJC1', I030, MXCOR)
C
         CALL GETALL(ICORE(I010), INSIZ, IRREPX, LSTSCR)
C
         CALL ALTSYMPCK0(ICORE(I010), ICORE(I000), INSIZ, TARSIZ, 
     &                   ICORE(I020), IRREPX, SSTSPN)
C
         CALL PUTALL(ICORE(I000), TARSIZ, IRREPX, LSTTAR)
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
