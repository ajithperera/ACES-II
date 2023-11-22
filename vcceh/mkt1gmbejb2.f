C
      Subroutine MKT1GMBEJB2(ICORE, MAXCOR, IUHF, IRREPX, IOFFSET)
C
C This subroutine computes two T1 contributions to the 
C W(MB,EJ) intermediate.  
C
C     G(mBEj) =  SUM T(j,f)*Hbar(mB,fE) - SUM T(N,B)*Hbar(mN,jE)
C     G(MbeJ) =  SUM T(J,F)*Hbar(Mb,Fe) - SUM T(n,b)*Hbar(Mn,Je) (UHF only)
C
C In the CC code F(EA) contibution from [SUM M F T(M,F)*Hbae(Ma,Fe)] was 
C also evaluated in this routine. That part is moved to a separate routine
C for convenience. The sign written here is opposite to the actual sign
C given in the JCP, 94, 4334, 1990. Here we compute the negtive of
C the actual expression given there. The reason will be clear when
C we evaluate the contribution of these ring terms to the doubles expression.
C
C This routine uses out of core algorithms.
C
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ONEM,ZILCH,ALPHA,BETA,FACTOR,SDOT
      CHARACTER*4 SSTSPN(2)
      CHARACTER*4 SPCASE(2)
      LOGICAL RHF, INCORE
C
      DIMENSION ICORE(MAXCOR),IOFFT1(8,2),SIZEVV(2),IOFFZLEJ(8,2),
     &          IOFFZLBM(8,2)
C
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
      DATA SSTSPN /'ABAB','BABA'/
      DATA SPCASE /'BAAB','ABBA'/
C
C First pick up T1 vector
C
      TLIST = IAPRT1AA
C
      CALL GETPERT1(ICORE, MAXCOR, MXCOR, IUHF, IRREPX, TLIST, 
     &              IOFFSET, IOFFT1)
      RHF = .FALSE.
      IF(IUHF .EQ. 0) RHF = .TRUE.
C
C Compute sizes of a full VV distribution which transform as the irrepx.
C     
      CALL IZERO(SIZEVV, 2)
      DO 3000 ISPIN = 1, 2
         DO 3001 IRREPL = 1, NIRREP
            IRREPR = DIRPRD(IRREPL, IRREPX)
            SIZEVV(ISPIN) = SIZEVV(ISPIN) + VRT(IRREPL,ISPIN)*
     &                      VRT(IRREPR,ISPIN)
 3001    CONTINUE
 3000 CONTINUE
C     
C Spin cases BAAB and ABBA respectively.
C     
      DO 10 ISPIN = 1, 1 + IUHF
         LSTTAR = (INGMCBAAB - 1) + ISPIN
         FACTOR = ONE
C
C Loop over irreps - In the first block of code this corresponds mB,
C while it is jE in the second block.
C
         DO 20 IRREPBM = 1, NIRREP
C
            IRREPWBM = IRREPBM
            IRREPWEF = IRREPWBM
            IRREPGBM = IRREPWBM
            IRREPGEJ = DIRPRD(IRREPGBM, IRREPX)
C
C Computes offfset into an i,A distribution for the target array.
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
C Computes dimensions of target matrix.
C
            IF (ISPIN .EQ. 1) THEN
               INDEX = 11
            ELSE IF (ISPIN .EQ. 2) THEN
               INDEX = 12
            ENDIF
C            
            DSZTAR = IRPDPD(IRREPGBM, INDEX)
            DISTAR = IRPDPD(IRREPGEJ, INDEX)
C
            TARSIZ = DISTAR*DSZTAR
C
C First do the Z(mB,Ej) = SUM T(j,f)*Hbar(mB,fE)
C  Z(Mb,eJ) = SUM T(J,F)*Hbar(Mb,Fe)
C This product is initially packed (jE,Bm) [(Je-Mb)].
C     
            LISTW1 = 31 - ISPIN
C
            DSZW = IRPDPD(IRREPWEF, ISYTYP(1, LISTW1))
            DISW = IRPDPD(IRREPWBM, ISYTYP(2, LISTW1))
C     
C I000 Holds the Z(jE,mB) target.
C I010 Holds an area eventually used as scratch in SYMTR1
C I011 and I012 scratch arrays for SYMTR1.
C I020 Holds the Hbar(Ef,Bm) integrals. 
C     
            I000 = 1
            I010 = I000 + IINTFP*DISTAR*DSZTAR
            I011 = I010 + IINTFP*MAX(DSZTAR, DISTAR, DSZW, 
     &                    SIZEVV(ISPIN))
            I012 = I011 + IINTFP*MAX(DSZTAR, DISTAR)  
            I020 = I012 + IINTFP*MAX(DSZTAR, DISTAR)
            I030 = I020 + IINTFP*DSZW*DISW
C
            CALL IZERO(ICORE, IINTFP*DISTAR*DSZTAR)
C
            IF (I030 .LE. MXCOR) THEN
               INCORE = .TRUE.
               CALL GETLST(ICORE(I020), 1, DISW, 2, IRREPWBM, LISTW1)
            ELSE
               INCORE = .FALSE.
               I030 = I020 + IINTFP*DSZW
            ENDIF
C
            DO 30 INUMBM = 1 , DISW
C
               IF(INCORE)THEN
                  IOFFW1R = (INUMBM - 1)*DSZW*IINTFP + I020
               ELSE
                  CALL GETLST(ICORE(I020), INUMBM, 1, 2, IRREPWBM, 
     &                        LISTW1)
                  IOFFW1R = I020
               ENDIF
C
               IOFFW1L = 0
               IOFFZR = (INUMBM - 1)*DISTAR*IINTFP + I000
C
               IF (ISPIN .EQ. 1) THEN
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
C                     
 40               CONTINUE
C
               ELSE
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
                     BETA = ZILCH
C
                     IF (MIN(NROW, NCOL, NSUM) .GT. 0) THEN
                        CALL XGEMM('T', 'N', NROW, NCOL, NSUM, ALPHA, 
     &                              ICORE(IOFFT), NSUM, ICORE(IOFFW1),
     &                              NSUM, BETA, ICORE(IOFFZ), NROW)
                     ENDIF
C
                     IOFFW1L = IOFFW1L + NCOL*NSUM*IINTFP
C 
 41               CONTINUE
               ENDIF
 30         CONTINUE
C     
C The code above produces mbej intermediates in the order.
C   jE,Bm (ISPIN = 1)
C   Je,Mb (ISPIN = 2)
C     
C The next piece will be evaluated as
C   Bm,Ej (ISPIN = 1, UHF)
C   bM,Je (ISPIN = 1, RHF)
C   bM,Je (ISPIN = 2)
C     
C We need to reorder what we have so that it will match up and 
C can be evaluated by XGEMM. Onec again note the change from
C DSZTAR to DISTAR.
C     
            IF (ISPIN .EQ. 1 .AND. .NOT. RHF) THEN
C
               CALL SYMTR3(IRREPGEJ, POP(1,2), VRT(1,1), DISTAR,
     &                     DSZTAR, ICORE(I000), ICORE(I010),  
     &                     ICORE(I011), ICORE(I012))
C
               I020 = I010 + DSZTAR*DISTAR*IINTFP
C     
               IF (I020 .GE. MXCOR) CALL INSMEM ('M3KT1GMBEJA2', I020,
     &         MXCOR)
C
               CALL TRANSP(ICORE(I000), ICORE(I010), DSZTAR, DISTAR)
c YAU : old
c              CALL ICOPY(DISTAR*DSZTAR*IINTFP,ICORE(I010),1,ICORE(I000),1)
c YAU : new
               CALL DCOPY(DISTAR*DSZTAR,ICORE(I010),1,ICORE(I000),1)
c YAU : end
C
            ELSE IF (ISPIN .EQ. 2) THEN
C
               CALL SYMTR1(IRREPGBM, POP(1,1), VRT(1,2), DISTAR,
     &                     ICORE(I000), ICORE(I010), ICORE(I011),
     &                     ICORE(I012))
C
               I020 = I010 + DSZTAR*DISTAR*IINTFP
C
               IF (I020 .GE. MXCOR) CALL INSMEM ('M3KT1GMBEJA2', I020,
     &         MXCOR)
C
               CALL TRANSP(ICORE(I000), ICORE(I010), DSZTAR, DISTAR)
c YAU : old
c              CALL ICOPY(DISTAR*DSZTAR*IINTFP,ICORE(I010),1,ICORE(I000),1)
c YAU : new
               CALL DCOPY(DISTAR*DSZTAR,ICORE(I010),1,ICORE(I000),1)
c YAU : end
C     
            ELSE IF (RHF) THEN
C
               I020 = I010 + DSZTAR*DISTAR*IINTFP
C
               IF (I020 .GE. MXCOR) CALL INSMEM ('M3KT1GMBEJA2', I020,
     &         MXCOR)
C
               CALL TRANSP(ICORE(I000), ICORE(I010), DSZTAR, DISTAR)
c YAU : old
c              CALL ICOPY(DISTAR*DSZTAR*IINTFP,ICORE(I010),1,ICORE(I000),1)
c YAU : new
               CALL DCOPY(DISTAR*DSZTAR,ICORE(I010),1,ICORE(I000),1)
c YAU : end
            ENDIF
C
C Now do the second contraction
C   Z(mB,Ej) = Z(mB,Ej) - SUM T(N,B)*Hbar(mN,jE) (ISPIN = 1)
C   Z(Mb,eJ) = Z(Mb,eJ) - SUM T(n,b)*Hbar(Mn,Je) (ISPIN = 2, or RHF)
C This product is initially packed [Bm,Ej] [(bm,Je)]
C     
            LISTW2 = 8 + ISPIN
            IF(RHF) LISTW2 = 10
C
            IRREPWEJ = IRREPGEJ
            IRREPWMN = IRREPWEJ
C            
            DSZW = IRPDPD(IRREPWMN, ISYTYP(1, LISTW2))
            DISW = IRPDPD(IRREPWEJ, ISYTYP(2, LISTW2))
C     
C I020 holds Hbar(Nm,Ej) [(Mn,Je)] integrals.
C     
            I020 = I010 + IINTFP*DSZW
            I030 = I020 + IINTFP*DSZW*DISW
C
            IF(I030 .LE. MXCOR) THEN
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
               IOFFZR = (INUMEJ - 1)*DSZTAR*IINTFP + I000
C
               IF(ISPIN .EQ. 1 .AND. .NOT. RHF) THEN
C
                  DO 60 IRREPM = 1, NIRREP
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
C
                     IOFFW2L = IOFFW2L + NCOL*NSUM*IINTFP
C
 60               CONTINUE
C
               ELSEIF (ISPIN .EQ. 2 .OR. RHF) THEN
C
                  DO 61 IRREPN = 1, NIRREP
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
                     IF (MIN(NROW, NCOL, NSUM) .GT. 0) THEN
                        CALL XGEMM('N', 'T', NROW, NCOL, NSUM, ALPHA,
     &                              ICORE(IOFFT), NROW, ICORE(IOFFW2),
     &                              NCOL, BETA, ICORE(IOFFZ), NROW)
                     ENDIF
C
                     IOFFW2L = IOFFW2L + NCOL*NSUM*IINTFP
C
 61               CONTINUE
               ENDIF
 50         CONTINUE
C
            IF (ISPIN .EQ. 2 .OR. RHF) THEN
               CALL SYMTR1(IRREPWEJ, POP(1,1), VRT(1,2), DSZTAR,
     &                     ICORE(I000), ICORE(I010), ICORE(I011),
     &                     ICORE(I012))
            ENDIF
C     
C Take the transpose and write these to disk for each irrep. These are ordered
C  (Bm,Ej) (ISPIN = 1, UHF)
C  (bM,eJ) (ISPIN = 2 or RHF) 
C     
            I020 = I010 + IINTFP*DISTAR*DSZTAR
C     
            CALL TRANSP(ICORE(I000), ICORE(I010), DISTAR, DSZTAR)
            CALL PUTLST(ICORE(I010), 1, DSZTAR, 1, IRREPGBM, LSTTAR)
C
 20      CONTINUE
C     
C Now switch ordering from (Bm,Ej) to (Bj,Em)
C     
         ISCSIZ = (NVRTO(1) + NVRTO(2))*(NOCCO(1) + NOCCO(2))
         INSIZ  = IDSYMSZ(IRREPX, INDEX, INDEX)
         TARSIZ = IDSYMSZ(IRREPX, ISYTYP(1, LSTTAR), ISYTYP(2, LSTTAR))
C
         I000 = 1
         I010 = I000 + TARSIZ*IINTFP
         I020 = I010 + INSIZ*IINTFP
         I030 = I020 + ISCSIZ
C
         IF (I030 .GT. MXCOR) CALL INSMEM('MKT1GMBEJB2', I030, MXCOR)
C
         CALL GETALL(ICORE(I010), INSIZ, IRREPX, LSTTAR)
         CALL ALTSYMPCK0(ICORE(I010), ICORE(I000), INSIZ, TARSIZ,
     &                   ICORE(I020), IRREPX, SSTSPN(ISPIN))
C     
C Now do transposition to get  E,m-B,j ordering irrep by irrep.
C     
         CALL PUTALL(ICORE(I000), TARSIZ, IRREPX, LSTTAR)
C
         IF (IFLAGS(1) .GE. 20) THEN
C
            NSIZE = TARSIZ
     
            CALL HEADER('Checksum @-QT1RING2', 0, LUOUT)
     
            WRITE(LUOUT, *) SPCASE(ISPIN), ' = ', SDOT(NSIZE,
     &                       ICORE(I000), 1, ICORE(I000), 1)
         ENDIF
C
 10   CONTINUE
C
      RETURN
      END
