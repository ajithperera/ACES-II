C
      SUBROUTINE MKT1GMBEJA2(ICORE, MAXCOR, IUHF, IRREPX, IOFFSET)
C
C This subroutine computes one of two T1*W contributions to the W(mbej)
C intermediate.  
C
C G(MB,EJ) = - SUM T(J,F)*Hbar(BM,FE) + SUM T(N,B)*Hbar(NM,JE)    (1)
C G(mb,ej) = - SUM T(j,f)*Hbar(bm,fe) + SUM T(n,b)*Hbar(nm,je)    (2)
C
C In the CC code F(EA) contibution from [SUM M F T(M,F)*Hbae(MA,FE)] was 
C also evaluated in this routine. That part is moved to a separate routine
C for convenience. The sign written here is opposite to the actual sign
C given in the JCP, 94, 4334, 1990. Here we compute the negtive of
C the actual expression given there. The reason will be clear when
C we evaluate the contribution of these ring terms to the doubles expression.
C
C This routine uses out of-core algorithms. 
C
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ONEM,ZILCH,ALPHA,BETA,FACTOR,SDOT
      DIMENSION ICORE(MAXCOR),IOFFT1(8,2),ABFULL(8),MNFULL(8)
C
      CHARACTER*4 SPCASE(2)
      LOGICAL INCORE
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYM/ POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF2AA,
     &             NF1BB,NF2BB
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /FILES/LUOUT,MOINTS
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
      DATA SPCASE /'AAAA','BBBB'/
C
C  First pick up T1 vector.
C
      TLIST = IAPRT1AA
C
      CALL GETPERT1(ICORE, MAXCOR, MXCOR, IUHF, IRREPX, TLIST, 
     &              IOFFSET, IOFFT1)
C
C  Spin cases AAAA and BBBB
C
      DO 10 ISPIN = 1, 1 + IUHF
         LSTTAR = (INGMCAA - 1) + ISPIN
         FACTOR = ONE
C
C Compute AB and MN DPD size for all irreps
C     
         DO 1000 IRREP = 1, NIRREP
            ICOUNT  = 0
            ICOUNT2 = 0
            DO 1001 IRREPE = 1, NIRREP
               IRREPF  = DIRPRD(IRREPE, IRREP)
               ICOUNT  = ICOUNT + VRT(IRREPF,ISPIN)*VRT(IRREPE,ISPIN)
               ICOUNT2 = ICOUNT2 + POP(IRREPF,ISPIN)*POP(IRREPE,ISPIN)
 1001       CONTINUE
            ABFULL(IRREP) = ICOUNT
            MNFULL(IRREP) = ICOUNT2
 1000    CONTINUE
C
C Loop over irreps : In the first block of code this corresponds to BM,
C while it is JE in the second block.
C First do the  G(MB,EJ) = SUM T(J,F)*HBar(BM,EF) contribution.
C The product is initialy packed as G(BM,EJ) matrix.
C
         DO 20 IRREPBM = 1, NIRREP
            IRREPWBM = IRREPBM
            IRREPWEF = IRREPWBM
            IRREPGBM = IRREPWBM
            IRREPGEJ = DIRPRD(IRREPGBM, IRREPX)
C
            ALPHA   = ONE*FACTOR   
            BETA    = ONE
            LSTWINT  = 26 + ISPIN
            IF(IUHF .EQ. 0) LSTWINT = 30
C     
            DSZTAR = IRPDPD(IRREPGBM, 8 + ISPIN)
            DISTAR = IRPDPD(IRREPGEJ, 8 + ISPIN)
            DSZW   = IRPDPD(IRREPWEF, ISYTYP(1, LSTWINT))
            DISW   = IRPDPD(IRREPWBM, ISYTYP(2, LSTWINT))
            DSZEXP = ABFULL(IRREPWEF)
C     
C I000 Holds the G(Je,bm) target for multiplication.
C I010 Holds area eventually used as scratch in symtr and the F(EA)
C      contribution computed by TRACER.
C I020 Holds the Hbar(EF,BM) integrals (enough space is allocated to hold
C      full list Hbar(EF,BM)).
C     
            I000 = 1
            I010 = I000 + IINTFP*DISTAR*DSZTAR
            I011 = I010 + IINTFP*DISTAR*DSZTAR
            I012 = I011 + IINTFP*DSZTAR
            I020 = I012 + IINTFP*DSZTAR
            I030 = I020 + IINTFP*DSZEXP*DISW
            I040 = I030 + IINTFP*MAX(DSZEXP, DISW)*(1-IUHF)
            IEND = I040 + IINTFP*MAX(DSZEXP, DISW)*(1-IUHF)
C     
            CALL IZERO(ICORE, DISTAR*DSZTAR*IINTFP)
C
            IF (IEND .LE. MXCOR) THEN
               INCORE = .TRUE.
               CALL GETLST(ICORE(I020), 1, DISW, 2, IRREPWBM, LSTWINT)
C
               IF(IUHF .NE. 0) THEN
                  CALL SYMEXP2(IRREPWEF, VRT(1, ISPIN), DSZEXP, DSZW,
     &                         DISW, ICORE(I020), ICORE(I020))
               ELSE
                  CALL ASSYM2A(IRREPWEF, VRT(1,1), DSZEXP, DISW, 
     &                         ICORE(I020), ICORE(I030), ICORE(I040))
               ENDIF
            ELSE
               INCORE = .FALSE.
C
               I030 = I020 + IINTFP*DSZEXP
               I030 = I020 + IINTFP*DSZW
               I040 = I030 + IINTFP*MAX(DSZEXP, DISW)*(1 - IUHF)
               IEND = I040 + IINTFP*MAX(DSZEXP, DISW)*(1 - IUHF)
               IF(IEND .GE. MXCOR) CALL INSMEM ('MKT1GMBEJA2', IEND,
     &         MXCOR)
C     
            ENDIF
C
            DO 30 INUMBM = 1, DISW
               IF(INCORE)THEN
                  IOFFW1R = (INUMBM - 1)*DSZEXP*IINTFP + I020
               ELSE
                  CALL GETLST(ICORE(I020), INUMBM, 1, 2, IRREPWBM,
     &                        LSTWINT)
                  IF(IUHF .NE. 0) THEN
                     CALL SYMEXP2(IRREPWEF, VRT(1,ISPIN), DSZEXP,
     &                             DSZW, 1, ICORE(I020), ICORE(I020))
                  ELSE
                     CALL ASSYM2A(IRREPWEF, VRT(1,ISPIN), DSZEXP, 1,
     &                            ICORE(I020), ICORE(I030), 
     &                            ICORE(I040))
                  ENDIF
                  IOFFW1R = I020
               ENDIF
C
               IOFFW1L = 0
               IOFFZL  = 0
               IOFFZR  = (INUMBM - 1)*DISTAR*IINTFP + I000
C
               DO 40 IRREPE = 1, NIRREP
C
                  IRREPF = DIRPRD(IRREPE, IRREPWEF)
                  IRREPJ = DIRPRD(IRREPF, IRREPX)
C     
                  IOFFT  = IOFFT1(IRREPF, ISPIN)
                  IOFFW1 = IOFFW1R + IOFFW1L
                  IOFFZ  = IOFFZR + IOFFZL 
C
                  NROW = POP(IRREPJ, ISPIN)
                  NCOL = VRT(IRREPE, ISPIN)
                  NSUM = VRT(IRREPF, ISPIN)
C
                  ALPHA = ONEM*FACTOR
                  BETA = ZILCH
C
                  IF (MIN (NROW, NCOL, NSUM) .GT. 0) THEN
                     CALL XGEMM('T', 'N', NROW, NCOL, NSUM, ALPHA, 
     &                           ICORE(IOFFT), NSUM, ICORE(IOFFW1),
     &                           NSUM, BETA, ICORE(IOFFZ), NROW)
                  ENDIF
C
                  IOFFW1L = IOFFW1L + NCOL*NSUM*IINTFP
                  IOFFZL  = IOFFZL + NROW*NCOL*IINTFP
C
 40            CONTINUE
 30         CONTINUE
C     
C Now transpose it to get it in the same ordering as the next piece.
C G(BM,JE)
C     
            IF (I020 .GE. MXCOR) CALL INSMEM ('MKT1GMBEJA2', I020,
     &         MXCOR)
C
            CALL TRANSP(ICORE(I000), ICORE(I010), DSZTAR, DISTAR)
c YAU : old
c           CALL ICOPY(IINTFP*DISTAR*DSZTAR,ICORE(I010),1,ICORE(I000),1)
c YAU : new
            CALL DCOPY(DISTAR*DSZTAR,ICORE(I010),1,ICORE(I000),1)
c YAU : end
C
C Now do the - SUM T(N,B)*Hbar(NM,JE) contribution - This is ordered
C G(BM,JE)
C     
            LSTINT = 6 + ISPIN
C
            IRREPWEJ = IRREPGEJ
            IRREPWMN = IRREPWEJ
C
            DSZW  = IRPDPD(IRREPWMN ,ISYTYP(1, LSTINT))
            DISW  = IRPDPD(IRREPWEJ, ISYTYP(2, LSTINT))
C
            DSZEXP = MNFULL(IRREPWMN)
C     
C I020 holds the Hbar(NM,JE) integrals (enough space is allocated to hold
C full list Hbar(NM,JE).
C     
            I020 = I010 + IINTFP*MAX(DSZEXP, MNFULL(1))
            I030 = I020 + IINTFP*DSZEXP*DISW
C
            IF (I030 .LE. MXCOR) THEN
               INCORE = .TRUE.
               CALL GETLST(ICORE(I020), 1, DISW, 2, IRREPWEJ, LSTINT)
               CALL SYMEXP2(IRREPWMN, POP(1,ISPIN), DSZEXP, DSZW, DISW,
     &                      ICORE(I020), ICORE(I020))
            ELSE
               INCORE = .FALSE.
               I030 = I020 + IINTFP*DSZW
            ENDIF
C
            DO 50 INUMJE = 1, DISW
               IF(INCORE) THEN
                  IOFFW2R = (INUMJE - 1)*DSZEXP*IINTFP + I020
               ELSE
                  CALL GETLST(ICORE(I020), INUMJE, 1, 2, IRREPWEJ,
     &                        LSTINT)
                  CALL SYMEXP2(IRREPWMN, POP(1,ISPIN), DSZEXP, DSZW, 1,
     &                         ICORE(I020), ICORE(I020))
                  IOFFW2R = I020
               ENDIF
C
               IOFFW2L = 0
               IOFFZR  = (INUMJE-1)*DSZTAR*IINTFP + I000
               IOFFZL  = 0
C
               DO 60 IRREPM = 1, NIRREP
C
                  IRREPN  = DIRPRD(IRREPM, IRREPWMN)
                  IRREPBO = DIRPRD(IRREPN, IRREPX) 
C
                  IOFFT  = IOFFT1(IRREPBO, ISPIN)
                  IOFFW2 = IOFFW2R + IOFFW2L
                  IOFFZ  = IOFFZR + IOFFZL
C
                  NROW = VRT(IRREPBO, ISPIN)
                  NCOL = POP(IRREPM, ISPIN)
                  NSUM = POP(IRREPN, ISPIN)
C
                  ALPHA = ONE*FACTOR
                  BETA  = ONE
C
                  IF (MIN (NROW, NCOL, NSUM) .GT. 0) THEN
                     CALL XGEMM('N', 'N', NROW, NCOL, NSUM, ALPHA,
     &                           ICORE(IOFFT), NROW, ICORE(IOFFW2),
     &                           NSUM, BETA, ICORE(IOFFZ), NROW)
                  ENDIF
C
                  IOFFZL  = IOFFZL + NROW*NCOL*IINTFP
                  IOFFW2L = IOFFW2L + NCOL*NSUM*IINTFP
C
 60            CONTINUE
 50         CONTINUE
C     
C Now switch ordering to (BM,EJ) and write these quantities to disk
C for each irrep.
C     
            CALL SYMTR1(IRREPGEJ, POP(1,ISPIN), VRT(1,ISPIN), DSZTAR,
     &                  ICORE(I000), ICORE(I010), ICORE(I011),
     &                  ICORE(I012))
C     
            I020 = I010 + DISTAR*DSZTAR*IINTFP
            CALL TRANSP(ICORE(I000), ICORE(I010), DISTAR, DSZTAR)
            CALL PUTLST(ICORE(I010), 1, DSZTAR, 1, IRREPGBM, LSTTAR)
 20      CONTINUE
C     
C Now switch ordering to (BJ,EM)
C     
         INDEX  = 8 + ISPIN
         ISCSIZ = NOCCO(ISPIN)*NVRTO(ISPIN)*2
         INSIZ  = IDSYMSZ(IRREPX, INDEX, INDEX)
         TARSIZ = IDSYMSZ(IRREPX, ISYTYP(1, LSTTAR), ISYTYP(2, LSTTAR))
C
         I000 = 1
         I010 = I000 + TARSIZ*IINTFP
         I020 = I010 + INSIZ*IINTFP
         I030 = I020 + ISCSIZ
C
         IF(I030 .GT. MXCOR) CALL INSMEM('T1QGMBAA', I030, MXCOR)
C
         CALL GETALL(ICORE(I010), INSIZ, IRREPX, LSTTAR)
         CALL ALTSYMPCK0(ICORE(I010), ICORE(I000), INSIZ, TARSIZ,
     &                   ICORE(I020), IRREPX, SPCASE(ISPIN))
C     
         CALL PUTALL(ICORE(I000), TARSIZ, IRREPX, LSTTAR)

C  Now do transposition to get G(EM,BJ) ordering by irrep
C     
            IF (IFLAGS(1) .GE. 20) THEN
C
               NSIZE = TARSIZ
C     
               CALL HEADER('Checksum @-QT1RING1', 0, LUOUT)
C     
               WRITE(LUOUT, *) SPCASE(ISPIN), ' = ', SDOT(NSIZE,
     &                      ICORE(I000), 1, ICORE(I000), 1)
            ENDIF
C
 10   CONTINUE
C
      RETURN
      END
