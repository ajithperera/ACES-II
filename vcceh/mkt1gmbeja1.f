C
      SUBROUTINE MKT1GMBEJA1(ICORE, MAXCOR, IUHF, IRREPX, IOFFSET)
C
C This subroutine computes one of two T1*W contributions to the W(mbej)
C intermediate.  
C
C  G(MB,EJ) = - SUM T(J,F)*Hbar(BM,FE) + SUM T(N,B)*Hbar(NM,JE)    (1)
C  G(mb,ej) = - SUM T(j,f)*Hbar(bm,fe) + SUM T(n,b)*Hbar(nm,je)    (2)
C
C In the CC code F(EA) contibution from [SUM M F T(M,F)*Hbae(MA,FE)] was 
C also evaluated in this routine. That part is moved to a separate routine
C for convenience. The sign written here is opposite to the actual sign
C given in the JCP, 94, 4334, 1990. Here we compute the negtive of
C the actual expression given there. The reason will be clear when
C we evaluate the contribution of these ring terms to the doubles expression.
C
C  This routine use in-core algorithms.
C
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ONEM,ZILCH,ALPHA,BETA,FACTOR,SDOT
      CHARACTER*4 SPCASE(2)
      DIMENSION ICORE(MAXCOR),IOFFT1(8,2),ABFULL(8),MNFULL(8),
     &          IZOFF1(8), IZOFF2(8)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /FILES/ LUOUT,MOINTS
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
      DATA ONE /1.0D+00/
      DATA ZILCH /0.0D+00/
      DATA ONEM /-1.0D+00/
      DATA SPCASE /'AAAA','BBBB'/
C
C First pick up T1 vector.
C
      TLIST = IAPRT1AA
C
      CALL GETPERT1(ICORE, MAXCOR, MXCOR, IUHF, IRREPX, TLIST, 
     &              IOFFSET, IOFFT1)
C
C Spin cases AAAA and BBBB
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
               ICOUNT  = ICOUNT + VRT(IRREPF, ISPIN)*VRT(IRREPE, ISPIN)
               ICOUNT2 = ICOUNT2 + POP(IRREPF,ISPIN)*POP(IRREPE, ISPIN)
 1001       CONTINUE
            ABFULL(IRREP) = ICOUNT
            MNFULL(IRREP) = ICOUNT2
 1000    CONTINUE
C
C Loop over irreps : In the first block of code this corresponds to BM,
C while it is JE in the second block.
C First do the  G(MB,EJ) = SUM T(J,F)*Hbar(BM,EF) contribution.
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
            LSTWINT = 26 + ISPIN
            IF(IUHF .EQ. 0) LSTWINT = 30
C     
            DSZTAR = IRPDPD(IRREPGBM, 8 + ISPIN)
            DISTAR = IRPDPD(IRREPGEJ, 8 + ISPIN)
            DSZW   = IRPDPD(IRREPWEF, ISYTYP(1, LSTWINT))
            DISW   = IRPDPD(IRREPWBM, ISYTYP(2, LSTWINT))
            DSZEXP = ABFULL(IRREPWEF)
C     
            I000 = 1
            I010 = I000 + IINTFP*DISTAR*DSZTAR
            I020 = I010 + DSZEXP*DISW*IINTFP
            I030 = I020 + DSZW*IINTFP
C     
            CALL IZERO(ICORE, IINTFP*DISTAR*DSZTAR)
C     
C Read Hbar(BM,FE) F < E or [Hbar(bm,fe) f < e].
C     
            IF(IUHF .NE. 0) THEN
               CALL GETTRN(ICORE(I010), ICORE(I020), DSZW, DISW, 2, 
     &                     IRREPWBM, LSTWINT)
            ELSE
               CALL GETTRN(ICORE(I010), ICORE(I020), DSZEXP, DISW, 2,
     &                     IRREPWBM, LSTWINT)
               CALL ASSYM2(IRREPWEF, VRT(1,1), DISW, ICORE(I010))
            ENDIF
C     
C Change meanings of DSZ and DIS
C     
            DSZW = DISW
            DISW = DSZEXP
C     
C Expand to I(BM,FE) [I(bm,fe)] and transpose ket indices.
C     
            CALL SYMEXP(IRREPWEF, VRT(1,ISPIN), DSZW, ICORE(I010))
C     
C Now loop over irreps of F and form the matrix product.
C     Z(BME;J) = I(BME;F)*T(F,J)  [ISPIN = 1]
C     Z(bme;j) = I(bme;f)*T(f,j)  [ISPIN = 2]
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
     &                           POP(IRREP, ISPIN)*IINTFP
 5000       CONTINUE
C
            DO 4101 IRREPF = 1, NIRREP
C     
C Figure out the dimensions of the matrices involved in the multiplication.
C     
               IRREPE = DIRPRD(IRREPF, IRREPWEF)
               IRREPJ = DIRPRD(IRREPF, IRREPX)
C     
               NROWI = DSZW*VRT(IRREPE, ISPIN)
               NROWZ = NROWI
               NCOLZ = POP(IRREPJ, ISPIN)
               NCOLI = VRT(IRREPF, ISPIN)
               NROWT = VRT(IRREPF, ISPIN)
               NCOLT = POP(IRREPJ, ISPIN)
C     
               IOFFT = IOFFT1(IRREPF, ISPIN)
               IOFFZ = IZOFF1(IRREPJ)
C     
               IF (MIN (NROWI, NCOLT, NCOLI) .GT. 0) THEN
                  CALL XGEMM('N', 'N', NROWZ, NCOLZ, NCOLI, ALPHA,
     &                        ICORE(IOFFI), NROWI, ICORE(IOFFT),
     &                        NROWT, BETA, ICORE(IOFFZ), NROWZ)
               ENDIF
C     
               IOFFI = IOFFI + NROWI*NCOLI*IINTFP
C
 4101       CONTINUE
C     
C Now reorder G(BM,EJ) to G(JE,MB).
C     
            I020 = I010 + IINTFP*MAX(DSZTAR, DISTAR)
            I030 = I020 + IINTFP*MAX(DSZTAR, DISTAR)
            I040 = I030 + IINTFP*MAX(DSZTAR, DISTAR)
C     
            CALL SYMTR1(IRREPGEJ, VRT(1,ISPIN), POP(1,ISPIN), DSZTAR,
     &                  ICORE(I000), ICORE(I020), ICORE(I030),
     &                  ICORE(I040))
            CALL SYMTR3(IRREPGBM, VRT(1,ISPIN), POP(1,ISPIN), DSZTAR,
     &                  DISTAR, ICORE(I000), ICORE(I020), ICORE(I030),
     &                  ICORE(I040))
C
C Allocate another array of size DISTAR*DSZTAR. This will handle
C general rectangular matrices.
C            
            I020 = I010 + IINTFP*DISTAR*DSZTAR
C
            CALL TRANSP(ICORE(I000), ICORE(I010), DISTAR, DSZTAR)
c YAU : old
c           CALL ICOPY(IINTFP*DISTAR*DSZTAR,ICORE(I010),1,ICORE(I000),1)
c YAU : new
            CALL DCOPY(DISTAR*DSZTAR,ICORE(I010),1,ICORE(I000),1)
c YAU : end
C     
C Now do the - SUM T(N,B)*Hbar(MN,EJ) contribution. This is ordered as
C (BM,JE).
C    
            ALPHA   = ONEM*FACTOR
            BETA    = ONE
            LSTWINT = 6 + ISPIN
C
            IRREPWEJ = IRREPGEJ
            IRREPWMN = IRREPWEJ
C
            DSZW = IRPDPD(IRREPWMN, ISYTYP(1, LSTWINT))
            DISW = IRPDPD(IRREPWEJ, ISYTYP(2, LSTWINT))
C
            DSZEXP = MNFULL(IRREPWMN)
C     
            I020 = I010 + DSZEXP*DISW*IINTFP
            I030 = I020 + DSZW*IINTFP
C     
C Read Hbar(JE,NM) N < M or [Hbar(je,nm) n < m]
C     
            CALL GETTRN(ICORE(I010), ICORE(I020), DSZW, DISW, 2, 
     &                  IRREPWEJ, LSTWINT)
C     
C Change meaning of DSZW and DISW.
C     
            DSZW = DISW
            DISW = DSZEXP
C     
C Expand to I(JE;NM) [I(je;nm)] and then negate, giving
C I(JE,MN) [I(je,mn)].
C   
            CALL SYMEXP(IRREPWMN, POP(1,ISPIN), DSZW, ICORE(I010))
C     
C Now loop over irreps of N and form the matrix product.
C Z(JEM;B) = I(JEM;N)*T(B,N) [ISPIN = 1]
C Z(jem;b) = I(jem;n)*T(b,n) [ISPIN = 2]
C First calculate Z matrix offsets for each irrep of B.
C     
            IZOFF2(1) = I000
C
            DO 5100 IRREPZB = 2, NIRREP
               IRREP = IRREPZB - 1
               IRREPZM = DIRPRD(IRREP, IRREPGBM)
               IZOFF2(IRREPZB) = IZOFF2(IRREPZB - 1) + 
     &                           DSZW*POP(IRREPZM, ISPIN)*
     &                           VRT(IRREP, ISPIN)*IINTFP
 5100       CONTINUE
C     
            IOFFI = I010
C     
            DO 5101 IRREPN = 1, NIRREP
C     
C Figure out the dimensions of the matrices involved in the multiplication.
C     
               IRREPM = DIRPRD(IRREPN, IRREPWMN)
               IRREPB = DIRPRD(IRREPN, IRREPX)
C     
               NROWI = DSZW*POP(IRREPM, ISPIN)
               NCOLI = POP(IRREPN, ISPIN)
               NROWT = VRT(IRREPB, ISPIN)
               NCOLT = POP(IRREPN, ISPIN)
C
               IOFFT = IOFFT1(IRREPB, ISPIN)
               IOFFZ = IZOFF2(IRREPB)
C     
               IF(MIN(NROWI, NCOLT, NCOLI) .GT. 0) THEN
                  CALL XGEMM('N', 'T', NROWI, NROWT, NCOLI, ALPHA,
     &                        ICORE(IOFFI), NROWI, ICORE(IOFFT), 
     &                        NROWT, BETA, ICORE(IOFFZ), NROWI)
               ENDIF
C     
               IOFFI = IOFFI + NROWI*NCOLI*IINTFP
C     
 5101       CONTINUE
C     
C Now switch ordering (JE,MB) to (EJ,BM) and write these quantities to disk for
C each irrep.
C     
            I020 = I010 + IINTFP*MAX(DSZTAR, DISTAR)
            I030 = I020 + IINTFP*MAX(DSZTAR, DISTAR)
            I040 = I030 + IINTFP*MAX(DSZTAR, DISTAR)
C     
C Note the changes in dimensions in these calls
C
            CALL SYMTR1(IRREPGBM, POP(1,ISPIN), VRT(1,ISPIN), DISTAR,
     &                  ICORE(I000), ICORE(I020), ICORE(I030), 
     &                  ICORE(I040))
            CALL SYMTR3(IRREPGEJ, POP(1,ISPIN), VRT(1,ISPIN), DISTAR,
     &                  DSZTAR, ICORE(I000), ICORE(I020), ICORE(I030), 
     &                  ICORE(I040))
C     
            CALL PUTLST(ICORE(I000), 1, DSZTAR, 2, IRREPGBM, LSTTAR)
C
 20      CONTINUE
C     
C Now switch ordering to E,M-B,J
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
         IF (I030 .GT. MXCOR) CALL INSMEM('MKT1GMBEJA1', I030, MXCOR)
C
         CALL GETALL(ICORE(I010), INSIZ, IRREPX, LSTTAR)
         CALL ALTSYMPCK0(ICORE(I010), ICORE(I000), INSIZ, TARSIZ,
     &                   ICORE(I020), IRREPX, SPCASE(ISPIN))
C     
         CALL PUTALL(ICORE(I000), TARSIZ, IRREPX, LSTTAR)
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
