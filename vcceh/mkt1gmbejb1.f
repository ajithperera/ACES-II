C
      SUBROUTINE MKT1GMBEJB1(ICORE, MAXCOR, IUHF, IRREPX, IOFFSET)
C
C This subroutine computes two T1 contributions to the 
C W(MB,EJ) intermediate.  
C
C     W(mBEj) =  SUM T(j,f)*Hbar(mB,fE) - SUM T(N,B)*Hbar(mN,jE)
C     W(MbeJ) =  SUM T(J,F)*Hbar(Mb,Fe) - SUM T(n,b)*Hbar(Mn,Je) (UHF only)
C
C In the CC code F(EA) contibution from [SUM M F T(M,F)*Hbae(Ma,Fe)] was 
C also evaluated in this routine. That part is moved to a separate routine
C for convenience. The sign written here is opposite to the actual sign
C given in the JCP, 94, 4334, 1990. Here we compute the negtive of
C the actual expression given there. The reason will be clear when
C we evaluate the contribution of these ring terms to the doubles expression.
C
C This routine uses in-core algorithms.
C
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ONEM,ZILCH,ALPHA,BETA,FACTOR,SDOT
      CHARACTER*4 SSTSPN(2)
      CHARACTER*4 SPCASE(2)
      LOGICAL RHF
C
      DIMENSION ICORE(MAXCOR),IOFFT1(8,2),SIZEVV(2),IZOFF1(8),
     &          IZOFF2(8)
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
      DATA SSTSPN /'ABAB','BABA'/
      DATA SPCASE /'BAAB','ABBA'/
C
C First pick up T1 vector
C
      TLIST = IAPRT1AA
C
      CALL GETPERT1(ICORE, MAXCOR, MXCOR, IUHF, IRREPX, TLIST, 
     &              IOFFSET, IOFFT1)
C
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
            IRREPWBM = IRREPBM
            IRREPWEF = IRREPWBM
            IRREPGBM = IRREPWBM
            IRREPGEJ = DIRPRD(IRREPGBM, IRREPX)
C     
C Compute dimension of target matrix.
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
C This product is initially packed in (Bm,Ej) [(Mb,eJ)].
C     
            ALPHA = ONE*FACTOR
            BETA = ONE
C
            LISTW1 = 31 - ISPIN
C
            DSZW  = IRPDPD(IRREPWEF, ISYTYP(1, LISTW1))
            DISW  = IRPDPD(IRREPWBM, ISYTYP(2, LISTW1))
C
            I000  = 1
            I010  = I000 + IINTFP*DISTAR*DSZTAR
            I020  = I010 + IINTFP*DISW*DSZW
            I030  = I020 + IINTFP*DSZW
            CALL IZERO(ICORE, IINTFP*DISTAR*DSZTAR)
C     
C Pick up the integrals an hold them in a matrix I(Bm,Ef) [I(Mb,Fe)].
C     
            CALL GETTRN(ICORE(I010), ICORE(I020), DSZW, DISW, 2, 
     &                  IRREPWBM, LISTW1)
C     
C Change meanings of DSZ and DIS and transpose ket indices of the 
C integrals if ISPIN = 2. This gives
C   I(Bm,Ef) [ISPIN = 1]
C   I(Mb,eF) [ISPIN = 2]
C     
            ITMP = DISW
            DISW = DSZW
            DSZW = ITMP
C
            IF (ISPIN .EQ. 2) THEN
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
C  Z(BmE,j) = I(BmE,f)*T(f,j)
C  Z(Mbe,J) = I(Mbe,F)*T(F,J)  
C     
            IOFFI = I010
C
C Set the offsets for the target array
C
            IZOFF1(1) = I000
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
               IOFFZ = IZOFF1(IRREPJ)
C
               IF (MIN(NROWI, NCOLT, NCOLI) .GT. 0) THEN
C
                  CALL XGEMM('N', 'N', NROWZ, NCOLZ, NCOLI, ALPHA, 
     &                        ICORE(IOFFI), NROWI, ICORE(IOFFT),
     &                        NROWT, BETA, ICORE(IOFFZ), NROWZ)
               ENDIF
C
               IOFFI = IOFFI + NROWI*NCOLI*IINTFP
C
 100        CONTINUE
C
C The code above produces mbej intermediates in the order.
C   Bm,Ej (ISPIN = 1)
C   Mb,eJ (ISPIN = 2)
C
C Now do the second contraction
C   Z(mB,Ej) = Z(mB,Ej) - SUM T(N,B)*Hbar(mN,jE) (ISPIN = 1)
C   Z(Mb,eJ) = Z(Mb,eJ) - SUM T(n,b)*Hbar(Mn,Je) (ISPIN = 2, or RHF)
C     
C We need to reorder what we have so that it will match up and 
C can be calculated by XGEMM.
C     
            IF (ISPIN .EQ. 1) THEN
C
               I020 = I010 + IINTFP*DSZTAR*DISTAR
C
               CALL TRANSP(ICORE(I000), ICORE(I010), DISTAR, DSZTAR)
c YAU : old
c              CALL ICOPY(IINTFP*DISTAR*DSZTAR,ICORE(I010),1,ICORE(I000),1)
c YAU : new
               CALL DCOPY(DISTAR*DSZTAR,ICORE(I010),1,ICORE(I000),1)
c YAU : end
C
               I020 = I010 + IINTFP*MAX(DSZTAR, DISTAR)
               I030 = I020 + IINTFP*MAX(DSZTAR, DISTAR)
               I040 = I030 + IINTFP*MAX(DSZTAR, DISTAR)
C
               CALL SYMTR1(IRREPGBM, VRT(1,1), POP(1,2), DISTAR,
     &                     ICORE(I000), ICORE(I010), ICORE(I020),
     &                     ICORE(I030))
            ELSE
C
               I020 = I010 + IINTFP*MAX(DSZTAR, DISTAR)
               I030 = I020 + IINTFP*MAX(DSZTAR, DISTAR)
               I040 = I030 + IINTFP*MAX(DSZTAR, DISTAR)
C
               CALL SYMTR1(IRREPGEJ, VRT(1,2), POP(1,1), DSZTAR,
     &                     ICORE(I000), ICORE(I010), ICORE(I020),
     &                     ICORE(I030))
C
               I020 = I010 + IINTFP*DSZTAR*DISTAR
C
               CALL TRANSP(ICORE(I000), ICORE(I010), DISTAR, DSZTAR)
c YAU : old
c              CALL ICOPY(IINTFP*DISTAR*DSZTAR,ICORE(I010),1,ICORE(I000),1)
c YAU : new
               CALL DCOPY(DISTAR*DSZTAR,ICORE(I010),1,ICORE(I000),1)
c YAU : end
            ENDIF
C
            LISTW2 = 8 + ISPIN
            IF (RHF) LISTW2 = 10
C
            ALPHA = ONEM*FACTOR
            BETA  = ONE
C
            IRREPWEJ = IRREPGEJ
            IRREPWMN = IRREPWEJ
C
            DSZW = IRPDPD(IRREPWMN, ISYTYP(1, LISTW2))
            DISW = IRPDPD(IRREPWEJ, ISYTYP(2, LISTW2))
C
            I020 = I010 + IINTFP*DSZW*DISW
            I030 = I020 + IINTFP*DSZW
C     
C Get the integrals and place into a matrix.
C  I(Ej,Nm)   [ISPIN = 1]
C  I(Je,Mn)   [ISPIN = 2 OR RHF]
C     
            CALL GETTRN(ICORE(I010), ICORE(I020), DSZW, DISW, 2, 
     &                  IRREPWEJ, LISTW2)
C
            IF (IUHF .EQ. 0) THEN
C
               I030 = I020 + MAX(DISW,DSZW)*IINTFP
               I040 = I030 + MAX(DISW,DSZW)*IINTFP
               I050 = I040 + MAX(DISW,DSZW)*IINTFP

               CALL SYMTR3(IRREPWEJ, POP(1, 1), VRT(1, 2), DISW, DSZW,
     &                     ICORE(I010), ICORE(I020), ICORE(I030),
     &                     ICORE(I040))
            ENDIF
C            
            ITMP = DISW
            DISW = DSZW
            DSZW = ITMP
C     
C Transpose ket indices for ISPIN = 1. This gives
C  I(Ej,mN) [ISPIN = 1]
C  I(Je,Mn) [ISPIN = 2 OR RHF]
C     
            IF(ISPIN .EQ. 1 .AND. .NOT. RHF) THEN
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
C Now do matrix multiplication as above, with the matrices
C  Z(Ej,mB) = I(Ejm,N)*T(B,N)
C  Z(Je,Mb) = I(JeM,n)*T(b,n)
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
     &                           DSZW*POP(IRREPZM, 3 - ISPIN)*
     &                           VRT(IRREP, ISPIN)*IINTFP
 5100       CONTINUE
C
            DO 200 IRREPN = 1, NIRREP
C
               IRREPM = DIRPRD(IRREPN, IRREPWMN)
               IRREPB = DIRPRD(IRREPN, IRREPX)
C
               NROWI = DSZW*POP(IRREPM, 3 - ISPIN)
               NCOLI = POP(IRREPN, ISPIN)
               NROWT = VRT(IRREPB, ISPIN)
               NCOLT = POP(IRREPN, ISPIN)
C
               NROWZ = NROWI
               NCOLZ = NROWT
C
               IOFFT = IOFFT1(IRREPB, ISPIN)
               IOFFZ = IZOFF2(IRREPB)
C
               IF (MIN(NROWI, NCOLT, NCOLI) .GT. 0) THEN
C
                  CALL XGEMM('N', 'T', NROWZ, NCOLZ, NCOLI, ALPHA,
     &                        ICORE(IOFFI), NROWI, ICORE(IOFFT), NROWT,
     &                        BETA, ICORE(IOFFZ), NROWZ)
               ENDIF
C
               IOFFI = IOFFI + NROWI*NCOLI*IINTFP
C
 200        CONTINUE
C     
C We now have
C  Z(Ej,mB)  [ISPIN = 1]
C  Z(Je,Mb)  [ISPIN = 2 OR RHF]
C     
C We want to resort these first to
C     
C     Z(Ej,Bm)  [ISPIN = 1]
C     Z(eJ,bM)  [ISPIN = 2]
C     
C Now transpose bra indices for target if ISPIN = 2, and ket indices in 
C either case. In the following calls DSZTAR should be DISTAR, but 
C it works because they are identical in size (serious assumptions
C made by original authors)
C    
            IF (ISPIN .EQ. 2) THEN
C
               I020 = I010 + IINTFP*MAX(DSZTAR, DISTAR)
               I030 = I020 + IINTFP*MAX(DSZTAR, DISTAR)
               I040 = I030 + IINTFP*MAX(DSZTAR, DISTAR)
C
               CALL SYMTR3(IRREPGEJ, POP(1,1), VRT(1,2), DISTAR, 
     &                     DSZTAR, ICORE(I000), ICORE(I010),
     &                     ICORE(I020), ICORE(I030))
            ENDIF
C
            I020 = I010 + IINTFP*MAX(DSZTAR, DISTAR)
            I030 = I020 + IINTFP*MAX(DSZTAR, DISTAR)
            I040 = I030 + IINTFP*MAX(DSZTAR, DISTAR)
C
            CALL SYMTR1(IRREPGBM, POP(1,3-ISPIN), VRT(1,ISPIN),
     &                  DISTAR, ICORE(I000), ICORE(I010), 
     &                  ICORE(I020), ICORE(I030))
C     
C Now dump them in to disk for each irrep
C     
            CALL PUTLST(ICORE(I000), 1, DSZTAR, 1, IRREPGBM, LSTTAR)
C
 20      CONTINUE
C     
C Now swith the ordering to and put back in the disk.
C Z(Ej, Bm) -> Z(Em,Bj) [ISPIN = 1]
C Z(eJ, bM) -> Z(eM,bJ) [ISPIN = 2 or RHF]
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
         IF(I030 .GT. MXCOR) CALL INSMEM('MKT1GMBEJB1', I030, MXCOR)
C         
         CALL GETALL(ICORE(I010), INSIZ, IRREPX, LSTTAR)
         CALL ALTSYMPCK0(ICORE(I010), ICORE(I000), INSIZ, TARSIZ,
     &                   ICORE(I020), IRREPX, SSTSPN(ISPIN))
C
         CALL PUTALL(ICORE(I000), TARSIZ, IRREPX, LSTTAR)
C
            IF (IFLAGS(1) .GE. 20) THEN

               NSIZE = TARSIZ
     
               CALL HEADER('Checksum @-QT1RING2', 0, LUOUT)
     
               WRITE(LUOUT, *) SPCASE(ISPIN), ' = ', SDOT(NSIZE,
     &                         ICORE(I000), 1, ICORE(I000), 1)
            ENDIF
C
 10   CONTINUE
C
      RETURN
      END
