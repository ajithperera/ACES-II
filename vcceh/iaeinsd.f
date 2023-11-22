C
      SUBROUTINE IAEINSD(ICORE, MAXCOR, IRREPX, IOFFSET, IUHF)
C
C This routine calculates the term
C
C  Z(IJ,AB) = P_(AB) SUM E T(IJ,AE) I(E,B)
C
C to the doubles equation using symmetry packed arrays.
C Spin-integrated formulas for this contribution 
C
C In UHF
C
C  SUM E T(IJ,AE) I(EB) - SUM E T(IJ,EB) I(EA)  [AAAA]
C
C  SUM e T(ij,ae) I(eb) - SUM e T(ij,eb) I(ea)  [BBBB]
C
C  SUM e T(Ij,Ae) I(eb) - SUM E T(Ij,Eb) I(EA)  [ABAB]
C
C In RHF 
C
C  SUM e T(Ij,Ae) I(eb) - SUM E  T(Ij,Eb) I(EA) [ABAB]
C
C In addition to those terms this subroutine also calculate 
C
C  Z(IA) = SUM E T(I,E) I(EA) [AA]
C
C  Z(ia) = SUM a T(i,e) I(ea) [BB]
C
C contribution to the singles equation.
C
C The term, - 1/2 SUM M Hbar(E,M) T(B,M) is added
C to the I(B,E) intermediate before contraction with T2
C amplitudes is carried out (here alpha and beta is implicit).
C
C Coded bu JG June/90 and July/90 and modified for quadratic term
C by Ajith 07/94
C
C     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIRPRD,DISSYT,DISSYZ,POP,VRT
      DIMENSION ICORE(MAXCOR), IIOFF1(8), IIOFF2(8)
      CHARACTER*6 SPCASE(2)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWO
      COMMON /SYM/ POP(8,2),VRT(8,2),NTAA,NTBB,
     &             NF1AA,NF1BB,NF2AA,NF2BB
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /FILES/ LUOUT, MOINTS
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
      DATA ZLICH, ONE, HALFM /0.00D+00, 1.0D0, -0.5D0/
      DATA SPCASE /'AA =  ', 'BB =  '/
C     
C Calculate size of I(A,E) array and load these arrays
C in to the memory. 
C
      NFVVAA = IRPDPD(IRREPX, 19)
      NFVVBB = IRPDPD(IRREPX, 20)
C
      I0AA  = MAXCOR + 1 - NFVVAA*IINTFP
      MXCOR = MAXCOR - NFVVAA*IINTFP
C
      IF (IUHF .EQ. 0) THEN
         I0BB = I0AA
      ELSE
         I0BB  = I0AA - NFVVBB*IINTFP
         MXCOR = MXCOR - NFVVBB*IINTFP
      ENDIF
C
      CALL GETLST(ICORE(I0AA), 1, 1, 1, 1, INTIAE)
C
      IF (IFLAGS(1) .GE. 40) THEN
         CALL HEADER('Quadratic I(AE) Int. @IAEINSD', 0, LUOUT)
      ENDIF
C
      IF (IUHF .EQ. 1) THEN
         CALL GETLST(ICORE(I0BB), 1, 1, 1, 2, INTIAE)
C
         IF (IFLAGS(1) .GE. 40) THEN
            CALL HEADER('Quadratic I(AE) Int. @IAEINSD', 0, LUOUT)
            CALL TAB(LUOUT, ICORE(I0BB), NFVVBB, 1, NFVVBB, 1)
         ENDIF
C
      ENDIF
C
      NFVOAAX = IRPDPD(IRREPX, 9)
      NFVOBBX = IRPDPD(IRREPX, 10)
      NFVOAA1 = NTAA
      NFVOBB1 = NTBB
C
      I0TA = I0BB - NFVOAAX*IINTFP
      I0ZA = I0TA - NFVOAA1*IINTFP
C     
      MXCOR = MXCOR - (NFVOAAX + NFVOAA1)*IINTFP
C
C Get the T1(BET) amplitudes and zero out the memory locations for 
C the product
C
      CALL GETLST(ICORE(I0TA), 1, 1, 1, 1, IBPRT1AA)
      CALL IZERO(ICORE(I0ZA), NFVOAA1*IINTFP)
C
      IF (IUHF .EQ. 0) THEN
         I0TB = I0TA
         I0ZB = I0ZA
      ELSE
         I0TB = I0ZA - NFVOBBX*IINTFP
         I0ZB = I0TB - NFVOBB1*IINTFP
C     
         MXCOR = MXCOR - (NFVOBBX + NFVOBB1)*IINTFP
C     
         CALL GETLST(ICORE(I0TB), 1, 1, 1, 2, IBPRT1AA)
         CALL IZERO(ICORE(I0ZB), NFVOBB1*IINTFP)
      ENDIF
C     
C Now we are in a position to evaluate I(A,E) contribution to the
C singles and doubles equation. First do the singles. 
C 
      DO 300 ISPIN = 1, IUHF + 1
C     
         IF (ISPIN .EQ. 1) THEN
            IIOFF1(1) = I0AA
            DO 5000 IRREPAO = 2, NIRREP
               IRREP = IRREPAO - 1
               IRREPE = DIRPRD(IRREP, IRREPX)
               IIOFF1(IRREPAO) = IIOFF1(IRREPAO - 1) + VRT(IRREP, 1)*
     &                          VRT(IRREPE, 1)*IINTFP
 5000       CONTINUE
            IOFFT = I0TA
            IOFFZ = I0ZA
            I0Z   = I0ZA
         ELSE
            IIOFF1(1) = I0BB
            DO 5100 IRREPAO = 2, NIRREP
               IRREP = IRREPAO - 1
               IRREPE = DIRPRD(IRREP, IRREPX)
               IIOFF1(IRREPAO) = IIOFF1(IRREPAO - 1) + VRT(IRREP, 2)*
     &                          VRT(IRREPE, 2)*IINTFP
 5100       CONTINUE
C
            IOFFT = I0TB
            IOFFZ = I0ZB
            I0Z   = I0ZB
C
         ENDIF
C     
         DO 250 IRREPI = 1, NIRREP
C
            IRREPTI = IRREPI
            IRREPTE = DIRPRD(IRREPTI, IRREPX)
            IRREPIE = IRREPTE
            IRREPIA = DIRPRD(IRREPIE, IRREPX)
            IRREPQI = IRREPTI
            IRREPQA = IRREPQI
C     
            NOCCTI = POP(IRREPTI, ISPIN)
            NVRTTE = VRT(IRREPTE, ISPIN)
            NVRTIE = VRT(IRREPIE, ISPIN)
            NVRTIA = VRT(IRREPIA, ISPIN)
C
            IOFFF = IIOFF1(IRREPIA)

            IF (NVRTTE .NE. 0) THEN
            CALL XGEMM('T', 'N', NVRTIA, NOCCTI, NVRTTE, ONE,
     &                  ICORE(IOFFF), NVRTIE, ICORE(IOFFT), NVRTTE,
     &                  ZLICH, ICORE(IOFFZ), NVRTIA)
            ELSE
            CALL ZERO(ICORE(IOFFZ),NVRTIA*NOCCTI)
            ENDIF 

            IOFFT = IOFFT + NVRTTE*NOCCTI*IINTFP
            IOFFZ = IOFFZ + NVRTIA*NOCCTI*IINTFP
C
 250     CONTINUE
C
C Write Q(AI) contribution to the singles list. Note that 
C here I can overwrite the CC lists since now I have a
C totally symmetric contribution.
C
         CALL PUTLST(ICORE(I0Z), 1, 1, 1, ISPIN + 2, 90)
C
         IF (IFLAGS(1) .GE. 20) THEN
C
            IF (ISPIN .EQ. 1) THEN
               NSIZE = NFVOAA1
            ELSE
               NSIZE = NFVOBB1
            ENDIF
C
            CALL HEADER('Checksum @-IAEINSD', 0, LUOUT)
            
            WRITE(LUOUT, *) SPCASE(ISPIN), SDOT(NSIZE, ICORE(I0Z),
     &                      1, ICORE(I0Z), 1)
         ENDIF
C            
 300  CONTINUE
C
C Before we calculate the I(B,E) contributions to the doubles equation we 
C need to add  - 1/2 SUM M Hbar(E,M)*T(B,M) to the I(B,E) intermediate.
C Note the difference from the CC equations. In CC equations F(E,M) 
C intermediate is used instead of Hbar(e,m) matrix elements.
C     
      IOFFEM = I0ZA
      IOFT   = I0TA 
C
C Overwrite the T1(beta) by T1(alpha) and get bare Hbar(E,M) and 
C overwrite the allocation for the Q(A,I).
C      
      CALL GETLST(ICORE(I0TA), 1, 1, 1, 1, IAPRT1AA)
      CALL GETLST(ICORE(I0ZA), 1, 1, 1, 1, 93)
C
      IIOFF2(1) = I0AA
C
      DO 5200 IRREPBO = 2, NIRREP
         IRREP = IRREPBO - 1 
         IRREPE = DIRPRD(IRREP, IRREPX)
         IIOFF2(IRREPBO) = IIOFF2(IRREPBO - 1) + VRT(IRREPE, 1)*
     &                     VRT(IRREP, 1)*IINTFP
 5200 CONTINUE
C
      DO 20 IRREPM = 1, NIRREP
C
         IRREPWM = IRREPM
         IRREPWE = IRREPWM
         IRREPTM = IRREPWM
         IRREPTB = DIRPRD(IRREPTM, IRREPX)
C
         NOCCM = POP(IRREPWM, 1)
         NVRTE = VRT(IRREPWE, 1)
         NVRTB = VRT(IRREPTB, 1)
C
         IOFFEB = IIOFF2(IRREPTB)
C
         CALL XGEMM('N', 'T', NVRTE, NVRTB, NOCCM, HALFM, 
     &               ICORE(IOFFEM), NVRTE, ICORE(IOFT), NVRTB, ONE,
     &               ICORE(IOFFEB), NVRTE)
C
         IOFFEM = IOFFEM + NOCCM*NVRTE*IINTFP
         IOFT   = IOFT   + NOCCM*NVRTB*IINTFP
C
 20   CONTINUE
C
      IF (IUHF .EQ. 1) THEN
C
         IOFFEM = I0ZB
         IOFT   = I0TB
C
         CALL GETLST(ICORE(I0TB), 1, 1, 1, 2, IAPRT1AA)
         CALL GETLST(ICORE(I0ZB), 1, 1, 1, 2, 93)
C
         IIOFF2(1) = I0BB
C
         DO 5300 IRREPBO = 2, NIRREP
            IRREP = IRREPBO - 1 
            IRREPE = DIRPRD(IRREP, IRREPX)
            IIOFF2(IRREPBO) = IIOFF2(IRREPBO - 1) + VRT(IRREPE, 2)*
     &                        VRT(IRREP, 2)*IINTFP
 5300    CONTINUE
C     
         DO 21 IRREPM = 1, NIRREP
C
            IRREPWM = IRREPM
            IRREPWE = IRREPWM
            IRREPTM = IRREPWM
            IRREPTB = DIRPRD(IRREPTM, IRREPX)
C
            NOCCM = POP(IRREPWM, 2)
            NVRTB = VRT(IRREPTB, 2)
            NVRTE = VRT(IRREPWE, 2)
C
            IOFFEB = IIOFF2(IRREPTB)
C
            CALL XGEMM('N', 'T', NVRTE, NVRTB, NOCCM, HALFM,
     &                  ICORE(IOFFEM), NVRTE, ICORE(IOFT), NVRTB, ONE,
     &                  ICORE(IOFFEB), NVRTE)
C
            IOFFEM = IOFFEM + NOCCM*NVRTE*IINTFP
            IOFT   = IOFT   + NOCCM*NVRTB*IINTFP
C
 21      CONTINUE
C
      ENDIF

C
C Lets do the contribution to the doubles equation 
C First AAAA and BBBB spin cases.
C     
      IF (IUHF .EQ. 1) THEN
C     
C These cases are only necessary in the UHF case and in the RHF
C AAAA amplitudes are calculated from ABAB.
C     
         DO 100 ISPIN = 1, 2
C     
            IF (ISPIN .EQ. 1) THEN 
               I000  = I0AA
               NFSIZ = NFVVAA
            ELSE
               I000 = I0BB
               NFSIZ = NFVVBB
            ENDIF
C
            LISTT = ISPIN + (IAPRT2AA2 - 1)
            LISTZ = ISPIN + 60
C     
            DO 50 IRREPIJ = 1, NIRREP
C    
               IRREPTIJ = IRREPIJ
               IRREPTAE = DIRPRD(IRREPTIJ, IRREPX)
               IRREPQIJ = IRREPTIJ
               IRREPQAB = IRREPQIJ
C     
               NVRTSQ1 = 0
               NVRTSQ2 = 0
C
               DO 45 IRREPAO = 1, NIRREP
                  IRREPBO = DIRPRD(IRREPAO, IRREPQAB)
                  IRREPTE = DIRPRD(IRREPAO, IRREPTAE)
                  NVRTSQ1 = NVRTSQ1 + VRT(IRREPAO, ISPIN)*
     &                     VRT(IRREPBO, ISPIN)
                  NVRTSQ2 = NVRTSQ2 + VRT(IRREPAO, ISPIN)*
     &                     VRT(IRREPTE, ISPIN)
 45            CONTINUE
C     
C Retrieve T2 amplitudes and calculate the corresponding contribution
C     
               DISSYT = IRPDPD(IRREPTAE, ISYTYP(1, LISTT))
               NUMSYT = IRPDPD(IRREPTIJ, ISYTYP(2, LISTT))
               DISSYZ = IRPDPD(IRREPQAB, ISYTYP(1, LISTZ))
               NUMSYZ = IRPDPD(IRREPQIJ, ISYTYP(2, LISTZ))
C
               I001 = 1
               I002 = I001 + IINTFP*MAX(NUMSYT*NVRTSQ1, NUMSYZ*NVRTSQ1,
     &                DISSYT*NUMSYT, NUMSYT*NVRTSQ2, DISSYZ*NUMSYZ)
               I003 = I002 + IINTFP*MAX(NUMSYZ*NVRTSQ1, NUMSYT*NVRTSQ1,
     &                DISSYT*NUMSYT, NUMSYT*NVRTSQ2, DISSYZ*NUMSYZ)
C
               IF (MIN (NUMSYT, NUMSYZ, DISSYT, DISSYZ) .NE. 0) THEN
                  I004 = I003 + IINTFP*MAX(NUMSYT, NUMSYZ)
                  IF (I004 .LT. MXCOR) THEN
C     
C In core version
C     
                     CALL MKDBLAA3(ICORE(I001), ICORE(I002), 
     &                             ICORE(I001), ICORE(I002),
     &                             ICORE(I000), VRT(1,ISPIN),
     &                             NVRTSQ1, NVRTSQ2, DISSYT, DISSYZ,
     &                             NUMSYT, NUMSYZ, NFSIZ, LISTT, LISTZ, 
     &                             IRREPTAE, IRREPTIJ, IRREPQAB,
     &                             IRREPQIJ, IRREPX, ISPIN,
     &                             ICORE(I003))
                  ELSE
                     CALL INSMEM('IAEINSD', I004, MXCOR)
                  ENDIF
C
               ENDIF
C
 50         CONTINUE
 100     CONTINUE
C
      ENDIF
       
C     
C AB spin case
C     
      LISTT = IAPRT2AB5
      LISTZ = 63
C     
C Loop over irreps
C     
      DO 200 IRREPIJ = 1, NIRREP
C
         IRREPTIJ = IRREPIJ
         IRREPTAE = DIRPRD(IRREPTIJ, IRREPX)
         IRREPQIJ = IRREPTIJ
         IRREPQAB = IRREPQIJ
C     
C Retrieve amplitudes calculate contribution to Z
C     
         DISSYT = IRPDPD(IRREPTAE, ISYTYP(1, LISTT))
         NUMSYT = IRPDPD(IRREPTIJ, ISYTYP(2, LISTT))
         DISSYZ = IRPDPD(IRREPQAB, ISYTYP(1, LISTZ))
         NUMSYZ = IRPDPD(IRREPQIJ, ISYTYP(2, LISTZ))
C         
         I001 = 1
         I002 = I001 + IINTFP*MAX(NUMSYT*DISSYT, NUMSYZ*DISSYZ)
         I003 = I002 + IINTFP*MAX(NUMSYZ*DISSYZ, NUMSYT*DISSYT)
C
         IF (MIN (NUMSYT, NUMSYZ, DISSYT, DISSYZ) .NE. 0) THEN
C
            I004 = I003 + IINTFP*MAX(NUMSYT, NUMSYZ, DISSYZ, 
     &             DISSYT)*3
C
            IF(I004 .LT. MXCOR) THEN
C     
C In core algorithm
C     
               CALL MKDBLAB3(ICORE(I001), ICORE(I002),
     &                       ICORE(I001), ICORE(I002),
     &                       ICORE(I0AA), ICORE(I0BB),
     &                       POP(1,1), POP(1,2), VRT(1,1), VRT(1,2),
     &                       DISSYT, DISSYZ, NUMSYT, NUMSYZ, NFVVAA, 
     &                       NFVVBB, LISTT, LISTZ, IRREPTAE, IRREPTIJ,
     &                       IRREPQAB, IRREPQIJ, IRREPX, IUHF,
     &                       ICORE(I003))
            ELSE
C     
C Out of core algorithm
C
               CALL INSMEM( 'IAEINSD', I004, MXCOR)
            ENDIF
C
         ENDIF
C
 200  CONTINUE
C
      RETURN
      END
