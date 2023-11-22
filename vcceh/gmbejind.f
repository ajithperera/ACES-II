C
      SUBROUTINE GMBEJIND(ICORE, MAXCOR, IRREPX, IOFFSET, ISPIN, IUHF)
C
C Driver for ring contributions to doubles equation
C
C   ab                                ae
C  Z (ALP,BET)   = P_(ij)P_(ab){SUM  t (BET)*G(ALP)}
C   ij                          m,e   im      (mb,ej) 
C
C                                ab
C                = P_(ij)P_(ab) S
C                                ij
C where G(ALP) is a two-particle intermediate, and spin
C orbital equation for G intermediate is as follows.
C             bf                    f                     b
C  (1/2) SUM t * Hbar(mn,ef) + SUM t * Hbar(mb,ef) - SUM t * Hbar(mn,ej) --> (1)
C        n,f  jn                f   j                 n   n
C
C Here and in the following equations alpha and beta is impilicit.
C
C The T2 contribution to G(mb,ej) is computed by t2qgmbej.f and 
C the T1 contributions to G(mb,ej) is computed by qt1ring.f. For the
C quadratic term there is no contribution from bare Hbar(mb,ej) 
C integrals. Actually, G(mb,ej) intermediate has the following
C structure. 
C 
C   AAAA - Negative of the equation (1)
C   BBBB - Negative of the equation (1)
C   ABAB - Identical to the equation (1)
C   BABA - Identical to the equation (1)
C   BAAB - Negative of the equation (1)
C   ABBA - Negative of the equation (1) 
C  
C The increments computed in various calls to MKDBLRNG1 and MKDBLRNG2 
C and added together and negated in routine SUMRNG. In addition to
C G(mb,ej) contributions this routine antisymmetrizes the contributions
C calculated in routine T1T1IND1 for AAAA and BBBB spin cases and add
C them to the final doubles list. The ABAB contributions from T1T1IND
C routine is appropriately rearrange and also add them to the final
C doubles list. 
C
C For various spin cases, Z is defined by
C
C   AB    AB   BA   BA   BA
C  Z   = S  + S  - S  + S     (AAAA)
C   IJ    IJ   JI   IJ   JI
C
C   Ab    Ab   bA   Ab   bA
C  Z   = S  + S  + S  + S     (ABAB)
C   Ij    Ij   jI   Ij   jI
C
C   aB    aB   Ba   aB   Ba
C  Z   = S  + S  + S  + S     (ABBA)
C   Ij    Ij   jI   Ij   jI
C
C The BBBB, BABA and BAAB spin cases can be obtained from these equations
C by changing all alpha labels to beta and all beta labels to alpha.
C
C The spin cases for S are given by,
C
C   AB             AE           Ae                
C  S  = S(AAAA) = t  G(MBEJ) + t  G(mBeJ)
C   IJ             IM           Im                   
C
C   ab             ae           aE                
C  S  = S(BBBB) = t  G(mbej) + t  G(MbEj)
C   ij             im           iM                   
C
C   Ab             AE           Ae                
C  S  = S(ABAB) = t  G(MbEj) + t  G(mbej)
C   Ij             IM           Im                   
C
C   aB             ae           aE
C  S  = S(BABA) = t  G(mBeJ) + t  G(MBEJ)
C   iJ             im           iM
C
C   Ab               eA
C  S  = S(BAAB) = - t  G(MbeJ)
C   iJ               iM      
C
C   aB               Ea
C  S  = S(ABBA) = - t  G(mBEj)
C   Ij               Im
C
C (All quantities are summed over m and e)
C
C Having the S quantities, the ring contribution to the doubles amplitudes
C (the Q(ij,ab)) can be computed from,
C
C                 AB             AB
C Spin case AA : Z   = P(IJ|AB) S
C                 IJ             IJ
C
C                 ab             ab
C Spin case BB : Z   = P(ij|ab) S
C                 ij             ij
C
C                 Ab    Ab    bA    bA    Ab
C Spin case AB : Z   = S   + S   - S   - S  
C                 Ij    Ij    jI    Ij    jI
C
C OR Z(AB) = S(ABAB) + S(BABA) - S(ABBA) - S(BAAB)
C
C
C The comments preceding each call to MKDBLRNG1 and MKDBLRNG2 include
C only the leading term of the G(mbej) intermediate which is used. It
C should be understood that the actual G is being used.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      INTEGER DIRPRD
      CHARACTER*4 SPCASE
      DIMENSION ICORE(MAXCOR)
C
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /FLAGS/ IFLAGS(100)
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
      NNP1O2(I)=I*(I+1)/2 
C
      NOCCA = NOCCO(1)
      NOCCB = NOCCO(2)
      NVRTA = NVRTO(1)
      NVRTB = NVRTO(2)
C
C For the quadratic term we update the disk.
C
      ITERM  = 1
      MAXSIZ = 0 
C     
      IF (ISPIN .LT. 3) THEN 
C
         IF(ISPIN .EQ. 1) SPCASE = 'AAAA'
         IF(ISPIN .EQ. 2) SPCASE = 'BBBB'
C     
C S(AAAA) and S(BBBB). Comments below refer to the AAAA case only.
C Flipping the sense of A and B will give the equations for the
C BBBB case. When the variable (ITERM = 1) we update the doubles
C increments already on the disk and when (ITERM = 0) we overwrite.
C
C             AB         AE
C Solve for  S   = SUM  T  * G(MB,EJ)  [First part of S(AAAA)]
C             IJ   M,E   IM
C
         LISTT = (IAPRT2AA1 - 1) + ISPIN
         LISTG = (INGMCAA -  1)  + ISPIN
C
         CALL MKDBLRNG1(ICORE, MAXCOR, LISTG, LISTT, 'TxG', MAXSIZ,
     &                 ITERM, 39 + ISPIN, IRREPX)
C
C            AB         Ae
C Solve for Q   = SUM  T  * G(mB,eJ)  [Second part of S(AAAA)]
C            IJ   m,e   Im
C
         LISTT = IAPRT2AB3 - ISPIN
         LISTG = INGMCBAAB - ISPIN
C     
         CALL MKDBLRNG1(ICORE, MAXCOR, LISTG, LISTT, 'TxG', MAXSIZ,
     &                  1, 39 + ISPIN, IRREPX)
C
C Resort and antisymmetrize T2 increment, and then augment
C the T2 amplitude list. Now the the total contribution is totally
C symmetric. So, I am free to use existing CC routines freely.
C
         NOCC  = NOCCO(ISPIN)
         NVRT  = NVRTO(ISPIN)
         NVSIZ = (NVRT*(NVRT - 1))/2
         NOSIZ = (NOCC*(NOCC - 1))/2
C
         NTOTSZ  = ISYMSZ(ISYTYP(1, 39 + ISPIN), ISYTYP(2, 39 + ISPIN))
         NTOTSZ2 = ISYMSZ(ISYTYP(1, 60 + ISPIN), ISYTYP(2, 60 + ISPIN))
         ISCRSZ  = NVSIZ + NOSIZ + NOCC*NVRT
C
         I000 = 1
         I010 = I000 + MAX(NTOTSZ, NTOTSZ2)*IINTFP
         I020 = I010 + MAX(NTOTSZ, NTOTSZ2)*IINTFP
         I030 = I020 + ISCRSZ
C
C Read the full list, and first symmetry pack and then anti-symmetrize the 
C (BJ,AI) quantity. The resulting quantity has the following structure
C (A<B,I<J).
C
         CALL GETALL(ICORE(I000), NTOTSZ, 1, 39 + ISPIN)
         CALL SST03I(ICORE(I000), ICORE(I010), NTOTSZ, NTOTSZ2,
     &               ICORE(I020), SPCASE)
C
C Update the appropriate lists   
C
         CALL GETALL(ICORE(I000), NTOTSZ2, 1, 60 + ISPIN)
C
         CALL SAXPY(NTOTSZ2, ONE, ICORE(I000), 1, ICORE(I010), 1)
         CALL PUTALL(ICORE(I010), NTOTSZ2, 1, 60 + ISPIN)
C
         IF (IFLAGS(1) .GE. 20) THEN
C
            NSIZE = NTOTSZ2
            CALL HEADER('Checksum @-GMBEJIND', 0, LUOUT)
C     
            WRITE(LUOUT, *) SPCASE, ' = ', SDOT(NSIZE, ICORE(I010),
     &                      1, ICORE(I010), 1)
         ENDIF
C
      ELSE IF (ISPIN .EQ. 3) THEN
C     
         IF (IUHF .NE. 0) THEN
C     
C UHF code. Spin case ABAB
C
            MAXSIZ = 0
C
C             Ab           AE
C Solve for  S   = - SUM  T  * G(Mb,Ej) [Part 1 of S(ABAB)]
C             Ij     M,E   IM
C
            LISTT = IAPRT2AA1
            LISTG = INGMCABAB
C
            CALL MKDBLRNG1(ICORE, MAXCOR, LISTG, LISTT, 'TxG', MAXSIZ,
     &                  ITERM, 42, IRREPX)
C
C             Ab         Ae
C Solve for  S   = SUM  T  * G(mb,ej) [Part 2 of S(ABAB)]
C             Ij   m,e   Im
C     
            LISTT = IAPRT2AB2
            LISTG = INGMCBB
C
            CALL MKDBLRNG1(ICORE, MAXCOR, LISTG, LISTT, 'TxG', MAXSIZ,
     &                     1, 42, IRREPX)
C
C Spin case BABA
C
C             aB         aE
C Solve for  S   = SUM  T  * G(MB,EJ) [part 2 of S(BABA)].
C             iJ   M,E   iM
C
            LISTT = IAPRT2AB1
            LISTG = INGMCAA
C
            CALL MKDBLRNG1(ICORE, MAXCOR, LISTG, LISTT, 'GxT', MAXSIZ,
     &                    ITERM, 42, IRREPX)
C
C             aB          ae
C Solve for  S   =  SUM  T  * G(Bm,Je) [part 1 of S(BABA)].
C             iJ    m,e   im
C
            LISTT = IAPRT2BB1
            LISTG = INGMCBABA
C
            CALL MKDBLRNG1(ICORE, MAXCOR, LISTG, LISTT, 'GxT', MAXSIZ,
     &                     1, 42, IRREPX)
C
C Spin case BAAB 
C
C               Ab         eA
C Solves for   S   = SUM  T  * G(Mb,Je)  [S(BAAB)]
C               iJ   M,e   iM
C
            LISTG = INGMCBAAB + IUHF
            LISTT = IAPRT2AB4
C
            CALL MKDBLRNG1(ICORE, MAXCOR, LISTG, LISTT, 'TxG', MAXSIZ,
     &                     1, 43, IRREPX)
C
C Spin case ABBA
C
C             aB         Ea
C Solves for S   = SUM  T  * G(mB,jE)  [S(ABBA)]
C             Ij   m,E   Im
C
            LISTT = IAPRT2AB3
            LISTG = INGMCBAAB
C            
            CALL MKDBLRNG1(ICORE, MAXCOR, LISTG, LISTT, 'GxT', MAXSIZ, 
     &                     1, 43, IRREPX)
C     
C Spin-adapted RHF code  
C     
         ELSE IF (IUHF .EQ. 0) THEN
C     
            NOCC = NOCCO(1)
            NVRT = NVRTO(1)
C     
            LISTG = INGMCABAB
            LISTT = IAPRT2AB2
C     
            CALL MKDBLRNG1(ICORE, MAXCOR, LISTG, LISTT, 'TxG', MAXSIZ,
     &                     ITERM, 42, IRREPX)
C
            LISTG = INGMCBAAB
            LISTT = IAPRT2AB4
C
            CALL MKDBLRNG2(ICORE, MAXCOR, LISTG, LISTT, 'TxG', MAXSIZ,
     &                     43, .TRUE., IRREPX)
         ENDIF
C     
C List 42 has the contribution (AI,bj) and list 43 has the contribution
C (Aj,bI). Add these together after appropriate resorting the list 43
C to get full S(ABAB) piece.
C The increments computed thus far (only Z(ABAB) here) are also negated in
C the routine SUMRNG. From now on I can use the same routines as in
C the original CC code since the now I have a totally symmetric 
C contribution.
C     
         NSZAB  = ISYMSZ(ISYTYP(1, 42), ISYTYP(2, 42))
         ISCRSZ = (NOCCA + NOCCB)*(NVRTA + NVRTB)
C
         I000 = 1
         I010 = I000 + NSZAB*IINTFP
         I020 = I010 + NSZAB*IINTFP
         I030 = I020 + ISCRSZ
C
         IF (I030 .GT. MAXCOR) CALL INSMEM('GMBEJIND', I030, MAXCOR)
C
         CALL SUMRNG(ICORE(I000), ICORE(I010), ICORE(I020), NSZAB,
     &               NSZAB, ISCRSZ)
C     
C Now convert ordering from (Ai,Bj) to (Ab,Ij) increments and augment the
C  T2 increment list.     
C
         I030 = I020 + ISCRSZ
         IF (I030 .GT. MAXCOR) CALL INSMEM('GMBEJIND', I030, MAXCOR)
C
         CALL SST02I(ICORE(I000), ICORE(I010), NSZAB, NSZAB, 
     &               ICORE(I020), 'AABB')
C     
C Now sum this with the existing (ABAB) increment and overwrite it.
C     
C
         IF (IFLAGS(1) .GE. 20) THEN
C
            NSIZE = NSZAB
C
            CALL HEADER('Checksum @-GMBEJIND ', 0, LUOUT)
C            
            WRITE(LUOUT, *) ' ABAB =  ', SDOT(NSIZE, ICORE(I010), 1,
     &                        ICORE(I010), 1)
C
         ENDIF
C
         CALL SUMSYM(ICORE(I010), ICORE(I000), NSZAB, 63)
C
      ENDIF
C
      RETURN
      END
