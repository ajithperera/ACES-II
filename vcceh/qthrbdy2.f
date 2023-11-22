C
      SUBROUTINE QTHRBDY2(ICORE, MAXCOR, IRREPX, IOFFSET, IUHF)
C
C Calculate three body contributions to the doubles equation.
C Spin orbital equations for this contribution is
C   
C  Z(IJ,AB) = - P_(AB) SUM (M,E;N,F) T(E,M)*T(A,N)*T(IJ,FB)
C             *Hbar(MN,EF)
C
C here alpha and beta is explicit. Spin integrated formulas for UHF
C and RHF cases are as given below.  
C 
C UHF
C  
C  Z(IJ,AB) = - P_(AB) {SUM(N,F) [ SUM (E,M) T(E,M)*Hbar(MN,EF) +
C               SUM(e,m) T(e,m)*Hbar(mN,eF)] T(N,A)*T(FB,IJ)}   [AAAA]
C  
C  Z(ij,ab) = - P_(ab) {SUM (n,f) [ SUM (e,m) T(e,m)*Hbar(mn,ef) + 
C               SUM(E,M) T(E,M)*Hbar(Mn,Ef)] T(n,a)*T(fb,ij}    [BBBB]
C 
C  Z(Ij,Ab) = - {{SUM(N,F) [ SUM (E,M) T(E,M)*Hbar(MN,EF) +
C               SUM(e,m) (T(e,m)T*Hbar(mN,eF))] T(N,A)*T(Fb,Ij)} +
C               {SUM (n,f) [ SUM(E,M) (T(E,M)*Hbar(Mn,Ef)) + 
C               SUM(e,m)*Hbar(mn,ef)] T(n,b)*T(Af,Ij)}}         [ABAB]
C
C RHF
C  
C  Z(Ij,Ab) = - {SUM(N,F) [ SUM (E,M) T(E,M)*Hbar(MN,EF) +
C               SUM(e,m) T(e,m)*Hbar(mN,eF)] T(N,A)*T(Fb,Ij)}  [ABAB]

C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C     
      INTEGER DISSYT, DISSYQ, DIRPRD, POP, VRT
      DIMENSION ICORE(MAXCOR)
C
      COMMON /SYM/ POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF1BB,
     &             NF2AA,NF2BB
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
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
      DATA ONEM /-1.0/
      DATA ZILCH /0.0/
C
      MXCOR = MAXCOR
C
C AAAA and BBBB spin cases (UHF only)
C
      IF (IUHF .NE. 0) THEN
C
         DO 100 ISPIN = 1, 2
C     
C First do the Z(F,N) = [T(E,M)*Hbar(MN,EF) + T(e,m)Hbar(Nm,Fe)] [ISPIN = 1]
C and Z(f,n) = [T(e,m)*Hbar(mn,ef) + T(E,M)Hbar(Mn,EF)] [ISPIN = 2]. 
C Allocate memory for Z(F,N) and leave the rest of the core to store
C T1 vector and Hbar(MN,EF) in the subroutine MKZNFTB. We need to consider
C only one DPD irrep; the irrepx.
C
            LISTT1A = IAPRT1AA
            LISTT1B = IBPRT1AA
            LISTH1  = 18 + ISPIN
            LISTH2  = 16 + ISPIN
C
            NFVOZ  = IRPDPD(IRREPX, 8 + ISPIN)
C
            I000 = 1
            I010 = MXCOR + 1 - NFVOZ*IINTFP
C
            MXCOR = MXCOR - NFVOZ*IINTFP
C
            CALL MAKZNFTB(ICORE(I010), ICORE(I000), MXCOR, NFVOZ, 
     &                    LISTT1A, LISTT1A, LISTH1, LISTH2, IUHF,
     &                    IRREPX, ISPIN)
C
C Now do the Y(F,A) = Z(F,N)*T(N,A) (ISPIN = 1) and Z(f,n)*T(n,a) (ISPIN = 2)
C The Z(F,N) is in array ICORE(I010). Allocate memory for the Y(F,A) and
C T(N,A). Now only three one body terms are held in memory simultaneously,
C which is not a problem.
C            
            NFVVY = IRPDPD(1, 18 + ISPIN)
            NFVOT = IRPDPD(IRREPX, 8 + ISPIN)
C
            I001 = I000 + NFVVY*IINTFP
            I002 = I001 + NFVOT*IINTFP
C         
            MXCOR = MXCOR - (NFVVY + NFVOT)*IINTFP
C
            IF (I002 .LT. MXCOR) THEN
C     
               CALL GETLST(ICORE(I001), 1, 1, 1, ISPIN, LISTT1B)
C
C Carry out the multiplication
C
               CALL MAKYFATB(ICORE(I010), ICORE(I000), ICORE(I001), 
     &                       ICORE(I002), NFVOZ, NFVVY, NFVOT, MXCOR,
     &                       IUHF, IRREPX, ISPIN)
            ELSE
               CALL INSMEM('QTHRBDY2', I002, MXCOR)
            ENDIF
C            
C Now do the final Q(AB,IJ) = T(FB,IJ)*Y(F,A) [ISPIN = 1] and
C Q(ab,ij) = T(fb,ij)*Y(f,a) [ISPIN = 2]. Allocate memory for
C the T2 vector and ICORE(I000) has the Y(F,A) array. Clean up
C the rest of the memory. In the next multiplication we have to
C keep three four-index arrays and one one-index array, resulting
C a requirement of an out-of core routine.
C     
            MXCOR = MAXCOR - NFVVY*IINTFP              
C
            LIST2 = 43 + ISPIN
            LISTQ = 60 + ISPIN
C            
            DO 10 IRREP = 1, NIRREP
C
C Now evrything is totally symmetric
C
               NVRT2SQ = 0
C
               DO 200 IRREPE = 1, NIRREP
                  IRREPF = DIRPRD(IRREPE, IRREP)
                  NVRT2SQ = NVRT2SQ + VRT(IRREPF, ISPIN)*VRT(IRREPE,
     &                      ISPIN)
 200           CONTINUE
C
               IRREPTBF = IRREP
               IRREPTIJ = IRREPTBF
               IRREPQIJ = IRREPTIJ
               IRREPQAB = IRREPQIJ
C
               DISSYT = IRPDPD(IRREPTBF, ISYTYP(1, LIST2))
               NUMSYT = IRPDPD(IRREPTIJ, ISYTYP(2, LIST2))
               DISSYQ = IRPDPD(IRREPQAB, ISYTYP(1, LISTQ))
               NUMSYQ = IRPDPD(IRREPQIJ, ISYTYP(2, LISTQ))
C
               I002 = I001 + NVRT2SQ*NUMSYT*IINTFP
               I003 = I002 + NVRT2SQ*NUMSYQ*IINTFP
C
               IF (MIN(DISSYT, NUMSYT, DISSYQ, NUMSYQ) .NE. 0) THEN
C
                  I004 = I003 + 3*IINTFP*MAX(DISSYT, NUMSYT, DISSYQ, 
     &                   NUMSYQ, NVRT2SQ)
C     
                  IF (I004 .LT. MXCOR) THEN
C
C In core version
C                  
                     CALL QTHRBAA2(ICORE(I000), ICORE(I001), 
     &                             ICORE(I002), NFVVY, DISSYT, NUMSYT,
     &                             DISSYQ, NUMSYQ, NVRT2SQ, IRREPTBF,
     &                             IRREPTIJ, IRREPQAB, IRREPQIJ, 
     &                             LIST2, LISTQ, POP(1, ISPIN),
     &                             VRT(1,ISPIN), ICORE(I003), IUHF,
     &                             ISPIN)
                  ELSE
                     CALL INSMEM('QTHRBDY2', I004, MXCOR)
                  ENDIF
C
               ENDIF
C
 10         CONTINUE
C
C ABAB spin case (UHF only). As shown in the expresion, Y(F,A)
C or Y(f,a) contribution for the two pieces identical to that of the 
C AAAA and BBBB spin cases. So we can do ABAB case inside the 
C same loop. 
C Now do the final Q(Ab,Ij) = T(Fb,Ij)*Y(F,A) [ISPIN = 1] and
C Q(Ab,Ij) = T(Af,Ij)*Y(f,b) [ISPIN = 2]. Allocate memory for
C the T2 vector and ICORE(I000) has the Y(N,I) arrays.
C the next multiplication we have to keep three four-index arrays and 
C one one-index array, resulting a requirement of an out-of core routine.
C
            LIST2 = 46
            LISTQ = 63 
C            
            DO 20 IRREP = 1, NIRREP
C               
C Now everything is totally symmetric
C
               IRREPTBF = IRREP
               IRREPTIJ = IRREPTBF
               IRREPQIJ = IRREPTIJ
               IRREPQAB = IRREPQIJ
C
               DISSYT = IRPDPD(IRREPTBF, ISYTYP(1, LIST2))
               NUMSYT = IRPDPD(IRREPTIJ, ISYTYP(2, LIST2))
               DISSYQ = IRPDPD(IRREPQAB, ISYTYP(1, LISTQ))
               NUMSYQ = IRPDPD(IRREPQIJ, ISYTYP(2, LISTQ))
C
               I002 = I001 + DISSYT*NUMSYT*IINTFP
               I003 = I002 + DISSYQ*NUMSYQ*IINTFP
C
               IF (MIN(DISSYT, NUMSYT, DISSYQ, NUMSYQ) .NE. 0) THEN
C
                  I004 = I003 + 3*IINTFP*MAX(DISSYT, NUMSYT, DISSYQ, 
     &                   NUMSYQ)
C     
                  IF (I004 .LT. MXCOR) THEN
C
C In core version
C                  
                     CALL QTHRBAB2(ICORE(I000), ICORE(I001),
     &                             ICORE(I002), NFVVY, DISSYT, NUMSYT,
     &                             DISSYQ, NUMSYQ, IRREPTBF, IRREPTIJ,
     &                             IRREPQAB, IRREPQIJ, LIST2, LISTQ,
     &                             POP(1, ISPIN), POP(1, 3-ISPIN),
     &                             VRT(1, ISPIN), VRT(1, 3-ISPIN),
     &                             ICORE(I003), IUHF, ISPIN)
                  ELSE
                     CALL INSMEM('QTHRBDY2', I004, MXCOR)
                  ENDIF
C
               ENDIF
C
 20         CONTINUE
C     
 100     CONTINUE
C     
      ENDIF
C
C Spin adapted code for RHF calculations (ABAB spin). The following
C piece of code is evaluated only for RHF calculations and (ISPIN = 1).
C     
C First do the Z(F,N) = [T(E,M)*Hbar(MN,EF) + T(E,M)Hbar(Mn,Ef)] [ISPIN = 1]
C Allocate memory for Z(F,N) and leave the rest of the core to store
C T1 vector and Hbar(MN,EF) in the subroutine MKZNF. We need to consider
C only one DPD irrep; the irrepx.
C
      IF (IUHF .EQ. 0) THEN
C
         LISTT1A = IAPRT1AA
         LISTT1B = IBPRT1AA
         LISTH1  = 19
         LISTH2  = 18
         ISPIN   = 1
C
         NFVOZ  = IRPDPD(IRREPX, 9)
C
         I000 = 1
         I010 = MXCOR + 1 - NFVOZ*IINTFP
C
         MXCOR = MXCOR - NFVOZ*IINTFP
C
         CALL MAKZNFTB(ICORE(I010), ICORE(I000), MXCOR, NFVOZ, LISTT1A,
     &                 LISTT1A, LISTH1, LISTH2, IUHF, IRREPX, 1)
C
C Now do the Y(F,A) = Z(F,N)*T(N,A) (ISPIN = 1). The Z(F,N) is in array ICORE(I010).
C Allocate memory for the Y(F,A) and T(N,A). Now only three one body terms are held
C in memory simultaneously, which is not a problem.
C            
         NFVVY = IRPDPD(1, 19)
         NFVOT = IRPDPD(IRREPX, 9)
C
         I001 = I000 + NFVVY*IINTFP
         I002 = I001 + NFVOT*IINTFP
C         
         MXCOR = MXCOR - (NFVVY + NFVOT)*IINTFP
C
         IF (I002 .LT. MXCOR) THEN
C     
            CALL GETLST(ICORE(I001), 1, 1, 1, 1, LISTT1B)
C
C Carry out the multiplication
C
            CALL MAKYFATB(ICORE(I010), ICORE(I000), ICORE(I001), 
     &                    ICORE(I002), NFVOZ, NFVVY, NFVOT, MXCOR,
     &                    IUHF, IRREPX, 1)
         ELSE
            CALL INSMEM('QTHRBDY2', I002, MXCOR)
         ENDIF
C
C Now do the final Q(Ab,Ij) = T(Fb,Ij)*Y(F,A) [ISPIN = 1] and
C Allocate memory for the T2 vector and ICORE(I000) has the Y(F,A) arrays.
C the next multiplication we have to keep three four-index arrays and 
C one one-index array, resulting a requirement of an out-of core routine.
C
         LIST2 = 46
         LISTQ = 63 
C            
         DO 30 IRREP = 1, NIRREP
C               
C Now everything is totally symmetric
C
            IRREPTBF = IRREP
            IRREPTIJ = IRREPTBF
            IRREPQIJ = IRREPTIJ
            IRREPQAB = IRREPQIJ
C
            DISSYT = IRPDPD(IRREPTBF, ISYTYP(1, LIST2))
            NUMSYT = IRPDPD(IRREPTIJ, ISYTYP(2, LIST2))
            DISSYQ = IRPDPD(IRREPQAB, ISYTYP(1, LISTQ))
            NUMSYQ = IRPDPD(IRREPQIJ, ISYTYP(2, LISTQ))
C     
            I002 = I001 + DISSYT*NUMSYT*IINTFP
            I003 = I002 + DISSYQ*NUMSYQ*IINTFP
C
            IF (MIN(DISSYT, NUMSYT, DISSYQ, NUMSYQ) .NE. 0) THEN
C     
               I004 = I003 + 3*IINTFP*MAX(DISSYT, NUMSYT, DISSYQ, 
     &                NUMSYQ)
C     
               IF (I004 .LT. MXCOR) THEN
C
C In core version
C                  
                  CALL QTHRBAB2(ICORE(I000), ICORE(I001), ICORE(I002),
     &                          NFVVY, DISSYT, NUMSYT, DISSYQ, NUMSYQ, 
     &                          IRREPTBF, IRREPTIJ, IRREPQAB, IRREPQIJ,
     &                          LIST2, LISTQ, POP(1, ISPIN),
     &                          POP(1, 3-ISPIN), VRT(1, ISPIN),
     &                          VRT(1, 3-ISPIN), ICORE(I003), IUHF,
     &                          ISPIN)
               ELSE
                  CALL INSMEM('QTHRBDY2', I004, MXCOR)
               ENDIF
C
            ENDIF
C
 30      CONTINUE
C
      ENDIF
C     
      RETURN
      END
