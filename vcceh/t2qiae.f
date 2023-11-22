C
      SUBROUTINE T2QIAE(ICORE, MAXCOR, IRREPX, IOFFSET, IUHF, FACT)
C
C This program computes,
C
C   Z(AE) = - 1/2 SUM F SUM M N  Hbar(MN,EF)*T(MN,AF)
C
C The multiplication is carried out using symmetry and the resulting
C product array is stored in a symmetry adapted way. Notice the
C difference from the CC code in which T(IN,EF) is in fact Tau(prime)(IN,EF).
C Spin integrated expressions are given as follows
C 
C UHF
C
C   Z(AE) = - SUM F SUM M<N Hbar(MN,EF)*T(MN,AF) [AAAA]
C
C   Z(ae) = - SUM f SUM m<n Hbar(mn,ef)*T(mn,af) [BBBB]
C
C   Z(AE) = - SUM f SUM M,n Hbar(Mn,Ef)*T(Mn,Af) [ABAB]
C
C   Z(ae) = - SUM F SUM m,N Hbar(mN,eF)*T(mN,aF) [BABA]
C
C RHF spin adapted expression
C
C  Z(AE) = - SUM f SUM M,n {2*Hbar(Mn,Ef) - Hbar(Mn,Fe)}*T(Mn,Af)
C
C Originaly coded by JG June/90 and modified by ajith 06/1994
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      INTEGER DIRPRD,DISSYW,DISSYT,POP,VRT
      CHARACTER*6 SPCASE(2)
      DIMENSION ICORE(MAXCOR)
C
      COMMON /SYM/ POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF1BB,NF2AA,NF2BB
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
      DATA SPCASE /'AA =  ', 'BB =  '/
C
      MXCOR = MAXCOR
C
      DO 1000 ISPIN = 1, IUHF + 1
C     
C AAAA and BBBB spin cases
C
         NVRTSQ = 0
         DO 50 IRREPE = 1, NIRREP
            IRREPF = DIRPRD(IRREPE, IRREPX)
            NVRTSQ = NVRTSQ + VRT(IRREPF, ISPIN)*
     &               VRT(IRREPE, ISPIN)
 50      CONTINUE
C
         I000 = MXCOR + 1 - NVRTSQ*IINTFP
C
         CALL GETLST(ICORE(I000), 1, 1, 1, ISPIN, INTIAE)
C
         IF (IUHF .EQ. 1) THEN
C     
            LISTT   = ISPIN + (IAPRT2AA2 - 1)
            LISTW   = ISPIN + 13
C     
C Loop over irreps of NM block
C     
            DO 100 IRREPMN = 1, NIRREP
               IRREPWMN = IRREPMN
               IRREPWEF = IRREPWMN
               IRREPTMN = IRREPWMN
               IRREPTAF = DIRPRD(IRREPTMN, IRREPX)
C     
C Retrieve Hbar integrals and T2 amplitudes
C     
               NVRT2SQ1 = 0
               NVRT2SQ2 = 0
C
               DO 110 IRREPE = 1, NIRREP
                  IRREPF = DIRPRD(IRREPE, IRREPWEF)
                  NVRT2SQ1 = NVRT2SQ1 + VRT(IRREPE, ISPIN)*
     &                       VRT(IRREPF, ISPIN)
 110           CONTINUE
C
               DO 120 IRREPF = 1, NIRREP
                  IRREPAO = DIRPRD(IRREPF, IRREPTAF)
                  NVRT2SQ2 = NVRT2SQ2 + VRT(IRREPAO, ISPIN)*
     &                       VRT(IRREPF, ISPIN)
 120           CONTINUE
C     
               DISSYW = IRPDPD(IRREPWEF, ISYTYP(1, LISTW))
               NUMSYW = IRPDPD(IRREPWMN, ISYTYP(2, LISTW))
               DISSYT = IRPDPD(IRREPTAF, ISYTYP(1, LISTT))
               NUMSYT = IRPDPD(IRREPTMN, ISYTYP(2, LISTT)) 
C               
               I001 = 1
               I002 = I001 + IINTFP*MAX(NUMSYW*NVRT2SQ1, 
     &                NUMSYT*NVRT2SQ2)
               I003 = I002 + IINTFP*MAX(NUMSYW*NVRT2SQ1, 
     &                NUMSYT*NVRT2SQ2)
               I004 = I003 + IINTFP*MAX(NUMSYT*DISSYT, NUMSYW*DISSYW)
C
               IF (MIN(NUMSYT, NUMSYW, DISSYT, DISSYW) .NE. 0) THEN
C
                  I005 = I004 + IINTFP*MAX(DISSYT, DISSYW, NUMSYW, 
     &                   NUMSYT)
C
                  IF(I005 .LT. MXCOR) THEN
C     
C In core version
C
                     CALL MKT2IAEAA(ICORE(I001), ICORE(I002), 
     &                              ICORE(I003), ICORE(I000), ISPIN,
     &                              POP(1,ISPIN), VRT(1,ISPIN),
     &                              NVRTSQ, NVRT2SQ1, NVRT2SQ2, DISSYW,
     &                              DISSYT, NUMSYW, NUMSYT, IRREPWEF,
     &                              IRREPTAF, IRREPWMN, IRREPTMN, 
     &                              IRREPX, ICORE(I004), IOFFSET, IUHF,
     &                              FACT)
                  ELSE
                     CALL INSMEM('T2QIAE', I005, MXCOR)
                  ENDIF
C
               ELSE
               ENDIF 
C     
 100        CONTINUE
         ENDIF
C     
C ABAB spin case
C     
         LISTT = IAPRT2AB5
         LISTW = 16
C     
C Loop over irreps of MN block
C     
         DO 200 IRREPMN = 1, NIRREP
C
            IRREPWMN = IRREPMN
            IRREPWEF = IRREPWMN
            IRREPTMN = IRREPWMN
            IRREPTAF = DIRPRD(IRREPTMN, IRREPX)
C     
C Retrive integrals and T2 amplitudes 
C     
            DISSYW = IRPDPD(IRREPWEF, ISYTYP(1, LISTW))
            NUMSYW = IRPDPD(IRREPWMN, ISYTYP(2, LISTW))
            DISSYT = IRPDPD(IRREPTAF, ISYTYP(1, LISTT))
            NUMSYT = IRPDPD(IRREPTMN, ISYTYP(2, LISTT))
C
            I001 = 1
            I002 = I001 + IINTFP*MAX(NUMSYW*DISSYW, NUMSYT*DISSYT)
            I003 = I002 + IINTFP*MAX(NUMSYW*DISSYW, NUMSYT*DISSYT)
            I004 = I003 + IINTFP*MAX(NUMSYW*DISSYW, NUMSYT*DISSYT)
C     
            IF (MIN(NUMSYT, NUMSYW, DISSYT, DISSYW) .NE. 0) THEN
C
               I005 = I004 + 3*IINTFP*MAX(DISSYT, DISSYW, NUMSYW,
     &                NUMSYT)
C
               IF(I005 .LE. MXCOR) THEN
C         
C In core version
C     
                  CALL MKT2IAEAB(ICORE(I001), ICORE(I002),
     &                           ICORE(I003), ICORE(I000), ISPIN, 
     &                           POP(1,ISPIN), POP(1,3-ISPIN), 
     &                           VRT(1,ISPIN), VRT(1,3-ISPIN), 
     &                           NVRTSQ, DISSYW, DISSYT, NUMSYW,
     &                           NUMSYT, IRREPWEF, IRREPTAF,
     &                           IRREPWMN, IRREPTMN, IRREPX,
     &                           ICORE(I004), IOFFSET, IUHF, FACT)
               ELSE
                  CALL INSMEM('T2QIAE', I005, MXCOR)
               ENDIF
C
            ELSE
            ENDIF
C     
 200     CONTINUE
C            
C Update the F(AE) contributions on the disk
C
         CALL PUTLST(ICORE(I000), 1, 1, 1, ISPIN, INTIAE)
C
         IF (IFLAGS(1) .GE. 20) THEN
C     
            NSIZE = NVRTSQ
C
            CALL HEADER('Checksum @-T2QIAE', 0, LUOUT)
            
            WRITE(LUOUT, *) SPCASE(ISPIN), SDOT(NSIZE, 
     &                      ICORE(I000), 1, ICORE(I000), 1)
         ENDIF
C
 1000 CONTINUE
C     
      RETURN
      END
