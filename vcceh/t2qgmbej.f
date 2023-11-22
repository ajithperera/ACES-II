C
      SUBROUTINE T2QGMBEJ(ICORE, MAXCOR, SPCASE, IUHF, IRREPX, 
     &                    IOFFSET)
C
C This routine and dependents compute the ring-type W(MBEJ)
C intermediate for all six possible spin cases. The algorithm
C assumes in-core storage of symmetry packed target, T2 and W
C vectors and uses the DPD symmetry approach to evaluate the
C contractions.
C
C Spin orbital equation for this intermediate (eventual modification
C is required to include T1 contributions)
C
C  W(mbej)= - (1/2)*SUM T(jn,[fb]) Hbar(mn,[ef])
C                   n,f
C
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION X,TWO,ONE,ONEM,HALF,ZILCH,TWOM,MHALF,SDOT
      CHARACTER*4 SPCASE 
      DIMENSION ICORE(MAXCOR)
C     
      COMMON /INFO/ NOCA,NOCB,NVRTA,NVRTB
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT1AA,NT1BB,NF1AA,NF1BB,NF2AA,
     &             NF2BB
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
C Common blocks used in the quadratic term
C
      DATA ZILCH /0.0D+00/
      DATA ONE /1.0D+00/
      DATA HALF /0.5D+00/
      DATA MHALF /-0.5D+00/
      DATA TWO /2.0D+00/
      DATA TWOM /-2.0D+00/
      DATA ONEM /-1.0D+00/
C
      NNP1O2(I) = I*(I+1)/2
      X = 1.0D+00 
C
C Do spin case AAAA and BBBB.
C
C  G(MBEJ)= (1/2)*SUM [T(JN,[FB])*Hbar(MN,[EF]) - T(Jn,Bf)*Hbar(Mn,Ef)]
C  G(mbej)= (1/2)*SUM [T(jn,[fb])*Hbar(mn,[ef]) - T(jN,bF)*Hbar(mN,eF)]
C
C Note the difference in sign from the sign given in the above spin 
C orbital expression. For the reasons become clear later on, here we 
C evaluate the negative of the AAAA and BBBB contributions. Also 
C note that the sign here is also consistent with the sign of T1
C ring AAAA and BBBB contributions, facilitating eventual addition.
C
      IF (SPCASE .EQ. 'AAAA' .OR. SPCASE .EQ. 'BBBB') THEN
         IF(SPCASE .EQ. 'AAAA') THEN
            ISPIN = 1
         ELSE IF (SPCASE .EQ. 'BBBB') THEN
            ISPIN = 2
         ENDIF
C
         MAXSIZ = 0
         LISTQ  = (INGMCAA - 1) + ISPIN
         LISTWA = 18 + ISPIN
         LISTWB = 19 - ISPIN  
         LISTTA = (IAPRT2AA1 - 1) + ISPIN
         LISTTB = (IAPRT2AB1 - 1) + ISPIN
C
         IF (IUHF .EQ. 0) LISTTB = IAPRT2AB2
C
         DO 100 IRREPFN = 1, NIRREP 
C
            IRREPWFN = IRREPFN
            IRREPWEM = IRREPWFN
            IRREPTFN = IRREPWFN
            IRREPTBJ = DIRPRD(IRREPTFN, IRREPX)
            IRREPGBJ = IRREPTBJ
            IRREPGEM = DIRPRD(IRREPGBJ, IRREPX)
C
            DSSYWA = IRPDPD(IRREPWEM, ISYTYP(1, LISTWA))
            NMSYWA = IRPDPD(IRREPWFN, ISYTYP(2, LISTWA))
            DSSYTA = IRPDPD(IRREPTFN, ISYTYP(1, LISTTA))
            NMSYTA = IRPDPD(IRREPTBJ, ISYTYP(2, LISTTA))
            DSSYWB = IRPDPD(IRREPWEM, ISYTYP(1, LISTWB))
            NMSYWB = IRPDPD(IRREPWFN, ISYTYP(2, LISTWB))
            DSSYTB = IRPDPD(IRREPTFN, ISYTYP(1, LISTTB))
            NMSYTB = IRPDPD(IRREPTBJ, ISYTYP(2, LISTTB))
            DISSYQ = IRPDPD(IRREPGEM, ISYTYP(1, LISTQ))
            NUMSYQ = IRPDPD(IRREPGBJ, ISYTYP(2, LISTQ))
C
            I000 = 1
            I010 = I000 + IINTFP*MAX(DISSYQ*NUMSYQ, NMSYWA*DSSYWA,
     &             NMSYTA*DSSYTA)
            I020 = I010 + IINTFP*MAX(DISSYQ*NUMSYQ, NMSYWA*DSSYWA,
     &             NMSYTA*DSSYTA)
            I030 = I020 + IINTFP*MAX(DISSYQ*NUMSYQ, NMSYWA*DSSYWA,
     &              NMSYTA*DSSYTA)
            I030 = MAX(I030, I010 + DISSYQ*NUMSYQ*IINTFP)
C
            NSIZE = DISSYQ*NUMSYQ
C
            CALL IZERO(ICORE(I000), DISSYQ*NUMSYQ*IINTFP)
C
            CALL GETLST(ICORE(I020), 1, NMSYTA, 1, IRREPTBJ, LISTTA)
C     
            CALL DOCNTRC(ICORE(I000), ICORE(I010), ICORE(I020),
     &                   DISSYQ, NUMSYQ, DSSYTA, NMSYTA, DSSYWA,
     &                   NMSYWA, IRREPWFN, LISTWA, LISTWA, LISTTA,
     &                  'WxT', X, HALF, 0)
C     
            I020 = I010 + IINTFP*MAX(NMSYWB*DSSYWB, NMSYTB*DSSYTB, 
     &             DISSYQ*NUMSYQ)
            I030 = I020 + IINTFP*MAX(NMSYWB*DSSYWB, NMSYTB*DSSYTB,
     &             DISSYQ*NUMSYQ)
C               
            CALL GETLST(ICORE(I020), 1, NMSYTB, 1, IRREPTBJ, LISTTB)
C
            CALL DOCNTRC(ICORE(I000), ICORE(I010), ICORE(I020),
     &                   DISSYQ, NUMSYQ, DSSYTB, NMSYTB, DSSYWB,
     &                   NMSYWB, IRREPWFN, LISTWB, LISTWB, LISTTB,
     &                   'WxT', X, MHALF, 0)
C
            CALL SUMSYM2(ICORE(I000), ICORE(I010), NSIZE, 1, IRREPGBJ,
     &                   LISTQ)
C
 100     CONTINUE
C
         IF (IFLAGS(1) .GE. 20) THEN
C
            IF (SPCASE .EQ. 'AAAA') THEN
               NTOTAL = IDSYMSZ(IRREPX, 9, 9)
            ELSE IF (SPCASE .EQ. 'BBBB') THEN
               NTOTAL = IDSYMSZ(IRREPX, 10, 10)
            ENDIF
C
            CALL GETALL(ICORE(I000), NTOTAL, IRREPX, LISTQ)
            CALL HEADER('Checksum @-T2QGMBEJ', 0, LUOUT)
            WRITE(LUOUT, *) SPCASE, ' = ', SDOT(NTOTAL,
     &                      ICORE(I000), 1, ICORE(I000), 1)
         ENDIF
C     
C Do spin case ABAB and  BABA
C
C  G(MbEj) = - (1/2)*SUM [T(jn,[fb])*Hbar(Mn,Ef) - T(jN,bF)*Hbar(MN,[EF])]
C  G(mBeJ) = - (1/2)*SUM [T(JN,[FB])*Hbar(mN,eF) - T(Jn,Bf)*Hbar(mn,[ef])]
C 
C Note that for ABAB and BABA the overall sign is consistent with that
C of given in above spin orbital expression. The sign is also consistent 
C with the the T1 ring ABAB and BABA contributions, facilitating eventual
C addition.
C
      ELSE IF (SPCASE .EQ. 'ABAB' .OR. SPCASE .EQ. 'BABA') THEN
C
         IF (IUHF .NE. 0) THEN
            IF (SPCASE .EQ. 'ABAB') THEN
C
               LISTQ  = INGMCABAB
               LISTQ1 = INGMCBABA
               LISTTB = IAPRT2AB2
               LISTWA = 18
               LISTWB = 19
C
               IF (IUHF .EQ. 0) THEN
                  LISTTA = IAPRT2AA1
               ELSE
                  LISTTA = IAPRT2BB1
               ENDIF
C               
            ELSE IF (SPCASE .EQ. 'BABA') THEN 
C
               LISTQ  = INGMCBABA
               LISTQ1 = INGMCABAB
               LISTTA = IAPRT2AA1
               LISTTB = IAPRT2AB1
               LISTWA = 17
               LISTWB = 20
            ENDIF
C
            DO 200 IRREPFN = 1 , NIRREP 
C
               IRREPWFN = IRREPFN
               IRREPWEM = IRREPWFN
               IRREPTFN = IRREPWFN
               IRREPTBJ = DIRPRD(IRREPTFN, IRREPX)
               IRREPGBJ = IRREPTBJ
               IRREPGEM = DIRPRD(IRREPGBJ, IRREPX)
C
               DSSYWA = IRPDPD(IRREPWEM, ISYTYP(1, LISTWA))
               NMSYWA = IRPDPD(IRREPWFN, ISYTYP(2, LISTWA))
               DSSYTA = IRPDPD(IRREPTFN, ISYTYP(1, LISTTA))
               NMSYTA = IRPDPD(IRREPTBJ, ISYTYP(2, LISTTA))
               DSSYWB = IRPDPD(IRREPWEM, ISYTYP(1, LISTWB))
               NMSYWB = IRPDPD(IRREPWFN, ISYTYP(2, LISTWB))
               DSSYTB = IRPDPD(IRREPTFN, ISYTYP(1, LISTTB))
               NMSYTB = IRPDPD(IRREPTBJ, ISYTYP(2, LISTTB))
               DISSYQ = IRPDPD(IRREPGEM, ISYTYP(1, LISTQ))
               NUMSYQ = IRPDPD(IRREPGBJ, ISYTYP(2, LISTQ))
C
               I000 = 1
               I010 = I000 + IINTFP*MAX(DISSYQ*NUMSYQ, NMSYWA*DSSYWA,
     &                NMSYTA*DSSYTA)
               I020 = I010 + IINTFP*MAX(DISSYQ*NUMSYQ, NMSYWA*DSSYWA,
     &                NMSYTA*DSSYTA)
               I030 = I020 + IINTFP*MAX(DISSYQ*NUMSYQ, NMSYWA*DSSYWA,
     &                NMSYTA*DSSYTA)
               I030 = MAX(I030, I010 + DISSYQ*NUMSYQ*IINTFP)
C
               NSIZE = DISSYQ*NUMSYQ
C
               CALL IZERO(ICORE(I000), DISSYQ*NUMSYQ*IINTFP)
C
               CALL GETLST(ICORE(I020), 1, NMSYTA, 1, IRREPTBJ,
     &                     LISTTA)
C
               CALL DOCNTRC(ICORE(I000), ICORE(I010), ICORE(I020),
     &                      DISSYQ, NUMSYQ, DSSYTA, NMSYTA, DSSYWA,
     &                      NMSYWA, IRREPWFN, LISTWA, LISTWA, LISTTA,
     &                     'WxT', X, MHALF, 0) 
C     
               I020 = I010 + IINTFP*MAX(NMSYWB*DSSYWB, NMSYTB*DSSYTB,
     &                DISSYQ*NUMSYQ)
               I030 = I020 + IINTFP*MAX(NMSYWB*DSSYWB, NMSYTB*DSSYTB,
     &                DISSYQ*NUMSYQ)
C
               I030 = MAX(I030, I010 + DISSYQ*NUMSYQ*IINTFP)
C
               CALL GETLST(ICORE(I020), 1, NMSYTB, 1, IRREPTBJ,
     &                     LISTTB)
C
               CALL DOCNTRC(ICORE(I000), ICORE(I010), ICORE(I020),
     &                      DISSYQ, NUMSYQ, DSSYTB, NMSYTB, DSSYWB,
     &                      NMSYWB, IRREPWFN, LISTWB, LISTWB, LISTTB,
     &                     'WxT', X, HALF, 0)
C
               CALL SUMSYM2(ICORE(I000), ICORE(I010), NSIZE, 1,
     &                      IRREPGBJ, LISTQ)
C
 200        CONTINUE   
C
            IF (IFLAGS(1) .GE. 20) THEN
C
               IF (SPCASE .EQ. 'ABAB') THEN
                  NTOTAL = IDSYMSZ(IRREPX, 9, 10)
               ELSE IF(SPCASE .EQ. 'BABA') THEN
                  NTOTAL = IDSYMSZ(IRREPX, 10, 9)
               ENDIF
C
               CALL GETALL(ICORE(I000), NTOTAL, IRREPX, LISTQ)
               CALL HEADER('Checksum @-T2QGMBEJ', 0, LUOUT)
               WRITE(LUOUT, *) SPCASE, ' = ', SDOT(NTOTAL,
     &                         ICORE(I000), 1, ICORE(I000), 1)
            ENDIF
C     
C Spin adapted code for RHF "ABAB" spin case.
C     
         ELSE
C
            LISTQ   = INGMCABAB
            LISTTA  = IAPRT2AB2
            LISTTB  = IAPRT2AB4
            LISTWA  = 18
            LISTWB  = 21
C
            DO 201 IRREPFN = 1, NIRREP
C
               IRREPWFN = IRREPFN
               IRREPWEM = IRREPWFN
               IRREPTFN = IRREPWFN
               IRREPTBJ = DIRPRD(IRREPTFN, IRREPX)
               IRREPGBJ = IRREPTBJ
               IRREPGEM = DIRPRD(IRREPGBJ, IRREPX)
C     
               DSSYWA = IRPDPD(IRREPWEM, ISYTYP(1, LISTWA))
               NMSYWA = IRPDPD(IRREPWFN, ISYTYP(2, LISTWA))
               DSSYTA = IRPDPD(IRREPTFN, ISYTYP(1, LISTTA))
               NMSYTA = IRPDPD(IRREPTBJ, ISYTYP(2, LISTTA))
               DSSYWB = IRPDPD(IRREPWEM, ISYTYP(1, LISTWB))
               NMSYWB = IRPDPD(IRREPWFN, ISYTYP(2, LISTWB))
               DSSYTB = IRPDPD(IRREPTFN, ISYTYP(1, LISTTB))
               NMSYTB = IRPDPD(IRREPTBJ, ISYTYP(2, LISTTB))
               DISSYQ = IRPDPD(IRREPGEM, ISYTYP(1, LISTQ))
               NUMSYQ = IRPDPD(IRREPGBJ, ISYTYP(2, LISTQ))
C
               I000 = 1
               I010 = I000 + IINTFP*MAX(DISSYQ*NUMSYQ, NMSYWA*DSSYWA,
     &                NMSYWB*DSSYWB, DSSYTB*NMSYTB, NMSYTA*DSSYTA)
               I020 = I010 + IINTFP*MAX(DISSYQ*NUMSYQ, NMSYWA*DSSYWA,
     &                NMSYWB*DSSYWB, DSSYTB*NMSYTB, NMSYTA*DSSYTA)
               I030 = I020 +  IINTFP*MAX(DISSYQ*NUMSYQ, NMSYWA*DSSYWA, 
     &                NMSYWB*DSSYWB, DSSYTB*NMSYTB, NMSYTA*DSSYTA)
C     
C Make {2 [T2(Nj,Fb)] - (Nj,Bf)}, Spin adapt the T2 amplitudes.
C Spin adaptation of Habr integrals is carried out in DOCNTRC.
C For RHF case NUMSYTA = NUMSYTB.
C     
               CALL GETLST(ICORE(I020), 1, NMSYTA, 1, IRREPTBJ,
     &                     LISTTA)
               CALL GETLST(ICORE(I010), 1, NMSYTA, 1, IRREPTBJ,
     &                     LISTTB)
C
               CALL SSCAL (NMSYTA*DSSYTA, TWO, ICORE(I020), 1)
               CALL SAXPY (NMSYTA*DSSYTA, ONEM, ICORE(I010), 1, 
     &                     ICORE(I020), 1)
C
               CALL DOCNTRC(ICORE(I000), ICORE(I010), ICORE(I020),
     &                      DISSYQ, NUMSYQ, DSSYTA, NMSYTA, DSSYWA,
     &                      NMSYWA, IRREPWFN, LISTWA, LISTWB, LISTTA,
     &                     'WxT', ZILCH, HALF, 1)
               NSIZE = DISSYQ*NUMSYQ
C     
C We also need spin-adapted T1RING contributions too.
C     
               CALL GETLST(ICORE(I010), 1, NUMSYQ, 1, IRREPGBJ,
     &                     INGMCABAB)
               CALL GETLST(ICORE(I020), 1, NUMSYQ, 1, IRREPGBJ,
     &                     INGMCBAAB)
C
               CALL SSCAL(NSIZE, TWO, ICORE(I010), 1)
               CALL SAXPY(NSIZE, ONEM, ICORE(I020), 1, ICORE(I010), 1)
               CALL SAXPY(NSIZE, ONE, ICORE(I010), 1, ICORE(I000), 1)
C
C Write RHF spin adapted [ABAB + BAAB] contribution to the disk. 
C
               CALL PUTLST(ICORE(I000), 1, NUMSYQ, 1, IRREPGBJ, LISTQ)
C
 201        CONTINUE
C
            IF (IFLAGS(1) .GE. 20) THEN
C     
               IF (SPCASE .EQ. 'ABAB') THEN
                  NTOTAL = IDSYMSZ(IRREPX, 9, 10)
               ENDIF
C
               CALL GETALL(ICORE(I000), NTOTAL, IRREPX, LISTQ)
               CALL HEADER('Checksum @-T2QGMBEJ', 0, LUOUT)
               WRITE(LUOUT, *) SPCASE, ' = ', SDOT(NTOTAL,
     &                         ICORE(I000), 1, ICORE(I000), 1)
            ENDIF
C
         ENDIF
C
      ELSE IF (SPCASE .EQ. 'ABBA' .OR. SPCASE .EQ. 'BAAB') THEN
C     
C Do spin cases ABBA and BAAB
C
C  W(Mb,eJ) = - (1/2)*SUM [T(Jn,Fb)*Hbar(Mn,Fe)]
C  W(mB,Ej) = - (1/2)*SUM [T(jN,fB)*Hbar(mN,fE)]
C
C Note the difference in sign from the sign given in the above spin 
C orbital expression. For the reasons become clear later on, here we 
C evaluate the negative of the AAAA and BBBB contributions. Also 
C note that the sign here is also consistent with the sign of T1
C ring AAAA and BBBB contributions, facilitating eventual addition.
C Also note that in the following piece of code when SPCASE is ABBA 
C we really calculating BAAB contribution. 
C
         IF (SPCASE .EQ. 'ABBA') THEN
C
            LISTQ  = INGMCBAAB
            LISTTA = IAPRT2AB3
            IF (IUHF .EQ. 0) LISTTA = IAPRT2AB4
            LISTWA = 21
C
         ELSE IF (SPCASE .EQ. 'BAAB') THEN
C
            LISTQ  = INGMCABBA
            LISTTA = IAPRT2AB4
            LISTWA = 22
C
         ENDIF
C
         DO 300 IRREPFN = 1, NIRREP 
C
            IRREPWFN = IRREPFN
            IRREPWEM = IRREPWFN
            IRREPTFN = IRREPWFN
            IRREPTBJ = DIRPRD(IRREPTFN, IRREPX)
            IRREPGBJ = IRREPTBJ
            IRREPGEM = DIRPRD(IRREPGBJ, IRREPX)
C            
            DSSYWA = IRPDPD(IRREPWFN, ISYTYP(1, LISTWA))
            NMSYWA = IRPDPD(IRREPWEM, ISYTYP(2, LISTWA))
            DSSYTA = IRPDPD(IRREPTFN, ISYTYP(1, LISTTA))
            NMSYTA = IRPDPD(IRREPTBJ, ISYTYP(2, LISTTA))
            DISSYQ = IRPDPD(IRREPGEM, ISYTYP(1, LISTQ))
            NUMSYQ = IRPDPD(IRREPGBJ, ISYTYP(2, LISTQ))
C
            I000 = 1
            I010 = I000 + IINTFP*MAX(DISSYQ*NUMSYQ, NMSYWA*DSSYWA,
     &             NMSYTA*DSSYTA)
            I020 = I010 + IINTFP*MAX(DISSYQ*NUMSYQ, NMSYWA*DSSYWA,
     &             NMSYTA*DSSYTA)
            I030 = I020 + IINTFP*MAX(DISSYQ*NUMSYQ, NMSYWA*DSSYWA,
     &             NMSYTA*DSSYTA)
            I030 = MAX(I030, I010 + DISSYQ*NUMSYQ*IINTFP)
C
            NSIZE = DISSYQ*NUMSYQ
C
            CALL IZERO(ICORE(I000), DISSYQ*NUMSYQ*IINTFP)
C
            CALL GETLST(ICORE(I020), 1, NMSYTA, 1, IRREPTBJ, LISTTA)
C
            CALL DOCNTRC(ICORE(I000), ICORE(I010), ICORE(I020),
     &                   DISSYQ, NUMSYQ, DSSYTA, NMSYTA, DSSYWA,
     &                   NMSYWA, IRREPWFN, LISTWA, LISTWA, LISTTA,
     &                  'WxT', X, MHALF, 0)
C
            CALL SUMSYM2(ICORE(I000), ICORE(I010), NSIZE, 1, IRREPGBJ,
     &                   LISTQ)
C
 300     CONTINUE
C
         IF (IFLAGS(1) .GE. 20) THEN
C     
            IF (SPCASE .EQ. 'ABBA') THEN
               NTOTAL = IDSYMSZ(IRREPX, 12, 12)
            ELSE IF(SPCASE .EQ. 'BAAB') THEN
               NTOTAL = IDSYMSZ(IRREPX, 11, 11)
            ENDIF
C
            CALL GETALL(ICORE(I000), NTOTAL, IRREPX, LISTQ)
            CALL HEADER('Checksum @-T2QGMBEJ', 0, LUOUT)
            WRITE(LUOUT, *) SPCASE, ' = ', SDOT(NTOTAL,
     &                      ICORE(I000), 1, ICORE(I000), 1)
            ENDIF
C
      ENDIF
C
      RETURN
      END
