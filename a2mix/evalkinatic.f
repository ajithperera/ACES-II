      SUBROUTINE EVALKINATIC(IPRIM, JPRIM, INCRF, JNCRF, ITYPE, JTYPE, 
     &                       IPRMCOUNT, JPRMCOUNT, ICFCOUNT, JCFCOUNT,
     &                       NTOTPRIM, NTOTCRF, MAXPRM, CNTMU, CNTNU,
     &                       EXPS, PCOEF, CENTER, KINTINT, 
     &                       TMP1, TMP2, TMP3, FAC, NATOMS)
C
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      DOUBLE PRECISION KINTINT(NTOTCRF, NTOTCRF), KINTIC
      LOGICAL OVERLAP_POTN, FILE_EXIST
C
      DIMENSION CNTMU(3), CNTNU(3), EXPS(NTOTPRIM), AAA(27), BBB(27),
     &          DISTN(3), FAC(9,9), LMN(27), JMN(27), CNTP(3), 
     &          PCOEF(NTOTPRIM, NTOTCRF),CENTER(3),TMP1(MAXPRM, MAXPRM),
     &          TMP2(MAXPRM, MAXPRM),TMP3(MAXPRM, MAXPRM), TEMP(4, 3)
C
      COMMON /HIGHL/ LMNVAL(3, 84), ANORM(84)
      DATA   DZERO /0.0D+00/, ONE /1.00D0/, TWO/2.00D+0/, FOUR/4.00D+0/
C
      IBTAND(I,J) = IAND(I,J)
      IBTOR(I,J)  = IOR(I,J)
      IBTXOR(I,J) = IEOR(I,J)
      IBTSHL(I,J) = ISHFT(I,J)
      IBTSHR(I,J) = ISHFT(I,-J)
      IBTNOT(I)   = NOT(I)
C
C Loop over number of conttracted functions
C
      ICFINDX = ICFCOUNT
      CALL ZERO(TMP1, MAXPRM*MAXPRM)
      CALL ZERO(TMP3, MAXPRM*MAXPRM)
C
      DO ICRF = 1, INCRF
C
         ICFINDX = ICFINDX + 1 
         JCFINDX = JCFCOUNT
C 
         DO JCRF = 1, JNCRF 
C
            JCFINDX = JCFINDX + 1 
C     
C Loop over primitive functions 
C
            IPRMINDX = IPRMCOUNT
            INDX = 0
            JNDX = 0

            DO LFTPRIM = 1, IPRIM
C
               IPRMINDX = IPRMINDX + 1 
               JPRMINDX = JPRMCOUNT
               JNDX = 0
               INDX = INDX + 1 
C                     
               DO RGTPRIM = 1, JPRIM
C                  
                  JPRMINDX = JPRMINDX + 1 
                  JNDX  = JNDX + 1 
C
                  TMP1(INDX, JNDX) = 0.0D0
                  CALL ZERO(TEMP, 12)
                  VAL_KINETIC = 0.0D00              
C
C Get the common center of the gaussians (LFTPRIM and RGTPRIM)
C
                  DO IXYZ = 1, 3 
                     
                     CALL GETIXYZ(IXYZ, IX, IY, IZ)
C
C This is the LB-1,MB-1,NB-1 contribution
C      
                     ILB = LMNVAL(1, ITYPE)
                     IMB = LMNVAL(2, ITYPE)
                     INB = LMNVAL(3, ITYPE)
                     JLB = LMNVAL(1, JTYPE)
                     JMB = LMNVAL(2, JTYPE)
                     JNB = LMNVAL(3, JTYPE)

                     CALL GETPRODUCT(ILB-IX, IMB-IY, INB-IZ,
     &                               JLB-IX, JMB-IY, JNB-IZ,
     &                               MAXPRM, EXPS(IPRMINDX),
     &                               EXPS(JPRMINDX), CENTER, CNTMU, 
     &                               CNTNU, PRDUCT)
C
                     TEMP(1, IXYZ) = (ILB*IX+IMB*IY+INB*IZ)*
     &                                (JLB*IX+JMB*IY+JNB*IZ)*PRDUCT

                     CALL GETPRODUCT(ILB-IX, IMB-IY, INB-IZ,
     &                               JLB+IX, JMB+IY, JNB+IZ,
     &                               MAXPRM, EXPS(IPRMINDX),
     &                               EXPS(JPRMINDX), CENTER, CNTMU, 
     &                               CNTNU, PRDUCT)
C
                     TEMP(2, IXYZ) = -TWO*(ILB*IX+IMB*IY+INB*IZ)*
     &                                 EXPS(JPRMINDX)*PRDUCT
C
                     CALL GETPRODUCT(ILB+IX, IMB+IY, INB+IZ,
     &                               JLB-IX, JMB-IY, JNB-IZ,
     &                               MAXPRM, EXPS(IPRMINDX),
     &                               EXPS(JPRMINDX), CENTER, CNTMU, 
     &                               CNTNU, PRDUCT)
C
                     TEMP(3, IXYZ) = -TWO*(JLB*IX+JMB*IY+JNB*IZ)*
     &                                 EXPS(IPRMINDX)*PRDUCT
C
                     CALL GETPRODUCT(ILB+IX, IMB+IY, INB+IZ,
     &                               JLB+IX, JMB+IY, JNB+IZ,
     &                               MAXPRM, EXPS(IPRMINDX),
     &                               EXPS(JPRMINDX), CENTER, CNTMU,
     &                               CNTNU, PRDUCT)
C
                     TEMP(4, IXYZ) = FOUR*EXPS(IPRMINDX)*EXPS(JPRMINDX)*
     &                               PRDUCT 
C
C Accumulate the contributions to form the full kinetic energy integral
C
                     VAL_KINETIC = TEMP(1, IXYZ) + TEMP(2, IXYZ) + 
     &                             TEMP(3, IXYZ) + TEMP(4, IXYZ) + 
     &                             VAL_KINETIC
C
C$$$                     Write(6,*) (TEMP1(1, I), I=1, 3)
C$$$                     Write(6,*) (TEMP2(1, J), J=1, 3)

                  ENDDO

C$$$                  Write(6,*) "INDX, JNDX =", INDX, JNDX
C$$$                  Write(6,*) "VAL, EFCTP =", VAL, EFACTP
C
                  TMP1(INDX, JNDX) = VAL_KINETIC 
C
C Loop over primitives end here!
C
               ENDDO
            ENDDO
C
C$$$            Write(6,*) "The Kinetic Energy integral"
C$$$            CALL OUTPUT(TMP1, 1, IPRIM, 1, JPRIM, IPRIM, JPRIM, 1)
C$$$            Write(6,*) "The product integral"
C$$$            CALL OUTPUT(TMP2, 1, IPRIM, 1, JPRIM, IPRIM, JPRIM, 1)
C$$$            Write(6,*) "CONTRACTION COEFICIENTS"
C$$$           CALL OUTPUT(PCOEF, 1, NTOTPRIM, 1, NTOTCRF, NTOTPRIM,
C$$$     &                  NTOTCRF, 1)
C 
C Built the contracted functions for this shell.
C
            IOFFC = IPRMCOUNT + 1
            JOFFC = JPRMCOUNT + 1
            CALL ZERO(TMP3, MAXPRM*MAXPRM)
C
C$$$            Write(6,*) "IOFFC, JOFFC =", IOFFC, JOFFC
            CALL XGEMM('N', 'N', IPRIM, 1, JPRIM, ONE, TMP1, MAXPRM, 
     &                  PCOEF(JOFFC, JCFINDX), NTOTPRIM, DZERO, TMP3,
     &                  MAXPRM)
            CALL XGEMM('T', 'N', 1, 1, IPRIM, ONE, 
     &                  PCOEF(IOFFC, ICFINDX), NTOTPRIM, TMP3, MAXPRM,
     &                  DZERO, KINTIC, 1)
C
C$$$            Write(6,*) 
C$$$            Write(6,*) "ICFINDX, JCFINDX =", ICFINDX, JCFINDX
C$$$            Write(6,*) "REPLC, REPLD =", REPLC, REPLD
            KINTINT(ICFINDX, JCFINDX) = KINTIC
C
C Loop over contracted functions end here!
C
         ENDDO
      ENDDO
C
      RETURN 
      END




