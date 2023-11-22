      SUBROUTINE EVALREPULI(IPRIM, JPRIM, INCRF, JNCRF, ITYPE, JTYPE, 
     &                      IPRMCOUNT, JPRMCOUNT, ICFCOUNT, JCFCOUNT,
     &                      NTOTPRIM, NTOTCRF, MAXPRM, CNTMU, CNTNU,
     &                      EXPS, PCOEF, CENTER, REPLINT, PRDTINT, 
     &                      TMP1, TMP2, TMP3, FAC, OVERLAP_POTN)
C
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      LOGICAL OVERLAP_POTN, FILE_EXIST
C
      DIMENSION CNTMU(3), CNTNU(3), EXPS(NTOTPRIM), AAA(27), BBB(27),
     &          DISTN(3), FAC(9,9), LMN(27), JMN(27), CNTP(3), 
     &          REPLINT(NTOTCRF, NTOTCRF), PRDTINT(NTOTCRF, NTOTCRF),
     &          PCOEF(NTOTPRIM, NTOTCRF),CENTER(3),TMP1(MAXPRM, MAXPRM),
     &          TMP2(MAXPRM, MAXPRM),TMP3(MAXPRM, MAXPRM)
C
      COMMON /HIGHL/ LMNVAL(3, 84), ANORM(84)
      DATA   DZERO /0.0D+00/, ONE /1.00D0/
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
C$$$      IPRMINDX = IPRMCOUNT  
      ICFINDX = ICFCOUNT
      CALL ZERO(TMP1, MAXPRM*MAXPRM)
      CALL ZERO(TMP2, MAXPRM*MAXPRM)

C
C$$$      LENGTH  = NTOTPRIM*NTOTCRF
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
C$$$               ICFINDX  = ICFINDX  + 1 
               JPRMINDX = JPRMCOUNT
               JNDX = 0
               INDX = INDX + 1 
C                     
               DO RGTPRIM = 1, JPRIM
C                  
                  JPRMINDX = JPRMINDX + 1 
C$$$                  JCFINDX  = JCFINDX  + 1
                  JNDX  = JNDX + 1 
C
                  TMP1(INDX, JNDX) = 0.D0
                  TMP2(INDX, JNDX) = 0.D0
C
C Get the common center of the gaussians (LFTPRIM and RGTPRIM)
C
                  CALL CIVPT(CNTMU, EXPS(IPRMINDX), CNTNU,
     &                       EXPS(JPRMINDX), CNTP, GAMA, EFACTP)
C
                  CALL RHFTCE(AAA, CNTMU, CNTP, FAC, ITYPE, ITM, LMN)
                  CALL RHFTCE(BBB, CNTNU, CNTP, FAC, JTYPE, JTM, JMN)
C
                  DO 9 K=1,3
                     DISTN(K) = CNTP(K) - CENTER(K)
 9                CONTINUE
                  VAL_REPLS = 0.0D0
                  VAL_OVLAP = 0.0D0
C
C$$$                  Write(6,*) "IPRMINDX, JPRMINDX", IPRMINDX, JPRMINDX
                  DO IIQ = 1, ITM

                     IF (.NOT. (AAA(IIQ) .EQ. 0.D0)) THEN 
C
                        IL = IBTAND(IBTSHR(LMN(IIQ),20),2**10 - 1)
                        IM = IBTAND(IBTSHR(LMN(IIQ),10),2**10 - 1)
                        IN = IBTAND(LMN(IIQ),2**10 - 1)
C
                        DO JJQ = 1, JTM
C
                           DD = AAA(IIQ)*BBB(JJQ)
C
                           IF (.NOT. (ABS(DD).LT.1.0D-19)) THEN
C
                              JL = IBTAND(IBTSHR(JMN(JJQ),20),2**10 - 1)
     &                           + IL
                              JM = IBTAND(IBTSHR(JMN(JJQ),10),2**10 - 1)
     &                           + IM
                              JN = IBTAND(JMN(JJQ),2**10 - 1) + IN
C
C For overlap potentials calculate the overlap integral.
C
                              IF (OVERLAP_POTN) THEN
C
                                 CALL OVERLAP(JL, JM, JN, GAMA, V)
                                 VAL_OVLAP = DD*V + VAL_OVLAP

C     
                              ELSE
C
C For Exchange/Correlation Potential calculate the repulsion integral.
C 
                                 CALL EVALREPLP(JL, JM, JN, GAMA, V,
     &                                          DISTN)
                                 VAL_REPLS = DD*V + VAL_REPLS
C
                              ENDIF
C
                           ENDIF
C                        
                        ENDDO
C                     
                     ENDIF 
C
                  ENDDO
C$$$                  Write(6,*) "INDX, JNDX =", INDX, JNDX
C$$$                  Write(6,*) "VAL, EFCTP =", VAL, EFACTP
C
                  IF (OVERLAP_POTN) THEN
                     TMP1(INDX, JNDX) = VAL_OVLAP*EFACTP + TMP1(INDX,
     &                                  JNDX)
                  ELSE
                     TMP1(INDX, JNDX) = VAL_REPLS*EFACTP + TMP1(INDX,
     &                                  JNDX)
                  ENDIF
C
C We can built the product here for Exchange/Correlation Potential
C
                  CALL BULTPRDUCT(INDX, JNDX, ITYPE, JTYPE, MAXPRM, 
     &                            EXPS(IPRMINDX), EXPS(JPRMINDX), 
     &                            CENTER, CNTMU, CNTNU, TMP2)
C
C Loop over primitives end here!
C
               ENDDO
            ENDDO
C
C$$$            Write(6,*) "The repulsion integral"
C$$$            CALL OUTPUT(TMP1, 1, IPRIM, 1, JPRIM, IPRIM, JPRIM, 1)
CSSS            Write(6,*) "The product integral"
CSSS            CALL OUTPUT(TMP2, 1, IPRIM, 1, JPRIM, IPRIM, JPRIM, 1)
CSSS            Write(6,*) "CONTRACTION COEFICIENTS"
CSSS           CALL OUTPUT(PCOEF, 1, NTOTPRIM, 1, NTOTCRF, NTOTPRIM,
CSSS     &                  NTOTCRF, 1)
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
     &                  DZERO, REPLC, 1)
C
C Built the contracted product functions for this shell.
C
            CALL ZERO(TMP3, MAXPRM*MAXPRM)
            CALL XGEMM('N', 'N', IPRIM, 1, JPRIM, ONE, TMP2, MAXPRM, 
     &                  PCOEF(JOFFC, JCFINDX), NTOTPRIM, DZERO, TMP3,
     &                  MAXPRM)
            CALL XGEMM('T', 'N', 1, 1, IPRIM, ONE, 
     &                  PCOEF(IOFFC, ICFINDX), NTOTPRIM, TMP3, MAXPRM,
     &                  DZERO, PRDCT, 1)
C
C$$$            Write(6,*) 
C$$$            Write(6,*) "ICFINDX, JCFINDX =", ICFINDX, JCFINDX
C$$$            Write(6,*) "REPLC, REPLD =", REPLC, REPLD
            REPLINT(ICFINDX, JCFINDX) = REPLC
            PRDTINT(ICFINDX, JCFINDX) = PRDCT            
C
C Loop over contracted functions end here!
C
         ENDDO
      ENDDO
C
      RETURN 
      END
