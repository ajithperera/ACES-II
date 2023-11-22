      SUBROUTINE EVALCORPOT(IANGMOM, JANGMOM, KANGMOM, LANGMOM, INCRF,
     &                      JNCRF, KNCRF, LNCRF, IBEGIN, JBEGIN, 
     &                      KBEGIN, LBEGIN, ICFNCOUNT, JCFNCOUNT, 
     &                      KCFNCOUNT, LCFNCOUNT, NTOTCRF, LENGTH, 
     &                      SHELIEQJ, SHELKEQL, SHELIJKL, ANGMOMIJ, 
     &                      ANGMOMKL, REPLINT, PRDUTINT, SCFDENT,
     &                      RELDENT, SCFDEND, RELDEND, GAMMATMP, 
     &                      GAMMA, POTNL_R, REPLSN_FACT, EXCHNG_FACT, 
     &                      CORLN_FACT_O, CORLN_FACT_T, IUHF)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      DIMENSION REPLINT(NTOTCRF, NTOTCRF), PRDUTINT(NTOTCRF, NTOTCRF),
     &          GAMMATMP(LENGTH), SCFDENT(NTOTCRF, NTOTCRF),
     &          RELDENT(NTOTCRF, NTOTCRF), SCFDEND(NTOTCRF, NTOTCRF),
     &          RELDEND(NTOTCRF, NTOTCRF), GAMMA(LENGTH)
C
      LOGICAL SHELIEQJ, SHELKEQL, ANGMOMIJ, ANGMOMKL, SHELIJKL
C
      DATA ONE, TWO, FOUR, FOURTH, HALF /1.0D0, 2.0D0, 4.0D0, 0.25D0, 
     &                                   0.50D0/
C
      INDEX(I,J,K,L) = INCRF*(JNCRF*(KNCRF*(L-1)+K-1)+J-1)+I
C
      FACIJKL = ONE
C
      IF (.NOT. SHELIJKL) FACIJKL=TWO*FACIJKL
      IF (.NOT. SHELIEQJ) FACIJKL=TWO*FACIJKL
      IF (.NOT. SHELKEQL) FACIJKL=TWO*FACIJKL
C  
C$$$      Write(6,*) "The gamma elements"
C$$$      WRITE(6,*) (GAMMA(I), I= 1, LENGTH)

C
C$$$      Write(6,*) "shell data and sym fac", shelijkl,shelieqj,shelkeql,
C$$$     &            angmomij,angmomkl, facijkl
C
      ICFNCOUNT = IBEGIN
C$$$      I         = 1 
C
      DO IMOMN = 1, IANGMOM
C
C$$$         IF (IMOMN .GT. 2) I = I + 1 
         JCFNCOUNT = JBEGIN
C$$$         ICFNCOUNT = (IMOMN - I)*INCRF + ICFNCOUNT 
C$$$         J         = 1 
C
         JMOMFINAL = JANGMOM
         IF (SHELIEQJ) JMOMFINAL = IMOMN
C
         DO JMOMN = 1, JMOMFINAL
C     
C$$$            IF (JMOMN .GT. 2) J = J + 1 
            KCFNCOUNT = KBEGIN
C$$$            JCFNCOUNT = (JMOMN - J)*JNCRF + JCFNCOUNT 
C$$$            k         = 1 
            FACIJMM   = ONE
            IF (SHELIEQJ .AND. IMOMN .NE. JMOMN) FACIJMM=TWO
C
            DO KMOMN = 1, KANGMOM
C
C$$$               IF (KMOMN .GT. 2) k = K + 1 
               LCFNCOUNT = LBEGIN
C$$$               KCFNCOUNT = (KMOMN - K)*KNCRF + KCFNCOUNT 
C$$$               L         = 1
C
               LMOMFINAL = LANGMOM
               IF (SHELKEQL) LMOMFINAL = KMOMN
     
               DO LMOMN = 1, LMOMFINAL
C
C$$$                  IF (LMOMN .GT. 2) L = L + 1 
C$$$                  LCFNCOUNT = (LMOMN - L)*LNCRF + LCFNCOUNT 
                  FACKLMM = ONE
                  IF (SHELKEQL .AND. KMOMN .NE. LMOMN) FACKLMM=TWO
C
C Loop over the number of contracted functions for angular momentum 
C component: Begin                                         

                  FACIJKLMM = FACIJMM*FACKLMM*FACIJKL
C
C The following subroutine put the appropriate gamma elements into
C the GAMMA array from GAMMATMP. Also take care of the symmetry 
C factors. 
C 
                  CALL MODFY_GAMMA(GAMMA, GAMMATMP, LENGTH, IMOMN,
     &                             JMOMN, KMOMN, LMOMN, INCRF, JNCRF,
     &                             KNCRF, LNCRF, IANGMOM, JANGMOM, 
     &                             KANGMOM, LANGMOM)

                  INDX = ICFNCOUNT
C
                  DO ICRF = 1, INCRF
C
                     JNDX = JCFNCOUNT
                     INDX = INDX + 1 
C
                     JCRFFINAL = JNCRF
                     IF (ANGMOMIJ) JCRFFINAL = ICRF
C     
                     DO JCRF = 1, JCRFFINAL
C     
                     KNDX = KCFNCOUNT
                     JNDX = JNDX + 1 
                     FACIJCR=ONE
                     IF (ANGMOMIJ .AND. ICRF .NE. JCRF) FACIJCR=TWO
C
                        DO KCRF = 1, KNCRF
C
                           LNDX = LCFNCOUNT
                           KNDX = KNDX + 1 
C
                           LCRFFINAL = LNCRF
                           IF (ANGMOMKL) LCRFFINAL = KCRF
C
                           DO LCRF = 1, LCRFFINAL
C
                              LNDX = LNDX + 1 
                              FACKLCR=ONE
                              IF (ANGMOMKL .AND. KCRF .NE. LCRF)
     &                            FACKLCR=TWO
C                              
C$$$       WRITE(6,*) "INDX, JNDX, KNDX, LNDX=", INDX, JNDX, KNDX, LNDX
C$$$       WRITE(6,*) "ICRF, JCRF, KCRF, LCRF=", ICRF, JCRF, KCRF, LCRF
C     
                              FACIJKLCR= FACIJCR*FACKLCR
C
C$$$           Write(6,*) "The overall factor = ",HALF*FACIJKLCR*FACIJKLMM
C$$$       Write(6,*)
C$$$
C$$$                              TERM_AB = HALF*FACIJKLMM*FACIJKLCR*
C$$$     &                                  REPLINT(JNDX, LNDX)*
C$$$     &                                  PRDUTINT(INDX, KNDX)
C
                              TERM_AB = FACIJKLMM*FACIJKLCR*
     &                                  REPLINT(INDX, JNDX)*
     &                                  PRDUTINT(KNDX, LNDX)
C                              
                              IGAMMA = INDEX(ICRF, JCRF, KCRF, LCRF)
C
                              DHFIJ  = SCFDENT(INDX, JNDX)
                              DHFIK  = SCFDENT(INDX, KNDX)
                              DHFJK  = SCFDENT(JNDX, KNDX)
                              DHFIL  = SCFDENT(INDX, LNDX)
                              DHFJL  = SCFDENT(JNDX, LNDX)
                              DHFKL  = SCFDENT(KNDX, LNDX)
C
                              DRELIJ = RELDENT(INDX, JNDX)-
     &                                 SCFDENT(INDX, JNDX)
                              DRELIK = RELDENT(INDX, KNDX)-
     &                                 SCFDENT(INDX, KNDX)
                              DRELJK = RELDENT(JNDX, KNDX)-
     &                                 SCFDENT(JNDX, KNDX)
                              DRELIL = RELDENT(INDX, LNDX)-
     &                                 SCFDENT(INDX, LNDX)
                              DRELJL = RELDENT(JNDX, LNDX)-
     &                                 SCFDENT(JNDX, LNDX)
                              DRELKL = RELDENT(KNDX, LNDX)-
     &                                 SCFDENT(KNDX, LNDX)
C
C Built the repulsion, exchange, correlation contributions to the potential
C
                              REPLSN = DHFIJ*DHFKL
                              EXCHNG = -FOURTH*(DHFIK*DHFJL +
     &                                  DHFIL*DHFJK)
                              CORLNO = DHFIJ*DRELKL + DRELIJ*DHFKL
     &                                 -FOURTH*(DHFIK*DRELJL +
     &                                  DHFIL*DRELJK + DRELIK*DHFJL +
     &                                  DRELIL*DHFJK)
                              CORLNT =  FOUR*GAMMA(IGAMMA)
C
                              POTNL = REPLSN_FACT*REPLSN  +
     &                                EXCHNG_FACT*EXCHNG  +
     &                                CORLN_FACT_O*CORLNO +  
     &                                CORLN_FACT_T*CORLNT
C
                              POTNL_R = POTNL*TERM_AB + POTNL_R
C
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO
C
C Loop over the number of contracted functions end here!
C
                  LCFNCOUNT = LCFNCOUNT + LNCRF
              ENDDO
              KCFNCOUNT = KCFNCOUNT + KNCRF
           ENDDO
           JCFNCOUNT = JCFNCOUNT + JNCRF
         ENDDO
         ICFNCOUNT = ICFNCOUNT + INCRF
      ENDDO
C
C Loop over the number of angular momentum components end here!
C
      RETURN
      END
