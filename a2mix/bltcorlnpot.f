      SUBROUTINE BLTCORLNPOT(IUHF, NTOTSHL, NTOTCRF, NANGMOMTSHL, 
     &                       NCONFUNTSHL, NOFFSETSHL, ILEFT,
     &                       REPLINT, PRDUTINT, SCFDENT, RELDENT, 
     &                       SCFDEND, RELDEND, GAMMA, POTNL_R, 
     &                       REPLSN_FACT, EXCHNG_FACT, CORLN_FACT_O,
     &                       CORLN_FACT_T)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      PARAMETER (LUPAO = 25)
      CHARACTER*80 FNAME
C
      LOGICAL SHELIEQJ, SHELIEQK, SHELJEQL, SHELKEQL, ANGMOMIJ, 
     &        ANGMOMKL, SHELIJKL
C
      DIMENSION NANGMOMTSHL(NTOTSHL), NCONFUNTSHL(NTOTSHL), 
     &          REPLINT(NTOTCRF, NTOTCRF), PRDUTINT(NTOTCRF,NTOTCRF),
     &          GAMMA(ILEFT), NOFFSETSHL(NTOTSHL),
     &          SCFDENT(NTOTCRF, NTOTCRF), RELDENT(NTOTCRF, NTOTCRF),
     &          SCFDEND(NTOTCRF, NTOTCRF), RELDEND(NTOTCRF, NTOTCRF)
C
      DATA IZERO, IONE /0, 1/
C
C Rewind the AO two particle density matrix file
C
      CALL GFNAME('AOGAM   ',FNAME, ILENGTH) 
      OPEN(UNIT=LUPAO,STATUS='UNKNOWN',ACCESS='SEQUENTIAL',
     &     FILE=FNAME(1:ILENGTH),FORM='UNFORMATTED')
      REWIND LUPAO
C
      INIT = IZERO
C      
C Loop over unique shell quadruplets
C
      DO ISHL = 1, NTOTSHL
C         
         DO JSHL = 1, ISHL

            SHELIEQJ  = (ISHL .EQ. JSHL)
C            
            DO KSHL = 1, ISHL
C
               SHELIEQK = (ISHL .EQ. KSHL)
               LFINAL = KSHL
               IF (SHELIEQK) LFINAL=JSHL
C
               DO LSHL = 1, LFINAL
C                                          
                  SHELJEQL = (JSHL .EQ. LSHL)
                  SHELKEQL = (KSHL .EQ. LSHL)
C
C Angular momentum components and number of contracted functions.
C
C$$$                  IF (ISHL .EQ. 7 .and. JSHL .EQ. 7 .and. KSHL .eq. 7
C$$$     &                .and. LSHL .eq. 4) CALL EMPTY

                  IANGMOM = NANGMOMTSHL(ISHL)
                  JANGMOM = NANGMOMTSHL(JSHL)
                  KANGMOM = NANGMOMTSHL(KSHL)
                  LANGMOM = NANGMOMTSHL(LSHL)
                  INCRF   = NCONFUNTSHL(ISHL)
                  JNCRF   = NCONFUNTSHL(JSHL)
                  KNCRF   = NCONFUNTSHL(KSHL)
                  LNCRF   = NCONFUNTSHL(LSHL)
C
                  ANGMOMIJ = SHELIEQJ.AND.((IANGMOM+JANGMOM-2) .EQ. 0)
                  ANGMOMKL = SHELKEQL.AND.((KANGMOM+LANGMOM-2) .EQ. 0)
                  SHELIJKL = SHELIEQK .AND. (JSHL .EQ. LSHL)
C
C$$$      WRITE(6,*) "ISHL, JSHL, KSHL, LSHL=", ISHL, JSHL, KSHL, LSHL
C$$$      WRITE(6,*) "IC, JC, KC, LC=", INCRF, JNCRF, KNCRF, LNCRF
C
C Read the AO two particle density matrix for this shell quadruplet.
C     
                  ISIZE = INIT + IANGMOM*INCRF*JANGMOM*JNCRF*
     &                    KANGMOM*KNCRF*LANGMOM*LNCRF
C
                  IF (ISIZE .GT. ILEFT) CALL INSMEM("BLTCORLNPOT", 
     &                                               ISIZE, ILEFT)
C
                  READ(LUPAO) LENGTH, (GAMMA(INIT + I), I=1, LENGTH)
C$$$                  WRITE(6,*) (GAMMA(INIT + I), I= 1, LENGTH)
C
                  IBEGIN = NOFFSETSHL(ISHL)
                  JBEGIN = NOFFSETSHL(JSHL)
                  KBEGIN = NOFFSETSHL(KSHL)
                  LBEGIN = NOFFSETSHL(LSHL)
C
                 CALL EVALCORPOT(IANGMOM, JANGMOM, KANGMOM, LANGMOM, 
     &                           INCRF, JNCRF, KNCRF, LNCRF, IBEGIN,
     &                           JBEGIN, KBEGIN, LBEGIN, IEND, 
     &                           JEND, KEND, LEND, NTOTCRF, LENGTH,
     &                           SHELIEQJ, SHELKEQL, SHELIJKL, 
     &                           ANGMOMIJ, ANGMOMKL, REPLINT,
     &                           PRDUTINT, SCFDENT, RELDENT, 
     &                           SCFDEND, RELDEND, GAMMA(INIT+1),
     &                           GAMMA(ISIZE + 1), POTNL_R, 
     &                           REPLSN_FACT, EXCHNG_FACT, 
     &                           CORLN_FACT_O, CORLN_FACT_T, IUHF)
C
C Loop over shell quadrupleuts:End!
C
C$$$                 WRITE(6,*) "The Potential =", POTNL_R
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
      CLOSE (UNIT=LUPAO, STATUS="KEEP")
C
      RETURN
      END
