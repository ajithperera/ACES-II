      SUBROUTINE REPPRDINT(MAXATMS, NATOMS, MAXSHELL, MAXPRM, IREDUN,
     &                     NSHL, NTOTPRIM, NTOTCRF, NOFFSETATMP,
     &                     NOFFSETATMC, NOFFSETPRM, NOFFSETCON, 
     &                     NANGMOMSHL, NPRIMFUNSHL, NCONFUNSHL, 
     &                     XPOINT, YPOINT, ZPOINT, ALPHA, PCOEF, 
     &                     COORD, REPLINT, PRDUTINT, KINTINT, TMP1,
     &                     TMP2, TMP3, OVERLAP_POTN)
C
C Evaluate integral of the form <Mu|1/(r1 - r2)|Nu> where r2 is 
C maintained as constant. Note the similarity to the standard 
C electron repulsion integral. Accordingly, this is built from
C routines that are used in the property integral package. 
C Ajith Perera 01/2000 
C     
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      DOUBLE PRECISION KINTINT(NTOTCRF, NTOTCRF)
      LOGICAL OVERLAP_POTN
C
      DIMENSION IP(20), FAC(9, 9),IOFFST(15), COORDI(3), COORDJ(3),
     &          COORD(NATOMS*3), ALPHA(NTOTPRIM), IREDUN(MAXATMS),
     &          PCOEF(NTOTPRIM*NTOTCRF), NANGMOMSHL(MAXATMS,MAXSHELL),
     &          NPRIMFUNSHL(MAXATMS,MAXSHELL),REPLINT(NTOTCRF,NTOTCRF),
     &          NCONFUNSHL(MAXATMS,MAXSHELL),PRDUTINT(NTOTCRF,NTOTCRF),
     &          CENTER(3), TMP1(MAXPRM, MAXPRM), TMP2(MAXPRM, MAXPRM),
     &          TMP3(MAXPRM, MAXPRM), NSHL(MAXATMS),
     &          NOFFSETATMP(MAXATMS), NOFFSETATMC(MAXATMS), 
     &          NOFFSETPRM(MAXATMS,MAXSHELL),
     &          NOFFSETCON(MAXATMS,MAXSHELL)
C      
C Do several important initilizations
C
      CALL SETRHF(FAC, IOFFST, IP)
C      
      IREDCNT = 0
      JREDCNT = 0
      IPRMCOUNT = 0
      JPRMCOUNT = 0
      ICFNCOUNT = 0
      JCFNCOUNT = 0
      CENTER(1) = XPOINT 
      CENTER(2) = YPOINT
      CENTER(3) = ZPOINT
C
C Loop over the symmetry unique atoms 
C      
      DO IATMS = 1, MAXATMS
C                  
         IPRMSTART =  NOFFSETATMP(IATMS)
         ICFNSTART =  NOFFSETATMC(IATMS)
         IREDATMS  =  IREDUN(IATMS)
C
C Loop over the redundent atoms 
C
         DO IRATMS = 1, IREDATMS         
C
C Loop over shells on symmetry unique atoms
C
            DO ISHL = 1, NSHL(IATMS)
C
               IREDCNT = IREDCNT + 1 
C
               DO JATMS = 1, MAXATMS
C                  
                  JPRMSTART =  NOFFSETATMP(JATMS)
                  JCFNSTART =  NOFFSETATMC(JATMS)
                  JREDATMS  =  IREDUN(JATMS)
C
                  DO JRATMS = 1, JREDATMS                 
C
                     DO JSHL = 1, NSHL(JATMS)
C
                        JREDCNT = JREDCNT + 1 
C                        
                        IANGMOM = NANGMOMSHL(IATMS, ISHL)
                        JANGMOM = NANGMOMSHL(JATMS, JSHL)
                        IPRIM   = NPRIMFUNSHL(IATMS, ISHL)
                        JPRIM   = NPRIMFUNSHL(JATMS, JSHL)
                        INCRF   = NCONFUNSHL(IATMS, ISHL)
                        JNCRF   = NCONFUNSHL(JATMS, JSHL)
                        IINTTYP = IOFFST(IANGMOM)
                        JINTTYP = IOFFST(JANGMOM)
                        
C$$$                        Write(6,*) "IREDCNT, JREDCNT", IREDCNT, JREDCNT
C$$$
C$$$                        Write(6,*) "IANGMOM, JANGMOM =", IANGMOM,
C$$$     &                               JANGMOM
C$$$                        Write(6,*) "IPRIM, JPRIM =", IPRIM, JPRIM
C$$$                        Write(6,*) "INCRF, JNCRF =", INCRF, JNCRF
C$$$                        Write(6,*) "ITYPE, JTYPE =", IINTTYP, JINTTYP
C$$$                        Write(6,*) "IATMS, JATMS =", IATMS, JATMS
C                        
                        ICOORD = (IATMS - 1)*3 
                        JCOORD = (JATMS - 1)*3 
C     
                        DO INDX = 1, 3
                           COORDI(INDX)  = COORD(ICOORD + INDX)
                           COORDJ(INDX)  = COORD(JCOORD + INDX)
                        ENDDO
C
C Loop over angular momentum components of each shell
C
                        IF (.NOT. ((IANGMOM*JANGMOM) .EQ. 0)) THEN
C
                           IPRMINTCNT = NOFFSETPRM(IATMS, ISHL) +
     &                                  IPRMSTART       
                           JPRMINTCNT = NOFFSETPRM(JATMS, JSHL) +
     &                                  JPRMSTART                           
                           ICFNINTCNT = NOFFSETCON(IATMS, ISHL) +
     &                                  ICFNSTART 
                           JCFNINTCNT = NOFFSETCON(JATMS, JSHL) +
     &                                  JCFNSTART 
C
                           IPRMCOUNT = IPRMINTCNT
                           ICFNCOUNT = ICFNINTCNT
                           ITYPE     = IINTTYP
C$$$                           I         = 1
C
                           DO IMOMN = 1, IANGMOM
                              ITYPE = ITYPE + 1 
C$$$                              IF (IMOMN .GT.2) I = I +  1
C$$$                              IPRMCOUNT = (IMOMN - I)*IPRIM +
C$$$     &                                     IPRMCOUNT 
C$$$                              ICFNCOUNT = (IMOMN - I)*INCRF +
C$$$     &                                     ICFNCOUNT
                              JPRMCOUNT = JPRMINTCNT
                              JCFNCOUNT = JCFNINTCNT
                              JTYPE     = JINTTYP
C$$$                              J = 1
C
                              DO JMOMN = 1, JANGMOM
                                 JTYPE = JTYPE + 1 
C$$$                                 IF (JMOMN .GT. 2) J = J +  1
C$$$                                 JPRMCOUNT = (JMOMN - J)*JPRIM + 
C$$$     &                                        JPRMCOUNT
C$$$                                 JCFNCOUNT = (JMOMN - J)*JNCRF +
C$$$     &                                        JCFNCOUNT
C
C Do the actual evaluation of the repulsion integral(A-Term in Eqn). Also bulit 
C the product function (B-Term in Eqn). 
C
C$$$                                 Write(6,*) "Loop over Ang-mom"
C$$$                                 Write(6,*) "IJPRMCUNT =", IPRMCOUNT, 
C$$$     &                                                     JPRMCOUNT
C$$$                                 Write(6,*) "IJCFNCUNT =", ICFNCOUNT,
C$$$     &                                                     JCFNCOUNT
 
                                 CALL EVALREPULI(IPRIM, JPRIM, INCRF, 
     &                                           JNCRF, ITYPE, JTYPE, 
     &                                           IPRMCOUNT, JPRMCOUNT,
     &                                           ICFNCOUNT, JCFNCOUNT,
     &                                           NTOTPRIM, NTOTCRF,
     &                                           MAXPRM, COORDI,COORDJ,
     &                                           ALPHA, PCOEF, CENTER,
     &                                           REPLINT,PRDUTINT,TMP1,
     &                                           TMP2, TMP3, FAC,
     &                                           OVERLAP_POTN)
C
C Do the actual evaluation of the kinetic energy and the nuclear repulsion 
C integarls. These two terms will contribute to the local representation 
C of the one-Hamiltonian contribution to the correlation energy. 
C 
                                 CALL EVALKINATIC(IPRIM, JPRIM, INCRF,
     &                                            JNCRF, ITYPE, JTYPE, 
     &                                            IPRMCOUNT, JPRMCOUNT,
     &                                            ICFNCOUNT, JCFNCOUNT,
     &                                            NTOTPRIM, NTOTCRF,
     &                                            MAXPRM, COORDI,COORDJ,
     &                                            ALPHA, PCOEF, CENTER,
     &                                            KINTINT,TMP1, TMP2,
     &                                            TMP3,FAC,NATOMS)
C
C Angular momentum loop end here!
C     
                                 JPRMCOUNT = JPRMCOUNT + JPRIM
                                 JCFNCOUNT = JCFNCOUNT + JNCRF
                              ENDDO
                              IPRMCOUNT = IPRMCOUNT + IPRIM
                              ICFNCOUNT = ICFNCOUNT + INCRF
                           ENDDO
C     
                        ENDIF
C
C Shell loop end here!
C

                     ENDDO
                  ENDDO
               ENDDO
               JREDCNT = 0
C
            ENDDO
         ENDDO
      ENDDO
C
C$$$      V_INT = 0.0D0
CSSS      Write(6,*) "The repulsion and product integral"
C$$$      DO I =1, NTOTCRF
C$$$         DO J= 1, I
C$$$            V_INT = V_INT + REPLINT(I,J)*REPLINT(I,J)
C$$$         ENDDO
C$$$      ENDDO
C$$$      Write(6,*) "The checksum of REPLINT =", V_INT
CSSS      CALL OUTPUT(REPLINT, 1, NTOTCRF, 1, NTOTCRF, NTOTCRF, NTOTCRF, 1)
C$$$      CALL OUTPUT(KINTINT, 1, NTOTCRF, 1, NTOTCRF, NTOTCRF, NTOTCRF, 1)
       
      RETURN
      END
 
