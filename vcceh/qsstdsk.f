C
      SUBROUTINE QSSTDSK(WIN, WOUT, NSIZE, POP1, POP2, POP3, POP4, 
     &                   IPOS, IRREPX, TYPE, LISTFROM, SYTYPL, SYTYPR, 
     &                   LHSEXP, LFTTYP)
C
C This routine resorts a W(PQ,RS) symmetry-packed quantity.
C The core requirements are somewhat less than SSTGEN, as only one
C copy of the full vector needs to be held in core.
C
C INPUT:
C         WIN - Scratch area holding individual distributions of the 
C               W(PQ,RS) quantity [must be dimensioned to maximum 
C               distribution size]
C        WOUT - The resorted quantity [Held in core]
C       NSIZE - The size of the W quantity 
C        POP1 - Population vector for P
C        POP2 - Population vector for Q
C        POP3 - Population vector for R
C        POP4 - Population vector for S
C      IRREPX - The overall symmetry of the list
C    LISTFROM - The list holding the PQRS quantity
C      SYTYPL - Left symmetry type of resorted quantity
C      SYTYPR - Right symmetry type of resorted quantity
C      LHSEXP - Logical flag set to .TRUE. if left-hand side needs to be expanded 
C               before resorting
C      LFTTYP - Left hand symmetry type [used only if LHSEXP = .TRUE. ]
C
C SCRATCH:
C
C        IPOS - A scratch vector of the maximum of the PS distribution length
C
C OUTPUT:
C        WOUT - The resorted P(PQ,RS) quantity.
C         
      IMPLICIT INTEGER (A-Z)
      CHARACTER*4 TYPE
      LOGICAL LHSEXP
      DOUBLE PRECISION WIN(NSIZE),WOUT(NSIZE)
      DIMENSION POP1(8),POP2(8),POP3(8),POP4(8),IPOS(100)
      DIMENSION INUMLFT(8),INUMRHT(8),IABSOFF(9)
      DIMENSION IPOSRHT(8,8),IPOSLFT(8,8)
C
      COMMON /SYMINF/ NSTART,NIRREP,IRRVEC(255,2),DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
C
      IF (TYPE .EQ. '1432') THEN
C
C            P Q R S -> P S R Q
C
C Compute size of PS and RQ distributions for each IRREP
C
         CALL IZERO(INUMLFT, 8)
         CALL IZERO(INUMRHT, 8)
C
         DO 100 IRREP = 1, NIRREP
            DO 110 IRR1 = 1, NIRREP
               IRR2 = DIRPRD(IRR1, IRREP)
               INUMLFT(IRREP) = INUMLFT(IRREP) + POP1(IRR2)*POP4(IRR1)
               INUMRHT(IRREP) = INUMRHT(IRREP) + POP3(IRR2)*POP2(IRR1)
 110        CONTINUE        
 100     CONTINUE
C
C Now compute absolute offsets to beginning of each irrep
C
         IABSOFF(1) = 1
C
         DO 120 IRRRHT = 1, NIRREP
            IRRLFT = DIRPRD(IRRRHT, IRREPX)
C
            IABSOFF(IRRRHT+1) = IABSOFF(IRRRHT) +
     &                          INUMLFT(IRRLFT)*INUMRHT(IRRRHT) 
            IPOSRHT(1, IRRRHT) = 0
            IPOSLFT(1, IRRLFT) = 0
C
            DO 125 IRREP1 = 1, NIRREP - 1
               IRREP2R = DIRPRD(IRREP1, IRRRHT)
               IRREP2L = DIRPRD(IRREP1, IRRLFT)
               IPOSRHT(IRREP1+1, IRRRHT) = IPOSRHT(IRREP1, IRRRHT) +
     &                                     POP3(IRREP2R)*POP2(IRREP1)*
     &                                     INUMLFT(IRRLFT)
               IPOSLFT(IRREP1+1, IRRLFT) = IPOSLFT(IRREP1, IRRLFT) +
     &                                     POP1(IRREP2L)*POP4(IRREP1)
 125        CONTINUE
 120     CONTINUE
C
         IOFFW = 1
C
C First loop over right-hand side irreps
C
         DO 10 IRRRHT = 1, NIRREP
            IDIS = 1
            IRRLFT = DIRPRD(IRRRHT, IRREPX) 
C
C Now loop over RS pairs
C
            DO 20 IRREP4 = 1, NIRREP
               IRREP3 = DIRPRD(IRREP4, IRRRHT)
               DO 30 INDEX4 =1, POP4(IRREP4)
                  DO 35 INDEX3 = 1, POP3(IRREP3)
C
C Read this distribution from the disk
C
                     CALL GETLST(WIN, IDIS, 1, 1, IRRRHT, LISTFROM)
C
                     IF (LHSEXP) THEN
                        DISSZ0 = IRPDPD(IRRRHT, ISYTYP(1, LISTFROM))
                        DISSZX = IRPDPD(IRRRHT, LFTTYP)
                        CALL SYMEXP2(IRRRHT, POP1, DISSZ0, DISSZX, 1,
     &                               WIN, WIN)
                     ENDIF
C
                     IDIS  = IDIS + 1
                     IOFFW = 1
C
C Now compute addresses for all distribution elements in target
C vector (pqrs -> psrq mapping) and scatter distribution into output
C array.
C
                     ITHRU = 1
C
                     DO 40 IRREP2 = 1, NIRREP
                        IRREP1  = DIRPRD(IRREP2, IRRLFT)
                        IRREPL0 = DIRPRD(IRREP1, IRREP4)
                        IRREPR0 = DIRPRD(IRREP3, IRREP2)
C
C Offset to this irrep of Q within this RQ irrep
C
                        IOFF0 = IABSOFF(IRREPR0) + IPOSRHT(IRREP2, 
     &                          IRREPR0)
C
                        DO 50 INDEX2 = 1, POP2(IRREP2)
C
C Offset to this Q within this irrep of Q within this RQ irrep
                           IOFF = IOFF0 + (INDEX2-1)*POP3(IRREP3)*
     &                            INUMLFT(IRREPL0)
C
C Offset to this R within
C
                           IOFF = IOFF + (INDEX3-1)*INUMLFT(IRREPL0)
C
C Offset to this irrep of S within this irrep of PS
C
                           IOFF = IOFF + IPOSLFT(IRREP4, IRREPL0)
C
C Offset to this particular value of S within this irrep of S
C within this irrep of S 
C
                           IOFF = IOFF + (INDEX4-1)*POP1(IRREP1) - 1 
C
                           DO 55 INDEX1 = 1, POP1(IRREP1)
                              IPOS(ITHRU) = IOFF + INDEX1
                              ITHRU = ITHRU + 1
 55                        CONTINUE
 50                     CONTINUE
 40                  CONTINUE
C
                     CALL SCATTER(ITHRU-1, WOUT, IPOS, WIN(IOFFW))
                     IOFFW = IOFFW + ITHRU - 1
C
 35               CONTINUE
 30            CONTINUE
 20         CONTINUE
 10      CONTINUE 
C
      ELSEIF (TYPE .EQ. '1324') THEN
C
C            P Q R S -> P S R Q
C
C Compute size of PS and RQ distributions for each irrep
C
         CALL IZERO(INUMLFT, 8)
         CALL IZERO(INUMRHT, 8)
C
         DO 1100 IRREP = 1, NIRREP
            DO 1110 IRR1 = 1, NIRREP
               IRR2 = DIRPRD(IRR1, IRREP)
               INUMLFT(IRREP) = INUMLFT(IRREP) + POP1(IRR2)*POP3(IRR1)
               INUMRHT(IRREP) = INUMRHT(IRREP) + POP2(IRR2)*POP4(IRR1)
 1110       CONTINUE        
 1100    CONTINUE
C
C Now compute absolute offsets to begining of each DPD irrep
C
         IABSOFF(1) = 1
C
         DO 1120 IRRRHT=1,NIRREP
            IRRLFT = DIRPRD(IRRRHT, IRREPX)
            IABSOFF(IRRRHT+1) = IABSOFF(IRRRHT)+
     &                          INUMLFT(IRRLFT)*INUMRHT(IRRRHT) 
            IPOSRHT(1, IRRRHT) = 0
            IPOSLFT(1, IRRLFT) = 0
C
            DO 1125 IRREP1 = 1, NIRREP - 1
               IRREP2R = DIRPRD(IRREP1, IRRRHT)
               IRREP2L = DIRPRD(IRREP1, IRRLFT)
               IPOSLFT(IRREP1+1, IRRLFT) = IPOSLFT(IRREP1, IRRLFT) +
     &                                     POP1(IRREP2L)*POP3(IRREP1)
               IPOSRHT(IRREP1+1, IRRRHT) = IPOSRHT(IRREP1,IRRRHT) +
     &                                     POP2(IRREP2R)*POP4(IRREP1)*
     &                                     INUMLFT(IRRLFT)
 1125       CONTINUE
 1120    CONTINUE
C
         IOFFW = 1
C
C First loop over right-hand side irreps
C
         DO 1010 IRRRHT = 1, NIRREP
            IDIS = 1
            IRRLFT = DIRPRD(IRRRHT, IRREPX) 
C
C Now loop over RS pairs
C
            DO 1020 IRREP4 =1, NIRREP
               IRREP3 = DIRPRD(IRREP4, IRRRHT)
               DO 1030 INDEX4 = 1, POP4(IRREP4)
                  DO 1035 INDEX3 = 1, POP3(IRREP3)
C
C Read this distribution from disk
C
                     CALL GETLST(WIN, IDIS, 1, 1, IRRRHT, LISTFROM)
C
                     IF(LHSEXP)THEN
                        DISSZ0 = IRPDPD(IRRRHT, ISYTYP(1, LISTFROM))
                        DISSZX = IRPDPD(IRRRHT, LFTTYP)
                        CALL SYMEXP2(IRRRHT, POP1, DISSZ0, DISSZX, 1, 
     &                               WIN, WIN)
                     ENDIF
C
                     IDIS  = IDIS + 1
                     IOFFW = 1
C
C Now compute adresses for all distribution elements in target
C vector (pqrs -> prqs mapping) and scatter distribution into
C output array
C
                     ITHRU = 1
                     DO 1040 IRREP2 = 1, NIRREP
                        IRREP1  = DIRPRD(IRREP2, IRRLFT)
                        IRREPL0 = DIRPRD(IRREP1, IRREP3)
                        IRREPR0 = DIRPRD(IRREP2, IRREP4)
C
C Offset to this irrep of S within this QS irrep
C
                        IOFF0  = IABSOFF(IRREPR0) + IPOSRHT(IRREP4,
     &                           IRREPR0)
                        DO 1050 INDEX2 = 1, POP2(IRREP2)
C
C Offset to this S within this irrep of S within this QS irrep
C
                           IOFF = IOFF0 + (INDEX4 - 1)*POP2(IRREP2)*
     &                            INUMLFT(IRREPL0)
C
C Offset to this Q within 
C
                           IOFF = IOFF + (INDEX2-1)*INUMLFT(IRREPL0)
C
C Offset to this irrep of R within this irrep of PR
C
                           IOFF = IOFF + IPOSLFT(IRREP3, IRREPL0)
C
C Offset to this particular value of R within this irrep of R 
C within this irrep of PR
C
                           IOFF = IOFF + (INDEX3-1)*POP1(IRREP1) - 1 
C
                           DO 1055 INDEX1 = 1, POP1(IRREP1)
                              IPOS(ITHRU) = IOFF + INDEX1
                              ITHRU = ITHRU + 1
 1055                      CONTINUE
 1050                   CONTINUE
 1040                CONTINUE
C
                     CALL SCATTER(ITHRU-1, WOUT, IPOS, WIN(IOFFW))
                     IOFFW = IOFFW + ITHRU - 1
C
 1035             CONTINUE
 1030          CONTINUE
 1020       CONTINUE
 1010    CONTINUE 
C     
      ELSEIF (TYPE .EQ. '1342') THEN
C
C            P Q R S -> P S Q R
C
C Compute size of PS and RQ distributions for each irrep
C
         CALL IZERO(INUMLFT, 8)
         CALL IZERO(INUMRHT, 8)
C
         DO 6100 IRREP = 1, NIRREP
            DO 6110 IRR1 = 1, NIRREP
               IRR2 = DIRPRD(IRR1, IRREP)
               INUMLFT(IRREP) = INUMLFT(IRREP) + POP1(IRR2)*POP3(IRR1)
               INUMRHT(IRREP) = INUMRHT(IRREP) + POP4(IRR2)*POP2(IRR1)
 6110       CONTINUE        
 6100    CONTINUE
C
C Now compute absolute offsets to beginning of each DPD irrep
C
         IABSOFF(1) = 1
C
         DO 6120 IRRRHT = 1, NIRREP
            IRRLFT = DIRPRD(IRRRHT, IRREPX)
            IABSOFF(IRRRHT + 1) = IABSOFF(IRRRHT) + 
     &                            INUMLFT(IRRLFT)*INUMRHT(IRRRHT) 
            IPOSRHT(1, IRRRHT) = 0
            IPOSLFT(1, IRRLFT) = 0
C
            DO 6125 IRREP1 = 1, NIRREP - 1
               IRREP2R = DIRPRD(IRREP1, IRRRHT)
               IRREP2L = DIRPRD(IRREP1, IRRLFT)
               IPOSLFT(IRREP1+1, IRRLFT) = IPOSLFT(IRREP1, IRRLFT) +
     &                                     POP1(IRREP2L)*POP3(IRREP1)
               IPOSRHT(IRREP1+1, IRRRHT) = IPOSRHT(IRREP1, IRRRHT) +
     &                                     POP4(IRREP2R)*POP2(IRREP1)*
     &                                     INUMLFT(IRRLFT)
 6125       CONTINUE
 6120    CONTINUE
C
         IOFFW = 1
C
C First loop over right-hand side irreps
C
         DO 6010 IRRRHT = 1, NIRREP
            IDIS = 1
            IRRLFT = DIRPRD(IRRRHT, IRREPX) 
C
C Now loop over RS pairs
C
            DO 6020 IRREP4 = 1, NIRREP
               IRREP3 = DIRPRD(IRREP4, IRRRHT)
               DO 6030 INDEX4 = 1, POP4(IRREP4)
                  DO 6035 INDEX3 = 1, POP3(IRREP3)
C
C Read this distribution from disk 
C
                     CALL GETLST(WIN, IDIS, 1, 1, IRRRHT, LISTFROM)
C
                     IF (LHSEXP) THEN
                        DISSZ0 = IRPDPD(IRRRHT, ISYTYP(1, LISTFROM))
                        DISSZX = IRPDPD(IRRRHT, LFTTYP)
                        CALL SYMEXP2(IRRRHT, POP1, DISSZ0, DISSZX, 1,
     &                               WIN, WIN)
                     ENDIF
C
                     IDIS = IDIS + 1
                     IOFFW = 1
C
C Now compute address for all distribution elements in target 
C vector (pqrs -> prsq mapping) and scatter distribution into 
C output array 
C
                     ITHRU = 1
                     DO 6040 IRREP2 = 1, NIRREP
                        IRREP1  = DIRPRD(IRREP2, IRRLFT)
                        IRREPL0 = DIRPRD(IRREP1, IRREP3)
                        IRREPR0 = DIRPRD(IRREP2, IRREP4)
C
C Offset to this irrep of Q within this SQ irrep
C
                        IOFF0 = IABSOFF(IRREPR0) +IPOSRHT(IRREP2, 
     &                          IRREPR0)
C
                        DO 6050 INDEX2 = 1, POP2(IRREP2)
C
C Offset to this Q within this irrep of Q within this SQ irrep
C
                           IOFF = IOFF0 + (INDEX2-1)*POP4(IRREP4)*
     &                            INUMLFT(IRREPL0)
C
C Offset to this S within 
C
                           IOFF = IOFF + (INDEX4-1)*INUMLFT(IRREPL0)
C
C Offset to this irrep of R within this irrep of R within this irrep of PR
C
                           IOFF = IOFF + IPOSLFT(IRREP3, IRREPL0)
C
C Offset to this particular value of R within this irrep of R
C within this irrep of PR
C
                           IOFF = IOFF + (INDEX3-1)*POP1(IRREP1) - 1 
C
                           DO 6055 INDEX1 = 1, POP1(IRREP1)
                              IPOS(ITHRU) = IOFF + INDEX1
                              ITHRU = ITHRU + 1
 6055                      CONTINUE
 6050                   CONTINUE
 6040                CONTINUE
C
                     CALL SCATTER(ITHRU-1, WOUT, IPOS, WIN(IOFFW))
                     IOFFW = IOFFW + ITHRU - 1
C
 6035             CONTINUE
 6030          CONTINUE
 6020       CONTINUE
 6010    CONTINUE 
C
      ELSEIF (TYPE .EQ. '1423') THEN
C
C            P Q R S -> P S R Q
C
C Compute size of PS and QR distributions for each irrep
C
         CALL IZERO(INUMLFT, 8)
         CALL IZERO(INUMRHT, 8)
C
         DO 2100 IRREP = 1, NIRREP
            DO 2110 IRR1 = 1, NIRREP
               IRR2 = DIRPRD(IRR1, IRREP)
               INUMLFT(IRREP) = INUMLFT(IRREP) + POP1(IRR2)*POP4(IRR1)
               INUMRHT(IRREP) = INUMRHT(IRREP) + POP2(IRR2)*POP3(IRR1)
 2110       CONTINUE        
 2100    CONTINUE
C
C Now compute absolute offsets to beginning of each DPD irrep
C
         IABSOFF(1) = 1
C
         DO 2120 IRRRHT = 1, NIRREP
            IRRLFT = DIRPRD(IRRRHT, IRREPX)
            IABSOFF(IRRRHT+1) = IABSOFF(IRRRHT) +
     &                          INUMLFT(IRRLFT)*INUMRHT(IRRRHT) 
            IPOSRHT(1, IRRRHT) = 0
            IPOSLFT(1, IRRLFT) = 0
C
            DO 2125 IRREP1 = 1, NIRREP - 1
               IRREP2R = DIRPRD(IRREP1, IRRRHT)
               IRREP2L = DIRPRD(IRREP1, IRRLFT)
C
               IPOSRHT(IRREP1+1, IRRRHT) = IPOSRHT(IRREP1, IRRRHT) +
     &                                     POP2(IRREP2R)*POP3(IRREP1)*
     &                                     INUMLFT(IRRLFT)
               IPOSLFT(IRREP1+1, IRRLFT) = IPOSLFT(IRREP1,IRRLFT)+
     &                                     POP1(IRREP2L)*POP4(IRREP1)
 2125       CONTINUE
 2120    CONTINUE
C
         IOFFW = 1
C
C First loop over right-hand side irreps  
C
         DO 2010 IRRRHT = 1, NIRREP
            IDIS = 1
            IRRLFT = DIRPRD(IRRRHT, IRREPX) 
C
C Now loop over RS pairs
C
            DO 2020 IRREP4 = 1, NIRREP
               IRREP3 = DIRPRD(IRREP4, IRRRHT)
               DO 2030 INDEX4 = 1, POP4(IRREP4)
                  DO 2035 INDEX3 = 1, POP3(IRREP3)
C
C Read this distribution from the disk
C
                     CALL GETLST(WIN, IDIS, 1, 1, IRRRHT, LISTFROM)
                     IF (LHSEXP) THEN
                        DISSZ0 = IRPDPD(IRRRHT, ISYTYP(1, LISTFROM))
                        DISSZX = IRPDPD(IRRRHT, LFTTYP)
                        CALL SYMEXP2(IRRRHT, POP1, DISSZ0, DISSZX, 1, 
     &                               WIN, WIN)
                     ENDIF
C
                     IDIS  = IDIS + 1
                     IOFFW = 1
C
C Now compute addresses for all distribution elements in target
C vector (pqrs -> psrq mapping) and scatter distribution into
C output array
C
                     ITHRU = 1
                     DO 2040 IRREP2 = 1, NIRREP
                        IRREP1  = DIRPRD(IRREP2, IRRLFT)
                        IRREPL0 = DIRPRD(IRREP1, IRREP4)
                        IRREPR0 = DIRPRD(IRREP3, IRREP2)
C
C Offsets to this irrep of R within this QR irrep
C
                        IOFF0 = IABSOFF(IRREPR0) + IPOSRHT(IRREP3,
     &                          IRREPR0)
C
C Offset to this R within this irrep of R within this QR irrep
C
                        IOFF0 = IOFF0 + (INDEX3-1)*POP2(IRREP2)*
     &                          INUMLFT(IRREPL0)
C
                        DO 2050 INDEX2 = 1, POP2(IRREP2)
C
C Offset to this Q within
C
                           IOFF = IOFF0 + (INDEX2-1)*INUMLFT(IRREPL0)
C
C Offset to this irrep of S within this irrep PS
C
                           IOFF = IOFF + IPOSLFT(IRREP4, IRREPL0)
C
C Offset to this particular value of S within this irrep of S
C within this irrep of PS
C
                           IOFF = IOFF + (INDEX4-1)*POP1(IRREP1) - 1 
                           DO 2055 INDEX1 = 1, POP1(IRREP1)
                              IPOS(ITHRU) = IOFF + INDEX1
                              ITHRU = ITHRU + 1
 2055                      CONTINUE
 2050                   CONTINUE
 2040                CONTINUE
C
                     CALL SCATTER(ITHRU-1, WOUT, IPOS, WIN(IOFFW))
                     IOFFW = IOFFW + ITHRU - 1
 2035             CONTINUE
 2030          CONTINUE
 2020       CONTINUE
 2010    CONTINUE 
C
      ELSEIF(TYPE .EQ. '3214') THEN
C
C            P Q R S -> R Q P S
C
C Compute size of PS and QR distributions for each irrep
C
         CALL IZERO(INUMLFT, 8)
         CALL IZERO(INUMRHT, 8)
C
         DO 3100 IRREP = 1, NIRREP
            DO 3110 IRR1 = 1, NIRREP
               IRR2 = DIRPRD(IRR1, IRREP)
               INUMRHT(IRREP) = INUMRHT(IRREP) + POP1(IRR2)*POP4(IRR1)
               INUMLFT(IRREP) = INUMLFT(IRREP) + POP3(IRR2)*POP2(IRR1)
 3110       CONTINUE        
 3100    CONTINUE
C
C Now compute absolute offsets beginning of each DPD irrep
C
         IABSOFF(1) = 1
C
       DO 3120 IRRRHT = 1, NIRREP
          IRRLFT = DIRPRD(IRRRHT, IRREPX)
          IABSOFF(IRRRHT + 1) = IABSOFF(IRRRHT)+
     &                          INUMLFT(IRRLFT)*INUMRHT(IRRRHT) 
          IPOSRHT(1, IRRRHT) = 0
          IPOSLFT(1, IRRLFT) = 0
C
          DO 3125 IRREP1 = 1, NIRREP - 1
             IRREP2R = DIRPRD(IRREP1, IRRRHT)
             IRREP2L = DIRPRD(IRREP1, IRRLFT)
C
             IPOSLFT(IRREP1+1, IRRLFT) = IPOSLFT(IRREP1, IRRLFT) +
     &                                   POP3(IRREP2L)*POP2(IRREP1)
C
             IPOSRHT(IRREP1+1, IRRRHT) = IPOSRHT(IRREP1,IRRRHT) +
     &                                   POP1(IRREP2R)*POP4(IRREP1)*
     &                                   INUMLFT(IRRLFT)
 3125     CONTINUE
 3120  CONTINUE
C
       IOFFW = 1
C
C First loop over right-hand side irreps  
C
       DO 3010 IRRRHT = 1, NIRREP
          IDIS   = 1
          IRRLFT = DIRPRD(IRRRHT, IRREPX) 
C
C Now loop over RS pairs
C
          DO 3020 IRREP4 = 1, NIRREP
             IRREP3 = DIRPRD(IRREP4, IRRRHT)
             DO 3030 INDEX4 = 1, POP4(IRREP4)
                DO 3035 INDEX3=1,POP3(IRREP3)
C
C Read this distribution from disk
C
                   CALL GETLST(WIN, IDIS, 1, 1, IRRRHT, LISTFROM)
                   IF(LHSEXP)THEN
                      DISSZ0 = IRPDPD(IRRRHT, ISYTYP(1, LISTFROM))
                      DISSZX = IRPDPD(IRRRHT, LFTTYP)
                      CALL SYMEXP2(IRRRHT, POP1, DISSZ0, DISSZX, 1,
     &                             WIN, WIN)
                   ENDIF
C
                   IDIS  = IDIS + 1
                   IOFFW = 1
C
C Now compute addrsses for all distribution elements in traget
C vector (pqrs -> rqps mapping) and scatter distribution into 
C output array
C
                   ITHRU = 1
                   DO 3040 IRREP2 = 1, NIRREP
                      IRREP1  = DIRPRD(IRREP2,IRRLFT)
                      IRREPL0 = DIRPRD(IRREP2,IRREP3)
                      IRREPR0 = DIRPRD(IRREP1, IRREP4)
C
C Offset to this irrep of S within this PS irrep
C
                      IOFF0 = IABSOFF(IRREPR0) + IPOSRHT(IRREP4,
     &                        IRREPR0)
C
C Offse to this S within this irrep of S within this RS irrep
C
                      IOFF0 = IOFF0 + (INDEX4-1)*POP1(IRREP1)*
     &                        INUMLFT(IRREPL0)
C
                      DO 3050 INDEX2 = 1, POP2(IRREP2)
                         DO 3055 INDEX1 = 1, POP1(IRREP1)
C
C Offset to this P within 
C
                            IOFF = IOFF0 + (INDEX1-1)*INUMLFT(IRREPL0)
C
C Offset to this irrep of Q within this irrep of RQ
C
                            IOFF = IOFF + IPOSLFT(IRREP2, IRREPL0)
C
C Offset to this particular value Q within this irrep of Q 
C within this irrep of RQ 
C
                            IOFF = IOFF+(INDEX2-1)*POP3(IRREP3) - 1 
                            IPOS(ITHRU) = IOFF + INDEX3
                            ITHRU = ITHRU + 1
C
 3055                    CONTINUE
 3050                 CONTINUE
 3040              CONTINUE
C
                   CALL SCATTER(ITHRU-1, WOUT, IPOS, WIN(IOFFW))
                   IOFFW = IOFFW + ITHRU - 1
C
 3035           CONTINUE
 3030        CONTINUE
 3020     CONTINUE
 3010  CONTINUE 
C     
      ELSEIF (TYPE .EQ. '2413') THEN
C
C            P Q R S -> Q S P R
C            1 2 3 4    2 4 1 3
C
C Compute size of QS and PR distributions for each irrep
C
         CALL IZERO(INUMLFT, 8)
         CALL IZERO(INUMRHT, 8)
C
         DO 4100 IRREP = 1, NIRREP
            DO 4110 IRR1 = 1, NIRREP
               IRR2 = DIRPRD(IRR1, IRREP)
               INUMRHT(IRREP) = INUMRHT(IRREP) + POP1(IRR2)*POP3(IRR1)
               INUMLFT(IRREP) = INUMLFT(IRREP) + POP2(IRR2)*POP4(IRR1)
 4110       CONTINUE        
 4100    CONTINUE
C
C Now compute absolute offsets to beginning of each DPD irrep
C
         IABSOFF(1) = 1
C
         DO 4120 IRRRHT = 1, NIRREP
            IRRLFT = DIRPRD(IRRRHT, IRREPX)
            IABSOFF(IRRRHT+1) = IABSOFF(IRRRHT)+
     &                          INUMLFT(IRRLFT)*INUMRHT(IRRRHT) 
            IPOSRHT(1, IRRRHT) = 0
            IPOSLFT(1, IRRLFT) = 0
C
            DO 4125 IRREP1 = 1, NIRREP - 1
               IRREP2R = DIRPRD(IRREP1, IRRRHT)
               IRREP2L = DIRPRD(IRREP1, IRRLFT)
               IPOSLFT(IRREP1+1, IRRLFT) = IPOSLFT(IRREP1, IRRLFT) +
     &                                     POP2(IRREP2L)*POP4(IRREP1)
               IPOSRHT(IRREP1+1, IRRRHT) = IPOSRHT(IRREP1, IRRRHT) +
     &                                     POP1(IRREP2R)*POP3(IRREP1)*
     &                                     INUMLFT(IRRLFT)
 4125       CONTINUE
 4120    CONTINUE
C
         IOFFW = 1
C
C First loop over right-hand side irreps
C
         DO 4010 IRRRHT = 1, NIRREP
            IDIS = 1
            IRRLFT = DIRPRD(IRRRHT, IRREPX) 
C
C Now loop over PR pairs
C
C            P Q R S -> Q S P R
C            1 2 3 4    2 4 1 3
C
            DO 4020 IRREP4 = 1, NIRREP
               IRREP3 = DIRPRD(IRREP4, IRRRHT)
               DO 4030 INDEX4 = 1, POP4(IRREP4)
                  DO 4035 INDEX3 = 1, POP3(IRREP3)
C
C Now compute addresses for all distribution elements in targets
C vector (pqrs -> qspr mapping) and scatter distribution in to 
C output array.
C
                     ITHRU = 1
                     DO 4040 IRREP2 = 1, NIRREP
                        IRREP1  = DIRPRD(IRREP2, IRRLFT)
                        IRREPL0 = DIRPRD(IRREP2, IRREP4)
                        IRREPR0 = DIRPRD(IRREP1, IRREP3)
C
C Offset to this irrep of R within this PR irrep
C
                        IOFF0 = IABSOFF(IRREPR0) + IPOSRHT(IRREP3, 
     &                          IRREPR0)
C
C Offset to this R within this irrep of R within this PR irrep
C
                        IOFF0 = IOFF0 + (INDEX3-1)*POP1(IRREP1)*
     &                          INUMLFT(IRREPL0)
C
                        DO 4050 INDEX2 = 1, POP2(IRREP2)
                           DO 4055 INDEX1 = 1, POP1(IRREP1)
C
C Offset to this P within 
C
                              IOFF = IOFF0 + (INDEX1-1)*
     &                               INUMLFT(IRREPL0)
C
C Offset to this irrep of S within this irrep  QS
C
                              IOFF = IOFF + IPOSLFT(IRREP4, IRREPL0)
C
C offset to this particular value of S within this irrep of S
C within this irrep of QS
C
                              IOFF = IOFF + (INDEX4-1)*POP2(IRREP2) - 1 
                              IPOS(ITHRU) = IOFF + INDEX2
                              ITHRU = ITHRU + 1
 4055                      CONTINUE
 4050                   CONTINUE
 4040                CONTINUE
C
                     CALL SCATTER(ITHRU-1, WOUT, IPOS, WIN(IOFFW))
                     IOFFW = IOFFW + ITHRU - 1

C
 4035             CONTINUE
 4030          CONTINUE
 4020       CONTINUE
 4010    CONTINUE 
C     
      ELSEIF (TYPE .EQ. '2314') THEN
C
C            P Q R S -> Q R P S
C
C Compute size of PS and QR distributions for each irrep
C
         CALL IZERO(INUMLFT, 8)
         CALL IZERO(INUMRHT, 8)
C
         DO 5100 IRREP = 1, NIRREP
            DO 5110 IRR1 = 1, NIRREP
               IRR2 = DIRPRD(IRR1,IRREP)
               INUMRHT(IRREP) = INUMRHT(IRREP) + POP1(IRR2)*POP4(IRR1)
               INUMLFT(IRREP) = INUMLFT(IRREP) + POP2(IRR2)*POP3(IRR1)
 5110       CONTINUE        
 5100    CONTINUE
C
C Now compute absolute offsets to beginning of each DPDP irrep
C
         IABSOFF(1) = 1
C
       DO 5120 IRRRHT = 1, NIRREP
          IRRLFT = DIRPRD(IRRRHT,IRREPX)
          IABSOFF(IRRRHT+1) = IABSOFF(IRRRHT) +
     &                        INUMLFT(IRRLFT)*INUMRHT(IRRRHT) 
          IPOSRHT(1, IRRRHT) = 0
          IPOSLFT(1, IRRLFT) = 0
C
          DO 5125 IRREP1 = 1, NIRREP - 1
             IRREP2R = DIRPRD(IRREP1, IRRRHT)
             IRREP2L = DIRPRD(IRREP1, IRRLFT)
             IPOSLFT(IRREP1+1, IRRRHT) = IPOSLFT(IRREP1, IRRRHT) +
     &                                   POP2(IRREP2R)*POP3(IRREP1)
             IPOSRHT(IRREP1+1, IRRLFT) = IPOSRHT(IRREP1,IRRLFT) +
     &                                   POP1(IRREP2L)*POP4(IRREP1)*
     &                                   INUMLFT(IRRLFT)
 5125     CONTINUE
 5120  CONTINUE
C
       IOFFW = 1
C
C First loop over right-hand side irreps
C
       DO 5010 IRRRHT = 1, NIRREP
          IDIS = 1
          IRRLFT = DIRPRD(IRRRHT, IRREPX) 
C
C Now loop over RS pairs
C
          DO 5020 IRREP4 = 1, NIRREP
             IRREP3 = DIRPRD(IRREP4, IRRRHT)
             DO 5030 INDEX4 = 1, POP4(IRREP4)
                DO 5035 INDEX3 = 1, POP3(IRREP3)
C
C Read this distribution from disk
C
                   CALL GETLST(WIN, IDIS, 1, 1, IRRRHT, LISTFROM)
                   IF(LHSEXP)THEN
                      DISSZ0 = IRPDPD(IRRRHT, ISYTYP(1, LISTFROM))
                      DISSZX = IRPDPD(IRRRHT, LFTTYP)
                      CALL SYMEXP2(IRRRHT, POP1, DISSZ0, DISSZX, 1,
     &                             WIN, WIN)
                   ENDIF
C
                   IDIS = IDIS + 1
                   IOFFW = 1
C
C Now compute addresses for all distribution elements in target
C vector (pqrs -> qrps mapping) and scatter distribution into output
C array 1234->2314
C
                   ITHRU = 1
C
                   DO 5040 IRREP2 = 1, NIRREP
                      IRREP1  = DIRPRD(IRREP2, IRRLFT)
                      IRREPL0 = DIRPRD(IRREP2, IRREP3)
                      IRREPR0 = DIRPRD(IRREP1, IRREP4)
C
C Offset to this irrep of S within this PS irrep
C
                      IOFF0 = IABSOFF(IRREPR0) + IPOSRHT(IRREP4,
     &                        IRREPR0)
C
C Offset to this S within this irrep of S within this RS irrep
C
                      IOFF0 = IOFF0 + (INDEX4-1)*POP1(IRREP1)*
     &                        INUMLFT(IRREPL0)
C
                      DO 5050 INDEX2 = 1, POP2(IRREP2)
                         DO 5055 INDEX1 = 1, POP1(IRREP1)
C
C Offset to this P within 
C
                            IOFF = IOFF0 + (INDEX1-1)*INUMLFT(IRREPL0)
C
C Offset to this irrep of R within this irrep of QR
C
                            IOFF = IOFF + IPOSLFT(IRREP3, IRREPL0)
C
C Offset to this particular value of R within this irrep of R
C within this irrep of QR
C
                            IOFF = IOFF + (INDEX3-1)*POP2(IRREP2) - 1 
                            IPOS(ITHRU) = IOFF + INDEX2
                            ITHRU = ITHRU + 1
 5055                    CONTINUE
 5050                 CONTINUE
 5040              CONTINUE
C
                   CALL SCATTER(ITHRU-1, WOUT, IPOS, WIN(IOFFW))
                   IOFFW = IOFFW + ITHRU - 1
 5035           CONTINUE
 5030        CONTINUE
 5020     CONTINUE
 5010  CONTINUE 
C
      ELSEIF (TYPE .EQ. '2431') THEN
C
C            P Q R S -> Q S R P
C            1 2 3 4    2 4 3 1
C
C Compute sizes of QS and PR distributions for each irrep
C
         CALL IZERO(INUMLFT, 8)
         CALL IZERO(INUMRHT, 8)
C
         DO 7100 IRREP = 1, NIRREP
            DO 7110 IRR1 = 1, NIRREP
               IRR2 = DIRPRD(IRR1, IRREP)
               INUMRHT(IRREP) = INUMRHT(IRREP) + POP3(IRR2)*POP1(IRR1)
               INUMLFT(IRREP) = INUMLFT(IRREP) + POP2(IRR2)*POP4(IRR1)
 7110       CONTINUE        
 7100    CONTINUE
C
C Now compute absolute offsets to begnning of each DPD irrep
C
         IABSOFF(1) = 1
C
         DO 7120 IRRRHT = 1, NIRREP
            IRRLFT = DIRPRD(IRRRHT, IRREPX)
            IABSOFF(IRRRHT+1) = IABSOFF(IRRRHT) + 
     &                          INUMLFT(IRRLFT)*INUMRHT(IRRRHT) 
C
            IPOSRHT(1,IRRRHT) = 0
            IPOSLFT(1,IRRLFT) = 0
C
            DO 7125 IRREP1 = 1, NIRREP - 1
               IRREP2R = DIRPRD(IRREP1, IRRRHT)
               IRREP2L = DIRPRD(IRREP1, IRRLFT)
               IPOSLFT(IRREP1+1, IRRRHT) = IPOSLFT(IRREP1, IRRLFT) +
     &                                     POP2(IRREP2L)*POP4(IRREP1)
               IPOSRHT(IRREP1+1, IRRLFT) = IPOSRHT(IRREP1, IRRRHT) +
     &                                     POP3(IRREP2R)*POP1(IRREP1)*
     &                                     INUMLFT(IRRLFT)
 7125       CONTINUE
 7120    CONTINUE
C
         IOFFW = 1
C
C First loop over right-hand side irreps
C
         DO 7010 IRRRHT = 1, NIRREP
            IDIS = 1
            IRRLFT = DIRPRD(IRRRHT, IRREPX) 
C
C Now loop over PR pairs 
C            P Q R S -> Q S R P
C            1 2 3 4    2 4 3 1
C
            DO 7020 IRREP4 = 1, NIRREP
               IRREP3 = DIRPRD(IRREP4, IRRRHT)
               DO 7030 INDEX4 = 1, POP4(IRREP4)
                  DO 7035 INDEX3=1,POP3(IRREP3)
C
C Read this distribution from disk
C
                     CALL GETLST(WIN, IDIS, 1, 1, IRRRHT, LISTFROM)
                     IF (LHSEXP) THEN
                        DISSZ0 = IRPDPD(IRRRHT, ISYTYP(1, LISTFROM))
                        DISSZX = IRPDPD(IRRRHT, LFTTYP)
                        CALL SYMEXP2(IRRRHT, POP1, DISSZ0, DISSZX, 1,
     &                               WIN, WIN)
                     ENDIF
C
                     IDIS = IDIS + 1
                     IOFFW = 1
C
C Now compute addresses for all distribution elements in target
C vector (pqrs -> qsrp mapping) and scatter distribution into
C output array
C
                     ITHRU = 1
C
                     DO 7040 IRREP2 = 1, NIRREP
                        IRREP1  = DIRPRD(IRREP2, IRRLFT)
                        IRREPL0 = DIRPRD(IRREP2, IRREP4)
                        IRREPR0 = DIRPRD(IRREP1, IRREP3)
C
C Offset to this irrep of P within this RP irrep
C
                        IOFF0 = IABSOFF(IRREPR0) + IPOSRHT(IRREP1, 
     &                          IRREPR0)
                        DO 7050 INDEX2 = 1, POP2(IRREP2)
                           DO 7055 INDEX1 = 1, POP1(IRREP1)
C
C Offset to this P within this irrep of P within this RP irrep
C
                              IOFF = IOFF0 + (INDEX1-1)*POP3(IRREP3)*
     &                               INUMLFT(IRREPL0)
C
C Offset to this R within 
C
                              IOFF = IOFF + (INDEX3-1)*INUMLFT(IRREPL0)
C
C Offset to this irrep of S within this irrep of QS
C
                              IOFF = IOFF + IPOSLFT(IRREP4, IRREPL0)
C
C Offset to this particular value of S within this irrep of S 
C within this irrep of QS
C
                              IOFF = IOFF + (INDEX4-1)*POP2(IRREP2) - 1 
                              IPOS(ITHRU) = IOFF + INDEX2
                              ITHRU = ITHRU + 1
 7055                      CONTINUE
 7050                   CONTINUE
 7040                CONTINUE
C     
                     CALL SCATTER(ITHRU-1, WOUT, IPOS, WIN(IOFFW))
                     IOFFW = IOFFW + ITHRU - 1
 7035             CONTINUE
 7030          CONTINUE
 7020       CONTINUE
 7010    CONTINUE 
C
      ENDIF
C
      RETURN
      END
