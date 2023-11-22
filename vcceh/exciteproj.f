C
      SUBROUTINE EXCITEPROJ(SCR, MAXCOR, IRREPX, SWITCH, NSIZEC, IUHF)
C
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      INTEGER POP, VRT
      DIMENSION SCR(MAXCOR), IOFFPOP(8, 2), IOFFVRT(8, 2),
     &          ICORINDX(2), IVALINDX(2), IINVALINDX(2),
     &          IOFFSET(8, 3)

      CHARACTER*8 SWITCH
C
      COMMON /FILES/ LUOUT, MOINTS
      COMMON /INFO/ NOCCO(2), NVRTO(2)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYM/ POP(8,2),VRT(8,2),NTAA,NTBB,ND1AA,ND1BB,
     &            ND2AA,ND2BB
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
C
C Read the SCF eigenvalues from the "JOBARC' file.      
C
      NBASA = NOCCO(1) + NVRTO(1)
C
      IF (IUHF .EQ. 1) THEN
         NBASB = NOCCO(2) + NVRTO(2)
      ELSE
         NBASB = 0
      ENDIF

      I000 = 1
      I010 = I000 + NBASA
      I020 = I010 + NBASB
C
      CALL ZERO(SCR(I000), NSIZEC)
      CALL UPDATES2(IRREPX, SCR(I000), 444, 0, 490, IUHF)
C
      CALL GETREC(20, 'JOBARC', 'SCFEVALA', IINTFP*NBASA, SCR(I000))
C      
      IF (IUHF. NE. 0) THEN
         CALL GETREC(20, 'JOBARC', 'SCFEVALB', IINTFP*NBASB, SCR(I010))
      ENDIF
C
C Calculate the offsets for energy eigenvalues.
C
      IOFFPOP(1, 1) = I000
      IOFFPOP(1, 2) = I010
C
      DO 10 IRREP = 2, NIRREP
         DO 20 ISPIN = 1, 1 + IUHF
            IOFFPOP(IRREP, ISPIN) = IOFFPOP(IRREP - 1, ISPIN) +
     $                              POP(IRREP - 1, ISPIN)
 20      CONTINUE
 10   CONTINUE
C
      IOFFVRT(1, 1) = IOFFPOP(NIRREP, 1) + POP(NIRREP, 1)
C
      IF (IUHF .NE. 0) THEN
         IOFFVRT(1, 2) = IOFFPOP(NIRREP, 2) + POP(NIRREP, 2)
      ENDIF
C
      DO 25 IRREP = 2, NIRREP
         DO 30 ISPIN = 1, 1 + IUHF
            IOFFVRT(IRREP, ISPIN) = IOFFVRT(IRREP - 1, ISPIN) +
     $                              VRT(IRREP - 1, ISPIN)
 30      CONTINUE
 25   CONTINUE
C     
C Now we need to divide the orbitals into three different groups - core,
C valence, and the inner valence. Do that based on orbital energies.
C Orbitals energies less than -4.0 a.u. are considered as core 
C orbitals and the orbital energies greater than 1.0 a.u. are considered
C as valence. The rest is considered as inner valence.
C
      ECORE = -4.0D+00
      EVAL  = -0.90D+00
      
      DO 40 ISPIN = 1, IUHF + 1
         DO 50 IRREP = 1, NIRREP
            DO 60 IPOP = 1, POP(IRREP, ISPIN)
C
               IF (SCR(IOFFPOP(IRREP, ISPIN) + IPOP - 1)
     &             .LE. ECORE .AND. SWITCH .EQ. 'CORE    ')
     &              THEN
C
C The current orbital is a core orbital
C     
C                  EDEBUG=SCR(IOFFPOP(IRREP, ISPIN) + IPOP - 1)
C                  Write(6, 300) 'The core orbital energy   ', EDEBUG

                  ICORINDX(ISPIN) = IPOP
C     
                  IF (IUHF .EQ. 0) THEN
                     ICORINDX(3 - ISPIN) = IPOP
                  ELSE
                     ICORINDX(3 - ISPIN) = 0
                  ENDIF
C
                  CALL BLTPROJS(IUHF, SCR(I020), MAXCOR - I020 + 1, 
     &                          IRREPX, IRREP, ICORINDX)
C
               ELSE IF (SCR(IOFFPOP(IRREP, ISPIN) + IPOP - 1)
     &                  .GE. EVAL .AND. SWITCH .EQ. 'VALENCE ') 
     &                  THEN
C 
C The current orbital is valence orbital
C
C                  EDEBUG=SCR(IOFFPOP(IRREP, ISPIN) + IPOP - 1)
C                  Write(6, 300) 'The valence orbital energy', EDEBUG

                  IVALINDX(ISPIN) = IPOP
C
                  IF (IUHF .EQ. 0) THEN
                     IVALINDX(3 - ISPIN) = IPOP
                  ELSE
                     IVALINDX(3 - ISPIN) = 0
                  ENDIF
C
                  CALL BLTPROJS(IUHF, SCR(I020), MAXCOR - I020 + 1, 
     &                             IRREPX, IRREP, IVALINDX)
C                  
               ELSE IF (SCR(IOFFPOP(IRREP, ISPIN) + IPOP - 1)
     &                  .GT. ECORE .AND. SCR(IOFFPOP(IRREP, ISPIN)
     &                  + IPOP - 1) .LT. EVAL .AND. SWITCH .EQ. 
     &                  'INNERVAL') THEN
C     
C The current orbital is inner-valence
C
C                  EDEBUG=SCR(IOFFPOP(IRREP, ISPIN) + IPOP - 1)
C                  Write(6, 300) 'The Inneval orbital energy', EDEBUG
                  
                  IINVALINDX(ISPIN) = IPOP
C     
                  IF (IUHF .EQ. 0) THEN
                     IINVALINDX(3 - ISPIN) = IPOP
                  ELSE
                     IINVALINDX(3 - ISPIN) = 0
                  ENDIF
C     
                  CALL BLTPROJS(IUHF, SCR(I020), MAXCOR - I020 + 1, 
     &                          IRREPX, IRREP, IINVALINDX)
                  
               ENDIF
C
 60         CONTINUE
 50      CONTINUE
 40   CONTINUE
C
C
C Now do the double excitation projection vectors.
C
      DO 70 ICASE = 1, 1 + 2*IUHF
C     
         IF (ICASE .EQ. 1) THEN
            ILOW  = 1
            IHIGH = 2
         ELSEIF (ICASE .EQ. 2) THEN
            ILOW  = 1
            IHIGH = 1
         ELSEIF (ICASE .EQ. 3) THEN
            ILOW  = 2 
            IHIGH = 2 
         ENDIF
C
         CALL ZERO(SCR(I020), MAXCOR - I020 + 1)            
C
         DO 80 ISPIN = ILOW,  IHIGH
C
            DO 90 IRREP = 1, NIRREP
C
               DO 100 IPOP = 1, POP(IRREP, ISPIN)
C
                  IF (SCR(IOFFPOP(IRREP, ISPIN) + IPOP - 1)
     &                .LE. ECORE .AND. SWITCH .EQ. 'CORE    ')
     &               THEN
C
C The current orbital is a core orbital
C     
C                     EDEBUG=SCR(IOFFPOP(IRREP, ISPIN) + IPOP - 1)
C                     Write(6, 300) 'The core orbital energy   ', EDEBUG

                     ICORINDX(ISPIN) = IPOP
C     
                     IF (IUHF .EQ. 0) THEN
                        ICORINDX(3 - ISPIN) = IPOP
                     ELSE
                        ICORINDX(3 - ISPIN) = 0
                     ENDIF
C
                     CALL BLTPROJD(IUHF, SCR(I020), MAXCOR - I020 + 1, 
     &                             IRREPX, IRREP, ICORINDX, ICASE, 
     &                             IOFFSET)
C
                  ELSE IF (SCR(IOFFPOP(IRREP, ISPIN) + IPOP - 1)
     &                    .GE. EVAL .AND. SWITCH .EQ. 'VALENCE ') 
     &                  THEN
C 
C The current orbital is valence orbital
C     
C                     EDEBUG=SCR(IOFFPOP(IRREP, ISPIN) + IPOP - 1)
C                     Write(6, 300) 'The valence orbital energy', EDEBUG

                     IVALINDX(ISPIN) = IPOP
C
                     IF (IUHF .EQ. 0) THEN
                        IVALINDX(3 - ISPIN) = IPOP
                     ELSE
                        IVALINDX(3 - ISPIN) = 0
                     ENDIF
C
                     CALL BLTPROJD(IUHF, SCR(I020), MAXCOR - I020 + 1, 
     &                             IRREPX, IRREP, IVALINDX, ICASE,
     &                             IOFFSET)
C                  
                  ELSE IF (SCR(IOFFPOP(IRREP, ISPIN) + IPOP - 1)
     &                    .GT. ECORE .AND. SCR(IOFFPOP(IRREP, ISPIN)
     &                    + IPOP - 1) .LT. EVAL .AND. SWITCH .EQ. 
     &                   'INNERVAL') THEN
C     
C The current orbital is inner-valence
C
C                     EDEBUG=SCR(IOFFPOP(IRREP, ISPIN) + IPOP - 1)
C                     Write(6, 300) 'The Inneval orbital energy', EDEBUG
                  
                     IINVALINDX(ISPIN) = IPOP
C     
                     IF (IUHF .EQ. 0) THEN
                        IINVALINDX(3 - ISPIN) = IPOP
                     ELSE
                        IINVALINDX(3 - ISPIN) = 0
                     ENDIF
C     
                     CALL BLTPROJD(IUHF, SCR(I020), MAXCOR - I020 + 1, 
     &                             IRREPX, IRREP, IINVALINDX, ICASE,
     &                             IOFFSET)
                  ENDIF
C
 100           CONTINUE
 90         CONTINUE
 80      CONTINUE
C
         CALL WRTDBLPRJ(SCR(I020), MAXCOR - I020 + 1,  IRREPX, ICASE,
     &                  IOFFSET)
C
 70   CONTINUE
C         
C Exciation projection patterns are determined and put on lists 444-446
C and 490. Now normalize the projection vectors.
C 
      CALL LOADVEC2(IRREPX, SCR, MAXCOR, IUHF, 490, 0, 443, NSIZEC,
     &              .FALSE.)
      
      DO 110 I = 1, NSIZEC
C
          IF (SCR(I) .GT. 0.75) THEN
            SCR(I) = 1.0D0
         ELSE
            SCR(I) = 0.0D0
         ENDIF
C
 110  CONTINUE
C
C Put normalized projection vectors back on the list
C
      CALL UPDATES2(IRREPX, SCR(I000), 444, 0, 490, IUHF)
C
 300  FORMAT (5X, A26, 3X, F7.4)
      RETURN
      END
