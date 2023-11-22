      SUBROUTINE DROPVC(NBAS,NDROP,IUHF,EVEC,IDROP,IPNTRA,IPNTRB,
     &                  IRPSVA,IRPSVB,EDROP)
C
C THIS ROUTINE GENERATES COMPRESSED MO POINTER LISTS, EIGENVECTORS AND
C  EIGENVALUES FOR POST-SCF CALCULATIONS IN WHICH CERTAIN MOLECULAR
C  ORBITALS WILL BE DROPPED.  THESE LISTS ARE CONSTRUCTED SUCH THAT THE
C  SYMMETRY REORDERING DONE IN ROUTINE GTNRB2 AND DESCENDENTS (PROGRAM
C  INTPROC) WILL PRODUCE THE CORRECT SYMMETRY-ORDERED LISTS AND
C  MO INTEGRAL LISTS.)  
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIRPRD
C
      DIMENSION EVEC(NBAS,NBAS),IDROP(2*NBAS),IPNTRA(NBAS)
      DIMENSION IPNTRB(NBAS),IRPSVA(NBAS),IRPSVB(NBAS),IRPDRP(8)
      DIMENSION EDROP(2*NBAS),NDRPOP(8),NDRVRT(8)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /FLAGS/  IFLAGS(100)
      COMMON /FLAGS2/ IFLAGS2(500)
      COMMON/SYMINF/NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
C
C FILL VECTOR IPNTR WITH ZEROS AND THEN PUT NONZERO VALUES IN
C  ALL POSITIONS WHICH CORRESPOND TO DROPPED MOs.  DO THIS FOR
C  BOTH A AND B BLOCKS.
C
      CALL IZERO(IDROP,2*NBAS)
      CALL IZERO(IPNTRA,NBAS)
      CALL IZERO(IPNTRB,NBAS)
      CALL ZERO(EVEC,NBAS*NBAS)
      NASCF=NOCCO(1)
      NBSCF=NOCCO(2)
      IONE=1
      ITWO=2
c     CALL PUTREC(20,'JOBARC','NOCCORB0',ITWO,NOCCO)
c     CALL PUTREC(20,'JOBARC','NVRTORB0',ITWO,NVRTO)
      NDROP=0
      CALL GETREC(20,'JOBARC','NUMDROPA',IONE,NDROP)
C
C PICK UP EIGENVALUES AND SORT THEM, SO THAT WE CAN DROP MOs
C  ACCORDING TO THEIR EIGENVALUE, RATHER THAN THEIR POSITION WITHIN
C  THE LIST OF SYMMETRY-BLOCKED SCF EIGENVALUES.
C
      CALL GETREC(20,'JOBARC','SCFEVLA0',NBAS*IINTFP,EVEC)
      CALL PUTREC(20,'JOBARC','SCFEVALA',NBAS*IINTFP,EVEC)
      if (iuhf.eq.0) 
     x   CALL PUTREC(20,'JOBARC','SCFEVALB',NBAS*IINTFP,EVEC)
      DO 303 I=1,NBAS
       IDROP(NBAS+I)=I 
303   CONTINUE
      CALL PIKSR2(NBAS,EVEC,IDROP(NBAS+1))
      IF (NDROP.NE.0) THEN
         CALL GETREC(20,'JOBARC','MODROPA ',NDROP,IDROP)
cmn
cmn In connection with fno: write original modropa0 record'
cmn
         CALL PUTREC(20,'JOBARC','MODROPA0',NDROP,IDROP)
      ENDIF

c-----------------------------------------------------------------------
c-----    test the input  ---- to trap a possible error ----- KB -------
c-----------------------------------------------------------------------
      if (ndrop.ne.0) then
c
       if (idrop(ndrop).gt.nbas) then
        write (6,9035) idrop(ndrop),nbas
 9035   format (////3x,'**** drop-mo indx',i4,' is larger than the',
     x          ' total number of mo',i4) 
        call errex 
       endif
c
CJDW 10/24/96.
C In analytical derivative calculations with dropped MOs, only lowest
C and/or highest can be dropped. This restriction does not apply to
C other calculations, for which any MOs can be dropped (in principle).
C Stop code from dying in ANY dropped MO calculation if the highest/lowest
C condition is not satisfied. And replace previous logic with something
C that may be easier to follow.
C
       IF(   (IFLAGS(47).NE.5 .AND. IFLAGS(47)  .NE.6
     &                        .AND. IFLAGS2(138).NE.1) .OR.
     &        IFLAGS(18).EQ.1                          .OR.
     &      ((IFLAGS(54).EQ.2  .OR. IFLAGS(54).EQ.3)   .AND.
     &       (IFLAGS2(138).EQ.0  .OR. 
     &        IFLAGS2(138).EQ.2)) )THEN
C
C Calculate number of dropped occupied orbitals and number of dropped
C virtual orbitals.
C
        NOCDRP=0
        NVRDRP=0
        DO 9040 I=1,NDROP
         IF(IDROP(I).LE.NOCCO(1))THEN
          NOCDRP = NOCDRP+1
         ELSE
          NVRDRP = NVRDRP+1
         ENDIF
 9040   CONTINUE
C
        IDRPERR = 0
C
C If there is one occupied orbital dropped, make sure it is the lowest one.
C
        IF(NOCDRP.EQ.1 .AND. IDROP(1).NE.1)THEN
         IDRPERR = 1
        ENDIF
C
C If there is one virtual orbital dropped, make sure it is the highest one.
C
        IF(NVRDRP.EQ.1 .AND. IDROP(NDROP).NE.NBAS)THEN
         IDRPERR = 2
        ENDIF
C
C If more than one occupied orbital is dropped, they must be 1,2,...,NOCDRP.
C
        IF(NOCDRP.GT.1)THEN
         IF(IDROP(1).NE.1)THEN
          IDRPERR = 3
         ENDIF
         DO 9050 I=2,NOCDRP
          IF( (IDROP(I) - IDROP(I-1)) .NE. 1 )THEN
           IDRPERR = 4
          ENDIF
 9050    CONTINUE
        ENDIF
C
C If more than one virtual orbital is dropped, they must be NBAS-NVRDRP+1,
C NBAS-NVRDRP+2,...,NBAS.
C
        IF(NVRDRP.GT.1)THEN
         IF(IDROP(NOCDRP+1).NE.NBAS-NVRDRP+1)THEN
          IDRPERR = 5
         ENDIF
         DO 9060 I=NOCDRP+2,NDROP
          IF( (IDROP(I) - IDROP(I-1)) .NE. 1 )THEN
           IDRPERR = 6
          ENDIF
 9060    CONTINUE
        ENDIF
C
        IF(IDRPERR.GT.0)THEN
         write (6,9036) 
         write (6,9031) nbas,ndrop,nocco(1)
         write (6,9032) (idrop(i),i=1,ndrop)
         call errex
        ENDIF
 9036  format (////3x,'**** only the lowest and/or highest',
     x                 ' mos can be dropped') 
 9031  format (4x,'# of MOs=',i5,10x,'# of drop-MOs=',i5,
     x        10x,'# of occupied alpha MOs=',i5)
 9032  format (/4x,'-- Your drop-MO indices =',20i3)
C
       ENDIF
      endif
c-----------------------------------------------------------------------
      CALL IZERO(IPNTRA,NBAS)
      DO 10 I=1,NDROP 
       NOFF=IDROP(I)
       INDXX=1
       IF(NOFF.LE.NASCF)INDXX=-1
       IDROP(I)=IDROP(NBAS+NOFF)
       IPNTRA(IDROP(NBAS+NOFF))=INDXX
10    CONTINUE
      IF(NDROP.NE.0)CALL PUTREC(20,'JOBARC','MODROPA ',NDROP,IDROP)
      CALL PUTREC(20,'JOBARC','IDROPA  ',NBAS,IPNTRA)
      IF(IUHF.EQ.1)THEN
       CALL GETREC(20,'JOBARC','SCFEVLB0',NBAS*IINTFP,EVEC)
       CALL PUTREC(20,'JOBARC','SCFEVALB',NBAS*IINTFP,EVEC)
      DO 304 I=1,NBAS
       IDROP(NBAS+I)=I 
304   CONTINUE
       CALL PIKSR2(NBAS,EVEC,IDROP(NBAS+1))
C$$$       CALL GETREC(20,'JOBARC','NUMDROPB',IONE,NDROP)
       IF (NDROP.NE.0) THEN 
          CALL GETREC(20,'JOBARC','MODROPB ',NDROP,IDROP)
cmn
cmn In connection with fno: write original modropa0 record'
cmn
          CALL PUTREC(20,'JOBARC','MODROPB0',NDROP,IDROP)
       ENDIF

       CALL IZERO(IPNTRB,NBAS)
       DO 20 I=1,NDROP 
        NOFF=IDROP(I)
        INDXX=1
        IF(NOFF.LE.NBSCF)INDXX=-1
        IDROP(I)=IDROP(NBAS+NOFF)
        IPNTRB(IDROP(I))=INDXX
20     CONTINUE
       IF(NDROP.NE.0)CALL PUTREC(20,'JOBARC','MODROPB ',NDROP,IDROP)
       CALL PUTREC(20,'JOBARC','IDROPB  ',NBAS,IPNTRB)
      ENDIF
C
C FIRST CLEAR UP THE OCCUPYA AND OCCUPYB VECTORS - THESE GIVE THE NUMBER
C  OF OCCUPIED ORBITALS IN A GIVEN SYMMETRY BLOCK.
C

c Nevin changed occupya0 to occupya for twodet to work
      CALL GETREC(20,'JOBARC','OCCUPYA ',NIRREP,IDROP)
      CALL GETREC(20,'JOBARC','IRREPALP',NBAS,IDROP(NBAS+1))
      CALL PUTREC(20,'JOBARC','OCCSCF  ',ITWO,NOCCO)
      CALL PUTREC(20,'JOBARC','VRTSCF  ',ITWO,NVRTO)
C-------------------------
      IF (NDROP.NE.0) then 
       CALL GETREC(20,'JOBARC','NAOBASFN',IONE,NNBAS) 
       CALL GETREC(20,'JOBARC','NUMBASIR',NIRREP,IDROP(NIRREP+1))
CJDW 7/1/96. I am resurrecting the next line, which was previously
C            commented out. This record will be read by GETAOINF.
C            Note that NUMBASIR gets modified in this routine and
C            is not sufficient in dropped core calculations (see
C            what games had to be played in PROPS).
       CALL PUTREC(20,'JOBARC','NUMBASI0',NIRREP,IDROP(NIRREP+1))
       IF (NBAS.EQ.NNBAS) THEN 
        CALL PUTREC(20,'JOBARC','FAOBASIR',NIRREP,IDROP(NIRREP+1))
       endif 
      endif 
C-------------------------
      ndrpopt = 0
      ndrvrtt = 0
      CALL PUTREC(20,'JOBARC','NDROTPOP',IONE,NDRPOPT) 
      CALL PUTREC(20,'JOBARC','NDROTVRT',IONE,NDRVRTT) 
      if (ndrop.ne.0) then 
      CALL IZERO(NDRPOP,NIRREP) 
      CALL IZERO(NDRVRT,NIRREP) 
      DO 50 I=1,NBAS
       IRROUT=IDROP(NBAS+I)
       IF(IPNTRA(I).NE.0.AND.I.LE.NASCF)THEN
        IDROP(IRROUT)=IDROP(IRROUT)-1
        NDRPOP(IRROUT) = NDRPOP(IRROUT) + 1
        NOCCO(1)=NOCCO(1)-1
       ENDIF
       IF(IPNTRA(I).NE.0.AND.I.GT.NASCF)THEN
        NDRVRT(IRROUT) = NDRVRT(IRROUT) + 1
       ENDIF
50    CONTINUE
      ndrpopt = 0
      ndrvrtt = 0
      do 3031 i=1,NIRREP
      ndrpopt = ndrpopt + ndrpop(i) 
      ndrvrtt = ndrvrtt + ndrvrt(i) 
 3031 continue 
      CALL PUTREC(20,'JOBARC','NDROPPOP',NIRREP,NDRPOP) 
      CALL PUTREC(20,'JOBARC','NDROPVRT',NIRREP,NDRVRT) 
      CALL PUTREC(20,'JOBARC','NDROTPOP',IONE,NDRPOPT) 
      CALL PUTREC(20,'JOBARC','NDROTVRT',IONE,NDRVRTT) 
C-------------------------
      CALL PUTREC(20,'JOBARC','OCCUPYA ',NIRREP,IDROP)
C-------------------------
       DO 501 I=1,NBAS
        IF(IPNTRA(I).NE.0) THEN
         IRROUT=IDROP(NBAS+I)
         IDROP(NIRREP+IRROUT)=IDROP(NIRREP+IRROUT)-1
        ENDIF
501    CONTINUE
       CALL PUTREC(20,'JOBARC','NUMBASIR',NIRREP,IDROP(NIRREP+1))
      endif 
c------------------
      NVRTO(1)=NBAS-NDROP-NOCCO(1)
      IF(IUHF.EQ.0)THEN
       NVRTO(2)=NBAS-NDROP-NOCCO(1)
       NOCCO(2)=NOCCO(1)
      ENDIF
      IF(IUHF.EQ.1)THEN

c Nevin changed occupyb0 to occupyb for twodet to work
       CALL GETREC(20,'JOBARC','OCCUPYB ',NIRREP,IDROP)
       CALL GETREC(20,'JOBARC','IRREPBET',NBAS,IDROP(NBAS+1))
       DO 51 I=1,NBAS
        IF(IPNTRB(I).NE.0.AND.I.LE.NBSCF)THEN
         IRROUT=IDROP(NBAS+I)
         IDROP(IRROUT)=IDROP(IRROUT)-1
         NOCCO(2)=NOCCO(2)-1
        ENDIF
51     CONTINUE
       CALL PUTREC(20,'JOBARC','OCCUPYB ',NIRREP,IDROP)
       NVRTO(2)=NBAS-NDROP-NOCCO(2)
      ENDIF
      CALL PUTREC(20,'JOBARC','NOCCORB ',ITWO,NOCCO)
      CALL PUTREC(20,'JOBARC','NVRTORB ',ITWO,NVRTO)
C
C PICK UP EIGENVECTORS AND STRIP OUT ALL COLUMNS CORRESPONDING
C  TO INACTIVE MOs.  THEN WRITE THEM BACK OUT.  
C
      CALL GETREC(20,'JOBARC','SCFEVCA0',NBAS*NBAS*IINTFP,
     &            EVEC)
      CALL PUTREC(20,'JOBARC','SCFEVECA',NBAS*NBAS*IINTFP,
     &            EVEC)
      ITHRU=0
      DO 15 I=1,NBAS
       IF(IPNTRA(I).EQ.0)THEN
        ITHRU=ITHRU+1
c YAU : old
c       CALL ICOPY(NBAS*IINTFP,EVEC(1,I),1,EVEC(1,ITHRU),1)
c YAU : new
        CALL DCOPY(NBAS,EVEC(1,I),1,EVEC(1,ITHRU),1)
c YAU : end
       ENDIF
15    CONTINUE
      CALL PUTREC(20,'JOBARC','SCFEVECA',ITHRU*NBAS*IINTFP,
     &            EVEC)
      NCOMPA=ITHRU
      CALL PUTREC(20,'JOBARC','NCOMPA  ',IONE,NCOMPA)
      IF(IUHF.EQ.1)THEN
       CALL GETREC(20,'JOBARC','SCFEVCB0',NBAS*NBAS*IINTFP,
     &             EVEC)
       CALL PUTREC(20,'JOBARC','SCFEVECB',NBAS*NBAS*IINTFP,
     &             EVEC)
       ITHRU=0
       DO 16 I=1,NBAS
        IF(IPNTRB(I).EQ.0)THEN
        ITHRU=ITHRU+1
c YAU : old
c        CALL ICOPY(NBAS*IINTFP,EVEC(1,I),1,EVEC(1,ITHRU),1)
c YAU : new
         CALL DCOPY(NBAS,EVEC(1,I),1,EVEC(1,ITHRU),1)
c YAU : end
        ENDIF
16     CONTINUE
       NCOMPB=ITHRU
       CALL PUTREC(20,'JOBARC','SCFEVECB',ITHRU*NBAS*IINTFP,
     &             EVEC)
       CALL PUTREC(20,'JOBARC','NCOMPB  ',IONE,NCOMPB)
      ENDIF
      CALL PUTREC(20,'JOBARC','NBASTOT ',IONE,NBAS)
      CALL PUTREC(20,'JOBARC','NBASCOMP',IONE,NCOMPA)
C
C NOW DO THE SAME THING FOR EIGENVALUES AND IRREP VECTORS.
C
      CALL GETREC(20,'JOBARC','SCFEVLA0',NBAS*IINTFP,EVEC)
      CALL GETREC(20,'JOBARC','IRREPALP',NBAS,IDROP)
c     CALL PUTREC(20,'JOBARC','IRREPAL0',NBAS,IDROP)
      DO 3210 I=1,NBAS
        IRPSVA(I)=IDROP(I)
 3210 CONTINUE
      IF(IFLAGS(1).GE.5)THEN
       IF(IUHF.EQ.0)WRITE(6,10000)(I,EVEC(I,1),I=1,NBAS)
       IF(IUHF.EQ.1)WRITE(6,10001)(I,EVEC(I,1),I=1,NBAS)
10000  FORMAT(T3,'SCF eigenvalues (a.u.):',/,(4('[',I3,'] ',F12.6)))
10001  FORMAT(T3,'SCF alpha eigenvalues (a.u.):',/,
     &             (4('[',I3,'] ',F12.6)))
      ENDIF
      ITHRU=0
      DO 25 I=1,NBAS
       IF(IPNTRA(I).EQ.0)THEN
        ITHRU=ITHRU+1
        EVEC(ITHRU,1)=EVEC(I,1)
        IDROP(ITHRU)=IDROP(I)
       ENDIF
25    CONTINUE
      NOUT=NCOMPA 
      CALL PUTREC(20,'JOBARC','SCFEVALA',NOUT*IINTFP,EVEC)
      CALL PUTREC(20,'JOBARC','IRREPALP',NOUT,IDROP)
      IF(IUHF.EQ.1)THEN
       CALL GETREC(20,'JOBARC','SCFEVLB0',NBAS*IINTFP,EVEC)
       CALL GETREC(20,'JOBARC','IRREPBET',NBAS,IDROP)
c      CALL PUTREC(20,'JOBARC','IRREPBE0',NBAS,IDROP)
       DO 3211 I=1,NBAS
         IRPSVB(I)=IDROP(I)
 3211  CONTINUE
       IF(IFLAGS(1).GE.5)THEN
        WRITE(6,10002)(I,EVEC(I,1),I=1,NBAS)
10002   FORMAT(T3,'SCF beta eigenvalues (a.u.):',/,
     &             (4('[',I3,'] ',F12.6)))
       ENDIF
       ITHRU=0    
       DO 35 I=1,NBAS
        IF(IPNTRB(I).EQ.0)THEN
         ITHRU=ITHRU+1
         EVEC(ITHRU,1)=EVEC(I,1)
         IDROP(ITHRU)=IDROP(I)
        ENDIF
35     CONTINUE
       NOUT=NCOMPB
c
      ENDIF
      CALL PUTREC(20,'JOBARC','SCFEVALB',NOUT*IINTFP,EVEC)
      CALL PUTREC(20,'JOBARC','IRREPBET',NOUT,IDROP)
C
C  NOW DO THE SAME FOR THE REORDERING VECTORS
C
      I8=8
      CALL IZERO(IDROP,2*NBAS) 
      CALL IZERO(IRPDRP,I8)
      CALL GETREC(20,'JOBARC','REORDERA',NBAS,IDROP)
      CALL PUTREC(20,'JOBARC','REORDEA0',NBAS,IDROP)
      NUMDRP=0
      DO 2010 I=1,NBAS
       DO 2011 J=1,NBAS
        IF(IDROP(J).EQ.I) THEN
         IF(IPNTRA(IDROP(J)).EQ.0) THEN
          IDROP(NBAS+J)=IDROP(J)-NUMDRP
          IDROP(J)=0
         ELSE
         DO 2012 K=IRPSVA(IDROP(J)),NIRREP
          IRPDRP(K)=IRPDRP(K)+1
 2012    CONTINUE
         NUMDRP=NUMDRP+1
         IDROP(J)=0
         ENDIF
        ENDIF
 2011  CONTINUE
 2010 CONTINUE
      CALL ICOPY(NBAS,IDROP(NBAS+1),1,IDROP,1)
      IPOS=0
      DO 2013 I=1,NBAS
       IF(IDROP(I).EQ.0)THEN
        IPOS=IPOS+1
       ELSE
        IDROP(NBAS+I-IPOS)=IDROP(I)
       ENDIF
2013  CONTINUE
      CALL PUTREC(20,'JOBARC','REORDERA',NBAS-NDROP,IDROP(NBAS+1))
c
      IF(IUHF.EQ.1) THEN
        CALL IZERO(IDROP,2*NBAS)
        CALL IZERO(IRPDRP,I8)
        CALL GETREC(20,'JOBARC','REORDERB',NBAS,IDROP)
c       CALL PUTREC(20,'JOBARC','REORDEB0',NBAS,IDROP)
        NUMDRP=0
        DO 2020 I=1,NBAS
          DO 2021 J=1,NBAS
            IF(IDROP(J).EQ.I) THEN
              IF(IPNTRB(IDROP(J)).EQ.0) THEN
                IDROP(NBAS+J)=IDROP(J)-NUMDRP
                IDROP(J)=0
              ELSE
                DO 2022 K=IRPSVB(IDROP(J)),NIRREP
                  IRPDRP(K)=IRPDRP(K)+1
 2022           CONTINUE
                NUMDRP=NUMDRP+1
                IDROP(J)=0
              ENDIF
            ENDIF
 2021     CONTINUE
 2020   CONTINUE
       CALL ICOPY(NBAS,IDROP(NBAS+1),1,IDROP,1)
       IPOS=0
       DO 2023 I=1,NBAS
        IF(IDROP(I).EQ.0)THEN
         IPOS=IPOS+1
        ELSE
         IDROP(NBAS+I-IPOS)=IDROP(I)
        ENDIF
2023   CONTINUE
       CALL PUTREC(20,'JOBARC','REORDERB',NBAS-NDROP,IDROP(NBAS+1))
      ENDIF
C
C
      CALL IZERO(IDROP,NBAS*2)
      CALL ZERO (EDROP,NBAS*2) 
      IF(IUHF.EQ.0)THEN
       IF(NCOMPA.NE.NBAS)THEN
        CALL GETREC (20,'JOBARC','SCFEVLA0',NBAS*IINTFP,EVEC)
        WRITE(6,1000)NBAS-NCOMPA
1000    FORMAT(/T3,'The following ',I3,' MOs will ',
     &            'be dropped:')
        ITHRU=0
        DO 100 I=1,NBAS
         IF(IPNTRA(I).NE.0)THEN
          ITHRU=ITHRU+1
          IDROP(ITHRU)=I
          EDROP(ITHRU)=EVEC(I,1) 
         ENDIF
100     CONTINUE
C---------------------------------------------
c       CALL PUTREC(20,'JOBARC','IIDROPA ',NBAS*2,IDROP) 
c       CALL PUTREC(20,'JOBARC','EEDROPA ',NBAS*2,EDROP) 
        DO 1015 I=1,ITHRU
         WRITE(6,1013) IDROP(I),EDROP(I)  
1013     FORMAT(3X,I3,5x,F15.7)
1015    CONTINUE
C---------------------------------------------
        WRITE(6,1002)NCOMPA
1002    FORMAT(T3,'There are ',I3,' active molecular ',
     &            'orbitals.'/)
       ENDIF
      ELSE
       IF(NCOMPA.NE.NBAS)THEN
        CALL GETREC (20,'JOBARC','SCFEVLA0',NBAS*IINTFP,EVEC)
        WRITE(6,2000)NBAS-NCOMPA
2000    FORMAT(/T3,'The following ',I3,' alpha MOs ',
     &            'will be dropped.')
        ITHRU=0
        DO 101 I=1,NBAS
         IF(IPNTRA(I).NE.0)THEN
          ITHRU=ITHRU+1
          IDROP(ITHRU)=I
          EDROP(ITHRU)=EVEC(I,1) 
         ENDIF
101     CONTINUE
c       CALL PUTREC(20,'JOBARC','IIDROPA ',NBAS*2,IDROP) 
c       CALL PUTREC(20,'JOBARC','EEDROPA ',NBAS*2,EDROP) 
        DO 1025 I=1,ITHRU
         WRITE(6,1013) IDROP(I),EDROP(I)  
1025    CONTINUE
        WRITE(6,1003)NCOMPA
1003    FORMAT(T3,'There are ',I3,' active alpha molecular ',
     &            'orbitals.')
       ENDIF
       IF(NCOMPB.NE.NBAS)THEN
        CALL GETREC (20,'JOBARC','SCFEVLB0',NBAS*IINTFP,EVEC)
        WRITE(6,6000)NBAS-NCOMPB
6000    FORMAT(/T3,'The following ',I3,' beta MOs ',
     &            'will be dropped.')
        ITHRU=0
        DO 102 I=1,NBAS
         IF(IPNTRB(I).NE.0)THEN
          ITHRU=ITHRU+1
          IDROP(ITHRU)=I
          EDROP(ITHRU)=EVEC(I,1) 
         ENDIF
102     CONTINUE
c       CALL PUTREC(20,'JOBARC','IIDROPB ',NBAS*2,IDROP) 
c       CALL PUTREC(20,'JOBARC','EEDROPB ',NBAS*2,EDROP) 
        DO 1035 I=1,ITHRU
         WRITE(6,1013) IDROP(I),EDROP(I)  
1035    CONTINUE
        WRITE(6,1004)NCOMPB
1004    FORMAT(T3,'There are ',I3,' active beta molecular ',
     &            'orbitals.'/)
       ENDIF
      ENDIF
c
      RETURN
      END
