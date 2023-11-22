      SUBROUTINE GETINV2U(INDXG,NTDIS,NRNUM,LISTW,IRREP,ISPIN)   
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER INDXG(NTDIS)
C
      COMMON /SYM2/NPOP(8,2),NVRT(8,2),NTF12(6)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),INDPD(8,8) 
      COMMON /SYMDROP/NDRPOP(8),NDRVRT(8)
      COMMON /SYMPOP/IRPDPD(8,22),ISYTYP(2,500),NTOT(18) 
C
      NBLOCK = IRPDPD(IRREP,ISYTYP(2,LISTW))  
      Nlengt = IRPDPD(IRREP,ISYTYP(1,LISTW))  
      NIJ = 0
      NIJ0 = 0
      NAB = 0
      NAB0 = 0
      DO 10 JRREP=1,NIRREP
       IJREP = INDPD(IRREP,JRREP)
       NBSI = NVRT(IJREP,ISPIN)
       NBSJ = NVRT(JRREP,ISPIN)
       NBSI0 = NBSI - NDRVRT(IJREP)
       NBSJ0 = NBSJ - NDRVRT(JRREP)
       if (jrrep.eq.ijrep) NIJ = NIJ + (NBSI*(NBSI-1))/2 
       if (jrrep.gt.ijrep) NIJ = NIJ + NBSI * NBSJ
       if (jrrep.eq.ijrep) NIJ0 = NIJ0 + (NBSI0*(NBSI0-1))/2 
       if (jrrep.gt.ijrep) NIJ0 = NIJ0 + NBSI0 * NBSJ0
       NBSA = NVRT(IJREP,ISPIN)
       NBSB = NVRT(JRREP,ISPIN)
       NBSA0 = NBSA - NDRVRT(IJREP)
       NBSB0 = NBSB - NDRVRT(JRREP)
       if (jrrep.eq.ijrep) NAB = NAB + (NBSA*(NBSA-1))/2 
       if (jrrep.gt.ijrep) NAB = NAB + NBSA * NBSB
       if (jrrep.eq.ijrep) NAB0 = NAB0 + (NBSA0*(NBSA0-1))/2 
       if (jrrep.gt.ijrep) NAB0 = NAB0 + NBSA0 * NBSB0
  10  CONTINUE
C
      IF (NIJ0.NE.NBLOCK .or. nab0.ne.nlengt .or. nij.ne.ntdis .or.
     x    nab.ne.nrnum)  THEN
       write(6,14) listw
  14   format (3x,'start of getinv2u with listw=',i5)
       WRITE (6,15) NIJ,NAB,nij0,nab0,nblock,nlengt
  15   FORMAT(//3X,'+++++  nij  nab  nij0  nab0  nblock nlengt',6i4) 
       WRITE (6,16) ntdis,nrnum
  16   FORMAT(3X,'+++++  ntdis  nrnum                         ',2i4) 
       call errex 
      ENDIF
C------------------------
      DO 200 I=1,NTDIS
       INDXG(I) = 0
 200  CONTINUE 
C------------------------
C---- MAKE DROP--INDX 
C------------------------
      NIJIN = NIJ 
      DO 700 KRREP=NIRREP,1,-1
       IKREP = INDPD(KRREP,IRREP)
       if (krrep.lt.ikrep) go to 678  
       NBSI = NVRT(IKREP,ISPIN) 
       NDI = NDRVRT(IKREP)
       NBSK = NVRT(KRREP,ISPIN) 
       NDK = NDRVRT(KRREP)
C
       if (nbsi.eq.0 .or. nbsk.eq.0) go to 678
C
       NIJF = NIJIN 
       if (NDK.EQ.0) GO TO 410
       NABIJ = NBSI * NDK 
       if (krrep.eq.ikrep) NABIJ = (NDK * (2*NBSI-NDK-1))/2 
       NIJI = NIJF
       NIJF = NIJF - NABIJ + 1 
       DO 400 I = NIJI,NIJF,-1
        INDXG(I) = 1
 400   CONTINUE
       NIJF = NIJF - 1
 410   CONTINUE
       IF (KRREP.EQ.IKREP) GO TO 610
       IF (NDI.EQ.0) GO TO 610
       NBSK0 = NBSK - NDK
       IF (NBSK0.GE.1) THEN
        NBSI0 = NBSI - NDI
        DO 600 I=NBSK0,1,-1
         DO 500 J=NDI,1,-1
          INDXG(NIJF) = 1
          NIJF = NIJF - 1
 500     CONTINUE
         NIJF = NIJF - NBSI0 
 600    CONTINUE
       ENDIF
 610   CONTINUE
       NIJAL = NBSI * NBSK
       IF (KRREP.EQ.IKREP) NIJAL = (NBSI*(NBSI-1))/2 
       NIJIN = NIJIN - NIJAL 
 678   CONTINUE
 700  CONTINUE
      IF (NIJIN.EQ.0) RETURN
       write (6,710) nijin   
 710   format (4x,'some error (getgv2u) = nijin=',i6,4x,'must be 0') 
       CALL ERREX 
      END
