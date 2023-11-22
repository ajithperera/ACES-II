      SUBROUTINE GETGOV5TU(G,TMP,LBL0,NBL0,I2,IRREP,LISTG,ISPIN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION G(2),TMP(2) 
      COMMON /SYM2/NPOP(8,2),NVRT(8,2),NTF12(6)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),INDPD(8,8) 
      COMMON /SYMDROP/NDRPOP(8),NDRVRT(8)
      COMMON /SYMPOP/IRPDPD(8,22),ISYTYP(2,500),NTOT(18) 
C
      NBLOCK = IRPDPD(IRREP,ISYTYP(2,LISTG))  
      LBLOCK = IRPDPD(IRREP,ISYTYP(1,LISTG))  
      NIJ = 0
      NIJ0 = 0
      NAB = 0
      NAB0 = 0
      DO 10 JRREP=1,NIRREP
       IJREP = INDPD(IRREP,JRREP)
       NBSI = NVRT(IJREP,ISPIN)
       NBSJ = NPOP(JRREP,ISPIN)
       NBSI0 = NBSI - NDRVRT(IJREP)
       NBSJ0 = NBSJ - NDRPOP(JRREP)
       NIJ = NIJ + NBSI * NBSJ
       NIJ0 = NIJ0 + NBSI0 * NBSJ0
       NBSA = NVRT(IJREP,ISPIN)
       NBSB = NVRT(JRREP,ISPIN)
       NBSA0 = NBSA - NDRVRT(IJREP)
       NBSB0 = NBSB - NDRVRT(JRREP)
       if (jrrep.eq.ijrep) NAB = NAB + (NBSA * (NBSA-1))/2 
       if (jrrep.gt.ijrep) NAB = NAB + NBSA * NBSB
       if (jrrep.eq.ijrep) NAB0 = NAB0 + (NBSA0 * (NBSA0-1))/2 
       if (jrrep.gt.ijrep) NAB0 = NAB0 + NBSA0 * NBSB0
  10  CONTINUE
C
      IF (NIJ0.NE.NBLOCK .or. nab0.ne.lblock .or. nij.ne.nbl0 .or.
     x    nab.ne.lbl0)  THEN
      write(6,111) listg
 111  format (3x,'start of getgov5tu with listg=',i5,5x) 
       WRITE (6,15) NIJ,NAB,nij0,nab0,nblock,LBLOCK,NBL0,Lbl0 
  15   FORMAT(//3X,'+++++  nij  nab  nij0  nab0  nblock LBLOCK',6i4
     x         /3x,'       nbl0   lbl0                        ',2i4) 
       call errex 
      ENDIF
C------------------------
      isize = nij*nab
      call zero(G,isize) 
      CALL GETTRN(G,TMP,NAB0,NIJ0,I2,IRREP,LISTG)
C------------------------
C---- EXPAND IJ BLOCK
C------------------------
       nnn = nij
       nij = nab
       nab = nnn
      ININEW = NIJ * NAB + 1
      INIOLD = NIJ0 * NAB0 + 1
      DO 900 KRREP=NIRREP,1,-1
       IKREP = INDPD(KRREP,IRREP)
       if (krrep.lt.ikrep) go to 900 
       NBSI = NVRT(IKREP,ISPIN) 
       NDI = NDRVRT(IKREP)
       NBSK = NVRT(KRREP,ISPIN) 
       NDK = NDRVRT(KRREP)
C
       if (nbsi.eq.0 .or. nbsk.eq.0) go to 900
C
       NABIJ = NAB * NBSI * NDK 
       if(krrep.eq.ikrep) NABIJ = NAB * (NDK*(2*NBSI-NDK-1))/2 
       ININEW = ININEW - NABIJ  
       NBSK0 = NDK + 1
       if (krrep.eq.ikrep) NBSK0 = NDK + 2
       DO 800 NK=NBSK,NBSK0,-1
        NABNDK = NAB*NDI 
        if(krrep.eq.ikrep) NABNDK = 0  
        ININEW = ININEW - NABNDK 
        NBSI0 = NDI + 1
        NBSII = NBSI
        if (krrep.eq.ikrep) nbsii = nk - 1
c       if (nbsii.le.0) go to 800 
        DO 700 NI=NBSII,NBSI0,-1
C-----------------------
C-----   EXPAND AB-LIST
C-----------------------
         DO 300 JRREP=NIRREP,1,-1
          IBREP = INDPD(IRREP,JRREP) 
          NABA = NVRT(IBREP,ISPIN)
          NDA  = NDRVRT(IBREP)
          NABA0 = NABA - NDA 
          NABB = NPOP(JRREP,ISPIN)
          NDB  = NDRPOP(JRREP)
          NABB0 = NABB - NDB 
          if (naba.eq.0 .or. nabb.eq.0) go to 300
c
          INIOLD = INIOLD - NABA0 * NABB0 
          ININEW = ININEW - NABA * NABB 
          CALL GCOPY (G(INIOLD),TMP,NABA0*NABB0) 
          CALL REPOTT(TMP,G(ININEW),
     X                    NABA0,NABB0,NABA,NABB,NDA,NDB,3)
 300     CONTINUE
C-----------------------
 700    CONTINUE
 800   CONTINUE
C-----------------------
 900  CONTINUE
       if (ININEW.EQ.1  .AND. INIOLD.EQ.1)  RETURN
       write(6,950) ininew,iniold
 950   format (4x,'Some error (getgov5tu) ; ininew=',i3,5x,'iniold=',i3)
       call errex 
      END
