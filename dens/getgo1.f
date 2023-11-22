      SUBROUTINE GETGO1 (G,GTMP,I1,NBL0,I2,IRREP,LISTG,ISPIN,listw,leng)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION G(2),GTMP(2)
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
       NBSI = NPOP(IJREP,1)
       NBSJ = NPOP(JRREP,2)
       NBSI0 = NBSI - NDRPOP(IJREP)
       NBSJ0 = NBSJ - NDRPOP(JRREP)
       NIJ = NIJ + NBSI * NBSJ
       NIJ0 = NIJ0 + NBSI0 * NBSJ0
       NBSA = NVRT(IJREP,1)
       NBSB = NVRT(JRREP,2)
       NBSA0 = NBSA - NDRVRT(IJREP)
       NBSB0 = NBSB - NDRVRT(JRREP)
       NAB = NAB + NBSA * NBSB
       NAB0 = NAB0 + NBSA0 * NBSB0
  10  CONTINUE
c
      IF (NIJ0.NE.NBLOCK .or. nab0.ne.nlengt .or. nij.ne.nbl0 .or.
     x    nab.ne.leng)  THEN
      write(6,111) listg,listw
 111  format (3x,'start of getgo1 with listg=',i5,5x,'listw=',i5)
       WRITE (6,15) NIJ,NAB,nij0,nab0,nblock,nlengt,nbl0,leng 
  15   FORMAT(//3X,'+++++  nij  nab  nij0  nab0  nblock nlengt',6i4,
     x         /3x,'       nbl0 leng ',2i4) 
       call errex 
      ENDIF
C------------------------
      isize = nij*nab
      call zero(G,isize) 
      CALL GETLST(G,I1,NIJ0,I2,IRREP,LISTG)
C------------------------
C---- EXPAND IJ BLOCK
C------------------------
      ININEW = NIJ * NAB + 1
      INIOLD = NIJ0 * NAB0 + 1
      DO 900 KRREP=NIRREP,1,-1
       IKREP = INDPD(KRREP,IRREP)
       NBSI = NPOP(IKREP,1) 
       NDI = NDRPOP(IKREP)
       NBSK = NPOP(KRREP,2) 
       NDK = NDRPOP(KRREP)
C
       if (nbsi.eq.0 .or. nbsk.eq.0) go to 900
C
       NBSK0 = NDK + 1
       DO 800 NK=NBSK,NBSK0,-1
        NBSI0 = NDI + 1
        DO 700 NI=NBSI,NBSI0,-1
C-----------------------
C-----   EXPAND AB-LIST
C-----------------------
         DO 300 JRREP=NIRREP,1,-1
          IBREP = INDPD(IRREP,JRREP) 
          NABA = NVRT(IBREP,1)
          NDA  = NDRVRT(IBREP)
          NABA0 = NABA - NDA 
          NABB = NVRT(JRREP,2)
          NDB  = NDRVRT(JRREP)
          NABB0 = NABB - NDB 
          if (naba.eq.0 .or. nabb.eq.0) go to 300
c
          INIOLD = INIOLD - NABA0 * NABB0 
          ININEW = ININEW - NABA * NABB 
          CALL GCOPY(G(INIOLD),GTMP,NABA0*NABB0) 
          CALL REPOTT(GTMP,G(ININEW),
     X                     NABA0,NABB0,NABA,NABB,NDA,NDB,1)
 300     CONTINUE
C-----------------------
 700    CONTINUE
        NABNDK = NAB*NDI 
        ININEW = ININEW - NABNDK 
 800   CONTINUE
       NABIJ = NAB * NBSI * NDK 
       ININEW = ININEW - NABIJ  
C-----------------------
 900  CONTINUE
       if (ININEW.EQ.1  .AND. INIOLD.EQ.1)  RETURN
       write(6,950) ininew,iniold
 950   format (4x,'Some error (getgo1) ; ininew=',i3,5x,'iniold=',i3)
       call errex 
      END
