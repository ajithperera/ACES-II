      SUBROUTINE GETINO4 (INDXG,NTDIS,NRNUM,LISTG,IRREP) 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER INDXG(NTDIS)
C
      COMMON /SYM2/NPOP(8,2),NVRT(8,2),NTF12(6)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),INDPD(8,8) 
      COMMON /SYMDROP/NDRPOP(8),NDRVRT(8)
      COMMON /SYMPOP/IRPDPD(8,22),ISYTYP(2,500),NTOT(18) 
C
      NBLOCK = IRPDPD(IRREP,ISYTYP(2,LISTG))  
      Nlengt = IRPDPD(IRREP,ISYTYP(1,LISTG))  
      NIJ = 0
      NIJ0 = 0
      NAB = 0
      NAB0 = 0
      DO 10 JRREP=1,NIRREP
       IJREP = INDPD(IRREP,JRREP)
       NBSI = NVRT(IJREP,1)
       NBSJ = NPOP(JRREP,2)
       NBSI0 = NBSI - NDRVRT(IJREP)
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
C
      IF (NIJ0.NE.NBLOCK .or. nab0.ne.nlengt .or. nij.ne.ntdis .or.
     x    nab.ne.nrnum)  THEN
       write(6,14) listg
  14   format (3x,'start of getino4 with listg=',i5)
       WRITE (6,15) NIJ,NAB,nij0,nab0,nblock,nlengt
  15   FORMAT(//3X,'+++++  nij  nab  nij0  nab0  nblock nlengt',6i4) 
       WRITE (6,16) ntdis,nrnum
  16   FORMAT(3X,'+++++  ntdis  nrnum                         ',2i4) 
       call errex 
      ENDIF
C------------------------
      DO 200 I=1,NTDIS
       INDXG(I)=0
 200  CONTINUE 
C------------------------
C---- MAKE DROP--INDX 
C------------------------
      NIJIN = NIJ  
      DO 700 KRREP=NIRREP,1,-1
       IKREP = INDPD(KRREP,IRREP)
       NBSI = NVRT(IKREP,1) 
       NDI = NDRVRT(IKREP)
       NBSK = NPOP(KRREP,2) 
       NDK = NDRPOP(KRREP)
C
       if (nbsi.eq.0 .or. nbsk.eq.0) go to 678 
C
       NIJF = NIJIN  
       NBSK0 = NBSK - NDK 
       if (ndi.ne.0) then 
        IF (NBSK0.GE.1) THEN
         NBSI0 = NBSI - NDI 
         DO 600 I=NBSK0,1,-1 
          DO 500 J=NDI,1,-1
           INDXG(NIJF) = 1
           NIJF = NIJF - 1 
 500      CONTINUE
          NIJF = NIJF - NBSI0 
 600     CONTINUE
        ENDIF 
       else
        NIJF = NIJF - NBSI * NBSK0
       endif  
       if (ndk.eq.0)  go to 410 
       NIJI = NIJF  
       NIJF = NIJF - NBSI * NDK + 1
       DO 400 I=NIJI,NIJF,-1
        INDXG(I) = 1
 400   CONTINUE 
       NIJF = NIJF - 1
 410   CONTINUE 
       NIJIN = NIJIN - NBSI * NBSK 
 678   CONTINUE 
 700  CONTINUE 
      IF (NIJIN.EQ.0) RETURN 
       write (6,710) nijin   
 710   format (4x,'some error (getino4) = nijin=',i6,4x,'must be 0') 
       CALL ERREX 
      END
