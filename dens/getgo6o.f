      SUBROUTINE GETGO6O(G,GTMP,I1,NRDIS,I2,IRREP,LISTG,ISPIN,LISTW,
     x                     NAB,INDXG,NIKI) 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION G(2),GTMP(2) 
      INTEGER INDXG(2)     
      COMMON /SYM2/NPOP(8,2),NVRT(8,2),NTF12(6)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),INDPD(8,8) 
      COMMON /SYMDROP/NDRPOP(8),NDRVRT(8)
      COMMON /SYMPOP/IRPDPD(8,22),ISYTYP(2,500),NTOT(18) 
C------------------------ 
      write (6,2913)
 2913 format (4x,'---- getgo6o  -----') 
c
      NAB0 = IRPDPD(IRREP,ISYTYP(1,LISTG))  
      isize = NRDIS*NAB 
      call zero (G,isize) 
      NIJDR = 0 
      NIKF = NIKF + NRDIS - 1 
      DO 100 I=NIKI,NIKF
       if (indxg(i).eq.0 .or. indxg(i).eq.1) then
        IF (INDXG(I).EQ.0) NIJDR = NIJDR + 1
       else
        write (6,9731) niki,nikf,i
 9731   format(4x,'niki nikf i ===',3i8)
        call errex
       endif 
 100  CONTINUE 
      CALL GETLST(G,I1,NIJ0,I2,IRREP,LISTG)
C------------------------
C---- EXPAND IJ BLOCK
C------------------------
      ININEW = NRDIS * NAB + 1
      INIOLD = NIJDR * NAB0 + 1
      DO 900 NIK=NIKF,NIKI,-1  
       IF (INDXG(NIK).EQ.0) THEN
C-----------------------
C-----   EXPAND AB-LIST
C-----------------------
        DO 300 JRREP=NIRREP,1,-1
          IBREP = INDPD(IRREP,JRREP) 
          NABA = NPOP(IBREP,ISPIN)
          NDA  = NDRPOP(IBREP)
          NABA0 = NABA - NDA 
          NABB = NPOP(JRREP,ISPIN)
          NDB  = NDRPOP(JRREP)
          NABB0 = NABB - NDB 
          if (naba.eq.0 .or. nabb.eq.0) go to 300
c
          INIOLD = INIOLD - NABA0 * NABB0 
          ININEW = ININEW - NABA * NABB 
          CALL GCOPY (G(INIOLD),GTMP,NABA0*NABB0) 
          CALL REPOTT(GTMP,G(ININEW),
     X                    NABA0,NABB0,NABA,NABB,NDA,NDB,4)
 300    CONTINUE
       ELSE
        ININEW = ININEW - NAB
       ENDIF 
C-----------------------
 900  CONTINUE
      NIKI = NIKI + NRDIS
      I1 = I1 + NIJDR
c
       if (ININEW.EQ.1  .AND. INIOLD.EQ.1)  RETURN
       write(6,950) ininew,iniold
 950   format (4x,'Some error (getgo6o) ; ininew=',i3,5x,'iniold=',i3)
       call errex 
      END
