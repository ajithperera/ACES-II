CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE EXPAND (DENS,TMP,NIRREP,ITYPE,ispin,LEN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION DENS(LEN),TMP(2) 
      COMMON /SYM2/NPOP(8,2),NVRT(8,2),NTF12(6)
      COMMON /SYMDROP/NDRPOP(8),NDRVRT(8)
C
      do 5 I=1,LEN
       TMP(I) = DENS(I) 
  5   CONTINUE
C
      NIJ1 = 0
      NIJ2 = 0
      DO 10 I=1,NIRREP
      if (itype.eq.4) then
C------------------------
C---- EXPAND O-O DENSITY BLOCK 
C------------------------
       ni = npop(i,ispin)
       nid= ndrpop(i)
       nj = npop(i,ispin)
       njd= ndrpop(i)
      endif 
      if (itype.eq.3) then 
C------------------------
C---- EXPAND V-O DENSITY BLOCK 
C------------------------
       ni = nvrt(i,ispin)
       nid= ndrvrt(i)
       nj = npop(i,ispin)
       njd= ndrpop(i)
      endif 
      if (itype.eq.2) then 
C------------------------
C---- EXPAND O-V DENSITY BLOCK 
C------------------------
       ni = npop(i,ispin)
       nid= ndrpop(i)
       nj = nvrt(i,ispin)
       njd= ndrvrt(i)
      endif 
      if (itype.eq.1) then 
C------------------------
C---- EXPAND V-V DENSITY BLOCK 
C------------------------
       ni = nvrt(i,ispin) 
       nid= ndrvrt(i) 
       nj = nvrt(i,ispin) 
       njd= ndrvrt(i) 
      endif 
      nij1 = nij1 + ni * nj
      nij2 = nij2 + (ni-nid)*(nj-njd)
  10  CONTINUE 
c----------------------------
      INI1 = NIJ1 + 1
      INI2 = NIJ2 + 1 
      inidif = nij1 - nij2 
c     if (inidif.gt.0) call zero (dens(ini2),inidif)
      call zero (dens,nij1)  
c 
      DO 100 I=NIRREP,1,-1
      if (itype.eq.4) then
       ni = npop(i,ispin)
       nid= ndrpop(i)
       nj = npop(i,ispin)
       njd= ndrpop(i)
      endif 
      if (itype.eq.3) then 
       ni = nvrt(i,ispin)
       nid= ndrvrt(i)
       nj = npop(i,ispin)
       njd= ndrpop(i)
      endif 
      if (itype.eq.2) then 
       ni = npop(i,ispin)
       nid= ndrpop(i)
       nj = nvrt(i,ispin)
       njd= ndrvrt(i)
      endif 
      if (itype.eq.1) then 
       ni = nvrt(i,ispin) 
       nid= ndrvrt(i) 
       nj = nvrt(i,ispin) 
       njd= ndrvrt(i) 
      endif
       ni2 = ni - nid
       nj2 = nj - njd
       INI1 = INI1 - NI*NJ 
       INI2 = INI2 - ni2 * nj2
       CALL REPOTT(TMP(INI2),DENS(INI1),ni2,nj2,ni,nj,nid,njd,itype) 
  100 CONTINUE
c----------------------------
      return
      end
