      SUBROUTINE UPDENS(DOO,DVV,DVO,DOV,CORE,ISPIN,FACT,NONHF)
C
C WRITES DENSITY MATRIX TO DISK
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER POP,VRT,DIRPRD
      LOGICAL NONHF
      DIMENSION DOO(*),DVV(*),DVO(*),DOV(*),CORE(*)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/SYMINF/NSTART,NIRREP,IRREPY(255,2),DIRPRD(8,8)
C
      DATA ONE,HALF/1.0D0,0.5d0/
C
C OCC-OCC
C
      CALL GETLST(CORE,1,1,1,ISPIN,160)
      ioff = 0
      ioff2 = 0
      do 1 irrep=1,nirrep
         len=pop(irrep,ispin)
         do 11 ii=1,len
            core(ioff2+ii)=core(ioff2+ii)+fact*doo(ioff+(ii-1)*len+ii)
 11      continue
         ioff=ioff+len*len
         ioff2 = ioff2 + len
 1    continue
      CALL PUTLST(CORE,1,1,1,ISPIN,160)
C
C VRT-VRT
C
      CALL GETLST(CORE,1,1,1,2+ISPIN,160)
      ioff = 0
      ioff2 = 0
      do 2 irrep=1,nirrep
         len=vrt(irrep,ispin)
         do 21 ii=1,len
            core(ioff2+ii)=core(ioff2+ii)+fact*dvv(ioff+(ii-1)*len+ii)
 21      continue
         ioff=ioff+len*len
         ioff2 = ioff2 + len
 2    continue
      CALL PUTLST(CORE,1,1,1,2+ISPIN,160)

      IF (NONHF) THEN
         call saxpy(nt(ispin),one,dvo,1,dov,1)
         call sscal(nt(ispin),half,dov,1)
         call getlst(core,1,1,1,4+ispin,160)
         CALL SAXPY(nt(ispin),fact,DVO,1,CORE,1)
         CALL PUTLST(CORE,1,1,1,4+ISPIN,160)
      ENDIF
      RETURN
      END
