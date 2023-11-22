      SUBROUTINE DINV11(VIA,DOO,ICORE,MAXCOR,IUHF)
C
C THIS SUBROUTINE CALCULATES THE V(AI) CONTRIBUTION
C DUE TO THE OCCUPIED-OCCUPIED BLOCK OF THE DOO MATRIX
C
C      V(AI) += SUM M,N <MI//NA> D(MN)
C
C THIS ROUTINE USES EXPLICITELY SYMMETRY
C
CEND
C
C  CODED SEPTEMBER/90
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION VIA
      INTEGER DIRPRD,DISSYW,POP,VRT,occ
      DIMENSION VIA(1),DOO(1),ICORE(MAXCOR)
      common/info/nocco(2),nvrto(2)
      COMMON /SYM/ POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF1BB,NF2AA,NF2BB
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
C
      MXCOR=MAXCOR

C Switch V(AI) to V(IA)
       call symtra(1,vrt(1,1),pop(1,1),1,via,icore)
       call dcopy(ntaa,icore,1,via,1)
       if (iuhf.ne.0) then
         call symtra(1,vrt(1,2),pop(1,2),1,via(1+ntaa),icore)
         call dcopy(ntbb,icore,1,via(ntaa+1),1)
       end if
C Add to V(IA)
       CALL DOOINVO(1,vIA,dOO,ICORE,MXCOR,IUHF,7)

C Switch back from V(IA) to V(AI)
       call symtra(1,pop(1,1),vrt(1,1),1,via,icore)
       call dcopy(ntaa,icore,1,via,1)
       if (iuhf.ne.0) then
         call symtra(1,pop(1,2),vrt(1,2),1,via(1+ntaa),icore)
         call dcopy(ntbb,icore,1,via(ntaa+1),1)
       end if

      return
      END
