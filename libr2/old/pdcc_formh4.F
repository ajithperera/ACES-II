      
      SUBROUTINE PDCC_FORMH4(ICORE,MAXCOR,IUHF)
C
C   THIS ROUTINE CALCULATE L(MN,EF) T(NJ,FB) => (MB,EJ) for
C   Distinguish cluster approximation. 
C
C  CCD
C   
C   H(ME,JB) = - SUM N,F  L(MN,EF) T(NJ,FB)
C
C  QCISD
C
C   H(ME,JB) = - SUM N,F L(MN,EF) T(NJ,FB)
C
C  CCSD
C
C   G(IA,JB) = - SUM N,F L(MN,EF) (T(NJ,FB) + ...)
C
C  THIS TERM IS VERY SIMILAR TO THE T1(IJ,AB) CONTRIBUTION TO
C  THE W-RING INTERMEDIATE
C
C  THE SPIN CASES ARE
C
C    AAAA : =  - SUM M,E T1(IM,BE) T1(MJ,EA) - SUM m,e T1(Im,Be) T1(Jm,Ae)
C
C    ABAB : =  - SUM m,E T1(Im,Eb) T1(Jm,Ea)
C
C    ABBA : =  - SUM M,E T1(IM,BE) T1(Mj,Ea) - SUM m,e T1(Im,Be) T1(mj,eb)
C
C FOR UHF IN ADDITION THE BBBB BABA AND BABA SPIN CASES HAS TO CALCULATES
C
C
C  AAAA : LIST 254  (CALCULATED FROM ABAB AND ABBA IN QUIKAA)
C  BBBB : LIST 255  (UHF ONLY)
C  ABBA : LIST 258
C  BAAB : LIST 259  (UHF ONLY)
C  ABAB : LIST 256
C  BABA : LIST 257  (UHF ONLY)
C
CEND
C 
C  CODED AUGUST/90  JG
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC
      LOGICAL DEBUG
      DIMENSION ICORE(MAXCOR)
      COMMON/METH/MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC
#include "sym.com"
#include "sympop.com"
#include "info.com"
#include "parcc.com"
C
#ifdef _DEBUG_LVL0
      write(6,"(a)") " ---Entered PDCC_FORMH4---"
#endif 

#ifdef _DEBUG_LVL0
      Imode = 0


      Call Inipck(1,14,14,253,Imode,0,1)
      Call Inipck(1,9,9,254,Imode,0,1)
      Call Inipck(1,9,10,256,Imode,0,1)
      Call Inipck(1,11,11,258,Imode,0,1)

      If (Iuhf .ne. 0) Then
         Call Inipck(1,10,10,255,Imode,0,1)
         Call Inipck(1,10,9,257,Imode,0,1)
         Call Inipck(1,12,12,259,Imode,0,1)
      Endif
#endif 
#ifdef _DEBUG_LVLM

      irrepx=1
      Nsize = Idsymsz(irrepx,isytyp(1,37),isytyp(2,37))
      nrow=irpdpd(Irrepx,isytyp(1,37))
      ncol=irpdpd(Ireepx,isytyp(2,37))

      i000 = 1
      i010 = i000 + Nsize
      i020 = i010 + Nsize
      i030 = i020 + Nsize

C      Call getall(Icore(i000),Nsize,irrepx,37)
C      call sstgen(icore(i000),icore(i010),Nsize,Vrt(1,1),Pop(1,1),
C     +            Vrt(1,2),Pop(1,2),icore(i020),irrepx,"1432")
C      write(6,*) "List-137-reordered"
C      call output(icore(i010),1,nrow,1,ncol,nrow,ncol,1)

      nrow=irpdpd(Irrepx,isytyp(1,137))
      ncol=irpdpd(Irrepx,isytyp(2,137))
      Nsize = Idsymsz(irrepx,isytyp(1,137),isytyp(2,137))
      Call getall(Icore(i000),Nsize,irrepx,137)
      write(6,*) "List-137"
      call output(icore(i000),1,nrow,1,ncol,nrow,ncol,1)

      nrow=irpdpd(Irrepx,isytyp(1,134))
      ncol=irpdpd(Irrepx,isytyp(2,134))
      Nsize = Idsymsz(irrepx,isytyp(1,134),isytyp(2,134))
      Call getall(Icore(i000),Nsize,Irrepx,134)
      write(6,*) "List-134"
      call output(icore(i000),1,nrow,1,ncol,nrow,ncol,1)

      nrow=irpdpd(Irrepx,isytyp(1,139))
      ncol=irpdpd(Irrepx,isytyp(2,139))
      Nsize = Idsymsz(irrepx,isytyp(1,139),isytyp(2,139))
      Call getall(Icore(i000),Nsize,irrepx,139)
      write(6,*) "List-139"
      call output(icore(i000),1,nrow,1,ncol,nrow,ncol,1)

#endif 

C Lets zero out the lists 54-59 and 254-259
#ifdef _DEBUG_LVLM
      IRREPX = 1
      NSIZE=IDSYMSZ(IRREPX,ISYTYP(1,54),ISYTYP(2,54))
      I000 = 1
      CALL DZERO(ICORE(I000),NSIZE)
      CALL PUTALL(ICORE(I000),NSIZE,IRREPX,54)
      CALL PUTALL(ICORE(I000),NSIZE,IRREPX,254)
      NSIZE=IDSYMSZ(IRREPX,ISYTYP(1,58),ISYTYP(2,58))
      CALL DZERO(ICORE(I000),NSIZE)
      CALL PUTALL(ICORE(I000),NSIZE,IRREPX,58)
      CALL PUTALL(ICORE(I000),NSIZE,IRREPX,258)
      NSIZE=IDSYMSZ(IRREPX,ISYTYP(1,56),ISYTYP(2,56))
      CALL DZERO(ICORE(I000),NSIZE)
      CALL PUTALL(ICORE(I000),NSIZE,IRREPX,56)
      CALL PUTALL(ICORE(I000),NSIZE,IRREPX,256)

      IF (IUHF .NE. 0) THEN
         NSIZE=IDSYMSZ(IRREPX,ISYTYP(1,55),ISYTYP(2,55))
         CALL DZERO(ICORE(I000),NSIZE)
         CALL PUTALL(ICORE(I000),NSIZE,IRREPX,55)
         CALL PUTALL(ICORE(I000),NSIZE,IRREPX,255)
         NSIZE=IDSYMSZ(IRREPX,ISYTYP(1,57),ISYTYP(2,57))
         CALL DZERO(ICORE(I000),NSIZE)
         CALL PUTALL(ICORE(I000),NSIZE,IRREPX,57)
         CALL PUTALL(ICORE(I000),NSIZE,IRREPX,257)
         NSIZE=IDSYMSZ(IRREPX,ISYTYP(1,59),ISYTYP(2,59))
         CALL DZERO(ICORE(I000),NSIZE)
         CALL PUTALL(ICORE(I000),NSIZE,IRREPX,59)
         CALL PUTALL(ICORE(I000),NSIZE,IRREPX,259)
      ENDIF 
#endif 

      CALL PDCC_H4ALL(ICORE,MAXCOR,'ABBA',IUHF)
      CALL PDCC_H4ALL(ICORE,MAXCOR,'ABAB',IUHF)
      CALL PDCC_H4ALL(ICORE,MAXCOR,'AAAA',IUHF)

CSSS      CALL PDCC_QUIKAA2(ICORE,MAXCOR,IUHF)

      IF (IUHF .NE. 0) THEN
         CALL PDCC_H4ALL(ICORE,MAXCOR,'BBBB',IUHF)
         CALL PDCC_H4ALL(ICORE,MAXCOR,'BABA',IUHF) 
      ENDIF 

      RETURN
      END
