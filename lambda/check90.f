      subroutine check90(icore, maxcor)
c      implicit DOUBLE PRECISION (A-H,O-Z)
      implicit integer (A-Z)
c      integer maxcor,icore, nt,i0l, numsyw, dissyw
      DIMENSION ICORE(MAXCOR)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NF1(2),NF2(2)

      I0L=1
      isize=iintfp*nt(1)
      isize2=iintfp*nt(2)

      DISSYW=IRPDPD(1,ISYTYP(1,110))
      NUMSYW=IRPDPD(1,ISYTYP(2,110))
c      DISSYW=IRPDPD(1,ISYTYP(1,90))
c      NUMSYW=IRPDPD(1,ISYTYP(2,90))

      print *,"<<<<<<<<<<< check110 >>>>>>>>>>>>"

c      print *,"nt(1) = ",isize
c      print *,"nt(2) = ",isize2
       print *,"numsyw, dissyw = ",numsyw, dissyw
      
c      CALL GETLST(ICORE(I0L),1,1,1,3,90)
      
       call getlst(icore(i0l),1,NUMSYW,2,1,110)
c       call getall(icore(i0l),NUMSYW*dissyw,1,130)
      call checksum("checksum   20 =",icore(i0L),20)
      call checksum("checksum  100 =",icore(i0L),100)
      call checksum("checksum full =",icore(i0L),numsyw*dissyw)
      print *

      return
      end

