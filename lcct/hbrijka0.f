      subroutine hbrijka0(icore, maxcor, iuhf)
C
C   Driver for hbarijka, to calculate the W<ij|ka> intermediates such as
C      W<IJ||KA> = <IJ||KA> - SUM(B) T1(B,K)* <IJ||AB>         ispin = 1
C
C  Written by Renee P Mattie,  Feb. 1991
CEND
C imp.inc should contain the appropriate implicit "undefined", "none" or
c  "logical" statement for this machine
C      include 'imp.inc'
      integer maxcor, iuhf, icore(maxcor)
      integer iintln, ifltln, iintfp, ialone, ibitwd
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      integer NTAA,NTBB,NF1AA,NF2AA, NF1BB,NF2BB, pop,vrt
      COMMON /SYM/ POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF2AA,
     &             NF1BB,NF2BB
      integer nstart, nirrep, irreps, dirprd
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      integer irpdpd, isytyp, id
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
C
      integer lstwin, lstwout, lstabij
      parameter (lstwin = 7-1, lstwout = 7-1, lstabij = 14-1)
C
      integer i, irrepb, irrepa, icount, ispin, lst, isiz, maxw,
     $     maxijab, abfull(8,3), istep,
     $     i0, i1, i2, newcor, js

C
      maxijab = 0
      maxw = 0
      istep1=3-2*iuhf
      istep2=2-1*iuhf
      do 100, i = 1, nirrep
 
         do 50, ispin = 1,4,istep1
            lst = lstwin + ispin 
            isiz = irpdpd(i,isytyp(1,lst))*irpdpd(i,isytyp(2,lst))
            maxw = max(isiz, maxw)
 50      continue
C
         do 70, ispin = 1,3,istep2
            lst = lstabij + ispin
            if ( ispin .lt. 3) then
               icount = 0
               do 60, irrepb = 1, nirrep
                  irrepa = dirprd(irrepb,i)
                  icount = icount + vrt(irrepa,ispin)*vrt(irrepb,ispin)
 60            continue
               abfull(i,ispin) = icount
            else
               abfull(i,ispin) = irpdpd(i,isytyp(1,lst))
            endif
            isiz = abfull(i,ispin) * irpdpd (i, isytyp(2,lst))
            maxijab = max(isiz, maxijab)
 70      continue

 100  continue
C
      i0 = 1
      i1 = i0 + maxw * iintfp
      i2 = i1 + maxijab * iintfp
      newcor = maxcor - i2
      if ( newcor .le. 0) then
         write (*,*) ' hbarijka0: need at least ',-newcor,' more icore',
     $        ' to call hbarijka'
         call errex
      endif
C
      istep = 1
      if (iuhf .eq. 0) istep = 3
      do 200, ispin = 1, 4, istep
         if (ispin .eq. 4 ) then
            js = 3
         else
            js = ispin
         endif
         write (*,*) ' '
         call hbarijka( icore(i0), lstwin+ispin, lstwout+ispin,
     $      icore(i1), lstabij +js , icore(i2), newcor, ispin, 
     $      iuhf, abfull(1,js) )
 200  continue
      return
      end
