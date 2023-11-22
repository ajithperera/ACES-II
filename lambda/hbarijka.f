      subroutine hbarijka ( W, lstwin, listw,  xint, listint,
     $     icore, maxcor,ispin, iuhf, abfull)
C                                               _   
C  Calculate the W<IJ|KA>  intermediate, = <IJ||W||KA>   (sum over b,B)
C      W<IJ||KA> = <IJ||KA> - SUM T1(B,K)* <IJ||AB>         ispin = 1
C      W<ij||ka> = <ij||ka> - SUM T1(b,k)* <ij||ab>         ispin = 2
C      W<Ij|Ak>  = <Ij|Ak>  + SUM T1(b,k)* <Ij|Ab>          ispin = 3
C      W<Ij|Ka>  = <Ij|Ka>  + SUM T1(B,K)* <Ij|Ba>          ispin = 4
C      
C   PARAM        USE                                                    CHANGED
C   -----        ------------------------------------------------       -------
C   ISPIN ...... spin case, 1-4, of W (see above)                           NO
C   IUHF ....... IUHF = 0/1 for RHF/ROHF case.                              NO
C   W .......... Scratch space used for W. large enough for largest         YES
C                 <IJ||KA> irrep.
C   XINT ....... Scratch space used for <IJ||AB>. large enough for          YES
C                largest irrep, AFTER A<B is expanded to all (A,B)
C   ICORE ...... Holds t1 amplitudes, as well as various scratch vectors.   YES
C                large enough to hold t1(ispin) as well as largest XINT 
C                irrep, unexpanded (same as expanded in ispin =3,4 cases)
C   MAXCOR ..... length of ICORE available, in integer words.                NO
C
C   LSTWIN ..... input list for W.                                           NO
C   LISTW ...... output list for W.  COULD be the same as W input list,      NO
C                lstwin
C   LISTINT .... input list for <ab|ij> integrals                            NO
C
C  Renee P Mattie Feb 1991
CEND
C imp.inc should contain the appropriate implicit "undefined", "none" or
c  "logical" statement for this machine
C      include 'imp.inc'
      integer lstwin, listint, abfull(8)
      integer listw, ispin, iuhf, maxcor, icore(maxcor)
      double precision w(*), xint(*)
C
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
      integer irrepij, irrepb, irrepa, iofft1(8),
     $     nblkw, lenw, nij, nab,  newab, nb, nija, nk, i0, i1, i2,i3,
     $     isb, isa, inext, newcor, ierr,
     $     ioffab, ioffak, iconv
      double precision one, zilch, fact,
     $     xnorm
      parameter (one = 1.d+00, zilch = 0.d+00)
C functions
      double precision  snrm2
C
      iconv = iintfp
      if (ispin .eq. 1 .or. ispin .eq. 2) then
         isb = ispin
         isa = isb
         fact = -one
      elseif (ispin .eq. 3 ) then
         isb = 2
         isa = 1
         fact = one
      elseif (ispin .eq. 4) then
         isb = 1
         isa = 2
         fact = one
      else
C        you have made a serious error
         write (*,*) ' hbarijka: inappropriate ispin value; no ispin =',
     $        ispin,' case defined'
         call errex
      endif
C
      call get1t1 (icore(1), maxcor, inext, isb, iuhf, iofft1, ierr)
      if (ierr .ne. 0) then
         write (*,*) 'Error in hbarijka: not enough room for t1(',ispin,
     $        ')'
         call errex
      endif
C
      xnorm =  zilch
      do 1000, irrepij = 1, nirrep
         nblkw = irpdpd(irrepij,isytyp(2,lstwin))
         lenw  = irpdpd(irrepij,isytyp(1,lstwin))
         call getlst(w, 1, nblkw,2,irrepij,lstwin)
         if (ispin .ne. 3 ) then
C         list is stored (ij;ka) : need (ij;ak) , unless ispin is 3
            i0 = inext
            i1 = i0 + lenw*iconv
            i2 = i1 + lenw*iconv
            i3 = i2 + lenw*iconv
            if (i3.gt. maxcor) then
               write (*,*) 'maxcor = ', maxcor
               write (*,*) '   @hbarijka-F: need at least ',i3 - maxcor,
     $            ' more core to call symtr1'
               call errex
            endif
            call symtr1(irrepij, pop(1,isb), vrt(1,isa), lenw, w,
     $         icore(i0), icore(i1), icore(i2) )
         endif
         nab = irpdpd( irrepij, isytyp(1,listint))
         nij = irpdpd( irrepij, isytyp(2,listint))
         if (nij .ne. lenw) then
            write (*,*) 'hbarijka: list ',listw,'<ij|ka> and list ',
     $         listint,' (<ab||ij>) do not conform in i,j indices'
            call errex
         endif
         i3 = inext + nij*nab*iconv
         if (i3 .gt. maxcor) then
            write (*,*)'hbarijka needs at least ', i3 - maxcor,
     $         ' more icore to read in <ab|ij> , list', listint
            call errex
         endif
         call getlst(icore(inext),1,nij, 1, irrepij, listint)
         call transp(icore(inext), xint, nij, nab)
         if (ispin .eq. 1 .or. ispin .eq. 2) then
            newab = abfull(irrepij)
            call symexp(irrepij, vrt(1,isa), nij, xint)
            nab = newab
         endif
         if ( ispin .eq. 4) then
C     <Ij|Ba> was transposed to (Ij;Ba), now  goes to (I,j; a,B)
            i0 = inext
            i1 = i0 + nij*iconv
            i2 = i1 + nij*iconv
            i3 = i2 + nab
            if (i3.gt. maxcor ) then
               write(*,*) ' hbarijka: need at least ', i3- maxcor,
     $            ' more icore to proceed with symtr1' ,
     $            ' of transposed list', listint
               call errex
            endif
            call symtr1( irrepij, vrt(1,isb), vrt(1,isa), nij, xint,
     $         icore(i0), icore(i1), icore(i2) )
         endif
         ioffab = 1
         ioffak = 1
         do 900, irrepb = 1, nirrep
            irrepa = dirprd(irrepij, irrepb)
            nb = vrt(irrepb,isb)
            nija = vrt(irrepa,isa)*nij
            nk = pop(irrepb,isb)
            if (nk .ne. 0 .and. nija .ne. 0 .and. nb .ne. 0)
     $         call xgemm('N', 'N', nija, nk, nb,    fact,
     $         xint(ioffab),nija, icore(iofft1(irrepb)),nb,
     $         one, W(ioffak), nija )
            ioffab = ioffab + nija*nb
            ioffak = ioffak + nija*nk
 900     continue
         if (ispin .ne. 3) then
            i0 = inext
            i1 = i0 + nij * iconv
            i2 = i1 + nij * iconv
            i3 = i2 + nij * iconv 
            if (i3.gt. maxcor) then
               write (*,*) 'hbarijka: need at least ',i3 - maxcor,
     $            ' more core to call symtr1'
               call errex
            endif
            call symtr1(irrepij, vrt(1,isa), pop(1,isb), 
     $         lenw, w, icore(i0), icore(i1), icore(i2) )
         endif
Cdebug 1
         xnorm = xnorm + snrm2( nblkw*lenw, w, 1)**2
         call putlst(w, 1, nblkw, 2, irrepij, listw)
 1000 continue
Cdebug 1
C      write (*,*) ' @hbarijka: xnorm(',ispin, ') = ', xnorm
      return
      end
