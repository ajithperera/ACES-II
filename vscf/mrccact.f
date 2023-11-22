      subroutine mrccact(scr, mxcor, ph2, nao, nbasx)
c     
c     this subroutine calculates the active character of the MO's 
c     given the zao_projectors on JOBARC
c
c     mxcor is the memory used in the routine = 3*nbasx*nbasx*iintfp
c     nao: number of MO's (number of symmetry adapted AO's, i.e. so's
c     nbasx: number of ZAO's (non-symmetry adapted AO's)
c     
      implicit none
c
      integer mxcor, nao, nbasx, ione, setscr, itodo,
     $    iscr1, iscr2, iscr3, i, ioff
      double precision scr(mxcor), zilch, one
      character*2 ph2
      character*8 string1, string2
c
      integer iintln,ifltln,iintfp,ialone,ibitwd
      common /machsp/ iintln,ifltln,iintfp,ialone,ibitwd
c     
      ione = 1
      one = 1.0d0
      zilch = 0.0d0
c
      iscr1 = 1
      iscr2 = iscr1 + nbasx*nbasx
      iscr3 = iscr2 + nbasx*nbasx
c     
      if (ph2 .eq. 'mm') then
         string1 = 'OCCACT_M'
         string2 = 'ZAO_MM  '
      elseif (ph2 .eq. 'jj') then
         string1 = 'OCCACT_J'
         string2 = 'ZAO_JJ  '
      elseif (ph2 .eq. 'ee') then
         string1 = 'VRTACT_E'
         string2 = 'ZAO_EE  '
      elseif (ph2 .eq. 'bb') then
         string1 = 'VRTACT_B'
         string2 = 'ZAO_BB  '
      endif
c     
      call getrec(-1, 'JOBARC', string1, ione, itodo)
      if (itodo .eq. 1) then
c     
c     projector exists on Jobarc. Pick it up and put it into iscr1
c     
      call getrec(20, 'JOBARC', string2, nbasx*nbasx*iintfp,scr(iscr1))
c     
c     transform to AO basis
c     
c     get cmp2zmat matrix.  A(zao, ao)
c     
         call getrec(20, 'JOBARC', 'CMP2ZMAT',nbasx*nao*iintfp,
     $    scr(iscr2))
c     
         call xgemm('N', 'N', nbasx, nao, nbasx, one,
     $        scr(iscr1), nbasx, scr(iscr2), nbasx,
     $        zilch, scr(iscr3), nbasx)
c
c     get zmat2cmp matrix. A-1(ao,zao)
c     
         call getrec(20, 'JOBARC', 'ZMAT2CMP',nbasx*nao*iintfp,
     $       scr(iscr2))
c     
         call xgemm('N', 'N', nao, nao, nbasx, one,
     $        scr(iscr2), nao, scr(iscr3), nbasx,
     $        zilch, scr(iscr1), nao)
c     
c     scr(iscr1) contains projector in AO basis
c     
c     form P * S * MO -> C (in scr2)
c     
c     get S in AO basis
c     
         call getrec(20, 'JOBARC', 'AOOVRLAP', nao*nao*iintfp,
     $       scr(iscr3))
c     
c     form P*S -> scr2
c     
         call xgemm('N', 'N', nao, nao, nao, one,
     $        scr(iscr1), nao, scr(iscr3), nao,
     $        zilch, scr(iscr2), nao)
c     
         call getrec(20, 'JOBARC', 'SCFEVCA0', nao*nao*iintfp,
     $       scr(iscr1))
c     
         call xgemm('N', 'N', nao, nao, nao, one,
     $        scr(iscr2), nao, scr(iscr1), nao,
     $        zilch, scr(iscr3), nao)
c     
c     C -> scr(3)
c     
         call getrec(20, 'JOBARC', 'AOOVRLAP', nao*nao*iintfp,
     $       scr(iscr2))
c     
c     finally form overlap: CT * S * C
c     
         call xgemm('N', 'N', nao, nao, nao, one, scr(iscr2),
     $        nao, scr(iscr3), nao, zilch, scr(iscr1), nao)
         call xgemm('T', 'N', nao, nao, nao, one, scr(iscr3),
     $        nao, scr(iscr1), nao, zilch, scr(iscr2), nao)
c     
c     finally put the diagonal elements in scr(1..nao)
c     
         ioff = iscr2
         do i = 1, nao
            scr(i) = scr(ioff)
            ioff = ioff+nao+1
         enddo
c     
         write(6,*) ' @mrccact: current active type ', ph2
         call output(scr, 1, nao, 1, 1, nao, 1, 1)
c     
       else
         call zero(scr, nao)
      endif
c     
      return
      end




c     
c Local Variables: c
c mode: fortran c
c fortran-do-indent: 2 c
c fortran-if-indent: 2 c
c fortran-continuation-indent: 4 c
c fortran-comment-indent-style: nil c
c End: c

