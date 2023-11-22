      subroutine mrcc_eval(scr, mxcor, ph2, nbas, iintfp)
c     
c     This subroutine calculates the active character of the orbitals 
c     given their approximate energy eigenvalues from a previous run.
c     
      implicit none
      integer mxcor, nactive, iactchar, iacteval, ieval, iintfp,
     $     nbas, iloc, i, ione, iprop, iactprop, ievalnew, ipropnew,
     $     ieval2
      double precision scr(mxcor), vmin
      character*8 string1, string2, string3
      character*2 ph2
      logical print, print2
c     
      ione = 1
      if (ph2 .eq. 'mm') then
         string1 = 'NUMACT_M'
         string2 = 'EVAL_M  '
         string3 = 'PROP_M  '
      elseif (ph2 .eq. 'jj') then
         string1 = 'NUMACT_J'
         string2 = 'EVAL_J  '
         string3 = 'PROP_J  '
      elseif (ph2 .eq. 'ee') then
         string1 = 'NUMACT_E'
         string2 = 'EVAL_E  '
         string3 = 'PROP_E  '
      elseif (ph2 .eq. 'bb') then
         string1 = 'NUMACT_B'
         string2 = 'EVAL_B  '
         string3 = 'PROP_B  '
      else
         write(6,*) ' unknown ph character in mrcc_eval'
         call errex
      endif
c     
      call getrec(-1, 'JOBARC', string1, ione, nactive)
      if (nactive .ne. 0) then
c     
c     Active orbital energies exists on Jobarc. Process.
c     
         iactchar = 1
         iacteval = iactchar + nbas
         iactprop = iacteval + nactive
         ieval = iactprop + nactive
         iprop = ieval + nbas
         ievalnew = iprop + nbas
         ipropnew = ievalnew + nactive
         ieval2 = ipropnew + nbas
c     
         call getrec(20, 'JOBARC', 'SCFEVLA0', nbas*iintfp,
     $        scr(ieval))
         call getrec(20, 'JOBARC', 'SCFPROP0', nbas*iintfp,
     $        scr(iprop))
         call getrec(20, 'JOBARC', string2, nactive*iintfp,
     $        scr(iacteval))
         call getrec(20, 'JOBARC', string3, nactive*iintfp,
     $        scr(iactprop))
c     
         print = .false.
         if (print) then
            write(6,*) ' @mrcc_eval: all eigenvalues '
            call output(scr(ieval), 1, 1, 1, nbas, 1, nbas, 1)
            write(6,*) ' active eigenvalues ', ph2
            call output(scr(iacteval),1, 1, 1, nactive,
     $           1, nactive, 1)
            write(6,*) ' @mrcc_eval: all props '
            call output(scr(iprop), 1, 1, 1, nbas, 1, nbas, 1)
            write(6,*) ' active eigenvalues ', ph2
            call output(scr(iactprop),1, 1, 1, nactive,
     $           1, nactive, 1)
         endif
c     
         call zero(scr(iactchar), nbas)
         do i = 1, nactive
            call scopy(nbas, scr(ieval), 1, scr(ieval2), 1)
            call fndclose2(nbas, scr(ieval2), scr(iacteval+i-1),
     $           scr(iprop), scr(iactprop+i-1), iloc, 0.3d0)
            if (iloc .eq. -1) then
               write(6,*) ' error in mrcc_eval'
               write(6,*) ' matching pair cannot be found'
               write(6,*) ' current pair '
               write(6,*) scr(iacteval+i-1), scr(iactprop+i-1)
               write(6,*) ' current eigenvalues'
               call output(scr(ieval), 1, 1, 1, nbas, 1, nbas, 1)
               write(6,*) ' current properties'
               call output(scr(iprop), 1, 1, 1, nbas, 1, nbas, 1)
               call errex
            else
               scr(ievalnew+i-1) = scr(ieval+iloc-1)
               scr(ipropnew+i-1) = scr(iprop+iloc-1)
               scr(ieval+iloc-1) = 1.0d30
               scr(iprop+iloc-1) = 1.0d30
               scr(iactchar+iloc-1) = 1.0d0
            endif
         enddo
c     
      else
         call zero(scr, nbas)
      endif
c     
      print2 = .true.
      if (print2 .and. nactive .gt. 0) then
         write(6,*) ' Energies old active orbitals'
         call output(scr(iacteval), 1, 1, 1, nactive, 1, nactive, 1)
         write(6,*) ' Energies new active orbitals'
         call output(scr(ievalnew), 1, 1, 1, nactive, 1, nactive, 1)
         write(6,*) ' properties old active orbitals'
         call output(scr(iactprop), 1, 1, 1, nactive, 1, nactive, 1)
         write(6,*) ' properties new active orbitals'
         call output(scr(ipropnew), 1, 1, 1, nactive, 1, nactive, 1)
      endif

      if (print) then
         write(6,*) ' Found active orbitals'
         call output(scr(iactchar), 1, 1, 1, nbas, 1, nbas, 1)
      endif
c     
      return
      end
