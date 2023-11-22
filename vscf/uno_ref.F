      subroutine uno_ref(pk,oneh,dens,fock,eval,evec,icore,maxmem,
     &     ldim1,ldim2,nbas,ipksiz,repuls,etot,iuhf, 
     &     luint,
     &     xform,densnew,
     &     densscr,scra,scrb,densfrac,
     &     naobasfn,scfks,scfksexact,scfkslastiter,
     &     V,z1,ksa,ksb,screxc,scr,scr2,valao,valgradao,totwt,
     &     max_angpts,natoms,intnumradpts,ncount,kshf)
c     
c     from the uhf or rohf density matrix the spatial density matrix is formed
c     dens_tot = dens_a + dens_b
c     this matrix is transformed to orthogonal basis and diagonalized to give 
c     new restricted mo's. the a new density matrix is formed specified by 
c     uno_charge and uno_mult. this defines the reference determinant in 
c     correlated calculations. 
c     the actual orbitals and orbital energies are defined in a semi-canonical 
c     fashion, but using a closed shell density with fractional occupations.
c     
      implicit double precision(a-h,o-z)
c     
      dimension oneh(ldim1),dens((iuhf+1)*ldim1),fock((iuhf+1)*ldim1)
      dimension evec((iuhf+1)*ldim2),eval((iuhf+1)*nbas)
      dimension icore(1),pk(ipksiz)
      dimension nocc(8,2)
c     
      dimension xform(ldim1), densnew(2*ldim2), densscr(2*ldim2),
     $     scra(2*ldim2), scrb(2*ldim2), densfrac(2*ldim2)
      dimension iadd(8,2),irem(8,2)
c---------------------------------------------------------------------
#include<aces.h>
#include<intgrtflags.com>
c
       integer iuhf,naobasfn,ncount,natoms,max_angpts,
     &         intnumradpts
       logical scfks,scfksexact,scfkslastiter,kshf
       integer z1(naobasfn,2)
       double precision V(naobasfn,naobasfn,iuhf+1),
     &                  ksa(ldim2),ksb(ldim2),
     &                  screxc(naobasfn,naobasfn),
     &                  scr(naobasfn,naobasfn,iuhf+1),
     &                  scr2(nbas,nbas),coef_nonloc,
     &                  valao(naobasfn,intnumradpts,
     &                        max_angpts,ncount),
     &                  valgradao(naobasfn,intnumradpts,
     &                            max_angpts,ncount,3),
     &                  totwt(ncount,intnumradpts,max_angpts)     
c----------------------------------------------------------------------------
c     
      common /machsp/ iintln,ifltln,iintfp,ialone,ibitwd
      common /files/ luout,moints
      common /flags/ iflags(100)
      common /flags2/ iflags2(500)
      common /popul/ nocc
c     nevin 6/5/95 modified so qrhf will work with fock=ao
      common /fock/ aofil
#include "symm2.com"
      logical aotran, print, aofil,noconv
      integer uno_charge, uno_open, uno_mult, uno_all, uno_closed,
     $     nclosed, n0_open
      double precision fraction
c
      if (iuhf .eq. 0) then
         if (iflags(22) .eq. 1) then
cmn
c   Everything is fine. Only use UNO_REF in first cycle of Brueckner calculation
c            
            return
         else
c
         write(6,*) ' @UNO_REF: should not be called with'
         write(6,*) ' closed shell reference state on input'
         call errex
      endif
      endif
      print = .true.
      uno_charge=iflags2(149)
      uno_mult = iflags2(150)
c     
c     calculate spatial density matrix.
c     
      call vadd(densnew, dens(itriof(1)), dens(ldim1+itriof(1)),
     $     ldim1, 1.0d0)
c     
c     get eigenvectors and eigenvalues for average density
c     
c     The density matrix is given in AO basis.
c     
      aotran = .true.
      ioff = 1
      do irp = 1, nirrep
         if (nbfirr(irp).eq. 0) go to 100
         nsize=nbfirr(irp)
         call expnd2(densnew(itriof(irp)), densscr, nsize)
         if (aotran) then
c     
            call expnd2(xform(itriof(irp)), scrb, nsize)
c     
c     calculate s^(1/2) by inverting xform. Does this work if linear dependencies?
c     
            call zero(scra,nsize*nsize)
            do i = 1, nsize
               scra((i-1)*nsize+i)= 1.0d0
            enddo
            info = 0
#ifdef _DEBUG_LVL0
       write(6,"(a)") "The overlap matirx"
       call output(scrb, 1, nsize, 1, nsize, nsize, nsize, 1)
CSSS       call eig(densscr, scra, nsize, nsize, 0)
CSSS       call output(densscr, 1, nsize, 1, nsize, nsize, nsize, 1)
#endif 
            call dgesv(nsize, nsize, scrb, nsize, icore, scra,
     $           nsize, info)
            if (info .ne. 0) then
               write(6,*) ' @uno_ref:something wrong inversion of xform'
               call errex
            endif
#ifdef _DEBUG_LVL0
       write(6,"(a)") "The inverse of overlap matirx"
       call output(scra, 1, nsize, 1, nsize, nsize, nsize, 1)
#endif 
c     
c     check inverse
c     
c                 write(6,*) ' checking inverse'
c                 call xgemm('N', 'N', nsize, nsize, nsize, 1.0d0,
c     $           scra, nsize, densnew, nsize, 0.0d0, scrb, nsize)
c                 call output(scrb, 1, nsize, 1, nsize, nsize, nsize, 1)
c     
            call xgemm('N', 'T', nsize, nsize, nsize, 1.0d0,
     $           densscr, nsize, scra, nsize, 0.0d0,
     $           scrb, nsize)
            call xgemm('N', 'N', nsize, nsize, nsize, 1.0d0, scra,
     $           nsize, scrb, nsize, 0.0d0, densscr, nsize)
            trace = 0.0
            do i = 1, nsize
               trace = trace + densscr((i-1)*nsize+i)
            enddo
c     write(6,*) ' @uno_ref: trace of transformed dens equals ' ,
c     $           trace
         endif
#ifdef _DEBUG_LVL0
       call output(densscr, 1, nsize, 1, nsize, nsize, nsize, 1)
#endif 
         call eig(densscr, scra, nsize, nsize, 0)
c     
c     copy eigenvalues (obtained on diagonal of denscr) in eval
c     and eigenvectors in evec
c     
         icount = 1
         do i = 1, nsize
            eval(i+ioff-1) = densscr(icount)
            icount = icount + nsize + 1
         enddo
c     
         call scopy(nsize*nsize, scra, 1, evec(isqrof(irp)), 1)
c     
         print = .true.
         if (print) then
            write(6,*) ' @uno_ref: symmetry block :', irp
            write(6,*) ' natural occupation numbers:'
            call output(eval(ioff), 1, 1, 1, nsize, 1, nsize, 1)
            trace = 0.0d0
            do i = 1, nsize
               trace = trace + eval(ioff+i-1)
            enddo
c     mn            write(6,*) ' trace of density matrix block', trace
c     mn            write(6,*) ' corresponding symmetry blocked natural',
c     mn     $           ' orbitals'
c     mn            call output(evec(isqrof(irp)), 1, nsize, 1, nsize,
c     mn     $           nsize, nsize, 1)
         endif
c     
         ioff = ioff + nsize
c     
 100  enddo
c     
c     current uno's and natural occupation numbers are there
c     
c     construct new density from uno's, also fill occupation vector
c     
      uno_open = uno_mult - 1
c     
c     determine # electrons in uno_ref
c     
      nelec = 0
      do ispin = 1, 2
         do irp = 1, nirrep
            nelec = nelec + nocc(irp,ispin)
         enddo
      enddo
c
c     
c     determine # of doubly occupied electrons in uno_ref
c     
      nclosed = 0
      do irp = 1, nirrep
cmn         nclosed = nclosed + min(nocc(irp,1), nocc(irp,2))
         nclosed = nclosed + nocc(irp,2)
      enddo
c     
c     get Hartree-Fock Charge
c     
      icharge = iflags(28)
      imult = iflags(29)
      nelec = nelec + icharge
c     
c     nelec denotes number of electrons in neutral system (i.e. same as nuclear charge)
c     
      uno_all = nelec - uno_charge
      uno_closed = uno_all - uno_open
c
      n0_open = nelec - 2*nclosed ! number of fractional electrons for neutral system.
      If (n0_open .NE. 0) Then
         fraction = n0_open / (2.0d0 * (imult-1))
      Else
         fraction = 1.0D0
      Endif 
c
      if (mod(uno_closed,2) .ne. 0) then
         write(6,*) ' impossible values of nelec, uno_charge, uno_mult'
         write(6,*) nelec, uno_charge, uno_mult
         call errex
      endif
      uno_closed = uno_closed/2
c     
c     construct density matrix (always construct both alpha and beta part).
c     
      factor = 1.0d0

      call zero(densscr,2*ldim2)
      call zero(densfrac,2*ldim2)
      call izero(nocc, 16)
c     
      do iorb = 1, nclosed
c     
c     find maximum occupation number
c     
         amax = -0.5
         icount = 0
         do irp = 1, nirrep
            do i = 1, nbfirr(irp)
               icount = icount + 1
               if (eval(icount) .gt. amax) then
                  irp0 = irp
                  icount0 = icount
                  i0 = i
                  amax = eval(icount)
               endif
            enddo
         enddo
c     
         eval(icount0) = -2.0d0
         nocc(irp0,1) = nocc(irp0,1) + 1
         nocc(irp0,2) = nocc(irp0,2) + 1
c     
c     add factor * c(p, icount) * c_dagger(q,icount) to density matrix
c     
         icount = isqrof(irp0) + (i0-1) * nbfirr(irp0)
         nsize  = nbfirr(irp0)
         Write(6,*) Icount
         call xgemm('N', 'T', nsize, nsize, 1, factor,
     $        evec(icount), nsize, evec(icount), nsize,
     $        1.0d0, densscr(isqrof(irp0)), nsize)
      enddo
c     
      do i = 1, ldim2
         densscr(ldim2+i) = densscr(i)
      enddo
c
c copy closed shell part to fractional density
c
      call scopy(2*ldim2, densscr, 1, densfrac, 1)
c
      do iorb = nclosed+1, uno_closed
c     
c     find maximum occupation number
c     
         amax = -0.5
         icount = 0
         do irp = 1, nirrep
            do i = 1, nbfirr(irp)
               icount = icount + 1
               if (eval(icount) .gt. amax) then
                  irp0 = irp
                  icount0 = icount
                  i0 = i
                  amax = eval(icount)
               endif
            enddo
         enddo
c     
         eval(icount0) = -2.0d0
         nocc(irp0,1) = nocc(irp0,1) + 1
         nocc(irp0,2) = nocc(irp0,2) + 1
c     
c     add factor * c(p, icount) * c_dagger(q,icount) to density matrix
c     
         icount = isqrof(irp0) + (i0-1) * nbfirr(irp0)
         nsize  = nbfirr(irp0)
         call xgemm('N', 'T', nsize, nsize, 1, factor,
     $        evec(icount), nsize, evec(icount), nsize,
     $        1.0d0, densscr(isqrof(irp0)), nsize)
c
c also add to densfrac
c
         call xgemm('N', 'T', nsize, nsize, 1, fraction,
     $        evec(icount), nsize, evec(icount), nsize,
     $        1.0d0, densfrac(isqrof(irp0)), nsize)
      enddo
c     
      do i = 1, ldim2
         densscr(ldim2+i) = densscr(i)
         densfrac(ldim2+i) = densfrac(i)
      enddo
c     
c     closed shell part is done, now open-shell part if it is there.
c     
      if (uno_open .ne. 0) then
         do iorb = 1, uno_open
c     
c     find maximum occupation number
c     
            amax = -0.5
            icount = 0
            do irp = 1, nirrep
               do i = 1, nbfirr(irp)
                  icount = icount + 1
                  if (eval(icount) .gt. amax) then
                     irp0 = irp
                     icount0 = icount
                     i0 = i
                     amax = eval(icount)
                  endif
               enddo
            enddo
c     
            eval(icount0) = -2.0d0
            nocc(irp0,1) = nocc(irp0,1) + 1
c     
c     add factor * c(p, icount) * ct(q,icount) to alpha density matrix
c     
            icount = isqrof(irp0) + (i0-1) * nbfirr(irp0)
            nsize  = nbfirr(irp0)
            call xgemm('N', 'T', nsize, nsize, 1, 1.0d0,
     $           evec(icount), nsize, evec(icount), nsize,
     $           1.0d0, densscr(isqrof(irp0)), nsize)
         enddo
      endif
c     
c     determine trace of current density matrix per symmetry block, alpha/beta
c     
      print = .true.
      if (print) then
         do irp = 1, nirrep
            nsize = nbfirr(irp)
            if (nsize .gt. 0) then
               tracea = 0.0d0
               traceb = 0.0d0
               do i = 1, nsize
                  tracea = tracea + densscr(isqrof(irp)+(i-1)*nsize+i-1)
               enddo
               Write(6,*) 
               do i = 1, nsize
                  traceb = traceb + densscr(ldim2+
     $                 isqrof(irp)+(i-1)*nsize+i-1)
               enddo
               write(6,*) ' partial traces new density matrix ', tracea,
     $              traceb
            endif
         enddo
      endif
c     
c     new density matrix is formed in orthogonal basis. transform to ao basis
c     
      aotran = .true.
      if (.not. aotran)
     $     write(6,*) ' @uno_ref: aotran is set to false!!!'
      if (aotran) then
         do irp = 1, nirrep
            if (nbfirr(irp).ne. 0) then
               nsize=nbfirr(irp)
               call expnd2(xform(itriof(irp)), scra, nsize)
c     alpha spin
               call xgemm('N', 'T', nsize, nsize, nsize, 1.0d0,
     $              densscr(isqrof(irp)), nsize, scra, nsize, 0.0d0,
     $              scrb, nsize)
               call xgemm('N', 'N', nsize, nsize, nsize, 1.0d0,
     $              scra, nsize,
     $              scrb, nsize, 0.0d0, densnew(isqrof(irp)), nsize)
c     beta spin
               call xgemm('N', 'T', nsize, nsize, nsize, 1.0d0,
     $              densscr(ldim2+isqrof(irp)), nsize, scra,
     $              nsize, 0.0d0, scrb, nsize)
               call xgemm('N', 'N', nsize, nsize, nsize, 1.0d0,
     $              scra, nsize, scrb, nsize,
     $              0.0d0, densnew(ldim2+isqrof(irp)), nsize)
c     
            endif
         enddo
      else
         call scopy(2*ldim2, densscr, 1, densnew, 1)
      endif
c     
c     fractional density matrix is formed in orthogonal basis. transform to ao basis
c     
      aotran = .true.
      if (.not. aotran)
     $     write(6,*) ' @uno_ref: aotran is set to false!!!'
      if (aotran) then
         do irp = 1, nirrep
            if (nbfirr(irp).ne. 0) then
               nsize=nbfirr(irp)
               call expnd2(xform(itriof(irp)), scra, nsize)
c     alpha spin
               call xgemm('N', 'T', nsize, nsize, nsize, 1.0d0,
     $              densfrac(isqrof(irp)), nsize, scra, nsize, 0.0d0,
     $              scrb, nsize)
               call xgemm('N', 'N', nsize, nsize, nsize, 1.0d0,
     $              scra, nsize,
     $              scrb, nsize, 0.0d0, densfrac(isqrof(irp)), nsize)
c     beta spin
               call xgemm('N', 'T', nsize, nsize, nsize, 1.0d0,
     $              densfrac(ldim2+isqrof(irp)), nsize, scra,
     $              nsize, 0.0d0, scrb, nsize)
               call xgemm('N', 'N', nsize, nsize, nsize, 1.0d0,
     $              scra, nsize, scrb, nsize,
     $              0.0d0, densfrac(ldim2+isqrof(irp)), nsize)
c     
            endif
         enddo
      endif
c
c we have overwritten densfrac: Now: neutral density in AO basis.
c     
c     sort eigenvectors according to occupied/virtual.
c     
      call scopy(ldim2, evec, 1, scra, 1)
      icount = 1
      ivec = 1
      iocca = 1
      ioccb = ldim2+1
      do irp = 1, nirrep
         nsize = nbfirr(irp)
         if (iocca .ne. isqrof(irp)) then
            write(6,*) ' @uno_ref: something wrong iocca'
            call errex
         endif
         if (ioccb .ne. ldim2+isqrof(irp)) then
            write(6,*) ' @uno_ref: something wrong ioccb'
            call errex
         endif
         ivrta = iocca + nocc(irp, 1) * nsize
         ivrtb = ioccb + nocc(irp, 2) * nsize
         do i =1, nsize
            if (eval(icount) .lt. -1.5) then
c     
c     doubly occupied orbital
c     
               call scopy(nsize, scra(ivec), 1, evec(iocca), 1)
               call scopy(nsize, scra(ivec), 1, evec(ioccb), 1)
               iocca = iocca + nsize
               ioccb = ioccb + nsize
            elseif (eval(icount) .lt. -0.5) then
c     
c     singly occupied orbital
c     
               call scopy(nsize, scra(ivec), 1, evec(iocca), 1)
               call scopy(nsize, scra(ivec), 1, evec(ivrtb), 1)
               iocca = iocca + nsize
               ivrtb = ivrtb + nsize
            else
c     
c     doubly virtual orbital
c     
               call scopy(nsize, scra(ivec), 1, evec(ivrta), 1)
               call scopy(nsize, scra(ivec), 1, evec(ivrtb), 1)
               ivrta = ivrta + nsize
               ivrtb = ivrtb + nsize
            endif
            ivec = ivec + nsize
            icount = icount + 1
         enddo
c     
         iocca = ivrta
         ioccb = ivrtb
      enddo
c     
c     transform eigvec to ao basis
c     
      if (aotran) then
         do irp = 1, nirrep
            if (nbfirr(irp).eq. 0) go to 300
            nsize=nbfirr(irp)
            call expnd2(xform(itriof(irp)), scra, nsize)
            call xgemm('N', 'N', nsize, nsize, nsize, 1.0d0,
     $           scra, nsize, evec(isqrof(irp)), nsize, 0.0d0,
     $           scrb(isqrof(irp)), nsize)
            call xgemm('N', 'N', nsize, nsize, nsize, 1.0d0,
     $           scra, nsize, evec(ldim2+isqrof(irp)),
     $           nsize, 0.0d0,
     $           scrb(ldim2+isqrof(irp)), nsize)
 300        continue
         enddo
c     
         if (print) then
            do irp = 1, nirrep
               nsize = nbfirr(irp)
               if (nsize .gt. 0) then
                  write(6,*) 'alpha eigenvectors in block ',
     $                 irp, nocc(irp,1)
                  call output(evec(isqrof(irp)), 1, nsize, 1,
     $                 nsize, nsize, nsize, 1)
                  write(6,*) 'beta eigenvectors in block ',
     $                 irp, nocc(irp,2)
                  call output(evec(ldim2+isqrof(irp)), 1, nsize, 1,
     $                 nsize, nsize, nsize, 1)
               endif
            enddo
         endif
c     
c     scrb contains all eigenvectors. copy to evec
c     
         call scopy(2*ldim2, scrb, 1, evec, 1)
      endif
c     
c     
c     construct density matrix in alternative way and subtract from
c     previous result. Check result to see all is fine.
c     
c     call scopy(2*ldim2, densnew, 1, scrb, 1)
      call zero(scrb, 2*ldim2)
      do irp = 1, nirrep
         if (nbfirr(irp).ne. 0) then
            nsize=nbfirr(irp)
            call xgemm('N', 'T', nsize, nsize, nocc(irp,1), 1.0d0,
     $           evec(isqrof(irp)), nsize, evec(isqrof(irp)),
     $           nsize, 1.0d0, scrb(isqrof(irp)), nsize)
            call xgemm('N', 'T', nsize, nsize, nocc(irp,2), 1.0d0,
     $           evec(isqrof(irp)), nsize,
     $           evec(isqrof(irp)), nsize,
     $           1.0d0, scrb(ldim2+isqrof(irp)), nsize)
         endif
      enddo
c     
c     determine trace of current density matrix per symmetry block, alpha/beta
c     
      if (print) then
      do irp = 1, nirrep
         nsize = nbfirr(irp)
         if (nsize .gt. 0) then
            tracea = 0.0d0
            traceb = 0.0d0
            do i = 1, nsize
               tracea = tracea + scrb(isqrof(irp)+(i-1)*nsize+i-1)
            enddo
            do i = 1, nsize
               traceb = traceb + scrb(ldim2+
     $              isqrof(irp)+(i-1)*nsize+i-1)
            enddo
            write(6,*) ' partial traces new density matrix ', tracea,
     $           traceb
         endif
      enddo
      endif
c     
      do i=1, 2*ldim2
         scrb(i) = scrb(i) - densnew(i)
      enddo
c     
      if (sdot(2*ldim2, scrb, 1, scrb, 1) .gt. 1.0d-3) then
         write(6,*) ' @uno_ref: something suspicious '
         write(6,*) ' Two ways of building density matrix do',
     $        ' not agree'
         write(6,*) ' difference matrix'
         do irp = 1, nirrep
            if (nbfirr(irp).ne. 0) then
               nsize=nbfirr(irp)
               call output(scrb(isqrof(irp)), 1, nsize, 1,
     $              nsize, nsize, nsize, 1)
               call output(scrb(ldim2+isqrof(irp)), 1, nsize,
     $              1, nsize, nsize, nsize, 1)
            endif
         enddo
         call errex
      endif
c     
      if (print .and.
     $     icharge .eq. uno_charge .and. imult .eq. uno_mult) then
c     
c     check new density with old density.
c     
         do irp = 1, nirrep
            if (nbfirr(irp).ne. 0) then
               nsize=nbfirr(irp)
               call squez2(densnew(isqrof(irp)),
     $              scrb(itriof(irp)), nsize)
               if (iuhf .ne. 0) then
                  call squez2(densnew(ldim2+isqrof(irp)),
     $                 scrb(ldim1+itriof(irp)), nsize)
               endif
            endif
         enddo
c     
         do i = 1, 2*ldim1
            scrb(i) = scrb(i) - dens(i)
         enddo
c     
         write(6,*) ' difference between new and old density: ',
     $        sdot(2*ldim1, scrb, 1, scrb, 1)
c     
      endif
c     
c     put density matrix in dens
c     
c     
      write(6,*)
      write(6,*) ' Final population in UNO_REF'
      WRITE(LUOUT,5302)(NOCC(I, 1),I=1,NIRREP)
      WRITE(LUOUT,5303)(NOCC(I, 2),I=1,NIRREP)
      CALL PUTREC(20,'JOBARC','OCCUPYA0',NIRREP,NOCC(1,1))
      IF(IUHF.EQ.1)THEN
       CALL PUTREC(20,'JOBARC','OCCUPYB0',NIRREP,NOCC(1,2))
      ENDIF
c     
 5302 FORMAT(T8,'   Alpha population by irrep: ',8(I3,2X))
 5303 FORMAT(T8,'    Beta population by irrep: ',8(I3,2X))
c
      if (.false.) then
c
c copy densnew to dens (original code)
      do irp = 1, nirrep
         if (nbfirr(irp).ne. 0) then
            nsize=nbfirr(irp)
            call squez2(densnew(isqrof(irp)),
     $           dens(itriof(irp)), nsize)
            if (iuhf .ne. 0) then
               call squez2(densnew(ldim2+isqrof(irp)),
     $              dens(ldim1+itriof(irp)), nsize)
            endif
         endif
      enddo
      else
      do irp = 1, nirrep
         if (nbfirr(irp).ne. 0) then
            nsize=nbfirr(irp)
            call squez2(densfrac(isqrof(irp)),
     $           dens(itriof(irp)), nsize)
            if (iuhf .ne. 0) then
               call squez2(densfrac(ldim2+isqrof(irp)),
     $              dens(ldim1+itriof(irp)), nsize)
            endif
         endif
      enddo
      endif
c     
c     construct fock matrix
c     
      if(aofil) then
         ilnbuf=600
         i000 = 1
         i010=i000+itriln(nirrep+1)*iintfp
         i020=i010+iintfp*ilnbuf
         i030=i020+ilnbuf
         i040=i030+nbas*nbas
         if(i040-i000.gt.maxmem) then
            call nomem('make fock matrix',
     &           '{mk.hff} <-- uno_ref <-- vscf',
     &           i040-i000,maxmem)
         endif
         if (iuhf .eq. 0) then
            call mkrhff(fock, dens, icore(i000),oneh,
     &           icore(i010),icore(i020),itriln(nirrep+1),
     &           nbas,nbfirr,icore(i030),ilnbuf,luint,.true.,
     &           naobasfn,0,scfks,scfksexact,scfkslastiter,
     &           V,z1,ksa,ksb,screxc,scr,scr2,valao,valgradao,
     &           totwt,max_angpts,natoms,intnumradpts,ncount,kshf)
         else
            call mkuhff(fock,fock(1+itriln(nirrep+1)),dens,
     &           dens(1+itriln(nirrep+1)),icore(i000),oneh,
     &           icore(i010),icore(i020),itriln(nirrep+1),
     &           nbas,nbfirr,icore(i030),ilnbuf,luint,.true.,
     &           naobasfn,iuhf,scfks,scfksexact,scfkslastiter,
     &           V,z1,ksa,ksb,screxc,scr,scr2,valao,valgradao,
     &           totwt,max_angpts,natoms,intnumradpts,ncount,kshf)
         endif
      else
         i000 = 1
         i010=i000+(iuhf+1)*itriln(nirrep+1)*iintfp
         i020=i010+mxirr2*iintfp
         i030=i020+nbas*nbas*iintfp
         if(i030-i000.gt.maxmem) then
            call nomem('make fock matrix',
     &           '{mkfock} <-- uno_ref <-- vscf',
     &           i030-i000,maxmem)
         endif
c     
         call mkfock(pk,oneh,fock,dens,icore(i000),icore(i010),
     &        icore(i020),ipksiz,itriln(nirrep+1),
     &        mxirr2,nbas,iuhf)
      endif
c     
c     now determine the qrhf energy.  the value for tol is irrelevant, since
c     it is used to determine whether to use a different format statement for
c     the energy.  this section in mkener is skipped for a qrhf case.
c     
      tol=0.0
      dmax=0.0
      iter=0
      noconv=.false.
      i010=i000+itriln(nirrep+1)*iintfp
      i020=i010+mxirr2*iintfp
      i030=i020+mxirr2*iintfp
      if(i030-i000.gt.maxmem) then
         call nomem('qrhf energy','{mkener} <-- uno_ref <-- vscf',
     &        i030-i000,maxmem)
      endif
c     
      call mkener(oneh,dens,fock,icore(i000),icore(i010),icore(i020),
     &     itriln(nirrep+1),mxirr2,repuls,dmax,iter,iuhf,etot,
     &     1,noconv,scfksiter,scfkslastiter)
c     
      write(luout,5000)etot
 5000 format(t3,'The UNO reference energy is ',f20.10,/)
c     
c     at this point, the values needed for the uno_ref calculation have
c     been placed in the appropriate arrays and will be dumped out
c     with everything else in dmpjob.
c     
c     we also need to calculate the new s**2 value, based on the new
c     occupancies.  this reference function is an eigenfunction of
c     spin, so we only determine the uncontaminated value.
c     
      nalpha=uno_closed+uno_open
      nbeta=uno_closed
c     
      savg=abs(nalpha-nbeta)/2.0d0
      s2=savg*(savg+1.0d0)
      amult=sqrt(1.0+4.0*s2)
      write(luout,5100)amult,s2
 5100 format(/,t3,'     the uno average multiplicity is ',f12.7,/,
     &     t3,'the uno expectation value of s**2 is ',f12.7,/)
      call putrec(20,'jobarc','S2SCF   ',iintfp,s2)
c
      if (uno_open .eq. 0) then
         write(6,*)
         write(6,*) ' @uno_ref: The reference state is changed',
     $        ' to closed shell'
cmn         call sscal(ldim2, 2.0d0, densnew, 1)
         iuhf = 1
      else
CSSS         iuhf = 1
      endif
c
c update flags !!
c
CSSS         iflags(28) = uno_charge
CSSS         iflags(29) = uno_mult
CSSS         iflags(11) = iuhf
CSSS         nflags=100
CSSS         call putrec(20, 'JOBARC','IFLAGS  ',nflags,iflags)
c
      return
      end
