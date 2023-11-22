
      subroutine numintint(dograd)

      implicit none
#include <aces.h>
#include <iuhf.com>
#include <machsp.com>
#include <intgrtflags.com>
#include <flags.com>
#include <kscore.com>
#include <mol.com>
#include <energ.com>
#include <grid.com>
#include <int.com>
#include <dftfunc.com>

      integer ncount,ncnt,nnn,ifctr,zirrtmp,zradpt,zrwt,
     &    zwtintr,zsize,zrelsiz,zrtmp,zrint,setptr,iradpt,iangpt,
     &    grid,zvalmo,zcdnt,zrsqrd,zvalprim,zxnat,zgradprim,zgradmo,
     &    zksvxc,znull,pnull,zvalao,zvalgradao,zscrnn,zfourind,
     &    zorder,rshift
      integer
     &    zksdint, zdens, zgradient, zgradscr, zgradord, ihfdftgrad

      M_REAL
     &    ftotele,ftottf,ftotw,
     &    factor,totwt,ddot

       M_REAL
     &    exch_energy,corr_energy,exch_corr_energy,energy

      character*2 atmnam,atomsymb

      logical evalpt,init,print_post_ks
      logical dograd,hfdftgrad 

      integer i, j, naonao

      M_REAL
     &    times(20)

      integer istat, fputc, mod
      character*1 twirl(8)
      data twirl /'|','/','-','\\','|','/','-','\\'/

#ifdef _DEBUG_INTGRT
      write(*,*) '@INTGRT: Entered ', __FILE__
#endif /* _DEBUG_INTGRT */

      call callstack_push('NUMINT')

c local pointers

      znull  =setptr(1,F_REAL,    1)
      pnull  =setptr(1,F_INTEGER, 1)

c zcdnt     : the coordinates of the current point with respect to each atom
c zrsqrd    : the radial distance sqared from each atom to the current point
c zvalprim  : the value of each primitive function at the current point
c zvalmo    : the value of each alpha/beta MO at the current point
c zvalao    : the value of each AO at the current point
c zvalgradao: the value of the x/y/z gradient component of each AO at
c           : the current point
c zgradprim : the value of the x/y/z gradient component of each primitive
c             function at the current point
c zgradmo   : the value of the alpha/beta MO x/y/z gradient component at the
c             current point

      zcdnt     = setptr(1,F_REAL, natoms*3)
      zrsqrd    = setptr(1,F_REAL, natoms)
      zvalprim  = setptr(1,F_REAL, totprim)
      zvalmo    = setptr(1,F_REAL, nbastot*2)
      zgradprim = setptr(1,F_REAL, totprim*3)
      zgradmo   = setptr(1,F_REAL, nbastot*3*2)

c zwtintr   : An intermediate in calculating the weight due to each atom.
c zradpt    : The array of radial points

      zxnat   = setptr(1,F_REAL, nbastot)
      zirrtmp = setptr(1,F_REAL, natoms)
      zradpt  = setptr(1,F_REAL, int_numradpts)
      zrwt    = setptr(1,F_REAL, int_numradpts*maxangpts)
      zwtintr = setptr(1,F_REAL, natoms)
      zsize   = setptr(1,F_REAL, natoms*natoms)
      zrelsiz = setptr(1,F_REAL, natoms*natoms)
      zrtmp   = setptr(1,F_REAL, natoms)
      zrint   = setptr(1,F_REAL, natoms*maxangpts)

c Pointers 
      zksvxc     = znull
      zfourind   = znull
      zvalao     = znull 
      zvalgradao = znull
      zscrnn     = znull
      zorder     = pnull
      if (int_ks) then
         zksvxc     = setptr(1,F_REAL,naobasfn*naobasfn*(iuhf+1))
         zvalao     = setptr(1,F_REAL, naobasfn)
         zvalgradao = setptr(1,F_REAL, 3*naobasfn)
      else
         if (int_tdks) then
            zvalao   = setptr(1, F_REAL, naobasfn)
            naonao   = rshift((nbastot*(nbastot+1)),1)
            zscrnn   = setptr(1, F_REAL, naonao)
            zfourind = setptr(1, F_REAL, naonao*(naonao+1)*(iuhf+1)/2)
            zksvxc   = setptr(1, F_REAL, naobasfn*naobasfn)
            zorder   = setptr(1, F_INTEGER, naobasfn)
         end if 
      end if

c Kohn-Sham gradient pointers
      if (dograd) then
         zksdint  =setptr(1,F_REAL   ,naobasfn*naobasfn*3*(iuhf+1))
         zdens    =setptr(1,F_REAL   ,naobasfn*naobasfn*(iuhf+1))
         zgradscr =setptr(1,F_REAL   ,naobasfn*naobasfn*3)
         zgradient=setptr(1,F_REAL   ,natoms*3)
         zgradord =setptr(1,F_INTEGER,naobasfn*2)
      else
         zksdint  =znull
         zdens    =znull
         zgradscr =znull
         zgradient=znull
         zgradord =pnull
      end if

      if (memknown.eq.0) goto 999

      if (int_ks) then
         call dzero(dcore(zksvxc),naobasfn*naobasfn*(iuhf+1))
         if (dograd) then
c         o zero the derivative-integral matrix
            call dzero(dcore(zksdint),naobasfn*naobasfn*3*(iuhf+1))
c         o set the HF exchange coefficient to zero since VDINT already
c           calculated its contribution
            coef_pot_exch(fun_exch_hf)=0.d0
         end if
      else
         if (int_tdks) then
            call dzero(dcore(zfourind),naonao*(naonao+1)*(iuhf+1)/2)
            call getrec(1,'JOBARC','KSSRTORD',naobasfn,kscore(zorder))
            call getrec(1,'JOBARC','CMP2ZMAT',naobasfn*nbastot*iintfp,
     &                  dcore(zksvxc))
         end if
      end if

c Determine relative atomic size for the atomic partitioning
      call intpath(kscore(pnshellatom),kscore(pnshellprim),
     &             dcore(zalpha),dcore(zprsqrd),dcore(zvalprim),
     &             dcore(zpcoeff(1)),dcore(zvalmo),
     &             dcore(zpthpt),dcore(zsize),dcore(zptdis),
     &             dcore(zrelsiz),dcore(zbslrd),kscore(pangfct))

c Skip integrating relaxed density if it doesn't exist, used for double
c integration involving CC amplitudes

      if (idns.eq.1.and.iflags(18).eq.0) goto 999

c Loop over integration centers (nuclei) and set the radial grid using
c Slater's rules.

c Set up ordered list of unique atoms
      ncount=compnorb
      if (c1symmet.eq.1) ncount=natoms

      do i=1,20
         times(i)=0.d0
      end do
C
      call getrec(20,'JOBARC','HFDFTGRA',1,ihfdftgrad)
C
      if (ihfdftgrad .eq. 1 ) then
         hfdftgrad=.true.
      else
         hfdftgrad=.false.
      end if
C
      if(hfdftgrad) then
        do i=1,fun_num_exch
            coef_pot_exch(i)=coef_exch(I)
        end do
         do i=1,fun_num_corr
            coef_pot_corr(i)=coef_corr(i)
         end do
        coef_pot_exch(fun_exch_hf)=0.d0
      end if

      ftotele = 0.d0
      ftottf  = 0.d0
      ftotw   = 0.d0
      do ncnt=1,ncount
         nnn=ncnt

c        Determine center and symmetry factor
c           (number of symmetry equivalent atoms)
         icntr=kscore(polist+ncnt-1)
         ifctr=kscore(pcomppopv+ncnt-1)
         factor=dble(ifctr)
         if (c1symmet.eq.1) then
            icntr=ncnt
            ifctr=1
            factor=1.d0
         end if

c        Set up grid for integrating only symmetry unique octants
         call oct(kscore(patomchrg),dcore(zatmvc),
     &            dcore(zrij),     kscore(pradgrid),
     &            dcore(zrelsiz),  dcore(zrtmp),
     &            dcore(zrint),    dcore(zradpt),
     &            dcore(zrwt),     dcore(zbslrd),
     &            dcore(zgridxyz), kscore(pgrdangpts))

         init=.true.
#ifdef _TWIRL
         j=rshift(int_numradpts,3)
         i=0
#endif /* _TWIRL */
         do iradpt=1,int_numradpts
#ifdef _TWIRL
            if (i.eq.j) then
               istat=fputc(0,'.')
               i=0
            end if
            i=i+1
#endif /* _TWIRL */
            grid=kscore(pradgrid+iradpt-1)
            do iangpt=1,kscore(pgrdangpts+grid-1)
#ifdef _TWIRL
               istat=fputc(0,twirl(1+mod(iangpt,8)))
#endif /* _TWIRL */
               evalpt=.true.

               call symoct(kscore(patomchrg), dcore(zatmvc),
     &                     dcore(zrij),      dcore(zaij),
     &                     dcore(zcdnt),     dcore(zrsqrd),
     &                     dcore(zirrtmp),   dcore(zwtintr),
     &                     totwt,            kscore(pradgrid),
     &                     dcore(zradpt),    dcore(zrwt),
     &                     kscore(pintegaxis),dcore(zgridxyz),
     &                     dcore(zgridwt),   kscore(pgrdangpts),
     &                     iradpt,           iangpt,
     &                     grid,             evalpt)

               if (evalpt) then
                  call integint(kscore(pnshellatom),kscore(pnshellprim),
     &                       dcore(zcdnt),      dcore(zalpha),
     &                       dcore(zrsqrd),     dcore(zvalprim),
     &                       dcore(zpcoeff(1)), dcore(zvalmo),
     &                       totwt,             dcore(zgradprim),
     &                       dcore(zgradmo),    dcore(zxnat),
     &                       dcore(zxocc),      dcore(zksvxc),
     &                       factor,            ifctr,
     &                       dcore(zvalao),     dcore(zvalgradao),
     &                       dcore(zfourind),   dcore(zscrnn),
     &                       dcore(zprimcoef),  dcore(zcoord),
     &                       kscore(pcompmemb),  kscore(pangfct),
     &                       kscore(zorder),     init,
     &                       iradpt,            iangpt,
     &                       times,             dcore(zksdint),dograd)
               end if
#ifdef _TWIRL
               istat=fputc(0,'\b')
#endif /* _TWIRL */
c           end do iangpt=1,kscore(pgrdangpts+grid-1)
            end do
c        end do iradpt=1,int_numradpts
         end do

c Print the value of the integrated density

         if (.not.int_ks) then
            if (int_printatom.gt.int_printlev) then
               atmnam=atomsymb(kscore(patomchrg+icntr-1))
               write(*,2000) icntr,atmnam,ifctr
               write(*,2001) totele
            end if
         end if

 2000    format(/'Atom',i3,' is ',a2,', it involves ',i3,
     &          ' equivalent atoms')
 2001    format('   The electron density                   =',
     &          f20.12,' electrons')

c        Add atomic result to total

         ftotele=ftotele+factor*totele
         do i=1,fun_num_exch-1
            tot_exch(i)=tot_exch(i)+factor*fun_exch(i)
         end do
         do i=1,fun_num_corr
            tot_corr(i)=tot_corr(i)+factor*fun_corr(i)
         end do

         ftottf=ftottf+factor*etottf
         ftotw=ftotw+factor*etotw

c     end do ncnt=1,ncount
      end do
#ifdef _TWIRL
      istat=fputc(0,'.')
      istat=fputc(0,'\n')
#endif /* _TWIRL */

c The energy corresponding to fun_exch_hf is ehfx

      tot_exch(fun_exch_hf)=ehfx

c Calculate hybrid values

      tot_hyb(fun_hyb_b3lyp)=
     &      (b3lypa*ehfx)+
     &      (1.d0-b3lypa)*tot_exch(fun_exch_lda)+
     &      b3lypb*(tot_exch(fun_exch_becke)-tot_exch(fun_exch_lda))+
     &      tot_corr(fun_corr_vwn)+
     &      b3lypc*(tot_corr(fun_corr_lyp)-tot_corr(fun_corr_vwn))

      if (int_ks) then
         if (dograd) then
            call dintproc(
     &              dcore(zksdint),
     &              dcore(zdens),
     &              dcore(zgradscr),
     &              dcore(zgradient),
     &              naobasfn,nbastot,
     &              dcore(zgradscr+naobasfn*naobasfn),
     &              kscore(zgradord),
     &              kscore(zgradord+naobasfn))
         else
            call putrec(1,'JOBARC','KSTOTELE',iintfp,ftotele)
            call putrec(1,'JOBARC','KSPOTENT',
     &          naobasfn*naobasfn*(iuhf+1)*iintfp,dcore(zksvxc))
         end if
      end if

c Print the value of the integrated density and functionals
      if (.not.int_ks) then

c         if (int_tdks) then
c We shall dump the integrals for needed for vtran
c        HERE!!!!!!!
c         end if

         write(*,2100) ftotele
c        write(*,2107) ftottf
c        write(*,2108) ftotw

         do i=1,fun_num_exch
            write(*,2121) nam_exch(i),tot_exch(i)
         end do
         do i=1,fun_num_corr
            write(*,2120) nam_corr(i),tot_corr(i)
         end do
         do i=1,fun_num_hyb
            write(*,2122) nam_hyb(i),tot_hyb(i)
         end do

c         print *, 'DEBUG: POTENTIAL/FUNCTIONAL COEFFICIENTS'
c         do i=1,fun_num_exch
c            write(*,2123) nam_exch(i),coef_pot_exch(i),coef_exch(i)
c         end do
c         do i=1,fun_num_corr
c            write(*,2123) nam_corr(i),coef_pot_corr(i),coef_corr(i)
c         end do

         write(*,2199)

c if a correlated density then read the Hartree energy from
c JOBARC (it was written to JOBARC from the first call to 
c intgrt; just after SCF)

CSSS         If (idns .eq. 1) call getrec(1,'JOBARC','HARTENRG',
CSSS     &                               iintfp,ehar)

         exch_energy=ddot(fun_num_exch,tot_exch,1,coef_pot_exch,1)
         corr_energy=ddot(fun_num_corr,tot_corr,1,coef_pot_corr,1)
         exch_corr_energy=exch_energy+corr_energy
         energy=ehar+exch_corr_energy

c Check if the density is from a KS calculation.
c If so, print the energies from the KS calculation.

         call getrec(1,'JOBARC','KSPRINT',iintfp,print_post_ks)

         if (print_post_ks) then
C#ifdef _DEBUG
            if (int_kspot.eq.fun_hyb_b3lyp) write(*,6664)

            write(*,6666) ehar
 
            write(*,6667) exch_energy
            write(*,6668) corr_energy
            write(*,6665) exch_corr_energy
            write(*,6669) energy
C#endif
         end if

         if (int_dft_fun.gt.fun_dft_none) then
            exch_energy=0.d0
            corr_energy=0.d0
            exch_corr_energy=0.d0
            energy=0.d0
   
            exch_energy=ddot(fun_num_exch,tot_exch,1,coef_exch,1)
            corr_energy=ddot(fun_num_corr,tot_corr,1,coef_corr,1)
   
            exch_corr_energy=exch_energy+corr_energy
            energy=ehar+exch_corr_energy 
   
            If (idns .eq. 0) Then
                write(*,6676) ehar
                write(*,6677) exch_energy
                write(*,6678) corr_energy
                write(*,6675) exch_corr_energy
                write(*,6679) energy  
            Else
                write(*,6686) ehar
                write(*,6687) exch_energy
                write(*,6688) corr_energy
                write(*,6685) exch_corr_energy
                write(*,6689) energy  
            Endif 

         endif

         call putrec(1,'JOBARC','TOTENERG',iintfp,energy)
         If (int_dft_fun.gt.fun_dft_none .and. idns .eq. 0)
     &       call putrec(1,'JOBARC','HARTENRG',iintfp,ehar)
c  
c      o write out the HF exchange coefficient for use in VDINT
         call putrec(1,'JOBARC','PCHF',
     &       iintfp,coef_pot_exch(fun_exch_hf))

 2100 format(/'   Total density integrates to             : ',
     &    f20.12,' electrons')
 2107 format('The total TF kinetic energy                  =',
     &    f20.12,' a.u.')
 2108 format('The total Weizsacker kinetic energy          =',
     &    f20.12,' a.u.')

 2120 format('   Corr tot contrib : ',A20,' : ',f20.12)
 2121 format('   Exch tot contrib : ',A20,' : ',f20.12)
 2122 format('   Hybr tot contrib : ',A20,' : ',f20.12)
 2123 format(A20,' : ',f20.12,', ',f20.12)

c Print some total energies
 2320 format('   Corr total energy: ',A20,' : ',f20.12)
 2321 format('   Exch total energy: ',A20,' : ',f20.12)
 2322 format('   Hybr total energy: ',A20,' : ',f20.12)

 6664 format('Gaussian94-like B3LYP with its SCF Density (Full VWN)'
     &,/,
     &t1,'The original B3LYP corresponds to func=b3lyp, kspot=lda,vwn'/) 
 6665 format(' Exc  from   KS   Calculation:     ',f20.12,' a.u.')
 6666 format(' Ehar from   KS   Calculation:     ',f20.12,' a.u.')
 6667 format(' Ex   from   KS   Calculation:     ',f20.12,' a.u.')
 6668 format(' Ec   from   KS   Calculation:     ',f20.12,' a.u.')
 6669 format(' Final Result from KS Calculation: ',f20.12,' a.u.'/)

 6675 format(' Exc  with the SCF Dens and Func:  ',f20.12,' a.u.')
 6676 format(' Ehar with the SCF Dens and Func:  ',f20.12,' a.u.')
 6677 format(' Ex   with the SCF Dens and Func:  ',f20.12,' a.u.')
 6678 format(' Ec   with the SCF Dens and Func:  ',f20.12,' a.u.')
 6679 format(' Final Result with the SCF Dens and Func:',
     &f20.12,' a.u.'/)

 6685 format(' Exc  with the CC Dens and Func:  ',f20.12,' a.u.')
 6686 format(' Ehar with the CC Dens and Func:  ',f20.12,' a.u.')
 6687 format(' Ex   with the CC Dens and Func:  ',f20.12,' a.u.')
 6688 format(' Ec   with the CC Dens and Func:  ',f20.12,' a.u.')
 6689 format(' Final Result with the SCF Dens and Func:',
     &f20.12,' a.u.'/)

      end if

 2199 format('     ')

  999 continue
      call relptr(1,F_REAL,   znull)
      call relptr(1,F_INTEGER,pnull)

      call callstack_pop
      return
      end

