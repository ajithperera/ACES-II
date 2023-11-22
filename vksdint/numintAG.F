
      subroutine numintAG(dograd)

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
     &    znull,pnull,zvalao,zvalgradao
      integer
     &   pmapatom,zgradS,zgradient,equ

      integer ppjj,imap2z,tgradient,NA
      integer basorder,zgradfinal,pchrgRatom

      integer zvalgrad2ao,zgrad2prim
      M_REAL
     &    ftotele,ftottf,ftotw,
     &    factor,totwt,ddot

       M_REAL
     &    exch_energy,corr_energy,exch_corr_energy,energy

      character*2 atmnam,atomsymb

      logical evalpt,init,print_post_ks
      logical dograd
      integer ihfdftgrad
      integer i, j,zdismap,zmapsymatom
      integer cosym,pk,xtm
      integer im2z,zgrad,zgrad2E,zcma,zcmB,zcc,zdd,zgradE

      M_REAL
     &    times(20)

      integer istat, fputc, mod
      character*1 twirl(8)

      integer zw,zdp,zdz,zp,zdmmuij,zdmmuji

      data twirl /'|','/','-','\\','|','/','-','\\'/


      call callstack_push('NUMINT')

c local pointers

      znull  =setptr(1,F_REAL,    1)
      pnull  =setptr(1,F_INTEGER, 1)

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
         zvalao     = setptr(1,F_REAL, naobasfn)
         zvalgradao = setptr(1,F_REAL, 3*naobasfn)

         zcmA=setptr(1,F_REAL,naobasfn*nbastot*iintfp)
         zcmB=setptr(1,F_REAL,naobasfn*nbastot*iintfp)
         zcc=setptr(1,F_REAL,naobasfn*naobasfn*(iuhf+1))
         zdd=setptr(1,F_REAL,naobasfn*nbastot*(iuhf+1))
         zgradE=setptr(1,F_REAL,natoms*3*natoms)
         zgradfinal=setptr(1,F_REAL,natoms*3*natoms)
         zgrad2E=setptr(1,F_REAL,natoms*3)
         zgrad  =setptr(1,F_REAL,natoms*3)
         im2z=setptr(1,F_INTEGER,natoms)
         pmapatom=setptr(1,F_INTEGER,natoms)
         zgradS=setptr(1,F_REAL,natoms*3)

         tgradient=setptr(1,F_REAL,natoms*3)

         zdismap=setptr(1,F_REAL,natoms*natoms)
         zmapsymatom=setptr(1,F_INTEGER,natoms*natoms*natoms)
         pchrgRatom=setptr(1,F_INTEGER,natoms)
         ppjj=setptr(1,F_INTEGER,natoms*natoms)
         imap2z=setptr(1,F_INTEGER,natoms)
         zw       =setptr(1,F_REAL,3*natoms)
         zdp      =setptr(1,F_REAL,3*natoms*natoms)
         zdz      =setptr(1,F_REAL,3*natoms)
         zp       =setptr(1,F_REAL,natoms)
         zdmmuij   =setptr(1,F_REAL,3)
         zdmmuji  =setptr(1,F_REAL,3)

         zgrad2prim=setptr(1,F_REAL,totprim*9)
         zvalgrad2ao=setptr(1,F_REAL,naobasfn*9) 

         basorder=setptr(1,F_INTEGER,naobasfn)
         equ=setptr(1,F_INTEGER,natoms*natoms) 
         zgradient=setptr(1,F_REAL   ,natoms*3)
        
          cosym=setptr(1,F_INTEGER,3)
          pk=setptr(1,F_INTEGER,natoms)
          xtm=setptr(1,F_INTEGER,3) 


       call getrec(1,'JOBARC','CNTERBF0',naobasfn,
     &              kscore(basorder))
       call getrec(1,'JOBARC','CCOEMOA',naobasfn*nbastot*iintfp,
     &              dcore(zcmA))
       call getrec(-1,'JOBARC','CCOEMOB',naobasfn*nbastot*iintfp,
     &              dcore(zcmB))
       call mkdenpra(dcore(zcmA),dcore(zcc),dcore(zxocc),dcore(zdd),
     &               dcore(zcmB),kscore(basorder))


c      call mkpra(dcore(zprimcoef))
       call getrec(20,'JOBARC','HFDFTGRA',1,ihfdftgrad)
       if(ihfdftgrad.eq.1) then   
        do i=1,fun_num_exch
            coef_pot_exch(i)=coef_exch(I)
         end do
         do i=1,fun_num_corr
            coef_pot_corr(i)=coef_corr(i) 
         end do
       end if
      coef_pot_exch(fun_exch_hf)=0.d0
      
      call intpath(kscore(pnshellatom),kscore(pnshellprim),
     &             dcore(zalpha),dcore(zprsqrd),dcore(zvalprim),
     &             dcore(zpcoeff(1)),dcore(zvalmo),
     &             dcore(zpthpt),dcore(zsize),dcore(zptdis),
     &             dcore(zrelsiz),dcore(zbslrd),kscore(pangfct))



       call dzero(dcore(zgradE),natoms*3*natoms)
       call dzero(dcore(zgradS),natoms*3)
       call izero(kscore(zmapsymatom),natoms*natoms*natoms)
       call izero(kscore(ppjj),natoms*natoms)

c Set up ordered list of unique atoms
      ncount=compnorb
      if (c1symmet.eq.1) ncount=natoms
c       ncount=natoms
      do i=1,20
         times(i)=0.d0
      end do
      ftotele = 0.d0
      ftottf  = 0.d0
      ftotw   = 0.d0
        if( natoms .ne. ncount) then
      call equi(kscore(pcompmemb),ncount,kscore(pcomppopv),kscore(equ))
      call dismap(kscore(pcomppopv),kscore(equ),ncount,dcore(zdismap),
     &            dcore(zcoord))
      call symatomap(dcore(zdismap),kscore(pcomppopv),kscore(equ),
     &               ncount,kscore(zmapsymatom),kscore(ppjj),
     &             dcore(zcoord),kscore(cosym),kscore(pk),
     &             kscore(xtm))
       end if
      do ncnt=1,ncount
         nnn=ncnt

c          if( ncnt .eq. 8 ) then
c        Determine center and symmetry factor
c           (number of symmetry equivalent atoms)
         icntr=kscore(polist+ncnt-1)
c         ifctr=kscore(pcomppopv+ncnt-1)
c         factor=dble(ifctr)
         if (c1symmet.eq.1) then
            icntr=ncnt
            ifctr=1
            factor=1.d0
         end if
c        icntr=ncnt
        factor=1.d0

c        Set up grid for integrating only symmetry unique octants
         call octag(kscore(patomchrg),dcore(zatmvc),
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
          call symoctag(kscore(patomchrg), dcore(zatmvc),
     &                     dcore(zrij),      dcore(zaij),
     &                     dcore(zcdnt),     dcore(zrsqrd),
     &                     dcore(zirrtmp),   dcore(zwtintr),
     &                     totwt,            kscore(pradgrid),
     &                     dcore(zradpt),    dcore(zrwt),
     &                     kscore(pintegaxis),dcore(zgridxyz),
     &                     dcore(zgridwt),   kscore(pgrdangpts),
     &                     iradpt,           iangpt,
     &                     grid,             evalpt,dcore(zw),
     &                    dcore(zdp),dcore(zdz),dcore(zp),
     &                    dcore(zdmmuji))

               if (evalpt) then
                  call integAG(kscore(pnshellatom),kscore(pnshellprim),
     &                       dcore(zcdnt),      dcore(zalpha),
     &                       dcore(zrsqrd),     dcore(zvalprim),
     &                       dcore(zpcoeff(1)), dcore(zvalmo),
     &                       totwt,             dcore(zgradprim),
     &                       dcore(zgradmo),    dcore(zxnat),
     &                       dcore(zxocc),      
     &                       factor,            ifctr,
     &                       dcore(zvalao),     dcore(zvalgradao),
     &                       dcore(zprimcoef),  dcore(zcoord),
     &                       kscore(pcompmemb),  kscore(pangfct),
     &                            init,
     &                       iradpt,            iangpt,
     &                       times,             dograd,
     &                       dcore(zcc),dcore(zgradE),kscore(pnaoatom),
     &                       dcore(zw),dcore(zgrad2E),
     &         dcore(zgrad2prim),dcore(zvalgrad2ao),kscore(basorder))
               end if
#ifdef _TWIRL
               istat=fputc(0,'\b')
#endif /* _TWIRL */
c           end do iangpt=1,kscore(pgrdangpts+grid-1)
            end do
c        end do iradpt=1,int_numradpts
         end do

c      end if
        if( natoms .ne. ncount) then
      call gradsymgen(ncnt,icntr,dcore(zgradE),dcore(zcoord),
     &     kscore(pcomppopv),ncount,kscore(equ),
     &  kscore(zmapsymatom))
      end if
c     end do ncnt=1,ncount
      end do

      call gradsymsum(ncount,dcore(zgradE),dcore(zgradfinal),
     &      dcore(zcoord),kscore(pcomppopv),kscore(equ),
     &      kscore(pmapatom),dcore(zgradS))
#ifdef _TWIRL
      istat=fputc(0,'.')
      istat=fputc(0,'\n')
#endif /* _TWIRL */
      CALL GETREC(20,'JOBARC','NATOMS',1,NA)

           call analgrad(dcore(zgradS),dcore(zgradient),
     &          dcore(zgrad),dcore(zgrad2E),kscore(im2z),
     &           kscore(patomchrg),kscore(imap2z),
     &           kscore(pchrgRatom),dcore(tgradient),NA)


      call relptr(1,F_REAL,   znull)
      call relptr(1,F_INTEGER,pnull)

      call callstack_pop
      return
      end

