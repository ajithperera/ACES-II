      subroutine gtao(spin,scr,scr2,cmp2c1,
     &    iordr,scr3,coef,nshellao,nangmom,
     &    scr4,xocc)

      implicit none
#include <aces.h>
#include <iuhf.com>
#include <intgrtflags.com>
#include <machsp.com>
#include <mol.com>
#include <int.com>
      integer
     &    spin,iordr(naobasfn),nangmom(natoms),
     &    nshellao(maxshell,natoms),ii,ij
      M_REAL
     &    scr(naobasfn,nbastot),scr2(naobasfn,nbastot),
     &    cmp2c1(naobasfn*naobasfn),
     &    scr3(nbastot),coef(naobasfn,nbastot),
     &    scr4(nbastot),xocc(nbastot,2)

      integer hh,ih
      integer
     &    iangmom,ifrom,imo,i,j,
     &    istart,iatom,shell,nao,nshell,ito,iao,ishell
      M_REAL
     &    fctr
      character*50 title

      character*1 type(2)
      data type/'A','B'/

      call callstack_push('GTAO')

c Get the AO to MO coefficients and the
c occupation. Convert occupation to real.
c Occupations are calculated in each iteration for KS
c calculation and in the dmpjob for any SCF calculation.
c With these occpuation numbers the density is properly
c symmetrized. 11/12/98 SI

      if(idns.eq.0)then
c SCF density
        call getrec(-1,'JOBARC','EVECAO_'//type(spin),
     &      iintfp*naobasfn*nbastot,coef)
        call getrec(-1,'JOBARC','SCFOCC'//spinc(spin),
     &      nbastot,iordr)
        do imo=1,nbastot
            xocc(imo,spin)=dble(iordr(imo))
        end do

      else
c Correlated density
        fctr=1.d0
        if (iuhf.eq.0) fctr=0.5d0
        call getrec(1,'JOBARC','AONTORB'//type(spin),
     &      iintfp*nbastot*nbastot,scr)
        call getrec(1,'JOBARC','OCCNUM_'//type(spin),
     &      iintfp*nbastot,scr4)
c       Transform to C1 cartesian AO's
        call getrec(1,'JOBARC','CMP2ZMAT',
     &      iintfp*naobasfn*nbastot,cmp2c1)
        call xgemm('n','n',naobasfn,nbastot,nbastot,1.d0,cmp2c1,
     &      naobasfn,scr,nbastot,0.d0,scr2,naobasfn)
c       Order vectors
        call getrec(1,'JOBARC','SCFEVL'//type(spin)//'0',
     &      iintfp*nbastot,scr3)
        do 10 imo=1,nbastot
          iordr(imo)=imo
   10   continue
        call piksr2(nbastot,scr3,iordr)
        do 20 imo=1,nbastot
          xocc(imo,spin)=fctr*scr4(iordr(imo))
          ifrom=iordr(imo)
          call B_COPY(naobasfn,scr2(1,ifrom),1,coef(1,imo),1)
   20   continue
      endif

c Put higher angular momentum vectors in correct order
c Coming in, the p orbitals are ordered (x,y,z,x,y,z,...)
c After reordering, they are ordered (x,x,...,y,y,...,z,z,...)
      istart=1
      do iatom=1,natoms
        shell=1
        do iangmom=1,nangmom(iatom)

          nao=nshellao(shell,iatom)
          if (nao.eq.0) goto 30
          nshell=iangmom*(iangmom+1)/2
          shell=shell+nshell
          ito=istart
          do iao=1,nao
            ifrom=istart+iao-1
            do ishell=1,nshell
              iordr(ito)=ifrom
              ifrom=ifrom+nao
              ito=ito+1
            end do
          end do
          istart=istart+nao*nshell
   30     continue
        end do
        end do

      call mat_reorder(1,0,naobasfn,nbastot,coef,scr2,scr3,iordr)
      call putrec(1,'JOBARC','KSSRTORD',naobasfn,iordr)

      if(spin .eq. 1) then
      call putrec(1,'JOBARC','CCOEMOA',naobasfn*nbastot*iintfp,
     &              coef)
       else  
      call putrec(1,'JOBARC','CCOEMOB',naobasfn*nbastot*iintfp,
     &              coef)
      end if

c Write out the MOs
      if (int_printmos.gt.int_printlev) then
        if (spin.eq.1) then
          title='*** MOs: ALPHA block'
        else
          title='*** MOs: BETA block'
        endif
        call mat_print(coef,naobasfn,nbastot,title,'MO','AO',10,5)
      endif

c Write out the occupancy.
       if (int_printocc.gt.int_printlev) then
        if (spin.eq.1) then
          title='*** Occupancy: ALPHA block'
        else
          title='*** Occupancy: BETA block'
        endif
        call mat_print(xocc(1,spin),1,nbastot,title,'MO','Occ',10,5)
      endif
      
      call callstack_pop
      return
      end
