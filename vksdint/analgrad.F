         subroutine analgrad(gradE,
     &     gradient,grad,grad2E,iscr,atomchrg,map2z,
     &     chrgRatom,tgradient,NA)

        implicit none
#include <aces.h>
#include <machsp.com>
#include <iuhf.com>
#include <mol.com>
#include<kscore.com>

        integer gg
      M_REAL
     &   gradE(natoms,3),
     &   gradient(natoms,3),tgradient(natoms,3),
     &   grad(natoms*3),
     &   grad2E(3,natoms),prak
        integer na,nR,dum
         integer atomchrg(natoms),map2z(na)
         integer iz,ia,tmpa,chrgRatom(na)
         integer iscr(natoms)
         integer igrid,ix
         integer i,j,iatom,atom,xyz,ixyz
         integer mm,m,offset,n
       integer zz,kk,jj
      call callstack_push('ANALGRAD')

      call getrec(1,'JOBARC','GRADIENT',natoms*3*iintfp,grad)

CSSS      write(6,*) 'The new gradient:'
CSSS      write(6,*)
CSSS     & 'atm    x               y               z'
      do M=1,natoms
         offset=3*m
         gradient(m,1)=grad(offset-2)
         gradient(m,2)=grad(offset-1)
         gradient(m,3)=grad(offset)
CSSS
CSSS         write(6,100)atomchrg(M), m, grad(offset-2),
CSSS     &                   grad(offset-1),
CSSS     &                   grad(offset)
      end do

      call getrec(1,'JOBARC','MAP2ZMAT',NA,map2z)
c      CALL GETREC(20,'JOBARC','NATOMS',1,NA)
c      CALL GETREC(20,'JOBARC','NREALATM',1,NR)
      call getrec(1,'JOBARC','ATOMCHRG',na,chrgRatom)
          ia=0
         dum=0
         do gg=1,NA
CSSS          write(*,*)gg, chrgRatom(gg),map2z(gg)
            if( chrgRatom(gg) .eq. 0) then
             dum=dum+1
             else 
              ia=ia+1
              map2z(ia)=map2z(ia)-dum
            end  if
         end do

c           write(*,*) ' n=',na,nr,natoms
c      call rotain(grade,dcore(zcoord))
           do ia=1,natoms
              tmpa=map2z(ia)
c               tmpa=ia 
CSSS           write(*,*) 't=',tmpa,ia
CSSS          write(*,*) gradE(ia,1),gradE(ia,2),gradE(ia,3)
               tgradient(ia,1)=gradE(tmpa,1)
               tgradient(ia,2)=gradE(tmpa,2)
               tgradient(ia,3)=gradE(tmpa,3)
           end do

      do j=1,natoms
            do xyz=1,3
               gradient(j,xyz)=
     &          gradient(j,xyz)+
     &         tgradient(j,xyz)
            end do
      end do
        call rotain(gradient,dcore(zcoord))
        
       Write(6,"(a)") " The final gradients in Cartesian frame"
       Write(6,"(a,a,5x,3a)") " Atom charge","     Atom no.", 
     &          " x", "               y", "              z"

       do MM=1,natoms
         offset=3*mm
         grad(offset-2)=gradient(mm,1)    
         grad(offset-1)=gradient(mm,2)
         grad(offset)=gradient(mm,3)

         write(6,100)atomchrg(mm), mm, grad(offset-2),
     &                   grad(offset-1),
     &                   grad(offset)
 100   format     (i3,15x,i3,3f16.10)
       end do
      
      Write(6,*)
      call putrec(1,'JOBARC','GRADIENT',natoms*3*iintfp,grad)

      call callstack_pop
      return
      end
 


         

           
