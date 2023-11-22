      subroutine intconst(atomchrg,coord,atmvc,rij,atmpath,ptdis,
     &    prsqrd,pathpt,bslrd,aij)
      
c This routine sets up factors for the numerical integration.
      
#include <aces.h>
#include <const.par>
#include <intgrtflags.com>
#include <mol.com>
#include <grid.com>
      integer atomchrg(natoms)
      M_REAL
     &    coord(3,natoms),atmvc(natoms,natoms,3),
     &    rij(natoms,natoms),atmpath(natoms,natoms,3,minpt),
     &    ptdis(natoms,natoms,minpt),
     &    prsqrd(natoms,natoms,minpt,natoms),
     &    pathpt(natoms,natoms,3,minpt,natoms),
     &    bslrd(natoms),aij(natoms,natoms)

      integer
     &    i,j,k,ipt,iatm
      M_REAL
     &    tdis(3),step(3),td,chi,chi2

      call callstack_push('INTCONST')

c Set up array of relative radii

      if(int_partpoly.eq.int_partpoly_bsrad)then
        do 10 i=1,natoms
          bslrd(i)=xbsl(atomchrg(i))
   10   continue
      else
        do 15 i=1,natoms
          bslrd(i)=one
   15   continue
      endif
      
c Set up atomic vectors
      
      do 30 i=1,natoms
        do 20 j=1,3
          atmvc(i,i,j)=zero
   20   continue
        rij(i,i)=zero
        aij(i,i)=zero
   30 continue
      
      if (natoms.eq.1) goto 900
      
      do 110 i=1,natoms
        do 100 j=1,natoms

          if(i.eq.j)goto 100
          rij(i,j)=zero
          do 40 k=1,3
            atmvc(i,j,k)=coord(k,i)-coord(k,j)
            atmvc(j,i,k)=-atmvc(i,j,k)
            rij(i,j)=rij(i,j)+atmvc(i,j,k)**2
   40     continue

          rij(i,j)=sqrt(rij(i,j))
          rij(j,i)=rij(i,j)

c If using fuzzy partitioning

          if(int_parttyp.eq.int_parttyp_fuzzy)then
            chi=bslrd(i)/bslrd(j)
            chi2=(chi-one)/(chi+one)
            aij(i,j)=chi2/(chi2*chi2-one)
            if(aij(i,j).gt.half)  aij(i,j)=half
            if(aij(i,j).lt.-half) aij(i,j)=-half
            
            chi=bslrd(j)/bslrd(i)
            chi2=(chi-one)/(chi+one)
            aij(j,i)=chi2/(chi2*chi2-one)
            if(aij(j,i).gt.half)  aij(j,i)=half
            if(aij(j,i).lt.-half) aij(j,i)=-half
            
          endif

c If using dynamic polyhedra partitioning

          if(int_partpoly.eq.int_partpoly_dynamic)then
c Set up integration for determining the point at which the electron density
c between two atoms is at a minimum.

            do 50 k=1,3
              tdis(k)=zero
              step(k)=atmvc(j,i,k)/minpt
   50       continue
            
            do 90 ipt=1,minpt
              td=zero
              do 60 k=1,3
                tdis(k)=tdis(k)+step(k)
                td=td+tdis(k)**2
                atmpath(i,j,k,ipt)=coord(k,i)+tdis(k)
   60         continue
              ptdis(i,j,ipt)=sqrt(td)
              do 80 iatm=1,natoms
                prsqrd(i,j,ipt,iatm)=zero
                do 70 k=1,3
                  pathpt(i,j,k,ipt,iatm)=
     &                atmpath(i,j,k,ipt)-coord(k,iatm)
                  prsqrd(i,j,ipt,iatm)=prsqrd(i,j,ipt,iatm)+
     &                pathpt(i,j,k,ipt,iatm)**2
   70           continue
   80         continue
   90       continue

          endif
          
  100   continue
  110 continue
      
  900 continue
      call callstack_pop
      return
      end
