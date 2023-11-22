      subroutine intpath(nshellatom,nshellprim,
     &    alpha,prsqrd,fctden,pcoeff,valmo,
     &    pathpt,size,ptdis,relsiz,bslrd,angfct)

c This routine determines the relative size of the atoms

#include <aces.h>
#include <const.par>
#include <io.par>
#include <intgrtflags.com>
#include <flags.com>
#include <sym.com>
#include <symm2_ks.com>
#include <mol.com>
#include <grid.com>
      integer
     &    nshellatom(natoms),nshellprim(natoms*maxshell),
     &    angfct(numangfct,3)
      M_REAL
     &    alpha(totprim),prsqrd(natoms,natoms,minpt,natoms),
     &    fctden(totprim),
     &    pcoeff(totprim,nbastot,2),
     &    valmo(totocca,2),
     &    pathpt(natoms,natoms,3,minpt,natoms),size(natoms,natoms),
     &    ptdis(natoms,natoms,minpt),relsiz(natoms,natoms),
     &    bslrd(natoms)

      integer
     &    iatom,ifct,iang,imomfct,iatm,jatm,i,ipt
      M_REAL
     &    angvl(8,3),angden,texp,roa,rob,dens,tdens(minpt)

      call callstack_push('INTPATH')

c If using dynamic polyhedra sizing
      if (int_partpoly.eq.int_partpoly_dynamic) then
        
c Loop over pairs of atoms
        do 90 iatm=1,natoms
          do 80 jatm=1,natoms
            
            if(iatm.eq.jatm) goto 80
            dens=1.d+99
            
c Loop over interatomic path
            do 70 ipt=1,minpt
              ifct=1
              do 60 iatom=1,natoms
                
c Angular momentum
                do 10 i=1,3
                  angvl(1,i)=one
                  angvl(2,i)=one
                  angvl(3,i)=pathpt(iatm,jatm,i,ipt,iatom)
   10           continue
                do 30 iang=3,maxangmom+2
                  do 20 i=1,3
                    angvl(iang+1,i)=angvl(iang,i)*
     &                  pathpt(iatm,jatm,i,ipt,iatom)
   20             continue
   30           continue
                
c Loop over angular momentum for each atom and find the density and
c value of certain factors for each function.
                do 50 iang=1,nshellatom(iatom)
                  angden=angvl(angfct(iang,1)+2,1)*
     &                angvl(angfct(iang,2)+2,2)*
     &                angvl(angfct(iang,3)+2,3)
                  
                  do 40 imomfct=1,nshellprim((iatom-1)*maxshell+iang)
                    texp=dexp(-alpha(ifct)*
     &                  prsqrd(iatm,jatm,ipt,iatom))
                    fctden(ifct)=angden*texp
                    ifct=ifct+1
   40             continue
   50           continue
   60         continue
              
c Contract functions with alpha coefficients and add to get alpha MO
c values at this point.  Then form the alpha density at this point.
c Then form the alpha exchange density.
              call xgemm('t','n',totocca,1,totprim,one,pcoeff(1,1,1),
     &            totprim,fctden,totprim,zero,valmo(1,1),totocca)
              call xgemm('t','n',1,1,totocca,one,valmo(1,1),totocca,
     &            valmo(1,1),totocca,zero,roa,1)
              
c Repeat for beta.
              if(iuhf.eq.1)then
                call xgemm('t','n',totoccb,1,totprim,one,pcoeff(1,1,2),
     &              totprim,fctden,totprim,zero,valmo(1,2),totoccb)
                call xgemm('t','n',1,1,totoccb,one,valmo(1,2),totoccb,
     &              valmo(1,2),totoccb,zero,rob,1)
                
              else
                rob=roa
              endif
              
c Find the point of minimum density
              tdens(ipt)=roa+rob
              if(tdens(ipt).lt.dens)then
                dens=tdens(ipt)
                size(iatm,jatm)=ptdis(iatm,jatm,ipt)
              endif
   70       continue
   80     continue
   90   continue
        
c Determine relative size of atoms
        do 110 iatm=1,natoms
          do 100 jatm=1,natoms
            if(iatm.eq.jatm)goto 100
            relsiz(iatm,jatm)=size(iatm,jatm)/
     &          (size(iatm,jatm)+size(jatm,iatm))
  100     continue
  110   continue
        
c Use Bragg-Slater radii for relative size of polyhedra

      elseif(int_partpoly.eq.int_partpoly_bsrad)then
        
        do 130 iatm=1,natoms
          do 120 jatm=1,natoms
            if(iatm.eq.jatm)goto 120
            relsiz(iatm,jatm)= bslrd(iatm)/
     &          (bslrd(iatm)+bslrd(jatm))
  120     continue
  130   continue
        
c Set relative size of atoms to 0.5
      else
        do 150 iatm=1,natoms
          do 140 jatm=1,natoms
            if(iatm.eq.jatm)goto 140
            relsiz(iatm,jatm)= half
  140     continue
  150   continue
        
      endif

  900 continue
      call callstack_pop
      return
      end
