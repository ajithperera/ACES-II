      subroutine crtgridAG(icntr,atomchrg,atmvc,rij,aij,cdnt,rsqrd,
     &    rrtmp,wtintr,atmwt,xx,gridwt,rwt,iangpt,iradpt,grid,
     &    dw,dp,dz,pp,dmmuji)

c This routine determines the cartesian coordinates for arbitrary grid
c point

      implicit none
#include <aces.h>
#include <intgrtflags.com>
#include <mol.com>
#include <grid.com>
#include <callstack.com>
      integer icntr,atomchrg(natoms)
      M_REAL
     &    cdnt(natoms,3),atmvc(natoms,natoms,3),
     &    rsqrd(natoms),rij(natoms,natoms),
     &    aij(natoms,natoms),wtintr(natoms),
     &    rrtmp(natoms),atmwt,xx(3)

      integer id, i,iatom,jatom
      M_REAL
     &  a,xmmuji, wttot,zmuij,xmuij,fuzzy,cutij


      M_REAL
     &  th,th1

      M_REAL pi,dw(3,natoms),W,DP(3,natoms,natoms),
     &     DZ(3,natoms),
     &       Z, PP(natoms),swt,rwt(int_numradpts,maxangpts),
     &       gridwt(maxangpts,numgrid),g4,ddcutij,ddcutji,
     &  mmuji,one,
     &       dmmuji(3), dcutij,dcutji,dcoeff,B,F2,F4   

       integer iangpt,iradpt,ia,ja,ijatom,jjatom,iiatom,k,grid
       one=1.0d0
       th=1.0d-12
       th1=1.0d-08
       pi=acos(-1.d0)
       callstack_curr='CRTGRID'

       call dzero(dw,natoms*3)
       call  dzero(dz,natoms*3)
       call dzero(dp,natoms*natoms*3)
       do 180 iatom=1,natoms
        rsqrd(iatom)=0.d0
        do 170 i=1,3
          cdnt(iatom,i)=atmvc(icntr,iatom,i)+xx(i)
          rsqrd(iatom)=rsqrd(iatom)+cdnt(iatom,i)**2
  170   continue
        rrtmp(iatom)=dsqrt(rsqrd(iatom))
  180  continue

       if(int_parttyp.eq.rigid)then
        atmwt=1.d0
       else
        wttot=0.d0
        do 200 iatom=1,natoms
          if (atomchrg(iatom).eq.110) then
            wtintr(iatom)=0.d0
            goto 200
          else
            wtintr(iatom)=1.d0
          endif
          do 190 jatom=1,natoms
            if (atomchrg(jatom).eq.110) goto 190
            if (iatom.ne.jatom) then 
            zmuij=(rrtmp(iatom)-rrtmp(jatom))/
     &          rij(iatom,jatom)
            xmuij=zmuij +aij(iatom,jatom)*(1.d0-zmuij**2)
            f4=fuzzy(int_fuzzyiter,xmuij)
            cutij=0.5d0*(1.d0-f4)
            wtintr(iatom)=wtintr(iatom)*cutij
            end if  
  190  continue
          wttot=wttot+wtintr(iatom)
  200  continue
        atmwt=wtintr(icntr)/wttot
       endif

c       if(icntr.eq.2)then
c        write(*,*) 'irad=',iradpt,'ia=',iangpT,'icntr=',icntr
       do 611 ja=1,natoms
         if(abs(aij(ja,1)+1.0d0) .lt. 1.0d-05) then
              pp(ja)=0.0d0
              goto 611
          else
              pp(ja)=1.0d0
         end if
       do 612 ia=1,natoms
c         write(*,*) 'ia=',ia
          if(ia .eq. ja) goto 612
           IF(ABS(AIJ(JA,IA)-1.0D0).LT.1.0D-05) GOTO 612
          mmuji=(rrtmp(ja)-rrtmp(ia))/rij(ja,ia)
          xmmuji=mmuji+aij(ja,ia)*(1.0d0-mmuji*mmuji) 
          F4=xmmuji
          G4=one
            do k=1,int_fuzzyiter
               g4=g4*(one-f4*f4)
               f4=f4*(1.5d0-0.5D0*f4*f4)
            end do
            f2=0.5d00*f4
            dcutij=0.5d00+f2
            dcutji=0.5d00-f2
           pp(ja)=pp(ja)*dcutji
            if(ja .eq. icntr) goto 612
            id=ja
            b=mmuji/(rij(ja,ia))**2
            a=1.0d0-2.0d0*aij(ja,ia)*mmuji
      dmmuji(1)=a*(-cdnt(id,1)/(rrtmp(id)*rij(ja,ia))-B*atmvc(ja,ia,1)) 
      dmmuji(2)=a*(-cdnt(id,2)/(rrtmp(id)*rij(ja,ia))-B*atmvc(ja,ia,2))
      dmmuji(3)=a*(-cdnt(id,3)/(rrtmp(id)*rij(ja,ia))-B*atmvc(ja,ia,3))

c      write(*,*)'d',dmmuji(3),B,atmvc(ja,ia,3),cdnt(ja,3)/rrtmp(ja),ja,ia

               dcoeff=-81.0d0/32.0d0*g4
              if(abs(dcutij).gt. th) then
                  ddcutij=dcoeff/dcutij
                  dp(1,ia,ja)=-ddcutij*dmmuji(1)
                  dp(2,ia,ja)=-ddcutij*dmmuji(2)
                  dp(3,ia,ja)=-ddcutij*dmmuji(3)
c              write(*,*) 'dp=',dp(3,ia,ja),ia,ja
c        write(*,*) 'ij=',dcoeff,ddcutij,'a=',a,'b=',b,ia,ja
              else
                 if(abs(g4) .gt. th1) then
                    call errex
                 endif
              endif
              
              if(abs(dcutji).gt. th) then
                  ddcutji=dcoeff/dcutji
                  dp(1,ja,ja)=dp(1,ja,ja)+ddcutji*dmmuji(1)
                  dp(2,ja,ja)=dp(2,ja,ja)+ddcutji*dmmuji(2) 
                  dp(3,ja,ja)=dp(3,ja,ja)+ddcutji*dmmuji(3)
c             write(*,*) 'dpp=',dp(3,ja,ja)

              else
                  if(abs(g4).gt. th1) then
                    call errex
                  end if
              end if
  612  continue
  611  continue
             Z=0.d0
             do iiatom=1,natoms
                z=z+pp(iiatom)
             end do 
       swt=rwt(iradpt,iangpt)*gridwt(iangpt,grid)
c       write(*,*) 'swt=',swt

            W=Pp(icntr)/Z*swt
c       write(*,*)'w=',W,Z

c      write(*,*) '1=',dp(3,1,1),'2=',dp(3,1,2),'3=',dp(3,1,3)
c      write(*,*) '4=',dp(3,2,1),'5=',dp(3,2,2),'6=',dp(3,2,3)
c      write(*,*) '7=',dp(3,3,1),'8=',dp(3,3,2),'9=',dp(3,3,3)

          
      do 613 jjatom=1,natoms
        if(jjatom .eq. icntr) goto 613
        do  ijatom=1,natoms
        dz(1,jjatom)=dz(1,jjatom)+DP(1,ijatom,jjatom)*Pp(ijatom)
        dz(2,jjatom)=dz(2,jjatom)+DP(2,ijatom,jjatom)*Pp(ijatom)
        dz(3,jjatom)=dz(3,jjatom)+DP(3,ijatom,jjatom)*Pp(ijatom)
        end do
        dw(1,jjatom)=w*(dp(1,icntr,jjatom)-dz(1,jjatom)/z)
        dw(2,jjatom)=w*(dp(2,icntr,jjatom)-dz(2,jjatom)/z)
        dw(3,jjatom)=w*(dp(3,icntr,jjatom)-dz(3,jjatom)/z)

        dw(1,icntr)=dw(1,icntr)-dw(1,jjatom)
        dw(2,icntr)=dw(2,icntr)-dw(2,jjatom)
        dw(3,icntr)=dw(3,icntr)-dw(3,jjatom)
c        write(*,*) 'dw=',dw(1,icntr),dw(2,icntr),dw(3,icntr),icntr
c      write(*,*) 'dww=',dw(1,jjatom),dw(2,jjatom),dw(3,jjatom),jjatom
  613  continue
c       endif
       return
       end
