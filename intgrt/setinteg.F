      subroutine setinteg(nshellatom,nshellprim,cdnt,alpha,rsqrd,
     &  valprim,
     &   gradprim,
     &  ifctr,
     &   coord,compmemb,angfct,iradpt,iangpt,
     &   times,ncount,valao,valgradao,
     &   dprimcoef,nprimao,nprimaol,ncnt)

      implicit none
#include <aces.h>
#include <maxangshell.par>
#include <iuhf.com>
#include <intgrtflags.com>
#include <mol.com>
#include <grid.com>
#include <int.com>
#include <dftfunc.com>

#include <timeinfo.com>

      integer
     &    nshellatom(natoms),nshellprim(maxshell,natoms),
     &    compmemb(natoms),ifctr,ncount,
     &    angfct(numangfct,3),iradpt,iangpt,ij,ik,ixyz
       
      integer ncnt

      M_REAL
     &    cdnt(natoms,3),alpha(totprim),rsqrd(natoms),
     &    coord(3,natoms)
      M_REAL pi,
     &    valprim(totprim),
     &    gradprim(totprim,3)


      M_REAL
     &    times(20),time

      integer
     &    iatom,ifct,i,iang,imomfct,ispin,xyz,l,m,n,naonao,j

      M_REAL
     &    xyzn(-2:maxangshell+2,3),denscomp,texp,roab(2),
     &    gradcomp(3,2),grad(3,2),twoa,xao(2),
     &    xxderao(2),xcderao(2),xcxc(2),
     &    comp(3),foura,func

      M_REAL
     &   valao(naobasfn,maxangpts,int_numradpts,ncount),
     &   valgradao(naobasfn,maxangpts,int_numradpts,ncount,3),
     &   dprimcoef(totprim,naobasfn)

      integer ii,iz,ip,nprimao(naobasfn),nprimaol(naobasfn),iao

      call callstack_push('SETINTEG')

      call timer(1)
      time=timenow

      pi = acos(-1.d0)


      ifct=1
      do iatom=1,natoms
c        Angular momentum
c                         xyzn(i,x) = x**i   if i>=0
c                                   = 0      if i<0
         do xyz=1,3
            xyzn(-2,xyz)=0.d0
            xyzn(-1,xyz)=0.d0
            xyzn( 0,xyz)=1.d0
         end do
         do iang=1,maxangmom+2
            do xyz=1,3
               xyzn(iang,xyz)=xyzn(iang-1,xyz)*cdnt(iatom,xyz)
            end do
         end do

c Loop over angular momentum for each atom and find the density and
c value of certain factors for each function.
c
c Each function is of the form:
c   x**l y**m z**n exp(-ar**2)
c
c The gradient is of the form:
c   [ lx**(l-1) - 2ax**(l+1) ] y**m z**n exp(-ar**2)
         do iang=1,nshellatom(iatom)
            l=angfct(iang,1)
            m=angfct(iang,2)
            n=angfct(iang,3)
            comp(1) = xyzn(m+0,2)*xyzn(n+0,3)
            comp(2) = xyzn(l+0,1)*xyzn(n+0,3)
            comp(3) = xyzn(l+0,1)*xyzn(m+0,2)

c           Density 0 component
            denscomp= xyzn(l+0,1)*comp(1)
c           Gradient -1 component
            gradcomp(1,1) = xyzn(l-1,1)*comp(1)*l
            gradcomp(2,1) = xyzn(m-1,2)*comp(2)*m
            gradcomp(3,1) = xyzn(n-1,3)*comp(3)*n

c           Gradient +1 component
            gradcomp(1,2) = xyzn(l+1,1)*comp(1)
            gradcomp(2,2) = xyzn(m+1,2)*comp(2)
            gradcomp(3,2) = xyzn(n+1,3)*comp(3)

            do imomfct=1,nshellprim(iang,iatom)
               texp=dexp(-alpha(ifct)*rsqrd(iatom))
               twoa=2.d0*alpha(ifct)
               foura=twoa*twoa

               valprim(ifct)=denscomp*texp
 
               gradprim(ifct,1)=(gradcomp(1,1)-twoa*gradcomp(1,2))*texp
               gradprim(ifct,2)=(gradcomp(2,1)-twoa*gradcomp(2,2))*texp
               gradprim(ifct,3)=(gradcomp(3,1)-twoa*gradcomp(3,2))*texp

               ifct=ifct+1
            end do

c        end do iang=1,nshellatom(iatom)
         end do
c     end do iatom=1,natoms
      end do
      call timer(1)
      times(1)=times(1)+(timenow-time)
      time=timenow
        


            
        
        do ii=1,naobasfn
          valao(ii,iangpt,iradpt,ncnt)=0.D0
             do iz=1,3
          valgradao(ii,iangpt,iradpt,ncnt,iz)=0.D0
             end do
        end do
   

       do iao=1,naobasfn
           do ip=nprimao(iao),nprimaol(iao)
                valao(iao,iangpt,iradpt,ncnt)=
     &          valao(iao,iangpt,iradpt,ncnt)+
     &          dprimcoef(ip,iao)*valprim(ip)
              do xyz=1,3
                valgradao(iao,iangpt,iradpt,ncnt,xyz)=
     &          valgradao(iao,iangpt,iradpt,ncnt,xyz)+
     &          dprimcoef(ip,iao)*gradprim(ip,xyz)
              end do
            end do
       end do


      call callstack_pop
      return

      end
   
