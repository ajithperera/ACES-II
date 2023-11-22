      subroutine integ(nshellatom,nshellprim,cdnt,alpha,rsqrd,valprim,
     &    pcoeff,totwt,xnat,
     &    valint,scr,aoovrlp,factor,ifctr,valao,xcoeff,
     &    coord,fullmemb,angfct,iradpt,iangpt,dimaux,
     &    times)


#include <aces.h>
#include <io.par>
#include <flags.com>
#include <intgrtflags.com>
#include <mol.com>
#include <const.par>
#include <grid.com>
#include <maxangshell.par>
#include <timeinfo.com>

      integer
     &    nshellatom(natoms),nshellprim(maxshell,natoms),
     &    fullmemb(natoms),ifctr,
     &    angfct(numangfct,3),iradpt,iangpt,dimaux

      M_REAL
     &    cdnt(natoms,3),alpha(totprim),rsqrd(natoms),
     &    valprim(totprim),pcoeff(totprim,nbastot,2),
     &    totwt,
     &    xnat(nbastot),
     &    valint((naobasfn*(naobasfn+1))/2,dimaux),
     &    scr((naobasfn*(naobasfn+1))/2,dimaux),
     &    aoovrlp(naobasfn,naobasfn),
     &    factor,
     &    valao(naobasfn),
     &    xcoeff(naobasfn*totprim),
     &    coord(3,natoms)

      M_REAL
     &    times(20),time

      integer
     &    iatom,ifct,iang,imomfct,xyz,l,m,n,i,j,naonao

      M_REAL
     &    xyzn(-2:maxangshell+2,3),denscomp,texp,
     &    twoa,xao,
     &    comp(3),foura

      call callstack_push('INTEG')

      naonao=(naobasfn*(naobasfn+1))/2
      call timer(1)
      time=timenow

      ifct=1
      do 160 iatom=1,natoms

c Angular momentum
c   xyzn(i,x) = x**i   if i>=0
c             = 0      if i<0
        do 110 xyz=1,3
          xyzn(-2,xyz)=zero
          xyzn(-1,xyz)=zero
          xyzn( 0,xyz)=one
  110   continue
        do 130 iang=1,maxangmom+2
          do 120 xyz=1,3
            xyzn(iang,xyz)=xyzn(iang-1,xyz)*cdnt(iatom,xyz)
  120     continue
  130   continue

c Loop over angular momentum for each atom and find the density and
c value of certain factors for each function.
c
c Each function is of the form:
c   x**l y**m z**n exp(-ar**2)
c
c The gradient is of the form:
c   [ lx**(l-1) - 2ax**(l+1) ] y**m z**n exp(-ar**2)

        do 150 iang=1,nshellatom(iatom)
          l=angfct(iang,1)
          m=angfct(iang,2)
          n=angfct(iang,3)
          comp(1) = xyzn(m+0,2)*xyzn(n+0,3)
          comp(2) = xyzn(l+0,1)*xyzn(n+0,3)
          comp(3) = xyzn(l+0,1)*xyzn(m+0,2)

c Density 0 component
          denscomp= xyzn(l+0,1)*comp(1)

            do 140 imomfct=1,nshellprim(iang,iatom)
            texp=dexp(-alpha(ifct)*rsqrd(iatom))
            twoa=two*alpha(ifct)
            foura=twoa*twoa

            valprim(ifct)=denscomp*texp
            ifct=ifct+1
  140     continue
  150   continue
  160 continue

      call timer(1)
      times(1)=times(1)+(timenow-time)
      time=timenow

c Contract functions with coefficients and add to get AO values at this point.

       call xgemm('t','n',naobasfn,1,totprim,one,xcoeff,totprim,
     &    valprim,totprim,zero,valao,naobasfn)

       xao=totwt*factor
       m=1
       do i=1,naobasfn
           do j=i,naobasfn
             scr(m,1)=valao(i)*valao(j)
             m=m+1
           end do
       end do

       if (int_overlp.eq.4) then
         call dsyr('l',dimaux,xao,scr(1,1),1,valint,dimaux)
       else
         call dsyr('l',naobasfn,xao,valao,1,aoovrlp,naobasfn)
         call xgemm('n','t',naonao,naobasfn,1,xao,scr,naonao,
     &         valao,naobasfn,one,valint,naonao)
       end if

      call timer(1)
      times(2)=times(2)+(timenow-time)
      time=timenow
      call callstack_pop
      return
      end
