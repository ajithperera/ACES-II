      subroutine integ(nshellatom,nshellprim,cdnt,alpha,rsqrd,valprim,
     &   pcoeff,valmo,totwt,gradprim,hessprim,gradmo,hessmo,xnat,stan,
     &   xocc,ksvxc,moovrlp,aoovrlp,factor,ifctr,valao,xcoeff,
     &   coord,xyzmat,fullmemb,angfct,init,iradpt,iangpt,
     &   times,ksdint,ngrdprim,valngrdao)

c This module determines the MO, density and energies at a given
c grid point

#include <aces.h>
#include <io.par>
#include <flags.com>
#include <intgrtflags.com>
#include <mol.com>
#include <const.par>
#include <grid.com>
#include <int.com>
#include <dftfunc.com>

#include <timeinfo.com>

      integer
     &   nshellatom(natoms),nshellprim(maxshell,natoms),
     &   fullmemb(natoms),ifctr,
     &   angfct(numangfct,3),iradpt,iangpt
      M_REAL
     &   cdnt(natoms,3),alpha(totprim),rsqrd(natoms),
     &   valprim(totprim),pcoeff(totprim,nbastot,2),
     &   valmo(nbastot,2),totwt,
     &   xnat(nbastot),xocc(nbastot,2),factor,
     &   moovrlp(naobasfn,naobasfn,iuhf+1),valao(naobasfn),
     &   xcoeff(naobasfn*totprim),aoovrlp(naobasfn,naobasfn),
     &   coord(3,natoms),xyzmat(naobasfn,naobasfn,3),stan(nbastot)
      M_REAL
     &   ksvxc(naobasfn,naobasfn,iuhf+1),
     &   gradprim(totprim,3),  gradmo(nbastot,3,2),
     &   hessprim(totprim,3,3),hessmo(nbastot,3,3,2)
      M_REAL
     &   ksdint(naobasfn,naobasfn,3,iuhf+1),
     &   ngrdprim(totprim,3),valngrdao(naobasfn,3)

      logical init
      M_REAL
     &   times(20),time

      integer
     &   iatom,ifct,i,iang,imomfct,j,ispin,xyz,xyz2,l,m,n

      M_REAL
     &   xyzn(-2:maxangshell+2,3),denscomp,texp,roab(2),grdaa,grdab,
     &   grdbb,gradcomp(3,2),grad(3,2),twoa,xao,
     &   hesscomp(3,3,3),comp(3),foura,hess(3,3,2),func

      M_REAL
     &   xtmp(fun_num_exch,2),ctmp(fun_num_corr,2),vtmp(2),ddot

      M_REAL
     &   tfkef,wkef,B_DOT

      call callstack_push('INTEG')

      call timer(1)
      time=timenow

      ifct=1
      do 160 iatom=1,natoms

c      o Angular momentum
c        xyzn(i,x) = x**i   if i>=0
c                  = 0      if i<0
         do 110 xyz=1,3
            xyzn(-2,xyz)=zero
            xyzn(-1,xyz)=zero
            xyzn( 0,xyz)=one
  110    continue
         do 130 iang=1,maxangmom+2
            do 120 xyz=1,3
               xyzn(iang,xyz)=xyzn(iang-1,xyz)*cdnt(iatom,xyz)
  120       continue
  130    continue

c        Loop over angular momentum for each atom and
c        find the density and certain factors for each function.
c
c        Each function is of the form:
c           x**l y**m z**n exp(-ar**2)
c
c        The gradient is of the form:
c           [ lx**(l-1) - 2ax**(l+1) ] y**m z**n exp(-ar**2)
c
c        The hessian d**2/dx**2 is of the form:
c           [ l(l-1)x**l-2 - 2a(2l+1)x**l + 4a**2 x**(l+2) ]
c           * y**m z**n exp(-ar**2)
c
c        The hessian d**2/dxdy is of the form:
c           [ lmx**(l-1)y**(m-1) - 2alx**(l-1)y**(m+1)
c             - 2amx**(l+1)y**(m-1) + 4a**2x**(l+1)y**(m+1) ]
c           z**n exp(-ar**2)

         do 150 iang=1,nshellatom(iatom)
            l=angfct(iang,1)
            m=angfct(iang,2)
            n=angfct(iang,3)
            comp(1) = xyzn(m+0,2)*xyzn(n+0,3)
            comp(2) = xyzn(l+0,1)*xyzn(n+0,3)
            comp(3) = xyzn(l+0,1)*xyzn(m+0,2)

c         o Density 0 component
            denscomp= xyzn(l+0,1)*comp(1)

c         o Gradient -1 component
            gradcomp(1,1) = xyzn(l-1,1)*comp(1)*l
            gradcomp(2,1) = xyzn(m-1,2)*comp(2)*m
            gradcomp(3,1) = xyzn(n-1,3)*comp(3)*n

c         o Gradient +1 component
            gradcomp(1,2) = xyzn(l+1,1)*comp(1)
            gradcomp(2,2) = xyzn(m+1,2)*comp(2)
            gradcomp(3,2) = xyzn(n+1,3)*comp(3)

c         o Hessian (xx) -2 component
            hesscomp(1,1,1) = xyzn(l-2,1)*comp(1)*l*(l-1)
            hesscomp(2,2,1) = xyzn(m-2,2)*comp(2)*m*(m-1)
            hesscomp(3,3,1) = xyzn(n-2,3)*comp(3)*n*(n-1)

c         o Hessian (xx) 0 component
            hesscomp(1,1,2) = denscomp*(2*l+1)
            hesscomp(2,2,2) = denscomp*(2*m+1)
            hesscomp(3,3,2) = denscomp*(2*n+1)

c         o Hessian (xx) +2 component
            hesscomp(1,1,3) = xyzn(l+2,1)*comp(1)
            hesscomp(2,2,3) = xyzn(m+2,2)*comp(2)
            hesscomp(3,3,3) = xyzn(n+2,3)*comp(3)

c         o Hessian (xy) -1 -1 component
            hesscomp(1,2,1) = xyzn(l-1,1)*xyzn(m-1,2)*xyzn(n+0,3)*l*m
            hesscomp(1,3,1) = xyzn(l-1,1)*xyzn(n-1,3)*xyzn(m+0,2)*l*n
            hesscomp(2,3,1) = xyzn(m-1,2)*xyzn(n-1,3)*xyzn(l+0,1)*m*n
            hesscomp(2,1,1) = hesscomp(1,2,1)
            hesscomp(3,1,1) = hesscomp(1,3,1)
            hesscomp(3,2,1) = hesscomp(2,3,1)

c         o Hessian (xy) -1 +1 component
            hesscomp(1,2,2) = xyzn(l-1,1)*xyzn(m+1,2)*xyzn(n+0,3)*l
            hesscomp(1,3,2) = xyzn(l-1,1)*xyzn(n+1,3)*xyzn(m+0,2)*l
            hesscomp(2,3,2) = xyzn(m-1,2)*xyzn(n+1,3)*xyzn(l+0,1)*m
            hesscomp(2,1,2) = xyzn(l+1,1)*xyzn(m-1,2)*xyzn(n+0,3)*m
            hesscomp(3,1,2) = xyzn(l+1,1)*xyzn(n-1,3)*xyzn(m+0,2)*n
            hesscomp(3,2,2) = xyzn(m+1,2)*xyzn(n-1,3)*xyzn(l+0,1)*n

c         o Hessian (xy) +1 +1 component
            hesscomp(1,2,3) = xyzn(l+1,1)*xyzn(m+1,2)*xyzn(n+0,3)
            hesscomp(1,3,3) = xyzn(l+1,1)*xyzn(n+1,3)*xyzn(m+0,2)
            hesscomp(2,3,3) = xyzn(m+1,2)*xyzn(n+1,3)*xyzn(l+0,1)
            hesscomp(2,1,3) = hesscomp(1,2,3)
            hesscomp(3,1,3) = hesscomp(1,3,3)
            hesscomp(3,2,3) = hesscomp(2,3,3)

            do 140 imomfct=1,nshellprim(iang,iatom)
               texp=dexp(-alpha(ifct)*rsqrd(iatom))
               twoa=two*alpha(ifct)
               foura=twoa*twoa

               valprim(ifct)=denscomp*texp

               gradprim(ifct,1)=(gradcomp(1,1)-twoa*gradcomp(1,2))*texp
               gradprim(ifct,2)=(gradcomp(2,1)-twoa*gradcomp(2,2))*texp
               gradprim(ifct,3)=(gradcomp(3,1)-twoa*gradcomp(3,2))*texp

c               ngrdprim(ifct,1)=-gradprim(ifct,1)
c               ngrdprim(ifct,2)=-gradprim(ifct,2)
c               ngrdprim(ifct,3)=-gradprim(ifct,3)

               hessprim(ifct,1,1)=(hesscomp(1,1,1)-twoa*hesscomp(1,1,2)+
     &            foura*hesscomp(1,1,3))*texp
               hessprim(ifct,2,2)=(hesscomp(2,2,1)-twoa*hesscomp(2,2,2)+
     &            foura*hesscomp(2,2,3))*texp
               hessprim(ifct,3,3)=(hesscomp(3,3,1)-twoa*hesscomp(3,3,2)+
     &            foura*hesscomp(3,3,3))*texp
               hessprim(ifct,1,2)=(hesscomp(1,2,1)-twoa*hesscomp(1,2,2)-
     &            twoa*hesscomp(2,1,2)+foura*hesscomp(1,2,3))*texp
               hessprim(ifct,1,3)=(hesscomp(1,3,1)-twoa*hesscomp(1,3,2)-
     &            twoa*hesscomp(3,1,2)+foura*hesscomp(1,3,3))*texp
               hessprim(ifct,2,3)=(hesscomp(2,3,1)-twoa*hesscomp(2,3,2)-
     &            twoa*hesscomp(3,2,2)+foura*hesscomp(2,3,3))*texp
               hessprim(ifct,2,1)=hessprim(ifct,1,2)
               hessprim(ifct,3,1)=hessprim(ifct,1,3)
               hessprim(ifct,3,2)=hessprim(ifct,2,3)

               ifct=ifct+1
  140       continue
  150    continue
  160 continue

      call timer(1)
      times(1)=times(1)+(timenow-time)
      time=timenow

      do ispin=1,iuhf+1

c        contract functions with alpha/beta MO coefficients
         call xgemm('t','n',nbastot,1,totprim,one,pcoeff(1,1,ispin),
     &      totprim,valprim,totprim,zero,valmo(1,ispin),nbastot)

c        add to get alpha/beta MO values
         call vecprd(valmo(1,ispin),valmo(1,ispin),xnat,nbastot)

c        form the alpha/beta density at this point
         roab(ispin)=B_DOT(nbastot,xnat,1,xocc(1,ispin),1)
        
         if (roab(ispin).ge.thresh) then

c           contract x/y/z gradients with alpha/beta coefficients and
c           add to get alpha/beta MO x/y/z gradient values at this point
            do xyz=1,3
               call xgemm('t','n',nbastot,1,totprim,
     &                    one,pcoeff(1,1,ispin),totprim,
     &                        gradprim(1,xyz),totprim,
     &                    zero,gradmo(1,xyz,ispin),nbastot)

c              form alpha/beta x/y/z density gradient at this point
               call vecprd(valmo(1,ispin),gradmo(1,xyz,ispin),
     &                     xnat,nbastot)
               grad(xyz,ispin)=two*B_DOT(nbastot,xnat,1,xocc(1,ispin),1)
            end do

c           Hessians are needed only for KS calcs.
            if (int_ks) then
               do xyz=1,3
                  do xyz2=xyz,3
                     call xgemm('t','n',nbastot,1,totprim,
     &                          one,pcoeff(1,1,ispin),totprim,
     &                              hessprim(1,xyz2,xyz),totprim,
     &                          zero,hessmo(1,xyz2,xyz,ispin),nbastot)
                     call vecprd(valmo(1,ispin),
     &                           hessmo(1,xyz2,xyz,ispin),xnat,nbastot)
                     call vecprd(gradmo(1,xyz,ispin),
     &                           gradmo(1,xyz2,ispin),stan,nbastot)
                     call daxpy(nbastot,one,stan,1,xnat,1)
                     hess(xyz2,xyz,ispin)=two*B_DOT(nbastot,
     &                                    xnat,1,xocc(1,ispin),1)
                     hess(xyz,xyz2,ispin)=hess(xyz2,xyz,ispin)
                  end do
               end do
            end if
         else
            roab(ispin)=zero
            do xyz=1,3
               grad(xyz,ispin)=zero
               do xyz2=1,3
                  hess(xyz,xyz2,ispin)=zero
               end do
            end do
         end if
      end do

      if (iuhf.eq.0) then
         roab(2)=roab(1)
         do i=1,3
            grad(i,2)=grad(i,1)
            do j=1,3
               hess(i,j,2)=hess(i,j,1)
            end do
         end do
      end if

      call timer(1)
      times(2)=times(2)+(timenow-time)
      time=timenow

c The common block

      ro     = roab(1)+roab(2)

      if (ro.lt.thresh) go to 600
      roinv  = one/ro
      rom    = roab(1)-roab(2)
      roa    = roab(1)
      rob    = roab(2)
      rs     = (four*pi*ro/three)**(-one/three)
      zeta   = (roa-rob)/ro
      gradx  = grad(1,1) + grad(1,2)
      grady  = grad(2,1) + grad(2,2)
      gradz  = grad(3,1) + grad(3,2)
      gradxm = grad(1,1) - grad(1,2)
      gradym = grad(2,1) - grad(2,2)
      gradzm = grad(3,1) - grad(3,2)
      gro2   = gradx**2 + grady**2 + gradz**2
      gro    = sqrt(gro2)
      gro2a  = grad(1,1)**2 + grad(2,1)**2 + grad(3,1)**2
      gro2b  = grad(1,2)**2 + grad(2,2)**2 + grad(3,2)**2

      hesxx  = hess(1,1,1) + hess(1,1,2)
      hesyy  = hess(2,2,1) + hess(2,2,2)
      heszz  = hess(3,3,1) + hess(3,3,2)
      hesxy  = hess(1,2,1) + hess(1,2,2)
      hesxz  = hess(1,3,1) + hess(1,3,2)
      hesyz  = hess(2,3,1) + hess(2,3,2)
      hesxxm = hess(1,1,1) - hess(1,1,2)
      hesyym = hess(2,2,1) - hess(2,2,2)
      heszzm = hess(3,3,1) - hess(3,3,2)
      hesxym = hess(1,2,1) - hess(1,2,2)
      hesxzm = hess(1,3,1) - hess(1,3,2)
      hesyzm = hess(2,3,1) - hess(2,3,2)
      xlap   = hesxx + hesyy + heszz
      xlapm  = hesxxm + hesyym + heszzm

      grdaa=grad(1,1)*grad(1,1)+grad(2,1)*grad(2,1)+grad(3,1)*grad(3,1)
      grdbb=grad(1,2)*grad(1,2)+grad(2,2)*grad(2,2)+grad(3,2)*grad(3,2)
      grdab=grad(1,1)*grad(1,2)+grad(2,1)*grad(2,2)+grad(3,1)*grad(3,2)

      trm1   = zero
      trm1m  = zero
      trm2   = zero
      trm2m  = zero

c Kohn-Sham potential

c     contract functions with coefficients and
c     add to get AO values at this point.

      call xgemm('t','n',naobasfn,1,totprim,one,xcoeff,totprim,
     &           valprim,totprim,zero,valao,naobasfn)
      do xyz=1,3
         call B_GEMM('t','n',naobasfn,1,totprim,
     &           -one,xcoeff,          totprim,
     &                gradprim(1,xyz), totprim,
     &           zero,valngrdao(1,xyz),naobasfn)
C         call B_GEMM('t','n',naobasfn,1,totprim,
C     &           one, xcoeff,          totprim,
C     &                ngrdprim(1,xyz), totprim,
C     &           zero,valngrdao(1,xyz),naobasfn)
      end do

      xtmp(fun_exch_hf,1)=zero
      xtmp(fun_exch_hf,2)=zero

      if (int_kspot.eq.fun_special) then

         if (coef_pot_exch(fun_exch_lda).ne.zero) then
            call pot_exch_lda(xtmp(fun_exch_lda,1),
     &                        xtmp(fun_exch_lda,2))
         end if

         if (coef_pot_exch(fun_exch_becke).ne.zero) then
            call pot_exch_becke(xtmp(fun_exch_becke,1),
     &                          xtmp(fun_exch_becke,2))
         end if

         if (coef_pot_exch(fun_exch_pbe).ne.zero) then
            call pot_exch_pbe(xtmp(fun_exch_pbe,1),
     &                        xtmp(fun_exch_pbe,2))
         end if

         if (coef_pot_exch(fun_exch_pw91).ne.zero) then
            call pot_exch_pw91(xtmp(fun_exch_pw91,1),
     &                         xtmp(fun_exch_pw91,2))
         end if

         if (coef_pot_corr(fun_corr_vwn).ne.zero) then
            call pot_corr_vwn(ctmp(fun_corr_vwn,1),
     &                        ctmp(fun_corr_vwn,2))
         end if
 
         if (coef_pot_corr(fun_corr_lyp).ne.zero) then
            call pot_corr_lyp(ctmp(fun_corr_lyp,1),
     &                        ctmp(fun_corr_lyp,2))
         end if 

         if (coef_pot_corr(fun_corr_pbe).ne.zero) then
            call pot_corr_pbe(ctmp(fun_corr_pbe,1),
     &                        ctmp(fun_corr_pbe,2))
         end if                                 

         if (coef_pot_corr(fun_corr_pw91).ne.zero) then
            call pot_corr_pw91(ctmp(fun_corr_pw91,1),
     &                         ctmp(fun_corr_pw91,2))
         end if

      else 
         if (int_kspot.eq.fun_hyb_b3lyp) then

            call pot_corr_vwn(ctmp(fun_corr_vwn,1),
     &                        ctmp(fun_corr_vwn,2))
            call pot_corr_lyp(ctmp(fun_corr_lyp,1),
     &                        ctmp(fun_corr_lyp,2))

            call pot_exch_lda(xtmp(fun_exch_lda,1),
     &                        xtmp(fun_exch_lda,2))
            call pot_exch_becke(xtmp(fun_exch_becke,1),
     &                          xtmp(fun_exch_becke,2))
         end if
      end if

      do 550 ispin=1,iuhf+1
         if (roab(ispin).lt.thresh) go to 550
         vtmp(ispin)=
     &        ddot(fun_num_exch,xtmp(1,ispin),1,coef_pot_exch,1)+
     &        ddot(fun_num_corr,ctmp(1,ispin),1,coef_pot_corr,1)
         xao=totwt*factor*vtmp(ispin)
         call dsyr('l',naobasfn,xao,valao,1,
     &             ksvxc(1,1,ispin),naobasfn)
         do xyz=1,3
            call B_GEMM('n','t',naobasfn,naobasfn,1,
     &                  xao,valao,                naobasfn,
     &                      valngrdao(1,xyz),     naobasfn,
     &                  one,ksdint(1,1,xyz,ispin),naobasfn)
         end do

c        Calculate Int vxc_alpha * n_alpha and Int vxc_beta * n_beta
c        needed for the ACTUAL KS SCF Energy. To print KS SCF energy
c        set "ener_ksscf"=on in ZMAT
         if (int_ks_scf) then 
            if (ispin.eq.1)  then 
               vxc_ksalpha=vxc_ksalpha+xao*roab(ispin)
            else
               vxc_ksbeta=vxc_ksbeta+xao*roab(ispin)
            end if
            if (iuhf.eq.0) then
               vxc_ksbeta=vxc_ksalpha
            end if
         end if  
 550  end do

 600  call timer(1)
      times(3)=times(3)+(timenow-time)
      time=timenow

      call callstack_pop
      return
      end

