      subroutine integint(nshellatom,nshellprim,cdnt,alpha,
     & rsqrd,valprim,
     &   pcoeff,valmo,totwt,gradprim,gradmo,xnat,
     &   xocc,ksvxc,factor,ifctr,valao,valgradao,fourind,scrnn,xcoeff,
     &   coord,compmemb,angfct,order,init,iradpt,iangpt,
     &   times,ksdint,dograd)

c This module determines the MO, density, energies and
c matrix elements involving the KS potential at a given
c grid point. It features nonstandard way of calculating <i|Vxc|j>.
c If F[n]=int f (n,|gradn|) dr, then <i|Vxc|j>=<i|df/dn|j> +
c int (df/d|gradn) (1/|grad n|) (grad n).grad(i*j) dr.
c 12/25/98 SI

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
     &    compmemb(natoms),ifctr,kol,ii,jj,
     &    angfct(numangfct,3),order(naobasfn),iradpt,iangpt

      M_REAL
     &    cdnt(natoms,3),alpha(totprim),rsqrd(natoms),
     &    valprim(totprim),pcoeff(totprim,nbastot,2),
     &    valmo(nbastot,2),
     &    totwt,gradprim(totprim,3),
     &    gradmo(nbastot,3,2),xnat(nbastot),
     &    xocc(nbastot,2),ksvxc(naobasfn,naobasfn,iuhf+1),factor,
     &    valao(naobasfn),valgradao(naobasfn,3),
     &    fourind((nbastot*(nbastot+1)/2)*
     &            (nbastot*(nbastot+1)/2+1)/2,iuhf+1),
     &    scrnn((nbastot*(nbastot+1))/2),
     &    xcoeff(naobasfn*totprim),
     &    coord(3,natoms)
      M_REAL
     &   ksdint(naobasfn,naobasfn,3,iuhf+1)
      M_REAL pi

      logical init
      logical dograd

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
     &   xtmp(fun_num_exch,2),ctmp(fun_num_corr,2),
     &   xdertmp(fun_num_exch,2),cdertmp(fun_num_corr,2),
     &   vtmp(2),vxdertmp(2),vcdertmp(2),
     &   abs_grcomp(2),xxdercomp(2),xcdercomp(2),
     &   lyptmp(2)

      M_REAL
     &    tfkef,wkef,B_DOT

#ifdef _DEBUG_INTGRT
c      write(*,*) '@INTGRT: Entered ', __FILE__
#endif /* _DEBUG_INTGRT */

      call callstack_push('INTEG')

      call timer(1)
      time=timenow

      pi = acos(-1.d0)

c Initialize integral and functional values

      if ((.not.int_ks).and.init) then
         do i=1,fun_num_exch
            fun_exch(i)=0.d0
         end do
         do i=1,fun_num_corr
            fun_corr(i)=0.d0
         end do
         do i=1,fun_num_hyb
            fun_hyb(i)=0.d0
         end do
         totele = 0.d0
         etottf = 0.d0
         etotw  = 0.d0
         init   = .false.
      end if

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
        
      do ispin=1,iuhf+1
c      o contract functions with alpha/beta MO coefficients
         call xgemm('t','n',nbastot,1,totprim,
     &              1.d0, pcoeff(1,1,ispin),totprim,
     &                   valprim,          totprim,
     &              0.d0,valmo(1,ispin),   nbastot)

c      o add to get alpha/beta MO values
         call vecprd(valmo(1,ispin),valmo(1,ispin),xnat,nbastot)
c      o form the alpha/beta density at this point
         roab(ispin)=B_DOT(nbastot,xnat,1,xocc(1,ispin),1)
         if (roab(ispin).ge.thresh) then
            do xyz=1,3
c            o contract x/y/z gradients with alpha/beta coefficients
               call xgemm('t','n',nbastot,1,totprim,
     &                    1.d0, pcoeff(1,1,ispin),  totprim,
     &                         gradprim(1,xyz),    totprim,
     &                    0.d0,gradmo(1,xyz,ispin),nbastot)
c            o add to get alpha/beta MO x/y/z gradients at this point
               call vecprd(valmo(1,ispin),gradmo(1,xyz,ispin),
     &                     xnat,nbastot)
c            o form alpha/beta x/y/z density gradient at this point
               grad(xyz,ispin)=2.d0*
     &         B_DOT(nbastot,xnat,1,xocc(1,ispin),1)
            end do
         else
            roab(ispin)=0.d0
            do xyz=1,3
               grad(xyz,ispin)=roab(ispin)
            end do
         end if
c     end do ispin=1,iuhf+1
      end do

      if (iuhf.eq.0) then
         roab(2)=roab(1)
         do i=1,3
            grad(i,2)=grad(i,1)
         end do
      end if

      call timer(1)
      times(2)=times(2)+(timenow-time)
      time=timenow

c The common block
         
      ro     = roab(1)+roab(2)
      if (ro.lt.thresh) then
         call timer(1)
         times(3)=times(3)+(timenow-time)
         time=timenow
         call callstack_pop
         return
      end if

      roinv  = 1.d0/ro
      rom    = roab(1)-roab(2)
      roa    = roab(1)
      rob    = roab(2)
      rs     = (4.d0*pi*ro/3.d0)**(-1.d0/3.d0)
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

      grdaa=grad(1,1)*grad(1,1)+grad(2,1)*grad(2,1)+grad(3,1)*grad(3,1)
      grdbb=grad(1,2)*grad(1,2)+grad(2,2)*grad(2,2)+grad(3,2)*grad(3,2)
      grdab=grad(1,1)*grad(1,2)+grad(2,1)*grad(2,2)+grad(3,1)*grad(3,2)

c Kohn-Sham potential
      if (int_ks) then

c Contract functions with coefficients and add to get AO values at this point.

       call xgemm('t','n',naobasfn,1,totprim,1.d0,xcoeff,totprim,
     &    valprim,totprim,0.d0,valao,naobasfn)

c Contract functions with coeff. and add to get the values of gradients of AO
c needed for computing matrix elements with the KS potential when integration
c by parts is used.

       do xyz=1,3
          call xgemm('t','n',naobasfn,1,totprim,1.d0,xcoeff,totprim,
     &    gradprim(1,xyz),totprim,0.d0,valgradao(1,xyz),naobasfn) 
       end do

       xtmp(fun_exch_hf,1)=0.d0
       xtmp(fun_exch_hf,2)=0.d0
       xdertmp(fun_exch_hf,1)=0.d0
       xdertmp(fun_exch_hf,2)=0.d0
       xdertmp(fun_exch_lda,1)=0.d0
       xdertmp(fun_exch_lda,2)=0.d0
       cdertmp(fun_corr_vwn,1)=0.d0
       cdertmp(fun_corr_vwn,2)=0.d0

       if (int_kspot.eq.fun_special) then

c Check what is needed and call it

          if (coef_pot_exch(fun_exch_lda).ne.0.d0) then
              call pot_exch_lda(xtmp(fun_exch_lda,1),
     &                          xtmp(fun_exch_lda,2))
          end if

          if (coef_pot_exch(fun_exch_becke).ne.0.d0) then
              call pot_exch_becke(xtmp(fun_exch_becke,1),
     &                            xdertmp(fun_exch_becke,1),
     &                            xtmp(fun_exch_becke,2),
     &                            xdertmp(fun_exch_becke,2))
          end if

          if (coef_pot_exch(fun_exch_pbe).ne.0.d0) then
              call pot_exch_pbe(xtmp(fun_exch_pbe,1),
     &                          xdertmp(fun_exch_pbe,1),
     &                          xtmp(fun_exch_pbe,2),
     &                          xdertmp(fun_exch_pbe,2))
          end if

          if (coef_pot_exch(fun_exch_pw91).ne.0.d0) then
              call pot_exch_pw91(xtmp(fun_exch_pw91,1),
     &                           xdertmp(fun_exch_pw91,1),
     &                           xtmp(fun_exch_pw91,2),
     &                           xdertmp(fun_exch_pw91,2))
          end if

          if (coef_pot_corr(fun_corr_vwn).ne.0.d0) then
               call pot_corr_vwn(ctmp(fun_corr_vwn,1),
     &                           ctmp(fun_corr_vwn,2))
          end if

          if (coef_pot_corr(fun_corr_lyp).ne.0.d0) then
                 call pot_corr_lyp(ctmp(fun_corr_lyp,1),
     &                             cdertmp(fun_corr_lyp,1),
     &                             lyptmp(1),
     &                             ctmp(fun_corr_lyp,2),
     &                             cdertmp(fun_corr_lyp,2),
     &                             lyptmp(2))
          end if 

          if (coef_pot_corr(fun_corr_pbe).ne.0.d0) then
                 call pot_corr_pbe(ctmp(fun_corr_pbe,1),
     &                             cdertmp(fun_corr_pbe,1),
     &                             ctmp(fun_corr_pbe,2),
     &                             cdertmp(fun_corr_pbe,2)) 
          end if                                 


          if (coef_pot_corr(fun_corr_pw91).ne.0.d0) then
               call pot_corr_pw91(ctmp(fun_corr_pw91,1),
     &                            cdertmp(fun_corr_pw91,1),
     &                            ctmp(fun_corr_pw91,2),
     &                            cdertmp(fun_corr_pw91,2))    
          end if

       else 

          if (int_kspot.eq.fun_hyb_b3lyp) then

c We know exactly what we need and call the components

          call pot_corr_vwn(ctmp(fun_corr_vwn,1),
     &                      ctmp(fun_corr_vwn,2))
          call pot_corr_lyp(ctmp(fun_corr_lyp,1),
     &                      cdertmp(fun_corr_lyp,1),
     &                      lyptmp(1), 
     &                      ctmp(fun_corr_lyp,2),
     &                      cdertmp(fun_corr_lyp,2),
     &                      lyptmp(2))

          call pot_exch_lda(xtmp(fun_exch_lda,1),
     &                      xtmp(fun_exch_lda,2))
           call pot_exch_becke(xtmp(fun_exch_becke,1),
     &                         xdertmp(fun_exch_becke,1),
     &                         xtmp(fun_exch_becke,2),
     &                         xdertmp(fun_exch_becke,2))
          end if
       end if

        abs_grcomp(1)=dsqrt(gro2a)
        abs_grcomp(2)=dsqrt(gro2b)

       do 500 ispin=1,iuhf+1
         if (roab(ispin).lt.thresh) go to 500

         vtmp(ispin)=
     &        B_DOT(fun_num_exch,xtmp(1,ispin),1,coef_pot_exch,1)+
     &        B_DOT(fun_num_corr,ctmp(1,ispin),1,coef_pot_corr,1)

         vxdertmp(ispin)=
     &        B_DOT(fun_num_exch,xdertmp(1,ispin),1,coef_pot_exch,1)
     &        + lyptmp(ispin)*coef_pot_corr(fun_corr_lyp)

         vcdertmp(ispin)=
     &        B_DOT(fun_num_corr,cdertmp(1,ispin),1,coef_pot_corr,1)

         xao(ispin)=totwt*factor*vtmp(ispin)
         xxderao(ispin)=totwt*factor*vxdertmp(ispin)
         xcderao(ispin)=totwt*factor*vcdertmp(ispin)
         call dsyr('l',naobasfn,xao(ispin),valao,1,
     &                                 ksvxc(1,1,ispin),naobasfn)

         do xyz=1,3
            xxdercomp(ispin)=
     &             xxderao(ispin)*grad(xyz,ispin)/abs_grcomp(ispin)

            xcdercomp(ispin)=
     &             xcderao(ispin)*(grad(xyz,1)+grad(xyz,2))/gro

            xcxc(ispin)=xxdercomp(ispin)+xcdercomp(ispin)

            call dsyr2('l',naobasfn,xcxc(ispin),valao,1,
     &             valgradao(1,xyz),1,ksvxc(1,1,ispin),naobasfn)

         end do 
         if (dograd) then
            do xyz=1,3
               call B_GEMM('n','t',naobasfn,naobasfn,1,
     &                     -xao(ispin),valao,                naobasfn,
     &                                 valgradao(1,xyz),     naobasfn,
     &                     1.d0,        ksdint(1,1,xyz,ispin),naobasfn)
            end do 
         end if

 500  continue
 
      else

       if (int_tdks) then
          call xgemm('t','n',naobasfn,1,totprim,1.d0,xcoeff,totprim,
     &         valprim,totprim,0.d0,valao,naobasfn)

c Reorder and transform to SA AO at each point of the grid. Slow but
c memory more efficient. Note that we have evaluated the MOS at each
c point of the grid in order to calculate the density. We could directly
c get the integrals in MO basis set and eliminate the transformations in
c vtran. It has to be checked and investigated further!!!
 
          do i=1,naobasfn
              scrnn(i)=valao(order(i))
          end do

          call xgemm('t','n',nbastot,1,naobasfn,1.d0,ksvxc,naobasfn,
     &         scrnn,naobasfn,0.d0,valao,nbastot)

c Transformations done. Calculate only unique combinations.

          naonao=(nbastot*(nbastot+1))/2

          m=1
          do i=1,nbastot
             do j=i,nbastot
               scrnn(m)=valao(i)*valao(j)
               m=m+1
             end do
          end do

          do ispin=1,iuhf+1
               call ker_exch_lda(vxdertmp(1),vxdertmp(2))
               call ker_corr_vwn(ro,vcdertmp(1),vcdertmp(2))
               xao(ispin)=totwt*factor*(vxdertmp(ispin)+vcdertmp(ispin))
c                xao(ispin)=totwt*factor

               m=1
               do i=1,naonao
                  do j=i,naonao
                     fourind(m,ispin)=xao(ispin)*scrnn(i)*scrnn(j)
                     m=m+1
                  end do
               end do
 
           end do
        end if

       totele= totele+totwt*ro

       call func_corr_vwn(func)
       fun_corr(fun_corr_vwn) = fun_corr(fun_corr_vwn) + totwt*func

       call func_corr_lyp(func)
       fun_corr(fun_corr_lyp) = fun_corr(fun_corr_lyp) + totwt*func

       call func_corr_pbe(func)
       fun_corr(fun_corr_pbe) = fun_corr(fun_corr_pbe) + totwt*func

       call func_corr_pw91(func)
       fun_corr(fun_corr_pw91) = fun_corr(fun_corr_pw91) + totwt*func

       call func_corr_wl(func)
       fun_corr(fun_corr_wl) = fun_corr(fun_corr_wl) + totwt*func

       call func_corr_wi(func)
       fun_corr(fun_corr_wi) = fun_corr(fun_corr_wi) + totwt*func

       call func_exch_lda(func)
       fun_exch(fun_exch_lda) = fun_exch(fun_exch_lda) + totwt*func

       call func_exch_becke(func)
       fun_exch(fun_exch_becke) = fun_exch(fun_exch_becke) + totwt*func

       call func_exch_pbe(func)
       fun_exch(fun_exch_pbe) = fun_exch(fun_exch_pbe) + totwt*func

       call func_exch_pw91(func)
       fun_exch(fun_exch_pw91) = fun_exch(fun_exch_pw91) + totwt*func

       etottf= etottf+totwt*tfkef(roab(1),roab(2),grdaa,grdab,grdbb)
       etotw = etotw +totwt*wkef(roab(1),roab(2),grdaa,grdab,grdbb)

      end if
      call timer(1)
      times(3)=times(3)+(timenow-time)
      time=timenow
      call callstack_pop
      return

      end
