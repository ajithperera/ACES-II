      Subroutine integAG(nshellatom,nshellprim,cdnt,alpha,
     & rsqrd,valprim,
     &   pcoeff,valmo,totwt,gradprim,gradmo,xnat,
     &   xocc,factor,ifctr,valao,valgradao,xcoeff,
     &   coord,compmemb,angfct,init,iradpt,iangpt,
     &   times,dograd,cc,gradE,naoatom,tgradwt,grad2E,
     &   grad2prim,valgrad2ao,basorder)


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

      integer bgrid, bx,pgrid,ncnt
  
      logical bada,ram 
      integer basorder(naobasfn) 
 
      integer iz,ikm,iao,ip,xxyz,izzz 
      integer
     &    nshellatom(natoms),nshellprim(maxshell,natoms),
     &    compmemb(natoms),ifctr,kol,ii,jj,
     &    angfct(numangfct,3),iradpt,iangpt

      M_REAL
     &    cdnt(natoms,3),alpha(totprim),rsqrd(natoms),
     &    valprim(totprim),pcoeff(totprim,nbastot,2),
     &    valmo(nbastot,2),
     &    totwt,gradprim(totprim,3),
     &    gradmo(nbastot,3,2),xnat(nbastot),
     &    xocc(nbastot,2),factor,
     &    valao(naobasfn),valgradao(naobasfn,3),
     &    fourind((nbastot*(nbastot+1)/2)*
     &            (nbastot*(nbastot+1)/2+1)/2,iuhf+1),
     &    scrnn((nbastot*(nbastot+1))/2),
     &    xcoeff(naobasfn*totprim),
     &    coord(3,natoms)
      M_REAL
     &   cc(naobasfn,naobasfn,iuhf+1),
     &   gradE(natoms,3,natoms),tgradwt(3,natoms),
     &   grad2E(natoms,3)

c      integer nprimao(naobasfn),nprimaol(naobasfn)

      M_REAL
     &   grad2prim(totprim,9),valgrad2ao(naobasfn,9)
c,dprimcoef(totprim,naobasfn),
c     &   dcoeff(naobasfn,nbastot)
      M_REAL pi

      logical init
      logical dograd

      M_REAL
     &    times(20),time

      integer naoatom(natoms),
     &    iatom,ifct,i,iang,imomfct,ispin,xyz,l,m,n,naonao,j

      M_REAL
     &    xyzn(-2:maxangshell+2,3),denscomp,texp,roab(2),
     &    gradcomp(3,2),grad(3,2),twoa,xao(2),
     &    xxderao(2),xcderao(2),xcxc(2,3),
     &    comp(3),foura,func

      M_REAL
     &   xtmp(fun_num_exch,2),ctmp(fun_num_corr,2),
     &   xdertmp(fun_num_exch,2),cdertmp(fun_num_corr,2),
     &   vtmp(2),vxdertmp(2),vcdertmp(2),
     &   abs_grcomp(2),xxdercomp(2),xcdercomp(2),
     &   lyptmp(2)

      M_REAL
     &    tfkef,wkef,B_DOT

      integer ax,agrid,nao,lll,iat,ibasis,jbasis

      M_REAL drxa,drya,drza
      M_REAL graddx,graddy,graddz,dra,dij,gaoxi,gaoyi,gaozi
#ifdef _DEBUG_INTGRT
c      write(*,*) '@INTGRT: Entered ', __FILE__
#endif /* _DEBUG_INTGRT */

      call callstack_push('INTEG')
        ncnt=icntr
      call timer(1)
      time=timenow
    
c      thresh=1.d0-19
      pi = acos(-1.d0)


              do bgrid=1,natoms
              do bx=1,3
               grad2E(bgrid,bx)=0.0d0
               end do
               end do


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

      bada=.true.
      ram= .TRUE.
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

               grad2prim(ifct,1)=(l*(l-1)*xyzn(l-2,1)
     &                          -twoa*(2.d0*l+1)*xyzn(l,1)
     &                          +4.d0*(alpha(ifct)**2.d0)*xyzn(l+2,1))
     &                      *comp(1)*texp

               grad2prim(ifct,4)=(m*(m-1)*xyzn(m-2,2)
     &                          -twoa*(2.d0*m+1)*xyzn(m,2)
     &                          +4.d0*(alpha(ifct)**2.d0)*xyzn(m+2,2))
     &                      *comp(2)*texp

               grad2prim(ifct,5)=(n*(n-1)*xyzn(n-2,3)
     &                          -twoa*(2.d0*n+1)*xyzn(n,3)
     &                          +4.d0*(alpha(ifct)**2.d0)*xyzn(n+2,3))
     &                      *comp(3)*texp

c grad2prim(ifct,2) is dp(y)/dx
      grad2prim(ifct,2)=((m*l*xyzn(l-1,1)-twoa*m*xyzn(l+1,1))*
     &                  xyzn(m-1,2)
     &               +(foura*xyzn(l+1,1)-twoa*l*xyzn(l-1,1))*
     &                 xyzn(m+1,2))*
     &                  xyzn(n,3)*texp

c grad2prim(ifct,3) is dp(z)/dx
      grad2prim(ifct,3)=((n*l*xyzn(l-1,1)-twoa*n*xyzn(l+1,1))*
     &                 xyzn(n-1,3)
     &               +(foura*xyzn(l+1,1)-twoa*l*xyzn(l-1,1))*
     &  xyzn(n+1,3))*
     &                  xyzn(m,2)*texp

C grad2prim(ifct,5) is dp(x)/dy
c         grad2prim(ifct,5)=((m*l*xyzn(m-1,2)-twoa*l*xyzn(m+1,2))*
c     &                 xyzn(l-1,1)
c     &               +(foura*xyzn(m+1,2)-twoa*m*xyzn(m-1,2))*
c     &     xyzn(l+1,1))*
c     &                  xyzn(n,3)*texp

cgrad2prim(ifct,6) is dp(z)/dY
       grad2prim(ifct,6)=((m*n*xyzn(m-1,2)-twoa*n*xyzn(m+1,2))*
     &              xyzn(n-1,3)
     &               +(foura*xyzn(m+1,2)-twoa*m*xyzn(m-1,2))*
     &              xyzn(n+1,3))*
     &                  xyzn(l,1)*texp
c grad2prim(ifct,8) is dp(x)/dz
c      grad2prim(ifct,8)=((n*l*xyzn(n-1,3)-twoa*l*xyzn(n+1,3))*
c     &               xyzn(l-1,1)
c     &               +(foura*xyzn(n+1,3)-twoa*n*xyzn(n-1,3))*
c     &              xyzn(l+1,1))*
c     &                  xyzn(m,2)*texp
c grad2prim(ifct,9) is dp(y)/dz
c          grad2prim(ifct,9)=((n*m*xyzn(n-1,3)-twoa*m*xyzn(n+1,3))*
c     &          xyzn(m-1,2)
c     &               +(foura*xyzn(n+1,3)-twoa*n*xyzn(n-1,3))*
c     &                 xyzn(m+1,2))*
c     &                  xyzn(l,1)*texp

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
         call xgemm('t','n',nbastot,1,totprim,
     &              1.d0, pcoeff(1,1,ispin),totprim,
     &                   valprim,          totprim,
     &              0.d0,valmo(1,ispin),   nbastot)

         call vecprd(valmo(1,ispin),valmo(1,ispin),xnat,nbastot)
         
          roab(ispin)=B_DOT(nbastot,xnat,1,xocc(1,ispin),1)
         
          if (roab(ispin).ge.thresh) then
            do xyz=1,3
               call xgemm('t','n',nbastot,1,totprim,
     &                    1.d0, pcoeff(1,1,ispin),  totprim,
     &                         gradprim(1,xyz),    totprim,
     &                    0.d0,gradmo(1,xyz,ispin),nbastot)
               call vecprd(valmo(1,ispin),gradmo(1,xyz,ispin),
     &                     xnat,nbastot)
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



       call xgemm('t','n',naobasfn,1,totprim,1.d0,xcoeff,totprim,
     &    valprim,totprim,0.d0,valao,naobasfn)

       do xyz=1,3
          call xgemm('t','n',naobasfn,1,totprim,1.d0,xcoeff,totprim,
     &    gradprim(1,xyz),totprim,0.d0,valgradao(1,xyz),naobasfn) 
       end do

       do izzz=1,6
          call xgemm('t','n',naobasfn,1,totprim,1.d0,xcoeff,totprim,
     &    grad2prim(1,izzz),totprim,0.d0,valgrad2ao(1,izzz),naobasfn)
       end do
 
       xtmp(fun_exch_hf,1)=0.d0
       xtmp(fun_exch_hf,2)=0.d0
       xdertmp(fun_exch_hf,1)=0.d0
       xdertmp(fun_exch_hf,2)=0.d0
       xdertmp(fun_exch_lda,1)=0.d0
       xdertmp(fun_exch_lda,2)=0.d0
       cdertmp(fun_corr_vwn,1)=0.d0
       cdertmp(fun_corr_vwn,2)=0.d0

c       if( bada) then
       if (int_kspot.eq.fun_special) then

          if (coef_pot_exch(fun_exch_lda).ne.0.d0) then
              call pot_exch_lda(xtmp(fun_exch_lda,1),
     &                          xtmp(fun_exch_lda,2))
           call func_exch_lda(func)
          do agrid=1,natoms
             grad2E(agrid,1)=
     &             tgradwt(1,agrid)*func
             grad2E(agrid,2)=tgradwt(2,agrid)*func
             grad2E(agrid,3)=tgradwt(3,agrid)*func
          end do
          end if

c prakash
c  below commented statement for HF-DFT gradient calculatio
          if (coef_pot_exch(fun_exch_becke).ne.0.d0) then
              call pot_exch_becke(xtmp(fun_exch_becke,1),
     &                            xdertmp(fun_exch_becke,1),
     &                            xtmp(fun_exch_becke,2),
     &                            xdertmp(fun_exch_becke,2))
             call func_exch_becke(func)
            do bgrid=1,natoms
              do bx=1,3
               grad2E(bgrid,bx)=
     &                tgradwt(bx,bgrid)*func
               end do
               end do
          end if

          if (coef_pot_exch(fun_exch_pbe).ne.0.d0) then
              call pot_exch_pbe(xtmp(fun_exch_pbe,1),
     &                          xdertmp(fun_exch_pbe,1),
     &                          xtmp(fun_exch_pbe,2),
     &                          xdertmp(fun_exch_pbe,2))
         call func_exch_pbe(func)
         do bgrid=1,natoms
              do bx=1,3
               grad2E(bgrid,bx)=
     &                tgradwt(bx,bgrid)*func
               end do
               end do 
          end if

          if (coef_pot_exch(fun_exch_pw91).ne.0.d0) then
              call pot_exch_pw91(xtmp(fun_exch_pw91,1),
     &                           xdertmp(fun_exch_pw91,1),
     &                           xtmp(fun_exch_pw91,2),
     &                           xdertmp(fun_exch_pw91,2))
         call func_exch_pw91(func)
          do bgrid=1,natoms
              do bx=1,3
               grad2E(bgrid,bx)=
     &                tgradwt(bx,bgrid)*func
               end do
               end do 
          end if

          if (coef_pot_corr(fun_corr_vwn).ne.0.d0) then
               call pot_corr_vwn(ctmp(fun_corr_vwn,1),
     &                           ctmp(fun_corr_vwn,2))
              call func_corr_vwn(func)
            do bgrid=1,natoms
              do bx=1,3
               grad2E(bgrid,bx)=
     &                tgradwt(bx,bgrid)*func
               end do
               end do 
          end if

          if (coef_pot_corr(fun_corr_lyp).ne.0.d0) then
                 call pot_corr_lyp(ctmp(fun_corr_lyp,1),
     &                             cdertmp(fun_corr_lyp,1),
     &                             lyptmp(1),
     &                             ctmp(fun_corr_lyp,2),
     &                             cdertmp(fun_corr_lyp,2),
     &                             lyptmp(2))
          call func_corr_lyp(func)
         do bgrid=1,natoms
              do bx=1,3
c              call func_corr_lyp(func)
               grad2E(bgrid,bx)=grad2E(bgrid,bx)
     &                +tgradwt(bx,bgrid)*func
               end do
               end do

          end if 

          if (coef_pot_corr(fun_corr_pbe).ne.0.d0) then
                 call pot_corr_pbe(ctmp(fun_corr_pbe,1),
     &                             cdertmp(fun_corr_pbe,1),
     &                             ctmp(fun_corr_pbe,2),
     &                             cdertmp(fun_corr_pbe,2)) 
          call func_corr_pbe(func)
           do bgrid=1,natoms
              do bx=1,3
               grad2E(bgrid,bx)=grad2E(bgrid,bx)
     &         +tgradwt(bx,bgrid)*func
             end do
          end do 
          end if                                 


          if (coef_pot_corr(fun_corr_pw91).ne.0.d0) then
               call pot_corr_pw91(ctmp(fun_corr_pw91,1),
     &                            cdertmp(fun_corr_pw91,1),
     &                            ctmp(fun_corr_pw91,2),
     &                            cdertmp(fun_corr_pw91,2))    
        call func_corr_pw91(func)
           do bgrid=1,natoms
              do bx=1,3
               grad2E(bgrid,bx)=grad2E(bgrid,bx)
     &                +tgradwt(bx,bgrid)*func
               end do
               end do   
          end if

       else 

          if (int_kspot.eq.fun_hyb_b3lyp) then
          call pot_corr_vwn(ctmp(fun_corr_vwn,1),
     &                      ctmp(fun_corr_vwn,2))
               do agrid=1,natoms
             call func_corr_Vwn(func)
             grad2E(agrid,1)=0.19d0*
     &             tgradwt(1,agrid)*func
             grad2E(agrid,2)=0.19d0*tgradwt(2,agrid)*func
             grad2E(agrid,3)=0.19d0*tgradwt(3,agrid)*func
          end do

          call pot_corr_lyp(ctmp(fun_corr_lyp,1),
     &                      cdertmp(fun_corr_lyp,1),
     &                      lyptmp(1), 
     &                      ctmp(fun_corr_lyp,2),
     &                      cdertmp(fun_corr_lyp,2),
     &                      lyptmp(2))
          do agrid=1,natoms
             call func_corr_lyp(func)
             grad2E(agrid,1)=grad2E(agrid,1)+0.81d0*
     &             tgradwt(1,agrid)*func
             grad2E(agrid,2)=grad2E(agrid,2)+0.81d0*
     &             tgradwt(2,agrid)*func
             grad2E(agrid,3)=grad2E(agrid,3)+0.81d0*
     &             tgradwt(3,agrid)*func
          end do

          call pot_exch_lda(xtmp(fun_exch_lda,1),
     &                      xtmp(fun_exch_lda,2))
                         do agrid=1,natoms
         call func_exch_lda(func)
         grad2E(agrid,1)=grad2E(agrid,1)+0.08d0*
     &             tgradwt(1,agrid)*func
        grad2E(agrid,2)=grad2E(agrid,2)+0.08d0*
     &             tgradwt(2,agrid)*func
        grad2E(agrid,3)=grad2E(agrid,3)+0.08d0*
     &             tgradwt(3,agrid)*func
 
            end do

           call pot_exch_becke(xtmp(fun_exch_becke,1),
     &                         xdertmp(fun_exch_becke,1),
     &                         xtmp(fun_exch_becke,2),
     &                         xdertmp(fun_exch_becke,2))
         do agrid=1,natoms
         call func_exch_becke(func)
         grad2E(agrid,1)=grad2E(agrid,1)+0.72d0*
     &             tgradwt(1,agrid)*func
        grad2E(agrid,2)=grad2E(agrid,2)+0.72d0*
     &             tgradwt(2,agrid)*func
        grad2E(agrid,3)=grad2E(agrid,3)+0.72d0*
     &             tgradwt(3,agrid)*func
       end do 

          end if
       end if

        abs_grcomp(1)=dsqrt(gro2a)
        abs_grcomp(2)=dsqrt(gro2b)

c         coef_pot_exch(fun_exch_becke)=1.d0
c         coef_pot_exch(5)=0.d0 
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

         xao(ispin)=factor*vtmp(ispin)
         xxderao(ispin)=factor*vxdertmp(ispin)
         xcderao(ispin)=factor*vcdertmp(ispin)

c       write(*,*) 'ir=',iradpt,iangpt
       do xyz=1,3
            xxdercomp(ispin)=
     &             xxderao(ispin)*grad(xyz,ispin)/abs_grcomp(ispin)
            xcdercomp(ispin)=
     &             xcderao(ispin)*(grad(xyz,1)+grad(xyz,2))/gro
            xcxc(ispin,xyz)=xxdercomp(ispin)+xcdercomp(ispin)
c            write(*,*)'xcx=',xyz,xcxc(ispin,xyz)
       end do 
         
         if(ram) then 
          nao=0
          lll=0
           do 4009 iat=1,natoms
             nao=naoatom(iat)
c           if(iradpt .eq. 5 .and. iangpt .eq.1)write(*,*) nao,iat 
             if(iat .eq. ncnt) goto 4008
               graddx=0.0d0
               graddy=0.0d0
               graddz=0.0d0
c                do ibasis=1,naobasfn
               do ibasis=1+lll,lll+nao
                 dra=0.0d0
                 drxa=0.0d0
                 drya=0.0d0
                 drza=0.0d0
                    do jbasis=1,naobasfn
                         dij=CC(ibasis,jbasis,ispin)
c               if(iradpt .eq. 5 .and. iangpt .eq.1)write(*,*) dIj
                         DRA=DRA+dij*valao(jbasis)
                         drxa=drxa+dij*valgradao(jbasis,1)
                         drya=drya+dij*valgradao(jbasis,2)
                         drza=drza+dij*valgradao(jbasis,3)
                     end do
                  gaoxi=valgradao(ibasis,1)
                  gaoyi=valgradao(ibasis,2)
                  gaozi=valgradao(ibasis,3)
  
       graddx=gaoxi*dra*xao(ispin)+
     &      (valgrad2ao(ibasis,1)*dra+gaoxi*drxa)*xcxc(ispin,1)+
     &      (valgrad2ao(ibasis,2)*dra+gaoxi*drya)*xcxc(ispin,2)+
     &      (valgrad2ao(ibasis,3)*dra+gaoxi*drza)*xcxc(ispin,3)
 
       graddy=gaoyi*dra*xao(ispin)+
     &      (valgrad2ao(ibasis,2)*dra+gaoyi*drxa)*xcxc(ispin,1)+
     &      (valgrad2ao(ibasis,4)*dra+gaoyi*drya)*xcxc(ispin,2)+
     &      (valgrad2ao(ibasis,6)*dra+gaoyi*drza)*xcxc(ispin,3)

      graddz=gaozi*dra*xao(ispin)+
     &      (valgrad2ao(ibasis,3)*dra+gaozi*drxa)*xcxc(ispin,1)+
     &      (valgrad2ao(ibasis,6)*dra+gaozi*drya)*xcxc(ispin,2)+
     &      (valgrad2ao(ibasis,5)*dra+gaozi*drza)*xcxc(ispin,3)
c        write(*,*)'z=',dra,drza,valgrad2ao(ibasis,7),graddz,ibasis

        gradE(iat,1,ncnt)=gradE(iat,1,ncnt)-2.0d0*totwt*graddx
        gradE(iat,2,ncnt)=gradE(iat,2,ncnt)-2.0d0*totwt*graddy
        gradE(iat,3,ncnt)=gradE(iat,3,ncnt)-2.0d0*totwt*graddz

      gradE(ncnt,1,ncnt)=gradE(ncnt,1,ncnt)+2.0d0*totwt*graddx
      gradE(ncnt,2,ncnt)=gradE(ncnt,2,ncnt)+2.0d0*totwt*graddy
        gradE(ncnt,3,ncnt)=gradE(ncnt,3,ncnt)+2.0d0*totwt*graddz
c       write(*,*)'icnt=',icntr,gradE(icntr,3),icntr
c       write(*,*)'iat=',iat,gradE(iat,3),totwt,graddz
c              end if
              end do
 4008     continue
              lll=lll+nao
 4009     continue
          end if
  500    continue
         call daxpy(natoms*3,1.d0,grad2E,1,gradE(1,1,ncnt),1)

c 500  continue
 

      call timer(1)
      times(3)=times(3)+(timenow-time)
      time=timenow

      call callstack_pop
      return

      end
