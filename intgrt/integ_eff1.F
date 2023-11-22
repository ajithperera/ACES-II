      subroutine integ_eff1( 
     &   totwt,
     &   ksvxc,factor,ifctr,valao,valgradao,fourind,scrnn,
     &   init,ncount,iradpt,iangpt,
     &   times,ksdint,dograd,dens,graddens,groupoint,
     &   intnumradpts,max_angpts,ncnt,kshf)


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
     &    ifctr,kol,ii,jj,qq,aa,ss,pp,
     &    iradpt,iangpt,ncnt,ncount

       M_REAL ero, temp,temp1,temp2
       logical kshf
       parameter (ero=0.0D+0)

      integer max_angpts,intnumradpts,i
      M_REAL
     &    totwt(ncount,intnumradpts,max_angpts),
     &    ksvxc(naobasfn,naobasfn,iuhf+1),factor,
     &    valao(naobasfn,max_angpts,intnumradpts,ncount),
     &    valgradao(naobasfn,max_angpts,intnumradpts,ncount,3),
     &    fourind((nbastot*(nbastot+1)/2)*
     &            (nbastot*(nbastot+1)/2+1)/2,iuhf+1),
     &    scrnn((nbastot*(nbastot+1))/2)
      M_REAL
     &   ksdint(naobasfn,naobasfn,3,iuhf+1)
      M_REAL pi

      logical init
      logical dograd

      M_REAL
     &    times(20),time,threshold

      integer
     &    ispin,xyz,naonao,j,m

      M_REAL
     &    roab(2),
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


      integer ikk,ill,ipp,groupoint

      M_REAL
     &     dens(groupoint,2),
     &     graddens(groupoint,3,2)

#ifdef _DEBUG_INTGRT
c      write(*,*) '@INTGRT: Entered ', __FILE__
#endif /* _DEBUG_INTGRT */

      call callstack_push('INTEG')

      call timer(1)
      time=timenow

      pi = acos(-1.d0)
      threshold=10.0**(-12.0)

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
     

       if(totwt(ncnt,iradpt,iangpt).gt.threshold)then 
      

      call timer(1)
      times(2)=times(2)+(timenow-time)
      time=timenow
c The common block
      roab(1)=dens(iangpt,1)
         do ikk=1,3
            grad(ikk,1)=graddens(iangpt,ikk,1)
         end do
    
 
       if(iuhf.eq.0)then
          roab(2)=roab(1)
            do ipp=1,3
              grad(ipp,2)=grad(ipp,1)
             end do
        else
           roab(2)=dens(iangpt,2)
             do ill=1,3
               grad(ill,2)=graddens(iangpt,ill,2)
             end do
        end if     
  
 
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


       xao(ispin)=totwt(ncnt,iradpt,iangpt)*factor*vtmp(ispin)
       xxderao(ispin)=totwt(ncnt,iradpt,iangpt)
     &                *factor*vxdertmp(ispin)
         xcderao(ispin)=totwt(ncnt,iradpt,iangpt)
     &   *factor*vcdertmp(ispin)


        do pp=1,naobasfn
          if(valao(pp,iangpt,iradpt,ncnt) .ne. ero) then
             temp=xao(ispin)*valao(pp,iangpt,iradpt,ncnt)
c             if(abs(temp) .gt. threshold) then 
             do qq=pp,naobasfn
               ksvxc(qq,pp,ispin)=ksvxc(qq,pp,ispin)+
     &          valao(qq,iangpt,iradpt,ncnt)*temp
             end do
c          end if
          end if
        end do





c         call dsyr('l',naobasfn,xao(ispin),valao(1,iangpt,
c     &     iradpt,icntr),1,
c     &                                 ksvxc(1,1,ispin),naobasfn)

         do xyz=1,3
            xxdercomp(ispin)=
     &             xxderao(ispin)*grad(xyz,ispin)/abs_grcomp(ispin)

            xcdercomp(ispin)=
     &             xcderao(ispin)*(grad(xyz,1)+grad(xyz,2))/gro

            xcxc(ispin)=xxdercomp(ispin)+xcdercomp(ispin)
          
         do ss=1,naobasfn
           if((valao(ss,iangpt,iradpt,ncnt) .ne. ero ) .OR.
     &     ( valgradao(ss,iangpt,iradpt,ncnt,xyz) .ne. ero)) then
            temp1=xcxc(ispin)*valgradao(ss,iangpt,iradpt,ncnt,xyz)
            temp2=xcxc(ispin)*valao(ss,iangpt,iradpt,ncnt)
               do aa=ss,naobasfn
                 ksvxc(aa,ss,ispin)=ksvxc(aa,ss,ispin)+
     &           valao(aa,iangpt,iradpt,ncnt)*temp1+
     &           valgradao(aa,iangpt,iradpt,ncnt,xyz)*temp2
               end do
           end if
         end do       


c            call dsyr2('l',naobasfn,xcxc(ispin),valao(1,iangpt,
c     & iradpt,icntr),
c     & 1,
c     & valgradao(1,iangpt,iradpt,icntr,xyz),1,
c     &   ksvxc(1,1,ispin),naobasfn)

         end do 

 500  continue
       else

       if (int_tdks) then

c Reorder and transform to SA AO at each point of the grid. Slow but
c memory more efficient. Note that we have evaluated the MOS at each
c point of the grid in order to calculate the density. We could directly
c get the integrals in MO basis set and eliminate the transformations in
c vtran. It has to be checked and investigated further!!!
 


          do ispin=1,iuhf+1
               call ker_exch_lda(vxdertmp(1),vxdertmp(2))
               call ker_corr_vwn(ro,vcdertmp(1),vcdertmp(2))
               xao(ispin)=totwt(ncnt,iradpt,iangpt)
     &            *factor*(vxdertmp(ispin)+vcdertmp(ispin))
c                xao(ispin)=totwt(ncnt,iradpt,iangpt)*factor

               m=1
               do i=1,naonao
                  do j=i,naonao
                     fourind(m,ispin)=xao(ispin)*scrnn(i)*scrnn(j)
                     m=m+1
                  end do
               end do
 
           end do
        end if

       totele= totele+totwt(ncnt,iradpt,iangpt)*ro

       call func_corr_vwn(func)
       fun_corr(fun_corr_vwn) = fun_corr(fun_corr_vwn) + 
     &    totwt(ncnt,iradpt,iangpt)*func

       call func_corr_lyp(func)
       fun_corr(fun_corr_lyp) = fun_corr(fun_corr_lyp) +
     & totwt(ncnt,iradpt,iangpt)*func

       call func_corr_pbe(func)
       fun_corr(fun_corr_pbe) = fun_corr(fun_corr_pbe) +
     &  totwt(ncnt,iradpt,iangpt)*func

       call func_corr_pw91(func)
       fun_corr(fun_corr_pw91) = fun_corr(fun_corr_pw91) +
     &  totwt(ncnt,iradpt,iangpt)*func

       call func_corr_wl(func)
       fun_corr(fun_corr_wl) = fun_corr(fun_corr_wl) +
     &  totwt(ncnt,iradpt,iangpt)*func

       call func_corr_wi(func)
       fun_corr(fun_corr_wi) = fun_corr(fun_corr_wi) + 
     & totwt(ncnt,iradpt,iangpt)*func

       call func_exch_lda(func)
       fun_exch(fun_exch_lda) = fun_exch(fun_exch_lda) +
     &  totwt(ncnt,iradpt,iangpt)*func

       call func_exch_becke(func)
       fun_exch(fun_exch_becke) = fun_exch(fun_exch_becke) +
     & totwt(ncnt,iradpt,iangpt)*func

       call func_exch_pbe(func)
       fun_exch(fun_exch_pbe) = fun_exch(fun_exch_pbe) +
     &   totwt(ncnt,iradpt,iangpt)*func

       call func_exch_pw91(func)
       fun_exch(fun_exch_pw91) = fun_exch(fun_exch_pw91) + 
     &  totwt(ncnt,iradpt,iangpt)*func

       etottf= etottf+totwt(ncnt,iradpt,iangpt)
     &   *tfkef(roab(1),roab(2),grdaa,grdab,grdbb)
       etotw = etotw +totwt(ncnt,iradpt,iangpt)
     & *wkef(roab(1),roab(2),grdaa,grdab,grdbb)
      end if

      end if
      call timer(1)
      times(3)=times(3)+(timenow-time)
      time=timenow
      call callstack_pop
      return

      end
