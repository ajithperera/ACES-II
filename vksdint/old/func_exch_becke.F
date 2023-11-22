c Becke exchange functional and potential 
c                                          
c       by Stan Ivanov, April 1998              


      subroutine func_exch_becke(func)

#include <aces.h>
#include <const.par>
#include <dftfunc.com>
#include <flags.com>
#include <callstack.com>

      M_REAL
     &    func

      M_REAL
     &    con13,con32,con43,con83,beta
      parameter (con13=one/three,con32=three/two,con43=four/three,
     &    con83=eight/three)
      parameter (beta=0.0042d0)

      M_REAL
     &    asinh,g,ro,x,x2,xa2,xb2,xa,xb,
     &    fbeckea,fbeckeb

      asinh(x)=dlog(x+dsqrt(x**two+one))
      g(x,x2)=beta*x2/(one+six*beta*x*asinh(x))

      callstack_curr='FUNC_EXCH_BECKE'

      func=zero
      fbeckea=zero
      fbeckeb=zero

      call func_exch_lda(func)

      if (roa.ge.thresh) then
        xa2=gro2a*roa**(-con83)
        xa=dsqrt(xa2)
        fbeckea=-((roa**con43)*g(xa,xa2))
      endif

      if (iuhf.eq.0) then
        fbeckeb=fbeckea
      else 
         if (rob.ge.thresh) then
            xb2=gro2b*rob**(-con83)
            xb=dsqrt(xb2)
            fbeckeb=-((rob**con43)*g(xb,xb2))
         endif
      end if

      func=func+fbeckea+fbeckeb

      return
      end


      subroutine pot_exch_becke(pota,potb)

#include <aces.h>
#include <flags.com>
#include <const.par>
#include <dftfunc.com>
#include <callstack.com>

      M_REAL
     &    pota,potb

      M_REAL
     &    potbeckea,potbeckeb,
     &    gradxa,gradya,gradza,hesxxa,hesyya,heszza,hesxya,hesxza,
     &    hesyza,groa,
     &    gradxb,gradyb,gradzb,hesxxb,hesyyb,heszzb,hesxyb,hesxzb,
     &    hesyzb,grob

      callstack_curr='POT_EXCH_BECKE'

      call pot_exch_lda(pota,potb)

      if (roa.ge.thresh) then
       gradxa=half*(gradx+gradxm)
       gradya=half*(grady+gradym)
       gradza=half*(gradz+gradzm)
       hesxxa=half*(hesxx+hesxxm)
       hesyya=half*(hesyy+hesyym)
       heszza=half*(heszz+heszzm)
       hesxya=half*(hesxy+hesxym)
       hesxza=half*(hesxz+hesxzm)
       hesyza=half*(hesyz+hesyzm)
       groa=dsqrt(gro2a)

c Since Becke's exchange functional has two pieces,
c we find the second piece that involves the 
c the gradient corrections in beckepot.
c The first piece is the LDA exchange. 

       call beckepot(potbeckea,roa,groa,gradxa,gradya,gradza,
     &    hesxxa,hesyya,heszza,hesxya,hesxza,hesyza)

       pota=pota+potbeckea
c       if (pota.gt.zero) pota=zero

      end if

      if (iuhf.eq.0) then
        potb=pota
      else        
        if (rob.ge.thresh) then

          gradxb=half*(gradx-gradxm)
          gradyb=half*(grady-gradym)
          gradzb=half*(gradz-gradzm)
          hesxxb=half*(hesxx-hesxxm)
          hesyyb=half*(hesyy-hesyym)
          heszzb=half*(heszz-heszzm)
          hesxyb=half*(hesxy-hesxym)
          hesxzb=half*(hesxz-hesxzm)
          hesyzb=half*(hesyz-hesyzm)
          grob=dsqrt(gro2b)

         call beckepot(potbeckeb,rob,grob,gradxb,gradyb,gradzb,
     &     hesxxb,hesyyb,heszzb,hesxyb,hesxzb,hesyzb)

          potb=potb+potbeckeb
c          if(potb.gt.zero) potb=zero
        end if

      endif

      return
      end


        subroutine beckepot(potbecke,ro,gro,gradx,grady,gradz,
     &    hesxx,hesyy,heszz,hesxy,hesxz,hesyz)

#include <aces.h>
#include <const.par>
#include <callstack.com>

      M_REAL
     &    ro,gro,gradx,grady,gradz,hesxx,hesyy,heszz,hesxz,hesyz,
     &    hesxy,potbecke

        M_REAL
     &    beta,x,x2,f,g,asinh,potbecke1,potbecke2,
     &    con43,con13,con53,con73,t1,t1a,
     &    t2x,t2y,t2z,t2,si1,si2,si,t3x,t3y,t3z,t3,dgx,dgy,dgz

      parameter (con43=four/three,con53=five/three,
     &   con73=seven/three,con13=one/three)
      parameter(beta=0.0042d0)

      asinh(x)=log(x+dsqrt(x**two+one))
      g(x)=(one+six*beta*x*asinh(x))
      f(x)=dsqrt(one+(x*x))

      callstack_curr='BECKEPOT'


          x=gro/(ro**con43)
          x2=x*x

          potbecke1=zero
          potbecke2=zero
          potbecke=zero


c We consider a functional of the form Int f(rho, |gradrho|) dr.
c potbecke1 is df/drho piece
c potbecke2 is grad.[(df/d|gradrho|)(gradrho/|gradrho|)]
c potbecke=potbecke1-potbecke2

       potbecke1=
     & (con43*(ro**con13)*x2*(beta*six*x2-f(x)))/
     & (g(x)*g(x)*f(x))


c     t1a is df/d|gradrho|
       t1a=(-two*x*(-f(x)-three*beta*x*f(x)*asinh(x)+
     & three*beta*x2))/(f(x)*g(x)*g(x))
       t1=(t1a/gro)*(hesxx+hesyy+heszz)

c     (d/dx + d/dy + d/dz) |grad n|
      dgx=(one/gro)*(gradx*hesxx+grady*hesxy+gradz*hesxz)
      dgy=(one/gro)*(gradx*hesxy+grady*hesyy+gradz*hesyz)
      dgz=(one/gro)*(gradx*hesxz+grady*hesyz+gradz*heszz)

      t2x=(-t1a/(gro*gro))*gradx*dgx
      t2y=(-t1a/(gro*gro))*grady*dgy
      t2z=(-t1a/(gro*gro))*gradz*dgz
      t2=t2x+t2y+t2z
 
c       si = d2f/dxd|gradn|
    
        si1=-two/((g(x)*f(x))**three)
        si2=-f(x)-x2*f(x)+18.0d0*beta*x2+
     & 15.0d0*beta*x2*x2+
     & 36.0d0*beta*beta*x2*x*asinh(x)+
     & 18.0d0*beta*beta*x2*x2*x*asinh(x)-
     & 36.0d0*beta*beta*x2*x2*f(x)
       si=si1*si2

       t3x=(si/gro)*gradx*
     &  ((dgx/(ro**con43))-(con43*x/ro)*gradx)
       t3y=(si/gro)*grady*
     &  ((dgy/(ro**con43))-(con43*x/ro)*grady)
       t3z=(si/gro)*gradz*
     &  ((dgz/(ro**con43))-(con43*x/ro)*gradz)
       t3=t3x+t3y+t3z

       potbecke2=t1+t2+t3

       potbecke=-beta*(potbecke1-potbecke2)

      return
      end

c Local Variables: c
c mode: fortran c
c fortran-do-indent: 2 c
c fortran-if-indent: 2 c
c fortran-continuation-indent: 4 c
c fortran-comment-indent-style: nil c
c End: c
