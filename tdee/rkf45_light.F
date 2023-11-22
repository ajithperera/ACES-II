*
* Numerical Recipes subroutine for fourth order Runge-Kutta integration:
*
* Given values for the variables y(1:n) and their derivatives dydx(i:n) 
* known at x, use the fourth-order Runge-Kutta method to advance the 
* solution over an interval h and return the incremented variables as 
* yout(1:n), which need not be a distinct array from y. The user supplies 
* the subroutine derivs(x,y,dydx), which returns derivatives dydx at x.
*
        subroutine rk45_light(derivs,n,yinout,T0,Tn,y)
        implicit none 
        double precision yinout(n),T0,Tn,y(5*n)
        integer n
        integer i000,i010,i020,i030,i040,i050
        external derivs
        
        i000 = 1
        i010 = i000 + n
        i020 = i010 + n
        i030 = i020 + n
        i040 = i030 + n
        i050 = i040 + n

        call dcopy(n,yinout,1,y(i000),1)
        Write(6,*) "The first"
        call prvecr(yinout,n)
        call derivs(t0,yinout,y(i010)) 

        call rk45_worker(y(i000),y(i010),n,t0,tn,yinout,y(i020),
     +                   y(i030),y(i040),derivs)
        
        return 
        end

	subroutine rk45_worker(y,dydx,n,x,h,yout,dym,dyt,yt,derivs)
	implicit none
	integer n,nmax,i
	double precision h,x,dydx(n),y(n),yout(n),s
	external derivs
	double precision h6,hh,xh,dym(n),dyt(n),yt(n)

	hh=h*0.5d0
	h6=h/6d0
	xh=x+hh
        call prvecr(dydx,n)
	do i=1,n
	   yt(i)=y(i)+hh*dydx(i)
	enddo

        call prvecr(yt,n)
        call checksum("YT",yt,n,s)
	call derivs(xh,yt,dyt)
	do i=1,n
	   yt(i)=y(i)+hh*dyt(i)
	enddo

        call prvecr(yt,n)
	call derivs(xh,yt,dym)
	do i=1,n
	   yt(i)=y(i)+h*dym(i)
	   dym(i)=dyt(i)+dym(i)
	enddo

        call prvecr(dyt,n)
	call derivs(x+h,yt,dyt)
	do i=1,n
	   yout(i)=y(i)+h6*(dydx(i)+dyt(i)+2d0*dym(i))
	enddo
        x = h

        call prvecr(yout,n)
	return
	end
