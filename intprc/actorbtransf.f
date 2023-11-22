        subroutine actorbtransf(hhhh,hhhl,hhll,hlhl,hlll,llll)
        implicit real*8(a-h,o-z)
        implicit integer*4(i-n)
        real*8 llll
        dimension a(5)
        pi=atan(1d0)*4d0
        write(6,*) '@GMOIAA-I orbital rotation analysis: ',
     +             hhhh,hhhl,hhll,hlhl,hlll,llll
        c=hlll-hhhl
        d=-hlhl+.25d0*hhhh+.25d0*llll-.5*hhll
        a(1)= -atan2(c,d)*.25d0
        a(2)=a(1)+.25d0*pi
        a(3)=a(1)+.5d0*pi
        a(4)=-a(1)
        a(5)=-a(2)
c test second derivative and the transformed integral for the two possibilities
        do i=1,5
         second=cos(4d0*a(i))*d - sin(4d0*a(i))*c
         ss=sin(a(i))
         cc=cos(a(i))
         transformed=(llll+hhhh)*cc*cc*ss*ss +
     +              (hlll-hhhl)*2*cc*ss*(cc*cc-ss*ss)+
     +               0-hhll*2*cc*cc*ss*ss+hlhl*(cc*cc-ss*ss)**2
        if(second.gt.0) then
        write(6,*)'@GMOIAA-I angle , second der, transf exchange: ',
     +       a(i),second,transformed
        goto 10
        endif
        enddo
10      continue
        return
        end





