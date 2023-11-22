subroutine normalize(contr,expon,lgto,ngau,ngto)
use constants
implicit double precision(a-h,o-z)
double precision,intent(inout),dimension(:,:)::contr,expon
integer,intent(in),dimension(:)::lgto,ngau
integer,intent(in)::ngto
pi3=pi*pi*pi
do i=1,ngto
       if(lgto(i).eq.1)then
           do j=1,ngau(i)
           contr(j,i)=contr(j,i)*(two*expon(j,i)/pi)**(.75d0)
           end do
       elseif(lgto(i).ge.2 .and. lgto(i).le.4)then
           do j=1,ngau(i)
           contr(j,i)=contr(j,i)*(128.0d0*expon(j,i)**5/pi3)**(.25d0)
           end do 
       elseif(lgto(i).ge.5 .and. lgto(i).le.7)then
           do j=1,ngau(i)
           contr(j,i)=contr(j,i)*(2048.0d0*expon(j,i)**7/9.0d0/pi3)**(.25d0)
           end do 
       elseif(lgto(i).ge.8 .and. lgto(i).le.10)then
          do j=1,ngau(i)
           contr(j,i)=contr(j,i)*(2048.0d0*expon(j,i)**7/pi3)**(.25d0)
           end do
           end if
end do
end subroutine normalize


























