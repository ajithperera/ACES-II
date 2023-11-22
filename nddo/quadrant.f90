subroutine quadrant(x,y,z,phi,sinth)
use constants
implicit double precision (a-h,o-z)
double precision, intent(inout)::phi,x,y,z
!print*,'in quadrant',phi*180.0d0/pi
!check to see if in x-y plane
!if(abs(sinth).gt.0.99999999999999999999999999d0)then







if(abs(sinth)==one)then
!   print*,'in xy plane'
        if((x.gt.zero).and.(y.gt.zero))then
        return
        elseif((x.lt.zero).and.(y.gt.zero))then
        return
        else
        phi=two*pi-phi
        end if
        


else





   if((x.gt.zero).and.(y.gt.zero).and.(z.gt.zero))then
   return
   elseif((x.lt.zero).and.(y.gt.zero).and.(z.gt.zero))then
   return
   elseif((x.lt.zero).and.(y.lt.zero).and.(z.gt.zero))then
   iquad=3
   elseif((x.gt.zero).and.(y.lt.zero).and.(z.gt.zero))then
   iquad=4
   elseif((x.gt.zero).and.(y.gt.zero).and.(z.lt.zero))then
   iquad=5
   elseif((x.lt.zero).and.(y.gt.zero).and.(z.lt.zero))then
   iquad=6
   elseif((x.lt.zero).and.(y.lt.zero).and.(z.lt.zero))then
   iquad=7
   elseif((x.gt.zero).and.(y.lt.zero).and.(z.lt.zero))then
   iquad=8
   end if


   if(iquad==3)then
   phi=two*pi-phi
   elseif(iquad==4)then
   phi=two*pi-phi
   elseif(iquad==7)then
   phi=two*pi-phi
   elseif(iquad==8)then
   phi=two*pi-phi
   end if


end if



end subroutine quadrant
