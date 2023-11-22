subroutine atoi(input,j)
character,intent(in)::input
integer,intent(inout)::j
if(input=='0')then
j=0
elseif(input=='1')then
j=1
elseif(input=='2')then
j=2
elseif(input=='3')then
j=3
elseif(input=='4')then
j=4
elseif(input=='5')then
j=5
elseif(input=='6')then
j=6
elseif(input=='7')then
j=7
elseif(input=='8')then
j=8
elseif(input=='9')then
j=9
end if
end subroutine atoi
