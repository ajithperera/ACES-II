subroutine initHess(nat,Hinit)
!use glovar
implicit none

integer,intent(in) :: nat
real*8,dimension(:),allocatable,intent(out) :: Hinit
integer :: ij,ios,i,l
real*8 :: bondl,angle,dihed
character(len=2) :: symbol
character(len=4) :: strangle,strdihed
integer,dimension(:),allocatable :: bond,conn,ang
real*8,dimension(:),allocatable :: rcov
real*8,dimension(:,:),allocatable :: dist
character(len=2),dimension(:),allocatable :: sym

        !allocate(cartcoord(3*natom))
        !call getrec(50, "JOBARC", "COORD", 1, cartcoord)
        
        open(unit=10,file="ZMAT",status="old",err=12)
        open(unit=13)

        allocate(Hinit((3*nat-6)))
        allocate(bond(nat),conn(nat),ang(nat))
        allocate(rcov(nat),sym(nat))
        allocate(dist(nat,nat))

        Hinit=0.0
        bond=0.0
        conn=0.0
        ang=0.0
        rcov=0.0
        dist=0.0
        
        ij=1
        
        read(10,*)

        do i=1,nat
                read(unit=10,fmt="(a2,1x)",iostat=ios,advance="no",eor=11) sym(i)
               
                read(unit=10,fmt="(i3,1x,a4,1x)",iostat=ios,advance="no",eor=11) bond(i), bondl
                
!               bondl=dist(i,bond(i))
!               dist(i,bond(i))=dist(bond(i),i)
!               conn(i)=conn(i)+1
!               conn(bond(i))=conn(bond(i))+1
!              
!               Hinit(ij)=0.3601*exp(-1.944*(bondl-rcov(i)-rcov(bond(i))))
!               ij=ij+1
               
                read(unit=10,fmt="(i3,1x,a4,1x,i3,1x,a4)",iostat=ios,err=14) ang(i), strangle, l, strdihed
                
11              do while (sym(i).ne.symbol)
                        read(13,iostat=ios) symbol, rcov(i)

                        if (ios.gt.0.OR.ios.lt.0) then
                                print *,"ERROR: Problem with rcov assignment!"
                                GOTO 14
                        else
                                continue
                        end if
                end do

                bondl=dist(i,bond(i))
                dist(i,bond(i))=dist(bond(i),i)
                conn(i)=conn(i)+1
                conn(bond(i))=conn(bond(i))+1
 
                Hinit(ij)=0.3601*exp(-1.944*(bondl-rcov(i)-rcov(bond(i))))
                ij=ij+1

        end do
        
        if (nat.gt.2) then
                do i=3,nat
                        Hinit(ij)=0.089+0.11*exp(-0.44*(dist(ang(i),bond(i))+dist(ang(i),i)-2*rcov(ang(i))-rcov(bond(i))-rcov(i)))
                        ij=ij+1
                        
                        if (i.gt.3) then
                                Hinit(ij)=0.0015+14.0*((conn(ang(i))+conn(bond(i))-2)**0.57)*exp(-2.85*(dist(ang(i),bond(i)) &
-rcov(bond(i))-rcov(ang(i))))/((dist(ang(i),bond(i))*(rcov(ang(i))+rcov(bond(i))))**4.00)
                                ij=ij+1
                        else
                                        continue
                        end if
                end do
        else
                        continue
        end if

!       deallocate(Hinit,bond,conn,ang,rcov,sym,dist)

        close (10)
        close (13)

14      return

        deallocate(Hinit,bond,conn,ang,rcov,sym,dist)
        
12      print *,'ERROR: File not found!'

end subroutine
