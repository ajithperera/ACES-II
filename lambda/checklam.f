       subroutine checklam(a)
C
       implicit double precision (a-h,o-z)
       integer dirprd
       dimension a(1)
       common/syminf/nirrep,nstart,irrepa(255,2),dirprd(8,8)
       COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),NTOT(18)

C
       sum=0.0
       sum1=0.0
       isize=isymsz(isytyp(1,63),isytyp(2,63))
       call getall(a,isize,1,63)
       do 1 i=1,isize
       sum=sum+a(i)*a(i)
       sum1=sum1+abs(a(i))
1      continue
       write(*,*) 'checksum',sum,sum1
       return
       end
