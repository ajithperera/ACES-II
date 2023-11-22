      subroutine printt1(string,t1,n)
      implicit double precision(a-h,o-z)
      dimension t1(n)
      character*7 string
c
      write(*,1000)  string
1000  format(a7)
c      write(*,2000) (t1(i),i=1,n)
      do i=1,n
        write(*,2000) (t1(i))
      enddo
2000  format(8f10.7)
      return
      end
