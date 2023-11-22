module indices
integer,dimension(:),allocatable::ifirst,ilast,ifirst2,ilast2
integer,dimension(:),allocatable::offset1,offset2
integer::ndim1,ndim2
double precision::scalar(10)=(/1.0d0,2.0d0,2.0d0,2.0d0,1.0d0,2.0d0,2.0d0,1.0d0,2.0d0,1.0d0/)

! these arrays are for indexing in the f2.f90 subroutine
!integer:: za(16)=(/0,1,2,3,1,4,5,6,2,5,7,8,3,6,8,9/)
!integer:: zb(16)=(/0,1,2,3,1,4,5,6,2,5,7,8,3,6,8,9/)
!integer:: zc(4)=(/0,4,8,12/)
!integer:: zd(4)=(/0,4,8,12/)

integer,dimension(9,9)::map=reshape( source =(/ &
1, 2, 3, 4, 11, 22, 37, 16, 29, &
2, 5, 6, 7, 12, 23, 38, 17, 30, &
3, 6, 8, 9, 13, 24, 39, 18, 31, &
4, 7, 9, 10, 14, 25, 40, 19, 32, &
11, 12, 13, 14, 15, 26, 41, 20, 33, &
22, 23, 24, 25, 26, 28, 43, 27, 35, &
37, 38, 39, 40, 41, 43, 45, 42, 44, &
16, 17, 18, 19, 20, 27, 42, 21, 34, &
29, 30, 31, 32, 33, 35, 44, 34, 36 /), shape=(/9,9/) )






end module indices
