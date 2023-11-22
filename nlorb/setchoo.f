      subroutine setchoo(nchoose,choose,mxchoose,bondsize)
      implicit none

      integer iii

      integer mxchoose, bondsize

      integer nchoose(bondsize),
     & choose(bondsize,mxchoose,bondsize)

CCCC disty bondsize 2 CCCC
CC      nchoose(1) = 0 
CC      nchoose(2) = 6

CC      choose(1,1,2) = 8
CC      choose(2,1,2) = 12
CC      choose(1,2,2) = 18
CC      choose(2,2,2) = 21
CC      choose(1,3,2) = 13
CC      choose(2,3,2) = 17
CC      choose(1,4,2) = 16
CC      choose(2,4,2) = 24
CC      choose(1,5,2) = 28
CC      choose(2,5,2) = 31
CC      choose(1,6,2) = 25
CC      choose(2,6,2) = 27

CCCC disty bondsize 6 CCCC
CC      nchoose(1) = 0 
CC      nchoose(2) = 0
CC      nchoose(3) = 9999
CC      nchoose(4) = 9999
CC      nchoose(5) = 9999
CC      nchoose(6) = -2

CC      choose(1,1,6) = 8
CC      choose(2,1,6) = 12
CC      choose(3,1,6) = 13
CC      choose(4,1,6) = 17
CC      choose(5,1,6) = 18
CC      choose(6,1,6) = 21

CC      choose(1,2,6) = 16 
CC      choose(2,2,6) = 24
CC      choose(3,2,6) = 25
CC      choose(4,2,6) = 27
CC      choose(5,2,6) = 28
CC      choose(6,2,6) = 31
      
CCCC 1-butylbenz bondsize 2 CCCC
CC      nchoose(1) = 0 
CC      nchoose(2) = 3
      
CC      choose(1,1,2) = 14
CC      choose(2,1,2) = 15
CC      choose(1,2,2) = 16
CC      choose(2,2,2) = 17
CC      choose(1,3,2) = 18
CC      choose(2,3,2) = 21

CCCC 1-butylbenz bondsize 6 CCCC
      if (.false.) then
      nchoose(1) = 0 
      nchoose(2) = 0
      nchoose(3) = 9999
      nchoose(4) = 9999
      nchoose(5) = 9999
      nchoose(6) = -1
      
      choose(1,1,6) = 14
      choose(2,1,6) = 15
      choose(3,1,6) = 16
      choose(4,1,6) = 17
      choose(5,1,6) = 18
      choose(6,1,6) = 21
      end if

CCCC 2-butylbenz bondsize 2 CCCC
CC      nchoose(1) = 0 
CC      nchoose(2) = 3
      
CC      choose(1,1,2) = 8
CC      choose(2,1,2) = 12
CC      choose(1,2,2) = 13
CC      choose(2,2,2) = 17
CC      choose(1,3,2) = 18
CC      choose(2,3,2) = 21

CCCC 2-butylbenz bondsize 6 CCCC
      if (.false.) then
      nchoose(1) = 0 
      nchoose(2) = 0
      nchoose(3) = 9999
      nchoose(4) = 9999
      nchoose(5) = 9999
      nchoose(6) = -1
      
      choose(1,1,6) = 8
      choose(2,1,6) = 12
      choose(3,1,6) = 13
      choose(4,1,6) = 17
      choose(5,1,6) = 18
      choose(6,1,6) = 21
      end if

CCCC disty bondsize 2 part 2 CCCC
      if (.false.) then
      nchoose(1) = 0 
      nchoose(2) = 6

      choose(1,1,2) = 8
      choose(2,1,2) = 13
      choose(1,2,2) = 12
      choose(2,2,2) = 18
      choose(1,3,2) = 17
      choose(2,3,2) = 21
      choose(1,4,2) = 16
      choose(2,4,2) = 24
      choose(1,5,2) = 28
      choose(2,5,2) = 31
      choose(1,6,2) = 25
      choose(2,6,2) = 27
      end if

CCC water trimer anion

      if (.false.) then
         nchoose(1) = 0
         nchoose(2) = 0
         nchoose(3) = 9999
         nchoose(4) = 9999
         nchoose(5) = 9999 
         nchoose(6) = 0
      end if

CCC water trimer linear anion

      if (.false.) then
         nchoose(1) = 0
         nchoose(2) = 0
         nchoose(3) = 9999
         nchoose(4) = 9999
         nchoose(5) = 0
      end if

      if (.false.) then
         nchoose(1) = 0
         nchoose(2) = 0
         nchoose(3) = 9999
         nchoose(4) = 9999
         nchoose(5) = 9999
         nchoose(6) = 9999
         nchoose(7) = 9999
         nchoose(8) = 9999
         nchoose(9) = 9999
         nchoose(10) = 0
      end if

CCC polyacetylene 

      if (.false.) then
         choose(1,1,8) = 14
         choose(2,1,8) = 16
         choose(3,1,8) = 18
         choose(4,1,8) = 20
         choose(5,1,8) = 22
         choose(6,1,8) = 24
         choose(7,1,8) = 26
         choose(8,1,8) = 28
C         choose(5,1,14) = 14 
C         choose(6,1,14) = 16 
C         choose(7,1,14) = 18 
C         choose(8,1,14) = 20 
C         choose(9,1,14) = 22 
C         choose(10,1,14) = 24 
C         choose(11,1,14) = 26 
C         choose(12,1,14) = 28 
C         choose(13,1,14) = 30 
C         choose(14,1,14) = 32 

         nchoose(1) = 0
         nchoose(2) = 0
         nchoose(3) = 9999
         nchoose(4) = 9999
         nchoose(5) = 9999
         nchoose(6) = 9999
         nchoose(7) = 9999
         nchoose(8) = -1
         go to 90
         nchoose(5) = 9999
         nchoose(6) = 9999
         nchoose(7) = 9999
         nchoose(8) = 9999
         nchoose(9) = 9999
         nchoose(10) = 9999
         nchoose(11) = 9999
         nchoose(12) = 9999
         nchoose(13) = 9999
         nchoose(14) = 9999
         nchoose(15) = 9999
         nchoose(16) = 9999
         nchoose(17) = 9999
         nchoose(18) = 0
 90      continue
      end if

CCC water dimer

      if (.true.) then
         nchoose(1) = 0
         nchoose(2) = 0
         nchoose(3) = -1

         choose(1,1,3) = 1
         choose(2,1,3) = 2
         choose(3,1,3) = 4
      end if

      return
      end

