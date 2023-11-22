
      subroutine FormPi(icore,maxcor,uhf,listl1,listl1off,listl2,
     &                  listl2rs,listr1,listr1off,listr2,listr2rs,
     &                  listz1,listz1off,listz2)
      implicit none
C     Input Variables
      integer uhf,maxcor,listl1,listl1off,listl2,listl2rs,listr1,
     &        listr1off,listr2,listr2rs,listz1,listz1off,listz2

C     Pre-allocated Local Variables
      integer icore(maxcor)

C     Local Variables
      integer i000
C     - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C     Memory Pointers
      i000=1
      call getdiag2(icore(i000),maxcor,uhf)
      call maklst2(listz1,listz1off,listz2,uhf)
      call calcxi(icore(i000),maxcor,uhf,listl1,listl1off,listl2,
     &            listl2rs,listr1,listr1off,listr2,listr2rs,listz1,
     &            listz1off,listz2,454,125)

      return
      end
