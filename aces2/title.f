
c This routine prints the ACES banner.

      SUBROUTINE TITLE
      WRITE(*,9000)
 9000 FORMAT(
     1/T15,'****************************************************'
     2/T15,'* ACES : Advanced Concepts in Electronic Structure *'
c Measure  '* 1 2 3 4 5 6 7 8 9 0 1 2  2 1 0 9 8 7 6 5 4 3 2 1 *'
     3/T15,'*                    Ver. 2.13.0                    *'
c     4/T15,'*               Release Candidate: 3               *'
c     5/T15,'*               exported DD MMM YYYY               *'
     6/T15,'****************************************************'
c This file was checked out with the tag $Name:  $
     &      )
      WRITE (*, 9100)
 9100 FORMAT(/T30,'Quantum Theory Project'
     1       /T30,'University of Florida'
     2       /T30,'Gainesville, FL  32611')
      RETURN
      END

