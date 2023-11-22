      subroutine addfdiag( eig, f, nirrep, ivec)
C
C This routine adds the occupied (virtual) values found in EIG to the
c  diagonal elements of the occupied(virtual) one-electron quantity F.
C Can be used to add the diagonal elements back to the the fock matrix or
c  an F intermediate.
C  
C   PARAM        USE                                                    CHANGED
C   -----        ------------------------------------------------       -------
C   EIG ........ Diagonal elements to be added to F                         NO
C   F .......... A symmetry-blocked one-electron list.                     YES
C   nirrep ..... Number of irreps.                                          NO
C   ivec ....... should be pop(1,ispin) or vrt(1,ispin), the symmetry       NO
C                 distribution of the occupied (virtual) orbitals
C 
C Written by Renee P Mattie, Feb 1991
CEND
C
C      include 'imp.inc'
      integer nirrep, ivec(nirrep)
      double precision eig(*), f(*)
C
      integer ieig, ifbar, irrep, na, a, inc
C
C      return
C      write (*,*) 'entering addfdiag'
      ieig = 0
      ifbar = 0
      do 100, irrep = 1, nirrep
         na = ivec(irrep)
         do 90, a = 1, na
            inc = a + (a-1)*na
            f(ifbar + inc) = f(ifbar + inc) + eig(ieig+a)
 90      continue
         ifbar = ifbar + na**2
         ieig = ieig + na
 100  continue
      return
      end
