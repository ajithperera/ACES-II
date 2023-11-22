
#ifndef _INT_COM_
#define _INT_COM_

c icntr     : The integration center
c idns      : a flag =0 for SCF orbitals and =1 for natural orbitals

      integer icntr,idns

      common /int/  icntr,idns
      save /int/

c array pointers

c zpcoeff(2): alpha/beta MO to primitive function transformation matrix
c zxocc     : alpha/beta orbital occupation

      integer
     &    zpcoeff(2),zxocc
      common /molecp/
     &    zpcoeff,zxocc
      save /molecp/

#endif /* _INT_COM_ */

