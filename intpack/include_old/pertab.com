
#ifndef _PERTAB_COM_
#define _PERTAB_COM_

c This common block contains information from the periodic table.

c maxpertab  : The current maximum atomic number in the table.
c atomicsym  : The 2 character atomic symbol for each element.
c atomicmass : The atomic mass of each element.

      M_REAL
     &    atomicmass(maxpertab)

      character*2 atomicsym(maxpertab)

      common /pertab/ atomicmass
      common /pertabc/atomicsym
      save /pertab/
      save /pertabc/

#endif /* _PERTAB_COM_ */

