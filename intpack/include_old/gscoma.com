
#ifndef _GSCOMA_COM_
#define _GSCOMA_COM_

c     gssopt - .true. if the guess file exists. if it is .true. then the
c              initial guess options are read from the guess file.
c     gssalw - if it is .true. then initial guess parameters are always
c              read from the guess file. if it is .false. and gssopt is
c              .true. then initial guess parameters are read from guess
c              for the first scf calculation only.
c     gssalt - if it is .true. then the initial guess orbitals will be
c              swapped according to parameters in the swap array.
c     gsslok - if it is .true. then attempts will be made in scfit to
c              maintain the character of the starting orbitals by monitoring
c              c(old)t * s c(new). which symmetry and spin blocks to try
c              to lock are specified by the lock array.
c     gssred - if it is .true. then the initial orbitals will be read from
c              the formatted file oldmos.
c     gsswrt - if it is .true. then the converged orbitals will be written
c              to the formatted file newmos (in exactly the same format
c              as oldmos).
c     gssufr - if it is .true. then a uhf guess is generated from a single
c              set of orbitals (the first set on oldmos). these may originate
c              from a multitude of sources, including a closed-shell rhf run
c              and the alpha orbitals of a previous  uhf run.

      logical         gssopt,gssalw,gssalt,gsslok,gssred,gsswrt,gssufr
      common /gscoma/ gssopt,gssalw,gssalt,gsslok,gssred,gsswrt,gssufr

#endif /* _GSCOMA_COM_ */

