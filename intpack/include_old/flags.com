
#ifndef _FLAGS_COM_
#define _FLAGS_COM_

c This common block contains values for all of the Aces3 keywords as
c specified by default or in the ZMAT file.  For historical reasons,
c the flags are broken up into two blocks (originally only 100 were
c used and later an additional 500 were added).  It is hoped that a
c more flexible way of doing this will be implemented.

c nflags,nflags2  : the number of flags in each block
c iflags,iflags2  : the values for all flags
c iprint          : a special flag containing print information
c iuhf            : 1 if a UHF calculation is done
c needbas         : set to 1 if basis set info is needed
c debug           : a logical flag which may be useful to some
c multipoint      : A flag that says this more than one single point calculation.
c getall          : logical flag that is used in mrcc
c get_nonsense    : if .not. get_nonsense aces3 stops calculation if it thinks it is wrong.
c development_version: indiacates that this is a development version rather than production.

      integer nflags,nflags2
      parameter (nflags=100)
      parameter (nflags2=500)
      
      integer iflags(nflags),iflags2(nflags2),iprint,iuhf,needbas
      logical debug, multipoint, get_nonsense, development_version

      common /flags/  iflags,iflags2,iprint,iuhf,needbas
      common /flagsl/ debug, multipoint, get_nonsense,
     &    development_version

#endif /* _FLAGS_COM_ */

