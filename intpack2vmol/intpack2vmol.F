      program intpack2vmol
c
c This module transforms the INTPACK's KSAOINT file into the VMOL
c standards. The transformed file can be read similarly to VMOL's
c two-electron integral files (IIII and the likes.)
c
c Igor Schweigert, Jan 2004
c $Id: intpack2vmol.FPP,v 1.2 2008/05/29 20:04:35 taube Exp $
c 
      implicit none
c
c     Common blocks
c
#include <icore.com>     
#include <flags.com>
#include <machsp.com>
c 
      logical exists, ok
c     
      character*8 label
c     
      integer
     &     i000, iuhf, maxmem, coretop,
     &     i010, i020, i030, itop,
     &     x, nbas, nbasx, naux, 
     &     luint, luaux, lnbuf, iauxaa, ntri
c     
      parameter (luint = 44)
      parameter (luaux = 45)
      parameter (lnbuf = 600)
      integer nut, ibuf (600)
      double precision buf (600)
c
c     Messages to the user
c     
 5900 format (
     &     /,
     &     'module INTPACK2VMOL',/,
     &     'Post-INTPACK integral file processing. ',/,
     &     'Development version. ',
     &     /)
 5905 format (
     &     t3,'@INTPACK2VMOL-F: The transformation of integrals with ',/
     &     t3,'a symmetry other than C1 has not been implemented yet',/)
 5910 format (
     &     t3,'@INTPACK2VMOL-F: Cannot find the KSAOINT file.',/)
 5920 format (
     &     t3,'@INTPACK2VMOL-I: The AO matrix elements of th auxiliary'/
     &     t3,'basis will be read from file KSAOINT and stored',/
     &     t3,'as a record in the OEPINT file.',/)
 5923 format (
     &     t3,'@INTPACK2VMOL-I: The AO matrix elements of the Slater',/
     &     t3,'potential will be read from file SLAT51AO and stored',/
     &     t3,'as a record in the OEPINT file.',/)
c
c     Print out the banner
c     
      write (6,5900)
c
c     Initialize the ACES environment: initialize the run-time
c     environment, set up the JOBARC subsystem, load ASV's, and allocate
c     memory according to IFLAGS(36).
c
c     Note that ACES_INIT also setups the ACES IO subsystem, which we do
c     not need, but I believe ACES_INIT will stay backward compatable
c     longer than all the lower level initializing routines.
c     
      call aces_init (icore, i000, maxmem, iuhf, .true.)
      coretop = i000
c
c     Check if the KSAOINT file exists.
c
      inquire (file = 'KSAOINT', exist = exists)
      if (.not.exists) then
         write (6, 5910)
         call errex
      endif
      write (6, 5920)
c
c     Get the AO, SO and auxiliary basis set sizes: <> Read AO and SO
c     basis set sizes from JOBARC. <> Read the auxiliary basis set size
c     from the INTPACK's file and put it in JOBARC.
c     
      call getrec(20, 'JOBARC', 'NBASTOT ', 1, nbas)
      call getrec(20, 'JOBARC', 'NAOBASFN', 1, nbasx)
      ntri = nbas * (nbas+1) / 2
c     
      open (
     &     unit = luint, file = 'NUMAUX', status = 'old',
     &     form = 'formatted')
      read (luint, '(i5)') naux
      close (luint, status = 'keep')
      call putrec (20, 'JOBARC', 'OEPNAUX ', 1, naux)
c
c     Transform the auxiliary integrals into the VMOL representations:
c     <> Allocate memory for the output and scratch arrays. <> Check if
c     there is enough memory available. <> Read the VMOL2ZMAT
c     transformation matrix. <> Read and transform the auxiliary
c     integrals.
c     
      iauxaa = coretop
      coretop = iauxaa + iintfp * naux * ntri
      i010 = coretop
      i020 = i010 + iintfp * nbasx * nbas
      i030 = i020 + iintfp * nbasx * nbasx
      itop = i030 + iintfp * nbasx * nbas
c     
      if (itop - i000 .gt. maxmem) call insmem (
     &     'Read the SCF eigenvalues and eigenvectors.',
     &     itop - i000, maxmem)
c     
      call getrec (
     &     1, 'JOBARC', 'ZMAT2CMP', nbas*nbasx*iintfp, icore (i010))
c     
      call oepreadaux (
     &     nbas, nbasx, naux, ntri,
     &     icore (i010), icore (i020),
     &     icore (i030), icore (iauxaa))
c
c     Write the transformed auxiliary integrals in a VMOL-like format:
c     <> Allocate memory for file writing buffers. <> Open the file
c     where the integrals will be stored <> Write the number of
c     auxiliary functions into the file.  <> Write the integrals into
c     the file. <> Check if the file with Slater AO MEs exists, if it
c     does, read it and store the integrals into OEPINT. <> Close the
c     file.
c     
      i010 = coretop
      i020 = i010 + iintfp * lnbuf
      itop = i020 + lnbuf
c     
      open (
     &     unit = luint, file = 'OEPINT',
     &     status = 'new', form = 'unformatted')
c     
      write (luint) '********'
      write (luint) 'NAUX    '
      write (luint) naux, '        '
c
      call oepsaveints (
     &     luint, 'AUX3CNTR', lnbuf, icore (i020), icore (i010),
     &     1.00d-14, naux * ntri, icore (iauxaa))
c
      inquire (file = 'SLAT51AO', exist = exists)
      if (exists) then
         write (6, 5923)
         open (
     &        unit = luaux, file = 'SLAT51AO',
     &        form = 'unformatted', status = 'old')
         call readarray (luaux, ntri, icore (iauxaa))
         close (unit=luaux, status='keep')
         call oepsaveints (
     &        luint, 'SLAT51AO', lnbuf, icore (i020), icore (i010),
     &        1.00d-14, ntri, icore (iauxaa))
      endif
c      
      close (unit = luint, status = 'keep')
c
c     Print the original integrals to test them.
c     
c$$$      write (6, *) 'INTPACK2VMOL: AUX3CNTR'
c$$$      call output (
c$$$     &     icore (iauxaa), 1, naux * ntri, 1, 1, naux * ntri, 1, 1)
c$$$      write (6, *) 'INTPACK2VMOL: AUX3CNTR END'
c
c     Finilize ACES environment and quit.
c     
      call aces_fin
      end


      subroutine readarray (unit, n, a)
      implicit none
      integer unit, n, i
      double precision a (n)
      read (unit) (a (i), i=1, n)
      return
      end
