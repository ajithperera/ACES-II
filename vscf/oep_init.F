      subroutine oep_init(petite_list, iuhf, iwork, imemleft,loslater)

      implicit double precision (a-h, o-z)

#include <oep.com>
#include <flags.com>
#include <flags2.com>
#include <machsp.com>
#include <symm2.com>

      Dimension iwork(imemleft)
      character*3 itertype
      logical petite_list, zmat_exsist, oepint_exsist
      logical loslater

      OEP1    = .FALSE.
      OEP2    = .FALSE.
      OEP1ON  = .FALSE.
      OEP2ON  = .FALSE.
      OEP2ITR = .FALSE.
      R12ON   = .FALSE.


      if (iflags2 (153).eq.3) then
         OEP1 = .TRUE.
         write (6, 5910)
      endif
      if (iflags2 (153).eq.4) then
         OEP2 = .TRUE.
         write (6, 5920)
      endif
      if (OEP1 .or. OEP2) then
c
         if (iuhf.ne.0) then
            write (6, 5901)
            call errex
         endif
         if (petite_list) then
            write (6, 5903)
            call errex
         endif

         i020 = i010 + iintfp * itriln(nirrep+1)
         i030 = i020 + 10000
         i040 = i030 + 500*20
         i050 = i040 + 500*20
         i060 = i050 + 500*1
         if ((i060-i000).gt.maxmem)
     &      call nomem (
     &        'Parse the OEP ASV list', '{OEPINI} <-- VSCF',
     &        i060-i000, maxmem)

         inquire (file = 'ZMAT', exist = zmat_exsist)
         if (.not. zmat_exsist) then
            write (6, 5930)
            call errex
         endif
         luoep = 100
         open (unit=luoep, file='ZMAT', form='formatted', status='old')
         call oepini (
     &        luoep, 10000, 500, 20,
     &        iwork (i020), iwork (i030), iwork (i040), iwork (i050),
     &        loslater)
         close (unit=luoep,status='keep')

         inquire (file = 'OEPINT', exist = oepint_exsist)
         if (.not. oepint_exsist) then
            write (6, 5933)
            call errex
         endif
         open (unit=luoep,file='OEPINT',form='unformatted',status='old')
         call locate (luoep, 'NAUX    ')
         read (unit=luoep) naux
         close (unit=luoep)
         write (6, 5936) naux
c
         if (OEP2) then
            OEPTOL = tol
            tol = 1.d-99
         endif
         itertype (1:3) = '   '
c
 5901    format(
     &        t3,'@VSCF-F: The OEP method for an unrestriced',/,
     &        t3,'reference is not implemented yet.')
 5903    format(
     &        t3,'@VSCF-F: The OEP method for any integral package',/,
     &        t3,'other than VMOL is not implemented yet.')
 5910    format(
     &        t3,'@VSCF-I: The first-order OEP method will be used.',/)
 5915    format(
     &        t3,'@VSCF-I: The exchange OEP-FIA method will be used.',/,
     &        t3,'The local exchange potential will be found from ',/,
     &        t3,'the condition that its occupied-virtual matrix ',/,
     &        t3,'elements are equal to those of non-local exchange',/,
     &        t3,'potential.',/)
5917    format(
     &        t3,'@VSCF-I: The exchange OEP-WFIA method will be used.',
     &        /,t3,'The local exchange potential will be found from ',/,
     &        t3,'the condition that its occupied-virtual matrix ',/,
     &        t3,'elements weighted with the energy denominators are',/,
     &        t3,'equal to those of non-local exchange potential.',/)
 5920    format(
     &        t3,'@VSCF-I: The second-order OEP method will be used.',/)
 5930    format (
     &        t3,'@VSCF-F: The ZMAT file is not found.',/)
 5933    format (
     &        t3,'@VSCF-F: The OEPINT file is not found.',/)
 5936    format (
     &        t3,'@VSCF-I: The auxiliary basis set consists of ',i4,/
     &        t3,'functions.',/)
 5950    format (
     &        t3,i6,1x,a3,1x,f20.10,10x,d20.10)
 5960    format(
     &        t3,1x,a10,f20.10,6x,a3,1x,d20.10/)
 5964    format(
     &        t3,'@VSCF-I: The HOMO condition on the OEP exchange',/,
     &        t3,'potential is violated by ',f10.5,' a.u. (',f10.5,
     &        ' eV).', /)
 5965    format(
     &        t3,'@VSCF-I: The OEP eigenvalues will be shifted by',/,
     &        t3,f10.5,' a.u. to satisfy the HOMO condition on',/,
     &        t3,'the OEP exchange potential.',/)
c
      endif
C
      return
      end

