      subroutine oepreadaux (
     &     nbas, nbasx, naux, ntri,
     &     zmtcmp, scrzxz, scraxz, auxao)
c
c This routine reads the KSAOINT file containing the auxiliary basis set
c ZMAT AO matrix elements and transform them to VMOL AO representation,
c and packs them into triangular matrix.
c
c This routine will be redundant when INTPACK is rewritten.
c
c In: nbas, nbasx, naux, ntri, zmtcmp
c Scr: srczxz, scraxz, scraxa,
c Out: auxao
c
c Igor Schweigert, Aug 2003
c $Id: oepreadaux.FPP,v 1.1.1.1 2008/05/29 19:35:39 taube Exp $
c
      implicit none
c
c     .. the input parameters                                         ..
c
      integer
     &     nbas, nbasx, naux, ntri
c
      double precision
     &     zmtcmp (nbas, nbasx), scrzxz (nbasx, nbasx),
     &     scraxz (nbas, nbasx), auxao (ntri, naux)
c
c     .. local variables                                              ..
c
      integer
     &     i, j, x
c
      logical exists
c
c     .. messages to the user                                         ..
c
 1005 format (
     &     t3, '@OEPREADAUX-F: File KSAOINT does not exist. ',
     &     'The INTPACK module should be run before performing ',
     &     'the OEP calculation.')
 1010 format (
     &     t3, '@OEPREADAUX-F: NBAS=', i4, ' is greater than NBASX=',
     &     i4, '.')
c
c     .. check if the KSAOINT file exists                             ..
c
      inquire (file = 'KSAOINT', exist = exists)
      if (.not.exists) then
         write (6, 1005)
         call errex
      endif
c
c     NBASX ge NBAS is assumed. Let's check it. (Stupid as it is , but I
c     want to make sure.)
c     
      if (nbas.gt.nbasx) then
         write (6, 1010)
         call errex
      endif
c
c     .. read the auxiliary ZMAT AO ME and transform them to VMOL AO  ..
c
      open (
     &     unit = 44, file = 'KSAOINT', status = 'old',
     &     form = 'unformatted')
c
      do x=1, naux
c
         do i=1, nbasx
            do j=1, nbasx
               read (44) scrzxz (i, j)
            enddo
         enddo
c
         call xgemm(
     &        'n', 'n', nbas, nbasx, nbasx,
     &        1.d0, zmtcmp, nbas,
     &        scrzxz, nbasx,
     &        0.d0, scraxz, nbas)
c
         call xgemm(
     &        'n', 't', nbas, nbas, nbasx,
     &        1.d0, scraxz, nbas,
     &        zmtcmp, nbas,
     &        0.d0, scrzxz, nbas)
c     
         call triang (scrzxz, auxao (1,x), nbas, 'y')
c
      enddo
c
c     .. done ..
c
      return
      end

         



      
      
