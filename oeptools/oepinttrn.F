      subroutine oepinttrn (
     &     nbas, nocc, luint, lnbuf, ibuf,
     &     buf, evec, paah, aaah, ppah, phah,
     &     ppph, pphh, phph, phhh)
c
c This routine transforms the two-electron integrals (read in from the
c VMOL file) from the AO to MO basis using a given set of SCF
c eigenvectors. See also documentation file OEPINTTRN.TEX for
c instructions on how the scratch and output arrays should be allocated
c and for the details of the implementation.
c
c Note that for simplicity the PHPH array is not packed!
c
c In: nbas, nocc, luint, lnbuf, ibuf, buf, evec
c Scr: paah, aaah, ppah, phah
c Out: ppph, pphh, phph, phhh
c
c Igor Schweigert, Dec 2003
c $Id: oepinttrn.FPP,v 1.1.1.1 2008/05/29 19:35:40 taube Exp $
c     
      implicit none
c
c     Arguments
c     
      integer
     &     nbas, nocc, luint, lnbuf, ibuf (lnbuf)
c 
      double precision
     &     buf (lnbuf), evec (nbas, nbas),
     &     paah (*), aaah (*), ppah (*), phah (*),
     &     ppph (*), pphh (*), phph (*), phhh (*)
c
c     Local variables
c
      logical
     &     i1i2, i3i4, i1i2i3i4
c           
      integer
     &     ind, nut, iand, ishft,
     &     n1, n2, i1, i2, i3, i4, i5
c     
      character*80
     &     fname
c
      double precision x
c
c     Define the two-electron integral offset operators and the
c     associated variables.
c
#define DEFINE_OPERATORS     
#include "oepintind.com"
#define DEFINE_VARIABLES
#define N_OCC nocc
#define N_BAS nbas
#include "oepintind.com"
c     
c
c     Read the 2e integrals in AO basis and transform one AO index into
c     a hole index: <> Initialize the array being calculated. <> Get the
c     name of the VMOL file from JOBARC, open the file, and shift the
c     pointer to the begining of the 2e integral record. <> Read the
c     integrals and transform one AO index to a hole index: <<>> Read
c     the value and indices of the current integral. <<>> Define the
c     symmetry of the integral. <<>> Loop over the hole index, check the
c     symmetry of the integral and update the corresponding integrals in
c     the AAAH array. <> Close the VMOL file.
c
      call zero (aaah, n_aaah)
c     
      call gfname('IIII    ',fname,i)
      open(unit=luint,file=fname(1:i),form='UNFORMATTED',
     &     access='SEQUENTIAL')
      call locate(luint,'TWOELSUP')
c     
      nut = lnbuf
      do while (nut.eq.lnbuf)
         read(luint) buf, ibuf, nut
         do ind=1,nut
c
            x=buf(ind)
            i1=iupki(ibuf(ind))
            i2=iupkj(ibuf(ind))
            i3=iupkk(ibuf(ind))
            i4=iupkl(ibuf(ind))
c
            i1i2 = i1.ne.i2
            i3i4 = i3.ne.i4
            i1i2i3i4 = i_aa (i1,i2) .ne. i_aa (i3,i4)
c     
            do i5=1, nocc
               aaah (i_aaah (i1, i2, i3, i5)) =
     &              aaah (i_aaah (i1, i2, i3, i5)) + evec (i4,i5)*x
               if (i3i4)
     &              aaah (i_aaah (i1, i2, i4, i5)) =
     &              aaah (i_aaah (i1, i2, i4, i5)) + evec (i3,i5)*x
               if (i1i2i3i4)
     &              aaah (i_aaah (i3, i4, i1, i5)) =
     &              aaah (i_aaah (i3, i4, i1, i5)) + evec (i2,i5)*x
               if (i1i2i3i4.and.i1i2)
     &              aaah (i_aaah (i3, i4, i2, i5)) =
     &              aaah (i_aaah (i3, i4, i2, i5)) + evec (i1,i5)*x
            enddo
         enddo
      enddo
c     
      close(unit=luint,status='KEEP')
c
c     AAAH->PAAH. Convert the first AO index to a particle index: <>
c     Initialize the target array. <> Loop over the new particle index
c     and remaining indices, and update the target array.
c     
      call zero (paah, n_paah)
c     
      do n1 = 1, n_aaah
         i4 = i4_aaah (n1)
         i3 = i3_aaah (n1, i4)
         i2 = i2_aaah (n1, i3, i4)
         i1 = i1_aaah (n1, i2, i3, i4)
         i1i2 = i1.ne.i2
         do i5 = 1, nbas-nocc
            n2 = i_paah (i5, i2, i3, i4)
            paah (n2) = paah (n2) + aaah (n1) * evec (i1, i5+nocc)
            if (i1.ne.i2) then
               n2 = i_paah (i5, i1, i3, i4)
               paah (n2) = paah (n2) + aaah (n1) * evec (i2, i5+nocc)
            endif
         enddo
      enddo
c
c     PAAH->PPAH,PHAH. Convert the second AO index to a particle and a
c     hole indices: <> Initialize the target arrays. <> Loop over the
c     source array and construct the target arrays.
c     
      call zero (ppah, n_ppah)
      call zero (phah, n_phah)
c     
      do n1 = 1, n_paah
         i4 = i4_paah (n1)
         i3 = i3_paah (n1, i4)
         i2 = i2_paah (n1, i3, i4)
         i1 = i1_paah (n1, i2, i3, i4)
         do i5 = i1, nbas-nocc
            n2 = i_ppah (i1, i5, i3, i4)
            ppah (n2) = ppah (n2) + paah (n1) * evec (i2, i5+nocc)
         enddo
         do i5 = 1, nocc
            n2 = i_phah (i1, i5, i3, i4)
            phah (n2) = phah (n2) + paah (n1) * evec (i2, i5)
         enddo
      enddo
c
c     PPAH->PPPH, PPHH. Convert the third AO index to a particle and a
c     hole indices: <> Initialize the target arrays. <> Loop over the
c     source array and construct the target arrays.
c     
      call zero (ppph, n_ppph)
      call zero (pphh, n_pphh)
c     
      do n1 = 1, n_ppah
         i4 = i4_ppah (n1)
         i3 = i3_ppah (n1, i4)
         i2 = i2_ppah (n1, i3, i4)
         i1 = i1_ppah (n1, i2, i3, i4)
         do i5 = 1, nbas-nocc
            n2 = i_ppph (i1, i2, i5, i4)
            ppph (n2) = ppph (n2) + ppah (n1) * evec (i3, i5+nocc)
         enddo
         do i5 = 1, i4
            n2 = i_pphh (i1, i2, i5, i4)
            pphh (n2) = pphh (n2) + ppah (n1) * evec (i3, i5)
         enddo
      enddo
c
c     PHAH->PHPH, PHHH. Convert the third AO index.: <> Initialize the
c     target arrays. <> Loop over the new indices and remaining indices,
c     and update the target arrays.
c     
      call zero (phph, n_phph)
      call zero (phhh, n_phhh)
c     
      do n1 = 1, n_phah
         i4 = i4_phah (n1)
         i3 = i3_phah (n1, i4)
         i2 = i2_phah (n1, i3, i4)
         i1 = i1_phah (n1, i2, i3, i4)
         do i5 = 1, nbas-nocc
            if (i_ph (i1, i2) .le. i_ph (i5, i4)) then
               n2 = i_phph (i1, i2, i5, i4)
               phph (n2) = phph (n2) + phah (n1) * evec (i3, i5+nocc)
            endif
         enddo
         do i5 = 1, i4
            n2 = i_phhh (i1, i2, i5, i4)
            phhh (n2) = phhh (n2) + phah (n1) * evec (i3, i5)
         enddo
      enddo
c     
      return
      end
