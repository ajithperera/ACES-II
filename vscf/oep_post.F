      subroutine oep_post(iuhf, nbas, ibufln, luint, nocc, icore, 
     &                    maxmem, scfevl, scfden, etot, dmax, 
     &                    itertype)
c
      implicit none
c
#include <oep.com>
#include <machsp.com>
#include <flags.com>
#include <flags2.com>
#include <symm2.com>
c
      integer
     &     nbas, iuhf, maxmem, nocc(16), ibufln, icore(*), i000,
     &     i010, i020, i030, i040, i050, i060, ione, luint
       integer i
      double precision au2ev, etot, dmax, homo, 
     &                 scfevl (nbas), scfden (itriln (nirrep+1))

      character*1 asv_null, itertype

      parameter (au2ev = 27.2113957D0)
c 
c     OEP post-SCF block: <> Print out the converged OEP energy. <> If
c     local exchange potential is used: Report by how much the HOMO
c     condition on the OEP exchange potential is violated. Shift the OEP
c     eigenvalues to satisfy the HOMO condition if no auxiliary exchange
c     potential is used. <> Allocate memory for scratch arrays necessary
c     to calculate the Fock matrix. <> Check if there is enough memory
c     available. <> Calculate the Fock matrix and write it over the OEP
c     Hamiltonian, in case a conventional correlation method will be
c     used with OEP orbitals. <> Calculate the components of the total
c     energy and report them. (WARNING: this requires reading the
c     integral file again!)
c
c     Note that orbital energies for ROHF are calculated somewhere
c     between the final Fock matrix construction and energy
c     calculation. So it is bypassed if the OEP method is invoked. If
c     one wants to do ROHF-like OEP method, one has to fix this.
c
      asv_null=achar(0)
C
#ifdef _DEBUG_LVL0
            do i=1,nbas
               write(6,*) i,scfevl(i)
            end do
#endif
C
      if (OEP1 .or. OEP2) then
c
         write (6, 5960) 'E (OEP) = ', etot, itertype
c
         call getrec(20, 'JOBARC', 'OEPXHOMO', iintfp, homo)
         write (6, 5964) homo, homo*au2ev
         if (OEPAXP (1:5) .eq. 'NONE'//asv_null) then
             call daxpy (nbas, 1.d0, homo, 0, scfevl, 1)
             write (6, 5965) homo
         endif
#ifdef _DEBUG_LVL0
            write(6,*) 'after hmo cont'
            do i=1,nbas
               write(6,*) i,scfevl(i)
            end do 
#endif
         call putrec(20, 'JOBARC', 'SCFEVLA0', nbas*iintfp, scfevl)
C
         ione = 1
         i010 = ione
         i020 = i010 + (iuhf+1) * isqrln (nirrep+1) * iintfp
         i030 = i020 + (iuhf+1) * isqrln (nirrep+1) * iintfp
         i040 = i030 + mod((i020-i000),iintfp)
         i050 = i040 + ibufln * iintfp
         i060 = i050 + ibufln 
         call oepenrcom (
     &        iflags2 (153), nbas, nocc (1), luint, ibufln,
     &        icore (i050), icore (i040),
     &        scfevl, scfden) 
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
C
       
C --- Make the standard from the main program vscf.F
c
CSSS         call mkrhff(
CSSS     &        icore(i040), icore(i030), icore(i010),
CSSS     &        icore(i055), icore(i056),
CSSS     &        itriln(nirrep+1), nbas, nbfirr,
CSSS     &        icore(i057), ibufln, luint, .true.)
c
cSSS         call oepfiachk (
cSSS     &        nbas, nocc (1), naux,
cSSS     &        icore (i030), maxmem-i080+i000, i080-i000,
cSSS     &        icore (oepauxao),
cSSS     &        scfevl, icore (i010), icore (i020))
c
CSSS         goto 5999
      endif

      return
      end

