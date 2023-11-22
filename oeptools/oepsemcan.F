      subroutine oepsemcan (
     &     nbas, nocc, luint, lnbuf, ibuf, 
     &     buf, evec, dens, oneh,
     &     scraa, scrpp, scrhh, scr1, scr2)
c     
c This routine build the Fock matrix and transforms given MO
c coefficients into the semicanonical basis.
c
c In: nbas, nocc, luint, lnbuf, ibuf, buf, evec, scfh
c Scr: scraa, scr1, scr2, scrpp, scrhh
c Out: evec is modified
c
c Igor Schweigert, Mar 2004
c $Id: oepsemcan.FPP,v 1.1.1.1 2008/05/29 19:35:40 taube Exp $
c     
      implicit none
c
c     Arguments
c     
      integer
     &     nbas, nocc, luint, lnbuf, ibuf (lnbuf)
c     
      double precision
     &     buf (lnbuf), evec (nbas, nbas), dens (*), oneh (*),
     &     scraa (*), scrpp (*), scrhh (*),
     &     scr1 (nbas, nbas), scr2 (nbas, nbas)
c
c     Local variables
c
      integer
     &     nut, n, ij, kl, ik, jl, il, jk, info
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
#include <oepintind.com>
#define DEFINE_VARIABLES
#define N_BAS nbas
#define N_OCC nocc
#include <oepintind.com>
c
c     Messages to the user
c     
 5905 format (
     &     t3,'@OEPSEMCAN-F: DSYEV failed to converge an eigenvalue',/)
 5910 format (
     &     t3,'@OEPSEMCAN-F: Illegal argument to DSYEV, argument #',
     &     i2,/)
c
c     Calculate the Fock matrix: <> Initialize the Fock matrix with
c     zeros. <> Get the name of the VMOL file from JOBARC, open the
c     file, and shift the pointer to the begining of the 2e integral
c     record. <> Read the integrals and calculate the two-electron
c     contribution to the Fock matrix: <<>> Read the value and indices
c     of the current integral. <<>> Define the symmetry of the
c     integral. <<>> Update the corresponding integrals. <> Close the
c     VMOL file. <> Scale the diagonal elements by factor of two. <> Add
c     the one-electron contribution to the Fock matrix.
c
      call zero (scraa, n_aa)
c     
      call gfname('IIII    ',fname, n)
      open(unit=luint,file=fname(1:n),form='UNFORMATTED',
     &     access='SEQUENTIAL')
      call locate(luint,'TWOELSUP')
c     
      nut = lnbuf
      do while (nut.eq.lnbuf)
         read(luint) buf, ibuf, nut
         do n=1,nut
c     
            x = buf(n)
            i = iupki ( ibuf(n))
            j = iupkj ( ibuf(n))
            k = iupkk ( ibuf(n))
            l = iupkl ( ibuf(n))
c     
            ij = i_aa_cr (j, i)
            kl = i_aa_cr (l, k)
            ik = i_aa_cr (k, i)
            jl = i_aa_cr (l, j)
            il = i_aa_cr (l, i)
            jk = i_aa_cr (j, k)
c     
            if (i.eq.j) x = .5d0 * x
            if (k.eq.l) x = .5d0 * x
            if (ij.eq.kl) x = .5d0 * x
c     
            scraa (ij) = scraa (ij) + 2.d0 * dens (kl) * x
            scraa (kl) = scraa (kl) + 2.d0 * dens (ij) * x
c     
            scraa (ik) = scraa (ik) - 0.5d0 * dens (jl) * x
            scraa (jl) = scraa (jl) - 0.5d0 * dens (ik) * x
            scraa (il) = scraa (il) - 0.5d0 * dens (jk) * x
            scraa (jk) = scraa (jk) - 0.5d0 * dens (il) * x
         enddo
      enddo
c     
      close(unit=luint,status='KEEP')
c     
      do i = 1, n_a
         n = i_aa (i, i)
         scraa (n) = 2.d0 * scraa (n)
      enddo
c     
      call daxpy (n_aa, 1.d0, oneh, 1, scraa, 1)
c     
c     Rotate orbitals to semicanonical: <> Transform the first AO index
c     into an MO index. <> Transform the second AO index into particle
c     and hole index to get the PP and HH blocks of the Fock matrix. <>
c     Diagonalize the PP and HH blocks of the Fock matrix. <> Check if
c     diagonalization was performed successfully. <> Form the
c     transformation matrix. <> Transform the MO coefficients.
c     
      call oepao2mo (
     &     nbas, nocc, evec,
     &     scraa, scr1, scrpp, scr2, scrhh)
c     
      call zero (scr1, n_ma)
      do n = 1, n_pp
         j = i2_pp (n) 
         i = i1_pp (n, j)
         scr1 (n_h+i, n_h+j) = scrpp (n)
      enddo
      do n = 1, n_hh
         j = i2_hh (n)
         i = i1_hh (n, j)
         scr1 (i, j) = scrhh (n)
      enddo
c     
      call dsyev (
     &     'v', 'u', n_a, scr1, n_a, scrpp, scr2, n_ma, info)
c     
      if (info.ne.0) then
         if (info.gt.0) then
            write (6, 5905)
         else
            write (6, 5910) abs (info)
         endif
         call errex
      endif
c     
      call dcopy (n_ma, evec, 1, scr2, 1)
      call xgemm (
     &     'n', 'n', n_a, n_a, n_a,
     &     1.d0, scr2, n_a,
     &     scr1, n_a,
     &     0.d0, evec, n_a)
c     
      return
      end
