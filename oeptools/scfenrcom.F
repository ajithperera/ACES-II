      subroutine scfenrcom (
     &     nbas, luint, lnbuf,
     &     ibuf, buf, dens, nucr, onee, kine, eext, coul, hfex)
c
c This routine calculates the components of the total SCF energy given
c the SCF density matrix.
c
c In: nbas, luint, lnbuf, dens
c Scr: ibuf, buf
c
c Igor Schweigert, Aug 2004
c $Id: scfenrcom.FPP,v 1.1.1.1 2008/05/29 19:35:40 taube Exp $
c 
      implicit none
c
c     Common blocks
c
#include <machsp.com>
c
c     Arguments
c     
      integer
     &     nbas, luint, lnbuf, ibuf (lnbuf)
c     
      double precision
     &     buf (lnbuf),  dens (*), 
     &     nucr, onee, kine, eext, coul, hfex
c
c     Local variables
c
      integer
     &     nut, iand, ishft, 
     &     n, i, j, k, l, ij, kl, ik, jl, il, jk,
     &     min, max, iupki, iupkj, iupkk, iupkl, i_aa_cr
c     
      character*80 fname
c
      double precision x
c
c     Define the VMOL index operators.
c     
      iupki (i) = iand (i,ialone)
      iupkj (i) = iand (ishft(i,-ibitwd),ialone)
      iupkk (i) = iand (ishft(i,-2*ibitwd),ialone)
      iupkl (i) = iand (ishft(i,-3*ibitwd),ialone)
      i_aa_cr (i, j) = min(i,j) + (max(i,j)*(max(i,j)-1))/2
c
c     Get the nuclear repulsion energy.
c
      call getrec (1, 'JOBARC', 'NUCREP  ', iintfp, nucr)
c
c     Calculate the one-electron, the Coulomb and non-local exchange
c     energies: <> Initialize the energies being calculated. <> Get the
c     name of the VMOL file from JOBARC, open the file. <> Scale the
c     offdiagonal elements of the density matrix by a factor of two. <>
c     Locate the one-electron integrals, read the integrals and
c     calculate the one-electron energy. <> Locate the kinetic-energy
c     intergrals, read them and calculate the kinetic energy. <>
c     Calculate the V external energy by subtracting the kinetic energy
c     from the one-electron energy. <> Restore the density matrix. <>
c     Locate the two-electrons integrals, read the integrals and
c     calculate the Coulomb and non-local exchange energies: <<>> Read
c     the value and indices of the current integral. <<>> Define the
c     symmetry of the integral. <<>> Update the energies. <> Close the
c     VMOL file.
c
      onee = 0.d0
      kine = 0.d0
      coul = 0.d0
      hfex = 0.d0
c 
      call gfname('IIII    ',fname,n)
      open(unit=luint,file=fname(1:n),form='UNFORMATTED',
     &     access='SEQUENTIAL')
c     
      call dscal (nbas*(nbas+1)/2, 2d0, dens, 1)
      do i = 1, nbas
         dens (i*(i+1)/2) = dens (i*(i+1)/2) * .5d0
      enddo
c     
      call locate(luint, 'ONEHAMIL')
      nut = lnbuf
      do while (nut.eq.lnbuf)
         read(luint) buf, ibuf, nut
         do n = 1, nut
            onee = onee + dens (ibuf(n)) * buf (n)
         end do
      end do
c
      call locate(luint, 'KINETINT')
      nut = lnbuf
      do while (nut.eq.lnbuf)
         read(luint) buf, ibuf, nut
         do n = 1, nut
            kine = kine + dens (ibuf(n)) * buf (n)
         end do
      end do
c
      eext = onee - kine
c
      call dscal (nbas*(nbas+1)/2, .5d0, dens, 1)
      do i = 1, nbas
         dens (i*(i+1)/2) = dens (i*(i+1)/2) * 2d0
      enddo
c      
      call locate(luint,'TWOELSUP')
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
            coul = coul + 4.0d0 * dens (ij) * x * dens (kl)
            hfex = hfex - 1.0d0 * dens (jl) * x * dens (ik)
            hfex = hfex - 1.0d0 * dens (jk) * x * dens (il)
         enddo
      enddo
c     
      close(unit=luint,status='KEEP')
c
      return
      end
