      subroutine oepenrcom (
     &     scftype, nbas, nocc, luint, lnbuf,
     &     ibuf, buf, eval, dens)
c
c This routine calculates and prints the components of the total
c energy. It also prints the order-by-order contribution to the total
c OEP energy if an OEP method has been used.
c
c In: nbas, nocc, luint, lnbuf, eval, dens
c Scr: ibuf, buf
c
c Igor Schweigert, Aug 2004
c $Id: oepenrcom.FPP,v 1.1.1.1 2008/05/29 19:35:39 taube Exp $
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
     &     scftype, nbas, nocc, luint, lnbuf, ibuf (lnbuf)
c     
      double precision
     &     buf (lnbuf), eval (nbas), dens (*)
c
c     Local variables
c
      integer
     &     ind, nut, iand, ishft,
     &     n, ij, kl, ik, jl, il, jk
c     
      character*80
     &     fname
c
      double precision
     &     x, nucr, onee, coul, hfex, exch, corr, enr0, enr1, enr2, etot
c
c     Define the two-electron integral offset operators and the
c     associated variables.
c
#define DEFINE_OPERATORS     
#include "oepintind.com"
#define DEFINE_VARIABLES
#define N_BAS nbas
#define N_OCC nocc
#include "oepintind.com"
c
c     Messages to the user.
c
 5900 format (
     &     t3,'@OEPENRCOM-I: The following interactions contribute',/,
     &     t3,'to the total energy',/)
 5903 format (
     &     t3,'@OEPENRCOM-I: KS DFT energy components are not yet',/,
     &     t3,'implemented.',/)
 5905 format (
     &     t3,'@OEPENRCOM-I: Unknown SCF type (SCF_TYPE = ',i2,')',/)
 5910 format (
     &     t3,a,1x,f20.10,' a.u.')
 5915 format (
     &     t3,a,1x,f20.10,' a.u.',1x,f20.10,' a.u.')
 5920 format (
     &     t3,'@OEPENRCOM-I: Order-by-order contributions to',/,
     &     t3,'the total energy are',/)
c
c     Get the nuclear repulsion energy.
c
      call getrec (1, 'JOBARC', 'NUCREP  ', iintfp, nucr)
c     Calculate the one-electron, the Coulomb and non-local exchange
c     energies: <> Initialize the energies being calculated. <> Get the
c     name of the VMOL file from JOBARC, open the file. <> Locate the
c     one-electron integrals, read the integrals and calculate the
c     one-electron energy. <> Locate the two-electrons integrals, read
c     the integrals and calculate the Coulomb and non-local exchange
c     energies: <<>> Read the value and indices of the current
c     integral. <<>> Define the symmetry of the integral. <<>> Update
c     the energies. <> Close the VMOL file.
c
      onee = 0.d0
      coul = 0.d0
      hfex = 0.d0
      call gfname('IIII    ',fname,n)
      open(unit=luint,file=fname(1:n),form='UNFORMATTED',
     &     access='SEQUENTIAL')
      onee = 0.d0
      call locate(luint, 'ONEHAMIL')
      nut = lnbuf
      do while (nut.eq.lnbuf)
         read(luint) buf, ibuf, nut
         do n = 1, nut
            onee = onee + dens (ibuf(n)) * buf (n)
         end do
      end do
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
c     Calculate the zero-order (or true SCF) energy.
c
      enr0 = 0.d0
      do n = 1, nocc
         enr0 = enr0 + 2.d0 * eval (n)
      enddo
c
c     Calculate the exchange, correlation, and total energies. They
c     depend on the SCF method used.
      if     (scftype.eq.0) then
         exch = hfex
         corr = 0.d0
         etot = nucr + onee + coul + hfex
         enr1 = etot - enr0
         enr2 = 0.d0
      elseif (scftype.eq.1) then
         write (6, 5903)
         call errex
      elseif (scftype.eq.2) then
         write (6, 5903)
         call errex 
      elseif (scftype.eq.3) then
         call getrec (1, 'JOBARC', 'OEP1ENER', iintfp, etot)
         exch = etot - nucr - onee - coul
         corr = 0.d0
         enr1 = etot - enr0
         enr2 = 0.d0
      elseif (scftype.eq.4) then
         call getrec (1, 'JOBARC', 'OEP2ENER', iintfp, etot)
         call getrec (1, 'JOBARC', 'OEP1ENER', iintfp, enr1)
         exch = enr1 - nucr - onee - coul
         corr = etot - enr1
         enr2 = etot - enr1
         enr1 = enr1 - enr0
      else
         write (6, 5905) scftype
         call errex
      endif
c     Report the energy components.
c
      write (6, 5900)
      write (6, 5915) 'Nuclear repulsion ',nucr,nucr
      write (6, 5915) 'One-electron      ',onee,nucr+onee
      write (6, 5915) 'Coulomb           ',coul,nucr+onee+coul
      write (6, 5915) 'Exchange          ',hfex,nucr+onee+coul+hfex
      write (6, 5915) 'Correlation       ',corr,nucr+onee+coul+hfex+corr
      write (6, 5910) 'Total             ',etot
      write (6, *)
c
      write (6, 5920)
      write (6, 5915) 'Zero-order        ',enr0,enr0
      write (6, 5915) 'First-order       ',enr1,enr0+enr1
      write (6, 5915) 'Second-order      ',enr2,enr0+enr1+enr2
      write (6, 5910) 'Total             ',etot
      write (6, *)
c
      return
      end
