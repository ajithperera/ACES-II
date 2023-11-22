      subroutine oep1ints (
     &     nbas, nocc, naux, luint, lnbuf, axptype, ibuf, 
     &     buf, evec, dens, scfh, scrah,
     &     auxaa, auxph, hfxaa, hfxph)
c
c This routine takes care of all the OEP1-related integrals. It reads
c (1) AO two-electron integrals from the VMOL integral file, (2) AO MEs
c of the auxiliary potential if requested, and (3) the AO MEs of the
c auxiliary functions. Then the Coulomb AOMEs are added to the SCF
c hamiltonian matrix (supplied via SCFH), the auxiliary and non-local
c exchange AOMEs are tranformed to the MO basis.
c
c If the auxiliary potential has been requested then its AOMEs are added
c to the SCF Hamiltonian matrix as well and subtracted from the
c non-local exchange AOMEs.
c
c Note that unlike OEP2INTS, here the auxiliary integrals are read in
c after two-electron integrals, because we use AUXAA as a scratch to
c store the Coulomb AOMEs (and auxiliary potential AOMEs if needed.)
c
c In: nbas, nocc, luint, lnbuf, ibuf, buf, evec, scfh
c Scr: scrah
c Out: auxaa, auxph, hfxaa, hfxhh
c
c Igor Schweigert, Jan 2004
c $Id: oep1ints.FPP,v 1.2 2008/06/03 15:57:36 taube Exp $
c     
      implicit none
c
c     Arguments
c     
      integer
     &     nbas, nocc, naux, luint, lnbuf, ibuf (lnbuf)
c
      character*(*)
     &     axptype
c     
      double precision
     &     buf (lnbuf), evec (nbas, nbas), dens (*),  scfh (*),
     &     scrah (*), auxaa (*), auxph (*), hfxaa (*), hfxph (*)
c
c     Local variables
      logical oepint_exsist
      integer
     &     ind, nut, iand, ishft,
     &     n1, n2, n3, i1, i2, i3, i4, 
     &     n, ij, kl, ik, jl, il, jk
c     
      character*80
     &     fname
      character*1 asv_null
c
      double precision x
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
c     Read the 2e integrals in AO basis, and calculate the Coulomb and
c     non-local exchange AOMEs: <> Initialize the arrays being
c     calculated. <> Get the name of the VMOL file from JOBARC, open the
c     file, and shift the pointer to the begining of the 2e integral
c     record. <> Read the integrals and calculate the AOMEs of the
c     Coulomb and non-local exchange operators (Coulomb AOMEs will be
c     stored in SCRMA): <<>> Read the value and indices of the current
c     integral. <<>> Define the symmetry of the integral. <<>> Update
c     the corresponding integrals. <> Close the VMOL file. <> Scale the
c     diagonal elements by factor of two. <> Add the Coulomb AOMEs to
c     the SCF Hamiltonian. 
c
      asv_null=achar(0)
      call zero (auxaa, n_aa)
      call zero (hfxaa, n_aa)
      call zero (hfxph, n_ph)
c     
      call gfname('IIII    ',fname,n1)
      open(unit=luint,file=fname(1:n1),form='UNFORMATTED',
     &     access='SEQUENTIAL')
      call locate(luint,'TWOELSUP')
c     
      nut = lnbuf
      do while (nut.eq.lnbuf)
         read(luint) buf, ibuf, nut
         do n1=1,nut
c
            x = buf(n1)
            i = iupki ( ibuf(n1))
            j = iupkj ( ibuf(n1))
            k = iupkk ( ibuf(n1))
            l = iupkl ( ibuf(n1))
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
            auxaa (ij) = auxaa (ij) + 2.d0 * dens (kl) * x
            auxaa (kl) = auxaa (kl) + 2.d0 * dens (ij) * x
c     
            hfxaa (ik) = hfxaa (ik) - 0.5d0 * dens (jl) * x
            hfxaa (jl) = hfxaa (jl) - 0.5d0 * dens (ik) * x
            hfxaa (il) = hfxaa (il) - 0.5d0 * dens (jk) * x
            hfxaa (jk) = hfxaa (jk) - 0.5d0 * dens (il) * x
         enddo
      enddo
c     
      close(unit=luint,status='KEEP')
c     
      do i = 1, n_a
         n = i_aa (i, i)
         auxaa (n) = 2.d0 * auxaa (n)
         hfxaa (n) = 2.d0 * hfxaa (n)
      enddo
c 
      call daxpy (n_aa, 1.d0, auxaa, 1, scfh, 1)
c
c     If requested, obtain the AOMEs of the auxiliary potential used to
c     correct the asymptotic behavior of the exchange OEP. <> If the
c     scaled Coulomb is used as the auxiliary potential, scale the
c     Coulomb AOMEs, update the SCF Hamiltonian, and non-local exchange
c     AOMEs. <> If the Slater potential is used as the auxiliary
c     potential, get the Slater AOMEs from the OEP integral file, update
c     the SCF Hamiltonian and non-local exchange AOMEs.
c
c
      if     (axptype (1:3) .eq. 'FA'//asv_null) then
c
         call dscal (n_aa, 1.d0/dble (2*nocc), auxaa, 1)
         call daxpy (n_aa,-1.d0, auxaa, 1, scfh, 1)
         call daxpy (n_aa, 1.d0, auxaa, 1, hfxaa, 1)
c     
      elseif (axptype (1:7) .eq. 'SLATER'//asv_null) then
       call getrec(1,'JOBARC','SLAT51AO',iintfp*n_aa,auxaa)

c         open(
c     &        unit = luint, file = 'SLAT51AB',
c     &        form = 'UNFORMATTED', access = 'SEQUENTIAL')
c         call locate (luint, 'SLAT51AB')
c     
c         nut = lnbuf
c         do while (nut.eq.lnbuf)
c            read(luint) buf, ibuf, nut
c            do ind = 1, nut
c               auxaa (ibuf (ind)) = buf(ind)
c            enddo
c         enddo
c     
         close (luint)
c
         call daxpy (n_aa, 1.d0, auxaa, 1, scfh, 1)
         call daxpy (n_aa,-1.d0, auxaa, 1, hfxaa, 1)
c
      endif
c
c     Calculate the PH part of non-local exchange MOME: <> Transform the
c     second AO index into a hole index. <> Transform the first AO index
c     into a particle index.
c     
      call zero (scrah, n_ah)
      do n1 = 1, n_ah
         i2 = i2_ah (n1)
         i1 = i1_ah (n1, i2)
         do i3 = 1, n_a
            n2 = i_aa_cr (i1, i3)
            scrah (n1) = scrah (n1) + hfxaa (n2) * evec (i3, i2)
         enddo
      enddo
c     
      do n1 = 1, n_ph
         i2 = i2_ph (n1)
         i1 = i1_ph (n1, i2) + n_h
         do i3 = 1, n_a
            n2 = i_ah (i3, i2)
            hfxph (n1) = hfxph (n1) + scrah (n2) * evec (i3, i1)
         enddo
      enddo
c
c     Read the auxiliary AOMEs from the OEPINT file and transform them
c     to MOs: <> Initialize the array being calculated.<> Open the
c     OEPINT file and locate the auxiliary integrals record. <> Read the
c     integrals. <> Close the OEPINT file. <> For each auxiliary
c     function, transform the AOMEs into PH part of the MOMEs.
c
      call zero (auxaa, n_aa * naux)
      call zero (auxph, n_ph * naux)
      inquire (file = 'OEPINT', exist = oepint_exsist)
         if (.not. oepint_exsist) then
            write (6,*)    'crappppp'
            call errex
         endif    
      open(
     &     unit = luint, file = 'OEPINT',
     &     form = 'UNFORMATTED', access = 'SEQUENTIAL')
      call locate (luint, 'AUX3CNTR')
c     
      nut = lnbuf
      do while (nut.eq.lnbuf)
         read(luint) buf, ibuf, nut
         do ind = 1, nut
            auxaa (ibuf (ind)) = buf(ind)
         enddo
      enddo
c     
      close (luint)
c  
      do i4 = 1, naux
         call zero (scrah, n_ah)
         do n1 = 1, n_ah
            i2 = i2_ah (n1)
            i1 = i1_ah (n1, i2)
            do i3 = 1, n_a
               n2 = i_aa_cr (i1, i3) + n_aa * (i4-1)
               scrah (n1) = scrah (n1) + auxaa (n2) * evec (i3, i2)
            enddo
         enddo
         do n1 = 1, n_ph
            i2 = i2_ph (n1)
            i1 = i1_ph (n1, i2) + n_h
            n3 = n1 + n_ph * (i4-1)
            do i3 = 1, n_a
               n2 = i_ah (i3, i2)
               auxph (n3) = auxph (n3) + scrah (n2) * evec (i3, i1)
            enddo
         enddo
      enddo
c     
      return
      end
