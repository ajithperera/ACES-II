      subroutine write_vibronic(freq, numvib, scr)
c     
c     Write an initial vibronic coupling file, which can be processed by VIBRON
c     and turned into a full Vibronic_Coupling file using info on PNTHEFF
c     which is written by write_pntheff
c     
      implicit none

#include <aces.h>
#include <vibron.com>
#include <machsp.com>
#include <mrcc_heff.com>
#include "flags.h"
#include "jodaflags.com"
c
      integer ione, i_sided, ijunk, irp, nmode2, npoint, irpi, irpj,
     $     i, j, k, nirpf, iunit, perline, ioff, ndim, ndim2,
     $     iref, ipnt, igrd, ifcm, ifreq, itop,
     $     ioff1, ioff2, ioff3, ioff4, numvib, icase
      character*4 doit, fnameirp
      character*8 labelf(14)
      character*80 string80
      character*6 nameheff0, ndimheff0
      M_REAL scr(*), stpsiz, freq(numvib), fact, e0, fred, two
      logical reduced

      two = 2.0d0

c     
c     Documentation on conversion factors:
c     
c     q_i = sqrt(M omega / hbar) * X_i = sqrt(omega * m_u / hbar) * a_0 * Q_i
c     = sqrt(2 pi c nu_i^tilde * m_u / hbar) * a_0 * Q_i
c     = f * sqrt(nu_i^tilde) * Q_i , where nu_i^tilde is the frequency of the
c     parent state normal mode in cm-1, Q_i is the mass weighted normal coordinate, 
c     while q_i = 1/sqrt(2) [ a_i^+ + a_i] is the dimensionless normal coordinate.
c     The conversion factor f evaluates to f=0.0911355
c     
      fred = 0.0911355d0
      fact = 27.2113956324672d0 / 0.0911355d0

 888  format(A50,A20)
c     
c     get some of the basic info from Jobarc / Flags
c     
      ione = 1
      i_sided = iflags2(h_IFLAGS2_grid_vibron)
      CALL GETREC(20,'JOBARC','NUMPOINT',IONE,npoint)
c     
      CALL GETREC(-1,'JOBARC','DANGERUS',IONE,IJUNK)
      if(iflags(79).eq.0.and.ijunk.eq.0)doit='FULL'
      if(iflags(79).eq.1.or.ijunk.eq.1)doit='COMP'
      CALL GETREC(20,'JOBARC',DOIT//'NIRX',IONE,nirrepf)
      call getrec(20, 'JOBARC', 'NMODES_F', nirrepf, numirpf)
c
      if (iflags(57) .lt. 0) then
         reduced = .true.
         stpsiz = dfloat(iflags(57)) /1000.0d0
      else
         reduced = .false.
         stpsiz=dfloat(iflags(57))*10.0d-5
      endif
c     
      write(6,*) ' @write_vibronic : normal mode information'
      write(6,800) nirrepf
 800  format('  Number of symmetries : ', I8)
      write(6,*) ' Number of normal modes per symmetry block '
      write(6,*) '   IRREP       # of Modes '
      do irp = 1, nirrepf
         write(6,801) irp, numirpf(irp)
      enddo
 801  format(2I8)
c     
      nmode = 0
      do irp = 1, nirrepf
         nmode = nmode + numirpf(irp)
      enddo
      write(6,802) nmode
 802  format( '  Total number of normal modes ', I8)
c     
      nmode2 = 0
      if (i_sided .ge. 3 .and. i_sided .ne. 5) then
         do irp = 1, nirrepf
            nmode2 = nmode2 + numirpf(irp)*(numirpf(irp)-1) / 2
         enddo
      endif
      if (i_sided .eq. 4 .or. i_sided .eq. 7) then
         do irpi = 1, nirrepf
            do irpj = 1, irpi - 1
               nmode2 = nmode2 + numirpf(irpi) * numirpf(irpj)
            enddo
         enddo
      endif
c     
      write(6,804) nmode2
 804  format(' Total number of quadratic displacements ', I8)
c     
c     Up to now the same as in write_pntheff 
c     


      do i = 1, nirrepf
         do k = 1, 8
            labelf(i)(k:k) = ' '
         enddo
      enddo
c     
      open(unit=70, file='IRPNAMES',status='OLD')
      rewind(70)
      read(70,*) string80
      read(70,*) nirpf
      do i = 1, nirpf
         read(70,*) labelf(i)
      enddo
      close(70)
c     
      if (nirpf .ne. nirrepf) then
         write(6,*) ' @write_vibronic, Something suspicious : '
         write(6,*) ' nirpf, nirrepf ', nirpf, nirrepf
c         call aces_exit(1)
      endif
c     
c     All results will be written to file 'Vibronic_Coupling'
c     
      iunit = 172
      perline = 8
      open(UNIT=iunit, FILE='Vibronic_Coupling',
     $     FORM='FORMATTED')
      rewind(iunit)
c     
      write(iunit,888) '  Vibronic Coupling Constants '
      write(iunit,888) '  All units are in eV or dimensionless '
      write(iunit,901)
      write(iunit, 888) ' '
      if (i_sided .eq. 1) then
         write(iunit, 888) ' Single-sided differentiation has been used'
         write(iunit, 888) ' Gradients (linear couplings) only'
      elseif (i_sided .eq. 2) then
         write(iunit,888) ' Double-sided differentiation has been used'
         write(iunit,888) ' Gradients and diagonal elements of Hessian.'
      elseif (i_sided .eq. 3) then
         write(iunit,888) ' Double-sided differentiation has been used'
         write(iunit,888) ' Gradients and totally symmetric Hessian '
      elseif (i_sided .eq. 4) then
         write(iunit,888) ' Double-sided differentiation has been used'
         write(iunit,888) ' Gradients and Complete Hessian '
      elseif (i_sided .eq. 5) then
         write(iunit,888) ' Double-sided differentiation has been used'
         write(iunit,888) ' Gradients and diagonal elements of Hessian',
     $        ' + Cubic/Quartic.'
      elseif (i_sided .eq. 6) then
         write(iunit,888) ' Double-sided differentiation has been used'
         write(iunit,888) ' Gradients and totally symmetric Hessian ',
     $        '+ Cubic/Quartic.'
      elseif (i_sided .eq. 7) then
         write(iunit,888) ' Double-sided differentiation has been used'
         write(iunit,888) ' Gradients and Complete Hessian ',
     $        '+ Cubic/Quartic.'
      endif
      write(iunit,888) ' '
c     
      call getrec(1, 'JOBARC', 'FREQS_0 ', numvib*iintfp,
     $     freq)
      write(iunit,888) ' Parent state Symmetric modes in eV '
      ioff = 1
      do irp = 1, 1
         write(iunit,704)
     $        (freq(ioff+i) /8065.5d0, i=0,numirpf(irp)-1)
         ioff = ioff + numirpf(irp)
      enddo
 704  format(20F12.8)
c     
      write(iunit,888) ' Parent state Non-Symmetric modes in eV '
      write(iunit,704) (freq(i) /8065.5d0, i=ioff, nmode)
      write(iunit,888) ' '
c     
c process ground state (a bit stupid for historical reasons ...
c
c     
         ndim = 1
         ndimheff0 = 'Neff_0'
         nameheff0 = 'Heff_0'
c         call getrec(-1, 'JOBARC', ndimheff(icase), ione, ndim)
         write(6,805) ndimheff0, ndim
c     
         write(iunit,888) ' Vibronic Coupling elements for case ',
     $        nameheff0
         write(iunit,901)
         write(iunit,901)
         write(iunit, 888)
c     
c     write common info for this calculation
c     
         write(iunit,811) i_sided
         write(iunit,812) nirrepf
         write(iunit,813) (numirpf(irp), irp=1,nirrepf)
c     
         call getrec(1, 'JOBARC', 'FREQS_0 ', numvib*iintfp,
     $        freq)
         ioff = 1
         do irp = 1, nirrepf
            if (numirpf(irp) .ne. 0) then
               do i = 1, 4
                  fnameirp(i:i) = labelf(irp)(i:i)
               enddo
               write(iunit,814) irp, numirpf(irp), fnameirp,
     $              (freq(ioff+i), i=0,numirpf(irp)-1)
            endif
            ioff = ioff + numirpf(irp)
         enddo
         write(iunit,888)
c     
c     process heff part
c     
         ndim2 = ndim*ndim
c     
         iref = 1
         ipnt = iref + 1
         igrd = ipnt + npoint
         ifcm = igrd + nmode
         ifreq = ifcm + nmode+nmode2
         itop = ifreq+nmode
c     
         do k = 1, nmode
            scr(ifreq+k-1) = fact / sqrt(abs(freq(k)))
         enddo
c     
c     first treat parent energy
c     
         call getrec(20,'JOBARC','REFPAR_E',
     $        iintfp, scr(iref))
         write(6,*) ' Ground state energy at reference point ',
     $        scr(iref)
c     
         call getrec(20,'JOBARC','PNTENERG',
     $        npoint*iintfp, scr(ipnt))
c     
         if (i_sided .eq. 1) then
            do k = 1, nmode
               scr(igrd+k-1) = 
     $              (scr(ipnt+k-1) - scr(iref)) / stpsiz
            enddo
         else
c     
c     two-sided gradients
c     
            ioff1 = ipnt-1
            ioff2 = ipnt + nmode -1
            do k = 1, nmode
               scr(igrd+k-1) = 
     $              (scr(ioff1+k) -
     $              scr(ioff2+k)) / (2.0d0*stpsiz)
            enddo
c     
c     diagonal elements of Hessian
c     
            do k = 1, nmode
               scr(ifcm+k-1) =
     $              (scr(ioff1+k) +
     $              scr(ioff2+k) - 2.0d0 * scr(iref)) / (stpsiz**2)
            enddo
c     
         endif
c     
         if (i_sided .ge. 3 .and. i_sided .ne. 5) then
c     
c     off-diagonal elements of Hessian (totally symmetric only)
c     
            ioff1 = ipnt - 1
            ioff2 = ioff1 + nmode
            ioff3 = ioff2 + nmode
            ioff4 = ioff3 + nmode2
            k=0
            do irp = 1, nirrepf
               do i = 1, numirpf(irp)
                  do j = i+1, numirpf(irp)
                     k=k+1
                     write(6,*)
     $                    scr(ioff3+k), scr(ioff4+k),
     $                    scr(ioff1+i), scr(ioff1+j),
     $                    scr(ioff2+i), scr(ioff2+j),
     $                    two*scr(iref)
c     
                     scr(ifcm+nmode+k-1) =
     $                    (scr(ioff3+k) + scr(ioff4+k)
     $                    - scr(ioff1+i) - scr(ioff1+j)
     $                    - scr(ioff2+i) - scr(ioff2+j)
     $                    + two*scr(iref) )
                     scr(ifcm+nmode+k-1) = scr(ifcm+nmode+k-1)
     $                    / (stpsiz**2)
                  enddo
               enddo
               ioff1 = ioff1 + numirpf(irp)
               ioff2 = ioff2 + numirpf(irp)
               ioff3 = ioff3 + numirpf(irp)*(numirpf(irp)-1) / 2
               ioff4 = ioff4 + numirpf(irp)*(numirpf(irp)-1) / 2
            enddo
c     
         endif
c     
         write(iunit,888)
         write(iunit,888) ' Gradients parent energy along normal modes'
         write(iunit, 901)
         write(iunit,888)
         do k = 1, nmode
            write(iunit,912) k, scr(igrd+k-1)
         enddo
c
         if (reduced) then
            fact = 27.2114 * 8065.5 
         if (i_sided .ne. 1) then
            write(iunit,888)
            write(iunit,888) ' Diagonal Hessian parent energy '
            write(iunit, 901)
            write(iunit,888)
            do k = 1, nmode
               write(iunit,914) k, 
     $              abs(scr(ifcm+k-1))*fact, 
     $              abs(scr(ifcm+k-1)) * fact / 8065.5
            enddo
         endif
c     
         if (.false. .and. (i_sided .ge. 3 .and. i_sided .ne. 5)) then
            write(6,*) ' off-diagonal elements hessian parent E'
c     
            k=0
            ioff = ifcm+nmode
            do irp = 1, nirrepf
               do i = 1, numirpf(irp)
                  do j = 1, i-1
                     write(iunit,916) irp, i, j,
     $                    abs(scr(ioff+k))*fact,
     $                    abs(scr(ioff+k)) * fact / 8065.5
                     k=k+1
                  enddo
               enddo
            enddo
         endif
         write(iunit,888)

         else
         if (i_sided .ne. 1) then
            write(iunit,888)
            write(iunit,888) ' Diagonal Hessian parent energy '
            write(iunit, 901)
            write(iunit,888)
            do k = 1, nmode
               write(iunit,914) k, 
     $              sqrt(abs(scr(ifcm+k-1)))*5.14048D03,
     $              sqrt(abs(scr(ifcm+k-1))) *5.14048D03 / 8065.5
            enddo
         endif
c     
         if (.false. .and. (i_sided .ge. 3 .and. i_sided .ne. 5)) then
            write(6,*) ' off-diagonal elements hessian parent E'
c     
            k=0
            ioff = ifcm+nmode
            do irp = 1, nirrepf
               do i = 1, numirpf(irp)
                  do j = 1, i-1
                     write(iunit,916) irp, i, j,
     $                    sqrt(abs(scr(ioff+k)))*5.14048D03,
     $                    sqrt(abs(scr(ioff+k))) *5.14048D03 / 8065.5
                     k=k+1
                  enddo
               enddo
            enddo
         endif
         write(iunit,888)
         endif
c     
c     This was known information concerning the Parent State.
c     
c         call getrec(20,'JOBARC',refheff(icase),
c     $        ndim2*iintfp, scr(iref))
c     
         call getrec(20,'JOBARC','REFPAR_E',
     $        iintfp, e0)
         scr(iref) = e0
c     
c     subract reference energy from diagonal and convert to EV's
c     
         j = iref
         do i = 1, ndim
            scr(j) = scr(j) - e0
            j = j + ndim+1
         enddo
         call SSCAL(ndim2, 27.2113956324672d0, scr(iref), 1)
c     
         write(iunit,888) ' Reference Hamiltonian '
         call put_matrix(scr(iref), ndim, ndim, perline, iunit,.false.)
c     
         write(iunit,888)
         write(iunit,888) ' All information concerning method ',
     $        nameheff0
         write(iunit,888)
         write(6,*) ' @write_vibronic: done with ', nameheff0
c     
c
      do 10 icase = 1, ncases
c     
         ndim = 0
         call getrec(-1, 'JOBARC', ndimheff(icase), ione, ndim)
         write(6,805) ndimheff(icase), ndim
 805     format(' Process vibronic_coupling ', A12, /,
     $        '   # of electronic states: ', I10)
         if (ndim .eq. 0) goto 10
c     
         write(iunit,888) ' Vibronic Coupling elements for case ',
     $        nameheff(icase)
         write(iunit,901)
         write(iunit,901)
 901     format(70('-'))
         write(iunit, 888)
c     
c     write common info for this calculation
c     
         write(iunit,811) i_sided
 811     format(t3,'Vibronic grid in ACESII : ', i4)
         write(iunit,812) nirrepf
 812     format(t3,'Total number of vibrational irreps', i4)
         write(iunit,813) (numirpf(irp), irp=1,nirrepf)
 813     format(t3,'Number of modes per symmetry', 20i4)
c     
         call getrec(1, 'JOBARC', 'FREQS_0 ', numvib*iintfp,
     $        freq)
         ioff = 1
         do irp = 1, nirrepf
            if (numirpf(irp) .ne. 0) then
               do i = 1, 4
                  fnameirp(i:i) = labelf(irp)(i:i)
               enddo
               write(iunit,814) irp, numirpf(irp), fnameirp,
     $              (freq(ioff+i), i=0,numirpf(irp)-1)
            endif
            ioff = ioff + numirpf(irp)
         enddo
 814     format(t3,'[',i2,']', t9,i4, t15, a4, t25, 20F8.2)
         write(iunit,888)
c     
c     process heff part
c     
         ndim2 = ndim*ndim
c     
         iref = 1
         ipnt = iref + 1
         igrd = ipnt + npoint
         ifcm = igrd + nmode
         ifreq = ifcm + nmode+nmode2
         itop = ifreq+nmode
c     
         do k = 1, nmode
            scr(ifreq+k-1) = fact / sqrt(abs(freq(k)))
         enddo
c     
c     first treat parent energy
c     
         call getrec(20,'JOBARC','REFPAR_E',
     $        iintfp, scr(iref))
         write(6,*) ' Ground state energy at reference point ',
     $        scr(iref)
c     
         call getrec(20,'JOBARC','PNTENERG',
     $        npoint*iintfp, scr(ipnt))
c     
         if (i_sided .eq. 1) then
            do k = 1, nmode
               scr(igrd+k-1) = 
     $              (scr(ipnt+k-1) - scr(iref)) / stpsiz
            enddo
         else
c     
c     two-sided gradients
c     
            ioff1 = ipnt-1
            ioff2 = ipnt + nmode -1
            do k = 1, nmode
               scr(igrd+k-1) = 
     $              (scr(ioff1+k) -
     $              scr(ioff2+k)) / (2.0d0*stpsiz)
            enddo
c     
c     diagonal elements of Hessian
c     
            do k = 1, nmode
               scr(ifcm+k-1) =
     $              (scr(ioff1+k) +
     $              scr(ioff2+k) - 2.0d0 * scr(iref)) / (stpsiz**2)
            enddo
c     
         endif
c     
         if (i_sided .ge. 3 .and. i_sided .ne. 5) then
c     
c     off-diagonal elements of Hessian (totally symmetric only)
c     
            ioff1 = ipnt - 1
            ioff2 = ioff1 + nmode
            ioff3 = ioff2 + nmode
            ioff4 = ioff3 + nmode2
            k=0
            do irp = 1, nirrepf
               do i = 1, numirpf(irp)
                  do j = i+1, numirpf(irp)
                     k=k+1
                     write(6,*)
     $                    scr(ioff3+k), scr(ioff4+k),
     $                    scr(ioff1+i), scr(ioff1+j),
     $                    scr(ioff2+i), scr(ioff2+j),
     $                    two*scr(iref)
c     
                     scr(ifcm+nmode+k-1) =
     $                    (scr(ioff3+k) + scr(ioff4+k)
     $                    - scr(ioff1+i) - scr(ioff1+j)
     $                    - scr(ioff2+i) - scr(ioff2+j)
     $                    + two*scr(iref) )
                     scr(ifcm+nmode+k-1) = scr(ifcm+nmode+k-1)
     $                    / (stpsiz**2)
                  enddo
               enddo
               ioff1 = ioff1 + numirpf(irp)
               ioff2 = ioff2 + numirpf(irp)
               ioff3 = ioff3 + numirpf(irp)*(numirpf(irp)-1) / 2
               ioff4 = ioff4 + numirpf(irp)*(numirpf(irp)-1) / 2
            enddo
c     
         endif
c     
         write(iunit,888)
         write(iunit,888) ' Gradients parent energy along normal modes'
         write(iunit, 901)
         write(iunit,888)
         do k = 1, nmode
            write(iunit,912) k, scr(igrd+k-1)
         enddo
 912     format(t3, 'Normal mode ', i4, t22, 'Gradient :', E16.10)
c
         if (reduced) then
            fact = 27.2114 * 8065.5 
         if (i_sided .ne. 1) then
            write(iunit,888)
            write(iunit,888) ' Diagonal Hessian parent energy '
            write(iunit, 901)
            write(iunit,888)
            do k = 1, nmode
               write(iunit,914) k, 
     $              abs(scr(ifcm+k-1))*fact, 
     $              abs(scr(ifcm+k-1)) * fact / 8065.5
            enddo
         endif
c     
         if (.false. .and. (i_sided .ge. 3 .and. i_sided .ne. 5)) then
            write(6,*) ' off-diagonal elements hessian parent E'
c     
            k=0
            ioff = ifcm+nmode
            do irp = 1, nirrepf
               do i = 1, numirpf(irp)
                  do j = 1, i-1
                     write(iunit,916) irp, i, j,
     $                    abs(scr(ioff+k))*fact,
     $                    abs(scr(ioff+k)) * fact / 8065.5
                     k=k+1
                  enddo
               enddo
            enddo
         endif
         write(iunit,888)

         else
         if (i_sided .ne. 1) then
            write(iunit,888)
            write(iunit,888) ' Diagonal Hessian parent energy '
            write(iunit, 901)
            write(iunit,888)
            do k = 1, nmode
               write(iunit,914) k, 
     $              sqrt(abs(scr(ifcm+k-1)))*5.14048D03,
     $              sqrt(abs(scr(ifcm+k-1))) *5.14048D03 / 8065.5
            enddo
         endif
 914     format(t3, 'Normal mode ', i4, t22, 
     $        ' Frequency : ', F8.2, ' cm-1', F12.6, ' eV')
c     
         if (.false. .and. (i_sided .ge. 3 .and. i_sided .ne. 5)) then
            write(6,*) ' off-diagonal elements hessian parent E'
c     
            k=0
            ioff = ifcm+nmode
            do irp = 1, nirrepf
               do i = 1, numirpf(irp)
                  do j = 1, i-1
                     write(iunit,916) irp, i, j,
     $                    sqrt(abs(scr(ioff+k)))*5.14048D03,
     $                    sqrt(abs(scr(ioff+k))) *5.14048D03 / 8065.5
                     k=k+1
                  enddo
               enddo
            enddo
         endif
 916     format(t3, '[', i1, ']', t8, i4, t15, i4, t23, F10.2,
     $        t35, 'cm-1', t40, F10.6, t52, 'eV')
         write(iunit,888)
         endif
c     
c     This was known information concerning the Parent State.
c     
         call getrec(20,'JOBARC',refheff(icase),
     $        ndim2*iintfp, scr(iref))
c     
         call getrec(20,'JOBARC','REFPAR_E',
     $        iintfp, e0)
         write(6,*) ' Do not subtract e0 '
         e0 = 0.0d0
c     
c     subract reference energy from diagonal and convert to EV's
c     
         j = iref
         do i = 1, ndim
            scr(j) = scr(j) - e0
            j = j + ndim+1
         enddo
         call SSCAL(ndim2, 27.2113956324672d0, scr(iref), 1)
c     
         write(iunit,888) ' Reference Hamiltonian '
         call put_matrix(scr(iref), ndim, ndim, perline, iunit,.false.)
c     
         write(iunit,888)
         write(iunit,888) ' All information concerning method ',
     $        nameheff(icase)
         write(iunit,888)
         write(6,*) ' @write_vibronic: done with ', nameheff(icase)
c     
 10   continue
c     
      close(unit=iunit, status='KEEP')
c     
      write(6,*) ' All done write_vibronic '
c     
      return
      end
