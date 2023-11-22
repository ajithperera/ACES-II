      subroutine write_vibron_input(freq, numvib, scr)
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
     $     ioff1, ioff2, ioff3, ioff4, numvib, icase,
     $     ncouple, nsymm
      character*4 doit, fnameirp
      character*8 labelf(14)
      character*50 string80, string
      M_REAL
     $     scr(*), stpsiz, freq(numvib), fact, e0, fred, two, x

      two = 2.0d0
c     

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
      stpsiz=dfloat(iflags(57))*10.0d-5
      call getrec(1, 'JOBARC', 'FREQS_0 ', numvib*iintfp,
     $     freq)
c     
      nmode = 0
      do irp = 1, nirrepf
         nmode = nmode + numirpf(irp)
      enddo
c     
c     All results will be written to file 'Vibronic_Input'
c     
      iunit = 172
      perline = 8
      open(UNIT=iunit, FILE='Vibron_Input',
     $     FORM='FORMATTED')
      rewind(iunit)
      nsymm = numirpf(1)
      ncouple = nmode - nsymm
c     
c     write common info for this calculation
c     
 922  format(i6, t20, ' ! units in eV')
 923  format('xxx', t20, ' ! Number of electronic states')
 924  format(i6, t20, ' ! Number of symmetric normal modes ')
 925  format(i6, t20, ' ! Number of non-A1 normal modes ')
 926  format(i6, t20, ' ! Number of Lanczos iterations ')
 927  format(50i4)
 928  format(F8.3)
 930  format(t2,A50)
c     
      i = 0
      write(iunit,922) i
      write(iunit, 923) 
      write(iunit,924) nsymm
      write(iunit,925) ncouple
      i = 500
      write(iunit,926) i
      do i = 1, nmode
         nquanta(i) = 0
      enddo
      write(iunit, 927) (nquanta(i), i = 1, nsymm)
      write(iunit,704)
     $     (freq(i) /8065.5d0, i=1, nsymm)
 704  format(50F12.6)
c     
      write(iunit,927) (nquanta(i), i = 1, ncouple)
      write(iunit,704) (freq(i) /8065.5d0, i=nsymm+1, nmode)
c
 660  format(20A)
      x = 0.1
      write(iunit,928) x
      x = 3.0
      write(iunit,928) x
      string='auto'
      write(iunit,*) string(1:4)
      write(iunit,*)
      string = '*vibron'
      write(iunit,660) string(1:7)
      string = 'short_title=cp'
      write(iunit,*) string(1:14)
      string = 'refine_coupling=on'
      write(iunit,*) string(1:18)
      string = 'cubic=on'
      write(iunit,*) string(1:8)
      string =  'quartic=on'
      write(iunit,*) string(1:10)
      string =  '*solve_coupling'
      write(iunit,660) string(1:16)
      string = 'purify_couplings=on'
      write(iunit,*) string(1:19)
      string = '*end'
      write(iunit,660) string(1:4)
      write(iunit,*)
c     
      close(unit=iunit, status='KEEP')
c     
      write(6,*) ' All done write_vibron_input '
c     
      return
      end
