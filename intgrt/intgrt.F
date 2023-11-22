
      program intgrt

c This program defines a grid of points and determines the value of
c the MO, the electron density, and several energy densities (such as
c the exchange density, the kinetic energy density, etc.) at each point
c (both at the SCF and correlated levels).  It also integrates over
c the gridpoints to determine total values.

c Written by Nevin Oliphant as hfdft, Quantum Theory Project,
c University of Florida, 1993.

c ReWritten by Nevin Oliphant and Sullivan Beck as potdens, 
c Quantum Theory Project, University of Florida, 1994.

c ReWritten by Nevin Oliphant and Sullivan Beck as intgrt, 
c Quantum Theory Project, University of Florida, 1995.

      implicit none
#include <aces.h>
#include <flags.com>
#include <flags2.com>
#include <iuhf.com>
#include <mol.com>
#include <kscore.com>
#include <symm2_ks.com>
#include <timeinfo.com>
#include <int.com>
#include <intgrtflags.com>

      integer intnumradpts,max_angpts,islateron
      integer maxmem, i0tmp,Length 
      logical kshf
      call callstack_init('INTGRT')

c INITIALIZATION

c Initializes constants used throughout the program.
      int_ks = .false.
      call initintgrt(kshf,intnumradpts,islateron)

      call timer(0)
      call aces_init_rte
      call aces_ja_init
      call getrec(1,'JOBARC','IFLAGS  ',100,iflags)
      call getrec(1,'JOBARC','IFLAGS2 ',500,iflags2)
      iuhf=0
      if (iflags(11).ge.1) iuhf=1
      call getrec(-1,'JOBARC','COMPNIRR',1,nirrep)
      call sb_com_symm2_ks
      maxmem=iflags(36)
      maxicore=500000

c      do memknown=0,1
      do memknown=-1,-1
#ifdef _DEBUG_INTGRT
         write(*,*) '@INTGRT: Starting main loop.'
         write(*,*) '         memknown is ', memknown
#endif /* _DEBUG_INTGRT */
         i0tmp = 1
         if (memknown.ne.0) then
            kscore(1) = 0
            do while ((kscore(1).eq.0).and.(maxmem.gt.0))
               call aces_malloc(maxmem,kscore,i0tmp)
               if (kscore(1).eq.0) maxmem = maxmem - 1024*1024
            end do
         end if
         call ks_mem_init(kscore(i0tmp),maxmem)
         call sb_com_mol

c Determine analytic values for the kinetic energy, nuclear 
c attraction energy, coulomb energy and exchange energy using
c the HF and correlated relaxed density.  This should be the first
c routine called (since it also determines the current print level).

c Start here for either SCF or Relaxed correlated density integrations

         Call getrec(0,"JOBARC","AONTORBA",Length,kscore)
         If (Length .gt.0) Then
            idns=1
         Else
            idns=0
         Endif

         call anlytc

c Set up grid for numerical integration

         call setupgrd(max_angpts)

c Start here again for natural orbitals after integrating over scf
c density
  900    continue

c Get the orbitals (SCF or natural orbitals) and determine the MO
c to primitive function transformation matrix

         call primvecint

c Set up the actual numerical grid for each symmetry unique atom and
c integrate over the symmetry unique octants of each symmetry unique
c atom

         call numintint(.false.)

c        if (pltflag.eq.1) call plotgrid

C This logic is obsolete and the control is tranfered to aces2
C driver. Ajith Perera, 06/2018.
c      o if there is a relaxed density do everything again
C         if (iflags(2).ne.0) then
C            if (iflags(2).ne.41) then
C               if (idns.eq.0) then
C                  idns=1
C                  goto 900
C               end if
C            end if
C         end if

         call sb_mem_fin

#ifdef _DEBUG_INTGRT
         call final_check
#endif /* _DEBUG_INTGRT */

c     end do memknown=0,1
      end do

      call aces_ja_fin
      call timer(1)
      write(*,9020) timetot
 9020 format(t3,'@CHECKOUT-I, Total execution time : ',f12.4,
     &       ' seconds.')

      call callstack_term

      stop
      end

