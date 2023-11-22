
      subroutine mrcc_branch

#include "aces2.h"
#include "flags.h"

      integer istat, ishell
      character*79 szHyper, szMolden

#ifdef _TRACE_ACES2
      write(*,*) '@ACES2: Entered ', __FILE__
#endif
      write(6,*) ' Entered mrcc_branch '
      write(6,*) ' geom_opt : ', geom_opt
      write(6,*) ' vib_specs: ', vib_specs
      write(6,*)
c
      if (geom_opt .and. vib_specs) then
         write(6,*) ' @MRCC_new: optimize geometries + frequencies'
         if (analytical_gradient) then
            istat=0
            do while (istat.eq.0)
               call mrcc_grad
               call a2putrec(1,'JOBARC','DIRTYFLG',1,0)
               call runit('xjoda')
               call a2getrec(1,'JOBARC','JODADONE',1,istat)
            end do
         else
            call opt_numrcl_4mrcc
         end if
         Call Prep4_post_opt_freq 
         call mrcc_vib
         return
      end if


      if (geom_opt) then
         if (analytical_gradient) then
            istat=0
            do while (istat.eq.0)
               call mrcc_grad
               call a2putrec(1,'JOBARC','DIRTYFLG',1,0)
               call runit('xjoda')
               call a2getrec(1,'JOBARC','JODADONE',1,istat)
            end do
         else
            call opt_numrcl_4mrcc
         end if
         return
      end if

      if (vib_specs) then
         call mrcc_vib
         return
      end if
c
      if (iflags2(3) .eq. 1) then
cmn
c Resonance Raman calculation
cmn
         write(6,*) ' @mrcc_branch: enter mrcc_grad'
          call mrcc_grad
C
c 
       else if (iflags2(h_IFLAGS2_prep_vibron).eq.1 .or.
     $         iflags2(h_IFLAGS2_geo_sample).eq.1 .or.
     $         iflags2(h_IFLAGS2_gen_quartic).eq.1 .or.
     $          iflags2(h_IFLAGS2_scan_normal).eq.1 ) then
c
         write(6,*) ' @mrcc_branch: enter mrcc_heff'
C
          call mrcc_heff

       else
         write(6,*) ' @mrcc_branch: enter energy_4mrcc', mrcc
          call energy_4mrcc
          szMolden  = 'xa2proc molden'
          call runit(szMolden)
      end if
c
      return
      end

