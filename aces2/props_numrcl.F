      Subroutine Props_numrcl(Calc_level, Energy, Raman_ints)

#include "aces2.h"
#include "parallel_aces.com"

      Character*4 Calc_level
      Logical Energy, Raman_ints 
   
#ifdef _TRACE_ACES2
      write(*,*) '@ACES2: Entered ', __FILE__
#endif
      istat = 0
      do while (istat.eq.0)
c
         If (Calc_level .Eq. "SCF ") Then
            If (Energy) Then
               Call Scf_eneg
            Else If (Raman_ints) Then
               Call Scf_eneg
               Call Runit('xvtran')
               Call Runit('xintprc')
               Call Runit(der_integral_package)
               Call Runit('xcphf')
            Endif
c
c Else for the Calc_level
c
         Else
c
           If (Energy) Then
              Call Scf_eneg 
              Call Change_orbital_space
              Call Prep4_post_scf
              Call Post_scf_eneg
           Else If (Raman_ints) Then
              Call Scf_eneg 
              Call Runit('xvprops')
              Call Change_orbital_space
              Call Prep4_post_scf 
              Call Post_scf_eneg
              Call Runit('xlambda')
              Call Runit('xvcceh')
              Call Runit('xdens')
              Call runit('xanti')
              Call runit('xbcktrn')
              Call runit(der_integral_package)
              Call Aces_ja_init
              Call Aces_io_init(icore,1,0,.false.)
              Call Aces_io_remove(54,'DERGAM')
              Call Aces_io_fin
              Call Aces_ja_fin
           Endif
         Endif
c
#ifdef _HAVE_MPI
         call a2getrec(1,'JOBARC','LASTGEOM',1,istat)
         if (nprocs.ne.1.and.istat.ne.0) call sync_fd
#endif
         call c_gtod(is,ius)
         print '(a,f10.1,a)', 'ACES2: Total elapsed time is ',
     &                        is+1.d-6*ius-dTimeStart,' seconds'
         Print*, "Writing the dirty flag"
         call a2putrec(1,'JOBARC','DIRTYFLG',1,0)
         call runit(xjoda)
         call a2getrec(1,'JOBARC','JODADONE',1,istat)
      end do

      Return
      End

