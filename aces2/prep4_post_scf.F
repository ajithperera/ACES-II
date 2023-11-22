#include "flags.h"
c
      Subroutine Prep4_post_scf
      implicit none
c
      integer ndrgeo, itest,ndropa
CSSS      logical Dropmo_analy_grads 
c
#include "flags.h"
#include "aces2.h"
c
c   o vtran/intprc were already run if HFSTAB=(ON|FOLLOW)
c     we have to re-run them if (HFSTAB=FOLLOW && 'SCFKICK'!=0) || DROPMO!=0
c
      if (iflags(h_IFLAGS_hfstab).eq.2) then
         call a2getrec(-1,'JOBARC','SCFKICK',1,itest)
      else
         itest=0
      end if
c
      Call A2getrec(20, 'JOBARC', 'NDROPGEO', 1, NDRGEO)
c
c   o use NUMDROPA since ndrgeo is 0 for single-point calcs.
      Call A2getrec(-1, 'JOBARC', 'NUMDROPA', 1, ndropa)
CSSS      dropmo_analy_grads = (Analytical_gradient .AND. 
CSSS     &                      Ndropa .ne. 0)
c
CSSS      Print*, ndropa, itest, iflags(h_IFLAGS_hfstab) 
      If (iflags(h_IFLAGS_hfstab).eq.0.or.ndropa.ne.0.or.
     &   (iflags(h_IFLAGS_hfstab).eq.2.and.itest.ne.0)) Then
c
         if (ndrgeo.ne.0)  then
            call runit('cp JOBARC JOBARC_AM')
            call runit('cp JAINDX JAINDX_AM')
         end if
c
         call runit('xvtran')
         call runit('xintprc')
c
C This block is commented to accomodate Monika's CCSDTQ programs.
c 05/21 Ajith Perera.
C
C        if (iflags2(h_IFLAGS2_kucharski).eq.1) then
C         o KJW 9/29/99 for Stan Kucharski's ccsdtq program
C            call runit('cp HF2 HF2N')
C            call runit('xintprc')
C            call runit('mv HF2N HF2')
C         else
C
C            call runit('xintprc')
C         end if

      Else 
C Commented on 06/09 by Ajith Perera in order to prevent problems
C arsing from rerunning these two MEs for the post-HF calcs with 
C HFSTABILITY=(ON|FOLLOW) but the analysis found no instabilities
C that can be followed. 
C Above comment is useless because I don't identify what the problem
C was and now I do not remember (could be symmetry is not lowered
C C1. Ajith Perera, 08/2012.
C
C Well, The problem is identified. If I let these two to rerun
C I get wrong correlation energies for the cases where user used
C HFSTAB=ON,CALC=POSTSCF, but found no instabilities. This was 
C notified to me once again by Brazilian Visting student 
C 11/06/2019. 

       If (Itest .Ne. 0 .OR. iflags(h_IFLAGS_hfstab).eq.2) Then
          Call runit('xvtran')
          Call runit('xintprc')
       Endif 

      Endif
c 
      Return
      End

