      Subroutine Post_scf_eneg_4brueck(Rot_grad)
c
#include "aces2.h"    
#include "flags.h"
c     
      Character*79 szEExtrap 
      Logical Rot_grad
      szEExtrap = 'xa2proc ener_extrp'
      Call a2getrec(20,'JOBARC','NDROPGEO',1,NDRGEO)
C
C The regular post_scr_energy routine has excited states logic. But
C we are now trying to develop Breuckner based excited state
C calculations. So, it become necessary first determine Brueckner
C and then proceed to exciated energies, 04/2015, Ajith Perera.
c
c This Saving of All MO JOBARC/JAINDX as JOBARC_AM/JAINDX_AM is only relevent
c for DROPMO gradients. It is done here even for energy only as a organizational
c convenince.
c
CSSS      If (Ndrgeo.ne.0)  then
CSSS         Call Runit('cp JOBARC JOBARC_AM')
CSSS         Call Runit('cp JAINDX JAINDX_AM')
CSSS      End if
C
c The Musial & Kucharski need to have HF2 file around which normaly
c get deleted after processing. I sure hope we can change this key-word!
c 
      
      If (Iflags2(h_IFLAGS2_kucharski) .EQ. 1)  Then 
c
         Call Runit('xvtran')
         Call runit('cp HF2 HF2N')
         Call Runit('xintprc')
         Call runit('mv HF2N HF2')
         Call Runit('xccsdtq < COM')
c
      Else If (Dirmp2) Then
c
         Call Runit('xdirmp2')
c
C Finite-order Brueckner work for COPT, 11/2016 (Ajith Perera for
C Y. Park). 
C
      Else if (Iflags2(h_IFLAGS2_copt) .EQ. 1) then
 
         Call Runit('xcopt')

      Else
C
         If (Iflags(h_IFLAGS_excite) .NE. 2 .AND. .NOT. Rot_grad) 
     &       Call Runit('xvcc')
         If (Iflags(h_IFLAGS_excite) .NE. 2 .AND. Rot_grad) 
     &       Call Runit('xpccd')
c
c Noniterative 5th-order triples and quadruples on top of CCSD
c
         If (Iflags(h_IFLAGS_calclevel) .EQ. 12 .OR.
     &      (Iflags(h_IFLAGS_calclevel) .GE. 26 .AND. 
     &       Iflags(h_IFLAGS_calclevel) .LE. 30)) Then
             Call Runit('xvcc5t')
             Call Runit('xvcc5q')
         Endif
C
         If (BlCCSDT) Then
           Call Runit('xlambda')
           Call Runit('xlcct')
         Endif
c
         If (Iflags(h_IFLAGS_calclevel) .EQ. 31) Call Runit('xvcc5t')
c
c Extrapolate energy (only for post-SCF methods), 01/2006, Ajith Perera
c
         If (bEextrap .OR. bCExtrap) Then
            Call Runit(szEExtrap)
         Endif
c
       Endif
c 
      Return
      End
