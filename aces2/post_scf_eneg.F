      Subroutine Post_scf_eneg
c
#include "aces2.h"    
#include "flags.h"
c     
      Character*79 szEExtrap 
      Logical parcc_exist
      szEExtrap = 'xa2proc ener_extrp'
      Call a2getrec(20,'JOBARC','NDROPGEO',1,NDRGEO)
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
      
      If (Iflags(h_IFLAGS_calc) .EQ. 39)  Then 
c
C         Call Runit('xvtran')
C         Call runit('cp HF2 HF2N')
C         Call Runit('xintprc')
C         Call runit('mv HF2N HF2')
C         Call Runit('xccsdtq < COM')

         Call Runit('xccsdtq')
c
      Else If (Dirmp2) Then
c
         Call Runit('xdirmp2')
c
      Else 
C
         If (Iflags(h_IFLAGS_excite) .NE. 2)  Call Runit('xvcc')
         If (Iflags(h_IFLAGS_calclevel) .EQ. 50 .OR.
     &       Iflags(h_IFLAGS_calclevel) .EQ. 51) then
             inquire(file="parcc",exist=parcc_exist)
             if (parcc_exist) Call Runit("rm parcc")
             Call Runit('xvcc') 
         Endif 
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
c EOM-CC excitation energies
C Do EOM-CC exciation energie with RPA/DRPA refeences. In principal
C MBPT, P-EOM must come here.
c
      If (Iflags(h_IFLAGS_excite) .GT. 2) Then
         Call Runit("xlambda")     
         Call Runit("xvee") 
      Else if (Iflags(h_IFLAGS_excite) .EQ. 2 .AND.
     &         Iflags2(h_IFLAGS2_eomref) .EQ. 0) Then
         Call Runit("xrpa") 
      Else if (Iflags(h_IFLAGS_excite) .EQ. 2 .AND. 
     &         Iflags2(h_IFLAGS2_eomref) .GE. 4) Then
         Call Runit("xrpa") 
         Call Runit("xvcc")
         Call Runit("xlambda")     
         Call Runit("xvee") 
      Endif 
C
c EA-EOM calculations with "xvea". Can also be done with "xmrcc"
c
      If (Iflags2(h_IFLAGS2_ea_calc) .EQ. 5) Then
         Call Runit("xlambda")
         Call Runit("xvea")
      Endif
C
C IP-EOM calculations with "xvip". Can also be done with "xmrcc"
c
      If (Iflags2(h_IFLAGS2_ip_calc) .EQ. 5) Then
         Call Runit("xlambda")
         Call Runit("xvip")
      Endif

      Return
      End
