      subroutine Tdee_reset_vars_4restart(Component,Irrepx,Naobfns,
     +                                    Nbfns,Ipert,Work_Label,
     +                                    D_pole,Q_pole)

      Implicit Integer (A-Z)
      Character*8 Label_D(3),Work_Label 
      Character*8 Label_Q(6)
      Logical D_pole,Q_pole

      Common /LISTDENS/LDENS

#include "syminf.com"
#include "sympop.com"
#include "flags.com"
#include "parcc.com"
#include "info.com"

      Data Label_D /'DIPOLE_X','DIPOLE_Y','DIPOLE_Z'/
      Data Label_Q /'QUAD_XX ','QUAD_YY ','QUAD_ZZ ','QUAD_XY ',
     +              'QUAD_XZ ','QUAD_YZ '/
      Data Ione /1/

      Nbfns   = Nocco(1) + Nvrto(1)
      Call Getrec(20, "JOBARC","NBASTOT ", Ione, Naobfns)

      Ldens = 160
      If (D_pole) Then
         Work_Label = Label_D(Component)
      Elseif (Q_pole) Then 
         Work_Label = Label_Q(Component)
      Endif 

      Ipert = Component 
     
      Return
      End
