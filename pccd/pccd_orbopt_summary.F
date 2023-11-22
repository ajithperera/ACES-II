      Subroutine Pccd_orbopt_summary(Iuhf,Ncycle)

      Implicit Double Precision(A-H,O-Z)

      Dimension Grda(6),Sgrda(6)
      Dimension Grdb(6),Sgrdb(6)
      Logical pCCD,Ccd,Lccd

      Character*4 Method

#include "flags.com"
#include "machsp.com"

      Common /Calc/Pccd,Ccd,Lccd

      Call Getrec(20,'JOBARC','SCFENEG ',Iintfp,Eref)
      Call Getrec(20,'JOBARC','TOTENERG',Iintfp,Ecor)

      Call Getrec(20,"JOBARC","GRDSTATA",6*Iintfp,Grda)
      Call Getrec(20,"JOBARC","SGRSTATA",6*Iintfp,Sgrda)
      If (Iuhf .Eq. 1) Then
         Call Getrec(20,"JOBARC","GRDSTATB",6*Iintfp,Grdb)
         Call Getrec(20,"JOBARC","SGRSTATB",6*Iintfp,Sgrdb)
      Endif 

      If (pCCD) Method = "pCCD"
      If (CCD)  Method = "CCD "
      If (LCCD) Method = "LCCD"

      write(6,*)
      Write(6,"(9x,2a,i2,a)")"Summary of the orbital optimization ",
     +                       "after ", Ncycle, " iterations"
      write(6,"(2a)") "  ---------------------------------------------",
     +               "---------------------"
      Write(6,"(a,F15.8)") "  The current reference energy: ", Eref
      Write(6,"(a,a,a,1x,F15.8)") "  The current ",method,
     +                            " energy     :", Ecor
      Write(6,*)
      Do Ispin = 1, Iuhf+1
      If (Ispin .Eq. 1) Then
      Write(6,"(3a,2(1x,E12.6))") "  The MAX and RMS alpha",
     +                            " gradients",
     +                            "        :",Grda(1),Grda(5)
      Write(6,"(2a,2(1x,E12.6))") "  The MAX and RMS scalled",
     +                            " alpha gradients:",Sgrda(1),Sgrda(5)
      Elseif (Ispin .Eq. 2) Then
      write(6,*)
      Write(6,"(3a,2(1x,E12.6))") "  The MAX and RMS beta",
     +                            "  gradients",
     +                            "        :",Grdb(1),Grdb(5)
      Write(6,"(2a,2(1x,E12.6))") "  The MAX and RMS scalled",
     +                            " beta  gradients:",Sgrdb(1),Sgrdb(5)
      Endif 
      Enddo 
      write(6,"(2a)") "  ---------------------------------------------",
     +               "---------------------"
     
      Return
      End
