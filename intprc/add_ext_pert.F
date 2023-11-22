#include "flags.h"

      Subroutine  Add_ext_pert(Work, Maxcor, X_field, Y_field, Z_field,
     &                         F, Nbas, Nbasx, Ispin)

      Implicit Double Precision (A-H, O-Z)
 
      Dimension Work(Maxcor), F(Nbas, Nbas)
      Logical X_field, Y_field, Z_field 
      Character*8 Label(3), DIRECTION

#include "machsp.com"
#include "jodaflags.com"

      Data Label /'DIPOLE_X','DIPOLE_Y','DIPOLE_Z'/
C      
      Length = Nbas*(Nbas+1)/2
C   
      IXYZt_efld = 1
      IXYZf_efld = IXYZt_efld + Length
      IAO_efld   = IXYZf_efld + Length
      IMO_efld   = IAO_efld   + Nbas*Nbas
      IEVEC      = IMO_efld   + Nbas*Nbas
      IScr       = IEVEC      + Nbas*Nbasx
      Inext      = Iscr       + 2*Nbax*Nbasx
C
      Call ZERO(Work(IXYZf_efld), Length)

      IF (X_field) Then
         IFLD_STREN  = Iflags(h_IOPPAR_xfield)
         DXFLD_STREN = DFLOAT(IFLD_STREN)*1.0D-06
         DIRECTION   = '     X  '
         WRITE(6,100) DIRECTION(6:6) , DXFLD_STREN

         Call Getrec(20, "JOBARC", Label(1), Length*Iintfp,
     &                         Work(IXYZt_efld))
         CALL DSCAL(Length, DXFLD_STREN, Work(IXYZt_efld), 1)
         Call Dcopy(Length, Work(IXYZt_efld), 1,  Work(IXYZf_efld), 1)
      Endif

      IF (Y_field) Then
         IFLD_STREN = Iflags(h_IOPPAR_xfield)
         DXFLD_sTREN = DFLOAT(IFLD_STREN)*1.0D-06
         DIRECTION   = '     Y  '
         WRITE(6,100) DIRECTION(6:6), DYFLD_STREN

         Call Getrec(20, "JOBARC", Label(2), Length*Iintfp,
     &               Work(IXYZt_efld))
         Call Daxpy(Length, DYFLD_sTREN, Work(IXYZt_efld), 1, 
     &              Work(IXYZf_efld), 1)

      Endif 

      IF (Z_field) Then
          IFLD_STREN = Iflags(h_IOPPAR_zfield)
          DZFLD_STREN = DFLOAT(IFLD_STREN)*1.0D-06
          DIRECTION   = '     Z  '
          WRITE(6,100) DIRECTION(6:6), DZFLD_STREN

         Call Getrec(20, "JOBARC", Label(3), Length*Iintfp,
     &               Work(IXYZt_efld))
         Call Daxpy(Length, DZFLD_sTREN, Work(IXYZt_efld), 1, 
     &              Work(IXYZf_efld), 1)

      Endif 

      Call Expnd2(Work(IXYZf_efld), Work(IAO_efld), Nbas)
#ifdef _DEBUG_LVLM
      write(6,"(a)") "The scaled AO basis perturbed vectors"
      call output(Work(IAO_efld), 1, Nbas, 1, Nbas, Nbas, Nbas, 1)
#endif
      Call AO2MO2(Work(IAO_efld), Work(IMO_efld), Work(Ievec),
     &                 Work(Iscr), Nbas, Nbasx, Ispin)

#ifdef _DEBUG_LVLM
      write(6,"(a)") "The scaled MO basis trans. perturbed vectors"
      call output(Work(IMO_efld), 1, Nbas, 1, Nbas, Nbas, Nbas, 1)
#endif 
      Call Daxpy(Nbas*Nbas, 1.0D0, Work(IMO_efld), 1, F, 1)

#ifdef _DEBUG_LVLM
      write(6,"(a,i2)") "The Fock matrix + perturbation and Ispin",
     &                    Ispin
      call output(Work(IMO_efld), 1, Nbas, 1, Nbas, Nbas, Nbas, 1)
#endif 

100   FORMAT(T3,'@FINFLD-I, Adding electric field perturbation to ',
     &          'Fock matrix.',/,T10,'Field direction',
     &          T30,': ',A1,/,T10,'Field strength',T30,': ',F10.7,
     &          ' a.u.')

       Return
       End 
