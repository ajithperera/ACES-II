      Subroutine Stpsize(Omega, Deltaq, Delx, Nvibs, B2Ang)
C 
      Implicit Double Precision (A-H, O-Z)
C
#include "machsp.com"
C
      Double Precision Meters_2Angs 
      Logical Step_length_present
C
      Dimension Deltaq(Nvibs), Omega(Nvibs)
  
      Parameter(H=6.62608D-34, C=2.99792458D+10, Amu=1.6605402D-27)
      Parameter(Meters_2Angs=1.0D-10)

      Data Two,Done /2.00D0,1.0D0/

C This is dimensionless reduced step-size, defined by Sqrt(omega/hbar)
C JCP, 121, 1383, 2004. The a.m.u. is 1.660538+D-27Kg. The following 
C is obtained from the first author of the cited paper. I can not 
C reconcile this with the Eqn. given in the paper. 
C
C 1/(c*Omega/h) = 1/(L*T^-1 * L^-1/M^2LT^-2 * T)  (omega is in wave numbers)
C               = 1/(ML^2T^-1/T^-1)
C               = ML^2
C 
C Sqrt(h/c*Omega) has units of Sqrt(M)*L 
C There is an external factor Sqrt(2*pi/(1/2*pi)) = 2*pi
C (Fred: Factor for reduced coordinates)

      Pi    = Acos(-1.0D+00)
      Fred  = Two*Pi*Dsqrt((c/h)*(Amu))*Meters_2angs*B2Ang

#ifdef _DEBUG_LVL0
      Write(6,"(a,F15.9)") " The factor for reduced coordinates: ", 
     +                        Fred
#endif 

      Do Imode = 1, Nvibs

         C1  = Fred*Dsqrt(Omega(Imode))
         Deltaq(Imode) = Delx/C1

      Enddo 
C
C Record the dimensional stepsize into a file "STEP_LENGTH". This is
C usefull to have in a file for when property deriavtives are 
C computed.

      Inquire(File="STEP_LENGT",Exist=Step_length_present)
      If (Step_length_present) Then
         Open(Unit=10,File="STEP_LENGT",Form="Formatted",
     +        Status="Old")
      Else 
         Open(Unit=10,File="STEP_LENGT",Form="Formatted",
     +        Status="New")
      Endif 
      Write(10,"(a)") "Dimensionless curveture independent step sizes"
      Write(10,"(5(1x,F15.9))") (Deltaq(I),I=1,Nvibs)
      Close(10)
C
      Return
      End



