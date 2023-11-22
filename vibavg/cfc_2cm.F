      Subroutine Cfc_2cm(Cfc,Omega,Nmodes,Nreals,Icol)

      Implicit Double Precision(A-H,O-Z)

      Dimension Cfc(Nmodes,Nmodes,Nmodes) 
      Dimension Omega(3*Nreals)

      Parameter(Hbar=1.054572266D-34)
      Parameter(au2cm=219474.63067D0)
      Parameter(C=2.99792458D+08)
      Parameter(B2ang=0.529177249D0)
      Parameter(Amu=1.6605402D-27)
      Parameter(Avgd=6.0221367D+23)
      Parameter(Thres1=1.0D-03,Thres2=1.0D-05)

      Data Two/2.0D0/

C Factor to convert to cubic force constants to cm-1.
C The naked units of cubic force constants are  E_h Kg^(-3/2)m^(-3)
C The dimensionless coordinates relates to mass-wegithed coordinate
C by  Fred* Omega*Q =  q
C
C The steps are taken in Q = (1/(Fred*Omega))*q. The cubic Force 
C constants in dimensionless normal coordinates is related to the
C standard coordinates by the relation F_ijk Q_iQ_jQ_k = f_ijkq_i
C q_jq_k.

      Pi          = Dacos(-1.0D0)
      Bohr2m      = 0.529177249*1.0D-10
      Factor_cfc1 = hbar/(Two*Pi*c*1.0D02)
      Factor_cfc2 = (Au2cm) * (1.0D0/Amu)**1.5 * (1.0D0/Bohr2m)**3

      Factor_cfc = Factor_cfc1**(1.5)*Factor_cfc2

#ifdef _DEBUG_LVL0
      Write(6,*) 
      Write(6,"(2a,ES24.15E2)") " The conversion factor for cubic",
     +                          " force constant: ", Factor_cfc
#endif
      Do I= 1, Nmodes 
         Do J = 1, Nmodes
           Do K = 1, Nmodes 
              If (Dabs(Min(Omega(K+Icol),Omega(J+Icol),
     +                     Omega(I+Icol))) .Gt. Thres1) Then
                  Frq = (Omega(K+Icol)*Omega(J+Icol)*Omega(I+Icol))
                  Z = Factor_cfc/Dsqrt(Frq)
                  If (Dabs(Cfc(K,J,I)) .Gt. Thres2) Then
                     Cfc(K,J,I) = Cfc(K,J,I)*Z   
                  Endif
              Endif 
           Enddo 
         Enddo 
      Enddo

      Return
      End 
