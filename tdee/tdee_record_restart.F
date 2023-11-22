      Subroutine Tdee_record_4restart(T0,Tn,T,Data_in_out,Nsize,Work,
     +                                Memleft,Irrepx,Iuhf,Iside,
     +                                First_block)

      Implicit Double Precision (A-H,O-Z)
      Logical First_block
      Character*80 Wrk

#include "tdee_intgrt_vars.com"

      Dimension Work(Memleft),Data_in_out(2*Nsize)

C Record the T0=Start time,Tn=Final time and the current time.
C After a restart the current time become T0 and Tn remains the
C same. In reallity we can modify Tn. That option can be added
C if required (by adding logic to read the new Tn).
C  
      If (First_block) Then
    
C Record the integration algorithm.

        If (Sg) Then
           Wrk(1:15) = "Shampine_Gordon"
           Write(10,"(15a)") Wrk(1:15)
        Else if (Rk) Then 
           Wrk(1:11) = "Runge_Kutta"
           Write(11,"(11a)") Wrk(1:11)
        Elseif (Lanczos) Then 
           Wrk(1:7)  = "Lanczos"
           Write(11,"(7a)") Wrk(1:7)
        Endif 

C Record the initial time and final time (in a.u)

        Write (10,"(E15.5,1x,E15.5)") T0, Tn

C Record the increment 

        Write (10,"(E15.5)") Increment 

C Record he relative and absolute error (change from default values)

        Write(10,"(E15.5,1X,E15.5)") Rel_error, Abs_error

        If (Iside .EQ. 1) Then
           Wrk(1:5) = "Right"
           Write(10,"(5a)") Wrk(1:5)
        Elseif (Iside .EQ. 2) Then
           Wrk(1:) = "Left"
           Write(10,"(4a)") Wrk(1:4)
        Endif 

        If (D_pole) Then
           If (Component .EQ. 1) Then
              Wrk(1:1) = "X"
              Write(10,"(1a)") Wrk(1:1)
           Elseif (Component .EQ. 2) Then
              Wrk(1:1) = "Y"
              Write(10,"(1a)") Wrk(1:1)
           Elseif (Component .EQ. 3) Then
              Wrk(1:1) = "Z"
              Write(10,"(1a)") Wrk(1:1)
           Endif 
        Else if (Q_pole) Then
           If (Component .EQ. 1) Then
              Wrk(1:2) = "XX"
              Write(10,"(2a)") Wrk(1:2)
           Elseif (Component .EQ. 2) Then
              Wrk(1:2) = "YY"
              Write(10,"(2a)") Wrk(1:2)
           Elseif (Component .EQ. 3) Then
              Wrk(1:2) = "ZZ"
              Write(10,"(2a)") Wrk(1:2)
           Elseif (Component .EQ. 4) Then
              Wrk(1:2) = "XY"
              Write(10,"(2a)") Wrk(1:2)
           Elseif (Component .EQ. 5) Then
              Wrk(1:2) = "XZ"
              Write(10,"(2a)") Wrk(1:2)
           Elseif (Component .EQ. 6) Then
              Wrk(1:2) = "YZ"
              Write(10,"(2a)") Wrk(1:2)
           Endif
         Endif 
        
C Record the symmetry of the current perturbation
    
         Write(10,"(i2)") Irrepx 

      Else

C Record the current time

        Write(10,"(E15.5)") T
        Backspace(10)

C Record the current dipole vector.
        If (D_pole) Then
        If (Iside .EQ. 1) Call Tdee_dmpgss(Data_in_out,2*Nsize,
     +                                     "RDIPOLEV")
        If (Iside .EQ. 2) Call Tdee_dmpgss(Data_in_out,2*Nsize,
     +                                     "LDIPOLEV")
        Else
        If (Iside .EQ. 1) Call Tdee_dmpgss(Data_in_out,2*Nsize,
     +                                     "RQUDPOLEV")
        If (Iside .EQ. 2) Call Tdee_dmpgss(Data_in_out,2*Nsize,
     +                                     "LQUDPOLEV")
        Endif 
 
#ifdef _DEBUG_LVL0
        If (Iside.eq.1) call checksum("current RDPV",Data_in_out,
     +                                 2*Nsize,s)
        If (Iside.eq.2) call checksum("current LDPV",Data_in_out,
     +                                 2*Nsize,s)
#endif 
      Endif 

      Return
      End 
