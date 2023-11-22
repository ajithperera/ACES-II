      Subroutine Tdee_read_tintgrt_vars(work,Maxcor,Nodiag_component)

      Implicit Double Precision (A-H,O-Z)

      Logical Zmat_Present
      Character*80 Wrk
      Integer Nondiag_component 
      Dimension Work(Maxcor)

#include "tdee_intgrt_vars.com"

      Inquire(File="ZMAT",Exist=Zmat_present)
      
      If (Zmat_present) Then
         Open(Unit=4,File="ZMAT",Form="FORMATTED",Status="OLD")
         Rewind(4)
      Else 
        Write(6,"(1x,a,a)") "@-tdee_read_tintgrt_vars; ZMAT file not",
     +                     " present"
        Call Errex
      Endif 

 300  Read(4,"(a)",End=900) Wrk
      If (Wrk(1:12) .NE. "*INTGRT_VARS") Goto 300 

C Read the integration algorithm. 

      Read(4,"(a)",End=900) Wrk
      SG = .FALSE.
      RK = .FALSE.
      LZ = .FALSE.

      If (Wrk(1:15) .EQ. "Shampine_Gordon") Sg =.TRUE.
      If (Wrk(1:11) .EQ. "Runge_Kutta")     Rk =.TRUE.
      If (Wrk(1:7)  .EQ. "Lanczos")         LZ =.TRUE.
   
C The initial time and final time (in a.u)

      Read(4,"(ES6.2E2,1x,ES6.2E2)",End=900) T_start, T_end

C The increment (or the size of the time step)

      Read(4,"(ES7.2E2)",End=900) Increment 

C The relative and absolute error (change from default values)

      Read(4,"(ES7.2E2,1X,ES7.2E2)",End=900) Rel_error, Abs_error 

      No_time_steps = Min(Max_internal,2*Int((T_end-T_start+Increment)
     +                   /Increment))

      Time_step_limit = (T_end-T_start)/Increment 

      Read(4,"(a)",End=900) Wrk

      If (Wrk(1:5) .EQ. "Right") Side  = 1
      If (Wrk(1:4) .EQ. "Left" ) Side  = 2

      Read(4,"(a)",End=900) Wrk

      D_pole = .False.
      Q_pole = .False.

      If (Wrk(1:6)  .EQ. "Dipole")     D_pole = .True.
      If (Wrk(1:10) .EQ. "Quadrupole") Q_pole = .True.

      Read(4,"(a)",End=900) Wrk

      If (D_pole) Then
         Nondiag_component = 0
         If (Wrk(1:1) .EQ. "X")  Component  = 1
         If (Wrk(1:1) .EQ. "Y")  Component  = 2
         If (Wrk(1:1) .EQ. "Z")  Component  = 3
         Read(4,"(a)",End=900) Wrk
         If (Wrk(1:6) .Eq. "Nodiag") Then
            Read(4,"(a)",End=900) Wrk
            If (Wrk(1:1) .EQ. "X")  Nondiag_component  = 1
            If (Wrk(1:1) .EQ. "Y")  Nondiag_component  = 2
            If (Wrk(1:1) .EQ. "Z")  Nondiag_component  = 3
         Endif 
      Elseif (Q_pole) Then 
         Nondiag_component = 0
         If (Wrk(1:2) .EQ. "XX") Component  = 1
         If (Wrk(1:2) .EQ. "YY") Component  = 2
         If (Wrk(1:2) .EQ. "ZZ") Component  = 3
         If (Wrk(1:2) .EQ. "XY") Component  = 4
         If (Wrk(1:2) .EQ. "XZ") Component  = 5
         If (Wrk(1:2) .EQ. "YZ") Component  = 6
         Read(4,"(a)",End=900) Wrk
         If (Wrk(1:6) .Eq. "Nodiag") Then
            Read(4,"(a)",End=900) Wrk
            If (Wrk(1:2) .EQ. "XX") Nodiag_component  = 1
            If (Wrk(1:2) .EQ. "YY") Nodiag_component  = 2
            If (Wrk(1:2) .EQ. "ZZ") Nodiag_component  = 3
            If (Wrk(1:2) .EQ. "XY") Nodiag_component  = 4
            If (Wrk(1:2) .EQ. "XZ") Nodiag_component  = 5
            If (Wrk(1:2) .EQ. "YZ") Nodiag_component  = 6
         Endif 
      Endif 
      Go to 999 

 900  Write(6,901) 
 901  Format(T3,"@-tdee_read_tintgrt_vars,INTGRT_VARS namelist not found',
     + or incomplete.")
      Write(6,"(1x,a)") " The default values are being used"
      
C Set the default values.

      SG        = .True.
      T_start   = 0.0E00
      T_end     = 1.0E03
      Increment = 5.0E-02
      Rel_error = 1.0D-06
      Abs_error = 1.0E-06
      No_time_steps = Min(Max_internal,
     +                2*Int((T_end-T_start+Increment)/Increment))

      Time_step_limit = (T_end-T_start)/Increment 

 999  Continue 

      Close(Unit=4,Status="KEEP")
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") " Entered - tdee_read_tintgrt_vars"
      Write(6,"(a)") " The settings for the time-integrator"
      Write(6,"(a,a)") " SG,RK,T_start,T_end,Rel_error,Abs_error,",
     +                 "increment,Max_time_steps,side,component"
      Write(6,"(1x,3l,2(1x,ES8.2E2),2(1x,ES8.2E2),1x,i8,1x,F5.2,
     +                                           1x,i1,1x,i1)") 
     +                              SG,RK,LZ,T_start,T_end,
     +         Rel_error,Abs_error,No_time_steps,increment,
     +         side,component 
      Write(6,*)
#endif
      
      Return
      End
