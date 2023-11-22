      Subroutine Tdee_read_from_restart(work,Maxcor,Irrepx,Done)

      Implicit Double Precision (A-H,O-Z)

      Dimension Work(Maxcor)
      Character*80 Wrk
      Logical Isopen,Done 

#include "tdee_intgrt_vars.com"

      Write(6,"(a,a)") "  This is a restart and initial data is read",
     +                 " from the restart file"
      Inquire(File="Restart_log",Opened=Isopen,Number=Iunit)
      If (Isopen) Then
         Rewind(10)
      Else
         Open(Unit=10,File="Restart_log",Form="formatted",
     +        Status="Old")
      Endif 

C Read the integration algorithm. 

      Read(10,"(a)",End=900) Wrk
      SG = .FALSE.
      RK = .FALSE.
      LZ = .FALSE.

      If (Wrk(1:15) .EQ. "Shampine_Gordon") Sg =.TRUE.
      If (Wrk(1:11) .EQ. "Runge_Kutta")     Rk =.TRUE.
      If (Wrk(1:7) .EQ. "Lanczos")          LZ =.TRUE.
   
C The initial time and final time (in a.u)

      Read(10,"(E15.5,1x,E15.5)",End=900) T_start, T_end

C The time increment 

      Read(10,"(E15.5)",End=900) Increment 

C The relative and absolute error (change from default values)

      Read(10,"(E15.5,1X,E15.5)",End=900) Rel_error, Abs_error 

C Left of right multiplication 

      Read(10,"(a)",End=900) Wrk

      If (Wrk(1:5) .EQ. "Right") Side  = 1
      If (Wrk(1:4) .EQ. "Left" ) Side  = 2

C Cartesian component of the Dipole vecor.

      Read(10,"(a)",End=900) Wrk

      If (Wrk(1:1) .EQ. "X") Component  = 1
      If (Wrk(1:1) .EQ. "Y") Component  = 2
      If (Wrk(1:1) .EQ. "Z") Component  = 3

C Read the symmetry of the current perturbation
 
      Read(10,"(I2)",End=900) Irrepx 

C Read last time recorded and compute the new start time.

      Read(10,"(E15.5)",End=900) T

C The increment (or the size of the time step)

CSSS      T_start = T + Increment 
      T_start = T 

      T_check = T + Increment

      If (T_check .eq. T_end) Then
          Write(6,*)
          Write(6,"(a,a)") "  Initial and final times are rqual. No",
     +                     " further time propagation is not needed!"
          Done = .True.
          Return 
      Else
          No_time_steps = 2*Int((T_end-T_start+Increment)/Increment)
          Time_step_limit = (T_end-T_start)/Increment
           Write(6,*)
           Write(6,"(a)") "  ---Restarting the time propgation---"
           Write(6,"(a,1x,E15.5)") "  The inital time: ", T_start
           Write(6,"(a,1x,E15.5)") "  The final  time: ", T_end 
      Endif 

      Go To 999

 900  Write(6,901) 
      Write(6,*)

 901  Format(T3,"@-tdee_read_from_restart, Read error occured",
     +          " while reading the",/,"  Restart_info file")
      Close(Unit=10,Status="KEEP")
      Write(6,*) 
      Call Errex 

 999  Continue 

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") " @-tdee_read_from_restart; Settings at exit!"
      Write(6,"(a)") " The settings for the time-integrator"
      Write(6,"(a,a)") " SG,RK,T_start,T_end,Rel_error,Abs_error",
     +                 " increment,Max_time_steps,side,component"
      Write(6,"(1x,3l,2(1x,ES8.2E2),2(1x,ES8.2E2),1x,i8,F5.2,
     +                                           1x,i1,1x,i1)") 
     +                              SG,RK,LZ,T_start,T_end,
     +         Rel_error,Abs_error,No_time_steps,increment,
     +         side,component 
      Write(6,*)
#endif
      
      Return
      End
