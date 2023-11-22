      Subroutine Tdee_time_propagate_SG(Scr,Data_in_out,Ac_RL_R,
     +                                  Ac_RL_I,Work,Memleft,
     +                                  Iuhf,Irrepx,Iside,Nsize,
     +                                  Restart)

      Implicit Double Precision(A-H,O-Z)

      Dimension Work(Memleft),Iwork(5),Scr(42*Nsize+100)
      Dimension Data_in_out(2*Nsize)
      External Tdee_pes_ode
      Integer t 
      Character*1 Str
      Logical Restart,First_block, Propgt_done

#include "tdee_intgrt_vars.com"

      Dimension Ac_RL_R(No_time_steps),Ac_RL_I(No_time_steps)

      Data Izero, Zero, Two /0, 0.0D0, 2.0D0/

      T0          = T_start
      Tn          = T_end 
      No_steps    = No_time_steps
      Relerr      = Rel_error 
      Abserr      = Abs_error 
      Step_size   = Increment 
      Iflag       = 1
      Icount      = 0
      First_block = .True. 

      If (Time_step_limit .Gt. Max_internal) Then
         Write(6,"(a,a,i6)") " The number of time-steps is more than",
     +                       " the currently allowed maximum of ",
     +                         Max_internal 
         Call Errex
      Endif 

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") "Entered - tdee_time_propagate"
      Write(6,"(a)") "The settings for the time-integrator"
      Write(6,"(a)") "SG,RK,T0,Tn,Relerr,Abserr,No_steps,Step_size"
      Write(6,"(1x,2l,2(1x,ES8.2E2),2(1x,ES8.2E2),1x,i8,1x,F7.4)")
     +                              SG,RK,T_start,T_end,
     +         Relerr,Abserr,No_steps,step_size
#endif

      Call Tdee_record_4restart(T0,Tn,T0,Data_in_out,Nsize,Work,
     +                          Memleft,Irrepx,Iuhf,Iside,
     +                          First_block)
      First_block = .False.

      Mubar_s_pq_t0 = 390
      Mubar_d_pq_t0 = 313
      Mubar_s_pq_tn = 394
      Mubar_d_pq_tn = 333

      If (Restart) Then

C Read the current starting vector from the restart file. This vector 
C correspond to the last time step that was successful. 

         If (D_Pole) TheN
         If (Iside .EQ. 1) Call Tdee_redgss(Data_in_out,2*Nsize,
     +                               "RDIPOLEV")
         If (Iside .EQ. 2) Call Tdee_redgss(Data_in_out,2*Nsize,
     +                               "LDIPOLEV")
         Elseif (Q_Pole) Then
         If (Iside .EQ. 1) Call Tdee_redgss(Data_in_out,2*Nsize,
     +                               "RQUDPOLEV")
         If (Iside .EQ. 2) Call Tdee_redgss(Data_in_out,2*Nsize,
     +                               "LQUDPOLEV")
         Endif 
#ifdef _DEBUG_LVL0
        If (Iside.eq.1) call checksum("current RDPV",Data_in_out,
     +                                 2*Nsize,s)
        If (Iside.eq.2) call checksum("current LDPV",Data_in_out,
     +                                 2*Nsize,s)
#endif
      Else

C Load the Mubar_R at t=0

         If (Iside .EQ. 1) Then
           Ioffr1 = 0
           Ioffr2 = 0
           Call Tdee_Load_vec(Irrepx,Data_in_out,Nsize,
     +                        Mubar_s_pq_t0,Ioffr1,Mubar_d_pq_t0,
     +                        Ioffr2,Iuhf,.False.)
           Call Dzero(Data_in_out(NSize+1),Nsize)
           Call Tdee_normaliz_init_vec(Data_in_out,2*Nsize)

         Elseif (Iside .EQ. 2) Then

           Ioffr1 = 2
           Ioffr2 = 10
           Call Tdee_Load_vec(Irrepx,Data_in_out,Nsize,
     +                        Mubar_s_pq_t0,Ioffr1,Mubar_d_pq_t0,
     +                        Ioffr2,Iuhf,.False.)
           Call Dzero(Data_in_out(NSize+1),Nsize)
           Call Tdee_normaliz_init_vec(Data_in_out,2*Nsize)
         Endif
      Endif 
        
      Do t = 1, No_steps 

         T_step = T0 + Sign(Step_size,Tn-T0)
         Icount = Icount + 1

C Check the boundness of the step

         If ((T_step-Tn)*(T_step-T0) .GT. 0.0D0) T_Step = Tn

#ifdef _DEBUG_LVL0
      Write(6,"(a,3(1x,2F10.5))") " The initial and final-time", T0,
     +                            T_step,Tn
#endif   
         If (T_step .LT. Tn) Then
            Call Ode(Tdee_pes_ode,2*Nsize,Data_in_out,T0,T_step,
     +               Relerr,Abserr,Iflag,Scr,Iwork)
         Else
             Write(6,*) 
             Write(6,"(a,a)") " The integration has reached the final",
     +                        " time"
             Go to 100
         Endif

#ifdef _DEBUG_LVL0
      Write(6,"(a,1x,F10.5)") "Time integration is completed", T0
      Write(6,"(a,1x,i2)") "Success??", Iflag 
#endif
       If (Iflag .NE. 2) Then
          Write(6,*)
          Write(Str,"(i1)") Iflag
          Write(6,"(a,a,a,2a)")" Error occured during integration step",
     +                         ". Use the value of the flag ",Str,
     +                         " to see"
          Write(6,"(a)") " whether integration can be restarted."
          Write(6,*)
          Call Errex
       Endif

       If (Iside .EQ. 1) Then

C Built <M(0)|M(t_n+1)> and store in a array dimensioned to number
C of time steps.
        
            Call Tdee_record_4restart(T0,Tn,T0,Data_in_out,Nsize,
     +                                Work,Memleft,Irrepx,Iuhf,
     +                                Iside,First_block)

            Call Tdee_built_ac_func(Data_in_out,Ac_RL_R,Ac_RL_I,
     +                              Work,Memleft,Nsize,Irrepx,
     +                              Iuhf,Iside,Icount,No_steps)
        
         Elseif (Iside .EQ. 2) Then

C Built <M(t_n+1)|M(0) > and store in a array dimensioned to number
C of time steps

            Call Tdee_record_4restart(T0,Tn,T,Data_in_out,Nsize,
     +                                Work,Memleft,Irrepx,Iuhf,
     +                                Iside,First_block)

            Call Tdee_built_ac_func(Data_in_out,Ac_RL_R,Ac_RL_I,
     +                              Work,Memleft,Nsize,Irrepx,
     +                              Iuhf,Iside,Icount,No_steps)
         Endif 

      Enddo 

 100  Continue 

#ifdef _DEBUG_LVL0
      If (Icount .ne. 0) Then
      If (iside .eq. 1) Then 
      Write(6,"(a,a)") " The real part of the right auto-correlation",
     +               " function."
      Write(6,"(6(1xF12.6))") (Ac_RL_R(i),i=1,Icount-1)
      Write(6,"(a,a)") " The imaginary part of the right auto-",
     +               "correlation function."
      Write(6,"(6(1xF12.6))") (Ac_RL_I(i),i=1,Icount-1)
      Else 
      Write(6,"(a,a)") " The real part of the left auto-correlation",
     +               " function."
      Write(6,"(6(1xF12.6))") (Ac_RL_R(i),i=1,Icount-1)
      Write(6,"(a,a)") " The imaginary part of the left auto",
     +               " function."
      Write(6,"(6(1xF12.6))") (Ac_RL_I(i),i=1,Icount-1)
      endif 
      endif
#endif 

      Return 
      End
