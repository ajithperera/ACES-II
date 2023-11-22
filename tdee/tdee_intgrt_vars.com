#ifndef _TDEE_INTGRT_VARS_COM_
#define _TDEE_INTGRT_VARS_COM_

      Logical Sg,Rk,Lz,Lanczos
      Integer No_time_steps,component,side,Time_step_limit
      Double Precision Increment 
      Logical D_pole,Q_pole 

      Parameter(Max_internal=100000)

      Common /Tdee_intgrt_vars/Sg,Rk,Lz,T_start,T_end,No_time_steps,
     +                         Rel_error,Abs_error,Increment,
     +                         component,side,Time_step_limit
      Common /Tdee_spect_type/D_pole,Q_pole

#endif  /* __TDEE_INTGRT_VARS_COM__ */

