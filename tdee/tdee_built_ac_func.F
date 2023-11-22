      Subroutine Tdee_built_ac_func(Mubar_tn,Ac_RL_R,Ac_RL_I,Work,
     +                              Memleft,Nsize,Irrepx,Iuhf,Iside,
     +                              Icount,No_steps)

      Implicit Double Precision(A-H,O-Z)

      Double Precision Mubar_tn,ML_t_MR_0

      Dimension Work(Memleft),Mubar_tn(2*Nsize)
      Dimension Ac_RL_R(No_steps),Ac_RL_I(No_steps)

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") "Entered - tdee_built_ac_func"
#endif 

      If (Iside .EQ. 1) Then
          Mubar_s_0    = 392
          Mubar_d_aa_0 = 326
          Ioffr1       = 0
          Ioffr2       = 0

      Else if (Iside .EQ. 2) Then
          Mubar_s_0    = 390
          Mubar_d_aa_0 = 316
          Ioffr1       = 0
          Ioffr2       = 0
      Endif

C Get Iside=1 Mbar_tilde(0) (see above)
C Get Iside=2 Mbar(0)

      Call Tdee_load_init_vec(Irrepx,Work,Nsize,Mubar_s_0,Ioffr1,
     +                        Mubar_d_aa_0,Ioffr2,Iuhf,.False.)

#ifdef _DEBUG_LVL0
      Call checksum("tdee_built_ac_func,M(0):",work,Nsize,s)
      Call checksum("tdee_built_ac_func,M(t):",Mubar_tn,2*Nsize,s)
#endif 

      ML_t_MR_0 = Ddot(Nsize,Mubar_tn,1,Work,1)
      Ac_RL_R(Icount) = ML_t_MR_0

      ML_t_MR_0 = Ddot(Nsize,Mubar_tn(Nsize+1),1,Work,1)
      Ac_RL_I(Icount) = ML_t_MR_0

      Write(6,"(2a,I8,2(1x,ES15.7E2))")"  The step and the autocorr.",
     +                                 " fn.:",
     +                                 Icount,Ac_RL_R(Icount),
     +                                 Ac_RL_I(Icount)

      Return 
      End
