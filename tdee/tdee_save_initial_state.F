      Subroutine Tdee_Save_initial_state(Work,Memleft,Irrepx,Iuhf,
     +                                   Iside,Nsize)

      Implicit Double Precision(A-H,O-Z)

      Dimension Work(Memleft)
      Logical Rhf

#include "machsp.com"
#include "syminf.com"
#include "sympop.com"
#include "flags.com"
#include "info.com"
#include "symloc.com"
#include "sym.com"
     
C Singles and doubles right and dipole functions list are saved. 

#ifdef _DEBUG_LVL0
        Write(6,*)
        Write(6,"(a)") "---Entered Tdee_Save_initial_state ---"
        Write(6,*)
#endif
      I000 = 1
      I010 = I000 + NSIZE 
      I020 = I010 + NSIZE

      Ioffr1 = 0
      Ioffr2 = 0
      Ioffsp = 2
      Mu_s   = 390
      Mu_d_o = 313
      Mu_d_d = 316

      Call Tdee_load_vec(Irrepx,Work(I000),Memleft,Mu_s,Ioffr1,Mu_d_o,
     +                   Ioffr2,Iuhf,.False.)
      Call Tdee_normaliz_init_vec(Work(I000),Nsize)

      Ioffr1  = 0
      Ioffr2  = 0
      Ioffsp  = 2
      Mut_s   = 392
      Mut_d_o = 323
      Mut_d_d = 326

      Call Tdee_load_vec(Irrepx,Work(I010),Memleft,Mut_s,Ioffr1,Mut_d_o,
     +                  Ioffr2,Iuhf,.False.)
   
      RL_NORM = Ddot(Nsize,Work(I000),1,Work(I010),1)
     
      Call Dscal(Nsize,1.0D0/RL_NORM,Work(I010),1)

C#ifdef _DEBUG_LVL0
      RL_NORM = Ddot(Nsize,Work(I000),1,Work(I010),1)
      Write(6,"(a,F15.7)") "@-Tdee_Save_initial_state,RL_NORM: ",
     +                      RL_NORM
C#endif 
      Call Tdee_dump_vec(Irrepx,Work(I000),Memleft,Mu_s,Ioffr1,Ioffsp,
     +                   Mu_d_d,Ioffr,Iuhf,.False.)
      Call Tdee_dump_vec(Irrepx,Work(I010),Memleft,Mut_s,Ioffr1,Ioffsp,
     +                   Mut_d_d,Ioffr,Iuhf,.False.)

#ifdef _NO_OBSOLETE
      If (Iside .EQ. 1) Then

         Ioffr1 = 0
         Ioffr2 = 0
         Ioffsp = 2
         Mu_s   = 390
         Mu_d_o = 313
         Mu_d_d = 316
       
         Call Tdee_load_vec(Irrepx,Work,Memleft,Mu_s,Ioffr1,Mu_d_o,
     +                      Ioffr2,Iuhf,.False.)
         Call Tdee_normaliz_init_vec(Work,Nsize)
         Call Tdee_dump_vec(Irrepx,Work,Memleft,Mu_s,Ioffr1,Ioffsp,
     +                      Mu_d_d,Ioffr,Iuhf,.False.)
  
      Elseif (Iside .EQ. 2) Then

         Ioffr1  = 0
         Ioffr2  = 0
         Ioffsp  = 2
         Mut_s   = 392 
         Mut_d_o = 323
         Mut_d_d = 326

         Call Tdee_load_vec(Irrepx,Work,Memleft,Mut_s,Ioffr1,Mut_d_o,
     +                      Ioffr2,Iuhf,.False.)
         Call Tdee_normaliz_init_vec(Work,Nsize)
         Call Tdee_dump_vec(Irrepx,Work,Memleft,Mut_s,Ioffr1,Ioffsp,
     +                      Mut_d_d,Ioffr,Iuhf,.False.)

      Endif 
#endif 

#ifdef _DEBUG_LVL0
        Write(6,*)
        Write(6,"(a)") "---Exit Tdee_Save_initial_state ---"
        Write(6,*)
#endif

      Return 
      End
