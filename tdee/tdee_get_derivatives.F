      Subroutine Tdee_get_derivatives(Irrepx,Iside,Iuhf,Memleft,Nsize,
     +                                Type,Dummy1,Dummy2)

      Implicit Double Precision (A-H,O-Z)
      Logical Type 
      Dimension Dummy1(Nsize),Dummy2(Nsize)
      Double Precision Mone

#include "icore.com"
#include "istart.com"
#include "machsp.com"

      Data Mone /-1.0D0/
    
C At t=0 (the initial state mubar(0) is already on the correct lists.
C Multiplication is entirly stand-alone and works from the data
C from the following lists.
C UHF   Right Left
C List: 390   392 AA(1),BB(2)
C List: 314   324 AAAA
C List: 315   325 BBBB
C List: 316   326 ABAB
C RHF
C List: 390   392 AA(1)
C List: 316   326 ABAB
C
C First do the -i Hbar * Mbar_R and put it in imaginary component

      Mubar_s_pq_t0 = 390
      Mubar_d_pq_t0 = 313

      If (Iside .EQ. 1) Then
         Ioffr1 = 0
         Ioffr2 = 0
         Ioffsp = 0

C This write Mubar(t_n),  iside=1

         Call Tdee_dump_vec(Irrepx,Dummy1,Nsize,Mubar_s_pq_t0,
     +                      Ioffr1,Ioffsp,Mubar_d_pq_t0,Ioffr2,
     +                      Iuhf,.False.)

      Elseif (Iside .EQ. 2) Then
         Ioffr1 = 2
         Ioffr2 = 10
         Ioffsp = 0

C This write Mubar(t_n),  iside=2

         Call Tdee_dump_vec(Irrepx,Dummy1,Nsize,Mubar_s_pq_t0,
     +                      Ioffr1,Ioffsp,Mubar_d_pq_t0,Ioffr2,
     +                      Iuhf,.False.)
      Endif

#ifdef _DEBUG_LVL0
      Write(6,*)
      If (type) Then
         write(6,"(a)") "The real part of the Mbar(t)"
         Call checksum("tdee_pes_ode_in:",Dummy1,Nsize,s)
      else
         write(6,"(a)") "The imginary part of the Mbar(t)"
         Call checksum("tdee_pes_ode_in:",Dummy1,Nsize,s)
      endif 
#endif

      Call Tdee_return_mubar_dot(Icore(I0),Memleft,Iuhf,Irrepx,
     +                           Iside)

C The data returned from the multiplication is stored as

C The data returned from the multiplication is stored as
C follows.
C UHF   Right Left
C List: 394   396 AA(1),BB(2)
C List: 334   344 AAAA
C List: 335   345 BBBB
C List: 336   346 ABAB
C RHF
C List: 394   396 AA(1)
C List: 336   346 ABAB

      Mubar_s_pq_tn = 394
      Mubar_d_pq_tn = 333

      If (Iside .EQ. 1) Then

        Ioffr1 = 0
        Ioffr2 = 0

C This Returns Mubar(t_n+1), iside=1

        Call Tdee_Load_vec(Irrepx,Dummy2,Nsize,Mubar_s_pq_tn,
     +                     Ioffr1,Mubar_d_pq_tn,Ioffr2,Iuhf,.False.)


      Elseif (Iside .EQ.2) Then

        Ioffr1 = 2
        Ioffr2 = 10

C This Returns Mubar(t_n+1), iside=2

        Call Tdee_Load_vec(Irrepx,Dummy2,Nsize,Mubar_s_pq_tn,
     +                     Ioffr1,Mubar_d_pq_tn,Ioffr2,Iuhf,.False.)

      Endif

      If (Type) Call Dscal(Nsize,Mone,Dummy2,1)

#ifdef _DEBUG_LVL0
      Write(6,*)
      If (type) Then
         write(6,"(a)") "The real part of the Mbar_dot(t)"
         Call checksum("tdee_pes_ode_in:",Dummy2,Nsize,s)
      else
         write(6,"(a)") "The imginary part of the Mbar_dot(t)"
         Call checksum("tdee_pes_ode_in:",Dummy2,Nsize,s)
      endif 
#endif

      Return
      End

