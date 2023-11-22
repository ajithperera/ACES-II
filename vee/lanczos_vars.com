#ifndef _LANCZOS_VARS_COM
#define _LANCZOS_VARS_COM

#ifdef _DULAL
C Nlim : The maximum allowed number of look-ahead Lanczos steps.
C Maxn : Lead dimension of the tridiagonal matrix
C Maxvw: Maximum allowed block size for the VW step.
C M    : Lead dimension of arrays iwk and dwk.
C Norm : Estimate of the norm of the matrix that is being
C        diagonalized. Typically set to 1. The algorithm
C        dynamically update this. Very large value effectively 
C        turn off look-ahead steps. 
C Info : Serve as both input/output. Only Info(1) is 
C        set by the caller. 

C The internally set values are the maximums allowed.

      Parameter(Maxn_int=1002,Maxvw_int=10,M_int=2*Maxvw_int+2)
      Parameter(No_max_vecs=10000)

      Integer Nlim,M,Maxn,MaxvW
      Double Precision Dwk(5*M_int*M_int+7*M_int)
      Double Precision HWk(4*Maxn_int*Maxn_int+8*Maxn_int),Norm 
      Integer Info(4),Iwk(4*M_int),Idx(3*Maxn_int)

      Common /lancsoz_vars/Nlim,Maxn,Maxvw,M,Norm,Iwk,Idx,Hwk,
     +                     Info 
#else 
 
C Maxn : Lead dimension of the tridiagonal matrix
C Info : Serve as both input/output. Only Info(1) is 
C        set by the caller. 
C Maxvw: This is set to 1. This allocate one extra vector 
C        of the diemension of Hbar (this space is not used).
C Convrg: Strict or light 

      Integer Maxn
      Parameter(Maxn_int=10000,Maxvw_int=1)
      Character*6 Convrg 

      Common /lanczos_vars/Maxn,Maxvw,Trd(Maxn_int,Maxn_int),
     +                     Convrg 
      
#endif 

#endif  /* _LANCZOS_VARS_COM__ */


