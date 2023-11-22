#include "flags.h"
      Subroutine Fix_qrhf_evecs(EvecA,EvecB,Dcore,Maxdcor,Nbas,Iuhf)

      Implicit Double Precision(A-H, O-Z)

#include "maxbasfn.par"
#include "machsp.com"
#include "symm2.com"
C
      Dimension EvecA(Nbas*Nbas),EvecB(Nbas*Nbas),Dcore(Maxdcor)
      Dimension Idummy(Maxbasfn),Idummy2(Maxbasfn)
      Dimension Idummy3(Maxbasfn),Ilocate(Maxbasfn)
      Dimension Irem(8,2),Iadd(8,2), Nocc_save(16)
C
      Common /Popul/ Nocc(8,2)
C    
      Ione = 1
      Call Izero(Irem,16)
      Call Icopy(16, Nocc, 1, Nocc_save, 1)
    
      Do I=1, Nbas
         Ilocate(I) = I
      Enddo

      CALL GETREC(-1,'JOBARC','QRHFTOT ',Ione,NMODIFY)
      CALL GETREC(-1,'JOBARC','QRHFIRR ',NMODIFY,IDUMMY)
      CALL GETREC(-1,'JOBARC','QRHFLOC ',NMODIFY,IDUMMY2)
      CALL GETREC(-1,'JOBARC','QRHFSPN ',NMODIFY,IDUMMY3)
      CALL GETREC(-1,'JOBARC','OSCALC  ',Ione,IOS)
      
      Call Getrec(20, "JOBARC", "SCFEVCA0", Nbas*Nbas*Iintfp,
     &            EvecA)
      If (Iuhf .Ne. 0) Call Getrec(20, "JOBARC", "SCFEVCB0", 
     &                      Nbas*Nbas*Iintfp, EvecB)

      Do I = 1, Nmodify

         Irrp = Idummy(I)
         ILoc = Idummy2(I)
         Ispn = Idummy3(I)
         
         If (Irrp .lt. 0) Then
             If (Ispn .ne. 1) Ispn = 2
                Irrp = -Irrp
                Iposorg = Irpoff(Irrp) + Nocc(Irrp,Ispn) + 1 -
     &                    Max(Iloc,1)  + Irem(Irrp,Ispn) - 
     &                    Iadd(Irrp,Ispn)
                Iposabs = Ilocate(Iposorg)
                Irem(Irrp,Ispn) = Irem(Irrp,Ispn) + 1
#ifdef _DEBUG_LVLM
      write(6,"(a,1x,7I4)") "Irrp,Iloc,Ispn,Iposorg,Iposabs:",
     &                       Irrp,Iloc,Ispn,Iposorg,Iposabs,
     &                       Inewvc(Iposabs),
     &                       Nocc(Irrp,Ispn)-Inewvc(Iposabs)
#endif 
                Ichk = 1 + Nocc(Irrp,Ispn) - Inewvc(Iposabs) 
                Ibgn = Nocc(Irrp,Ispn)
                Iend = Nocc(Irrp,Ispn) + 1  - Ichk
#ifdef _DEBUG_LVLM
       write(6,"(a,1x,3I4)") "Ichk,Ibgn,Iend", Ichk,Ibgn, Iend
#endif 
C                If (Ichk .gt. 1) Then
                    If (Ios .eq. 0) then
                        Ioff_bgn  = Nbas*(Ibgn-1) + 1
                        Ioff_end  = Nbas*(Iend-1) + 1
                        Call Dswap(Nbas, EvecA(Ioff_end), 1, 
     &                             EvecA(Ioff_bgn), 1)
                        Call Dswap(Nbas, EvecB(Ioff_end), 1, 
     &                             EvecB(Ioff_bgn), 1)

                        Itmp = Ilocate(Ibgn)
                        Ilocate(Ibgn) = Ilocate(Iend) 
                        Ilocate(Iend) = Itmp
                    Else
                        if (ispn .eq. 1) then
                           Ioff_bgn  = Nbas*(Ibgn-1) + 1
                           Ioff_end  = Nbas*(Iend-1) + 1
                           Call Dswap(Nbas, EvecA(Ioff_end), 1,
     &                                EvecA(Ioff_bgn), 1)

                           Itmp = Ilocate(Ibgn)
                           Ilocate(Ibgn) = Ilocate(Iend)
                           Ilocate(Iend) = Itmp
                        elseif (ispn .eq. 2) then
                           Ioff_bgn  = Nbas*(Ibgn-1) + 1
                           Ioff_end  = Nbas*(Iend-1) + 1
                           Call Dswap(Nbas, EvecA(Ioff_end), 1,
     &                                EvecB(Ioff_bgn), 1)

                           Itmp = Ilocate(Ibgn)
                           Ilocate(Ibgn) = Ilocate(Iend)
                           Ilocate(Iend) = Itmp
                        endif 
                    Endif 
C                Endif 

                 Nocc(Irrp, Ispn) = Nocc(Irrp,Ispn) - 1
         Else
                if (Ispn .ne. 2) Ispn = 1
                   Iposorg = Irpoff(Irrp) + Nocc(Irrp,Ispn) +
     &                       Max(Iloc,1)  + Irem(Irrp,Ispn) -
     &                       Iadd(Irrp,Ispn)
                Iposabs = Ilocate(Iposorg)
                Iadd(Irrp,Ispn) = Iadd(Irrp,Ispn) + 1
#ifdef _DEBUG_LVLM
      write(6,"(a,1x,5I4)") "Irrp,Iloc,Ispn,Iposorg,Iposabs:",
     &                       Irrp,Iloc,Ispn,Iposorg,Iposabs
#endif 
                Ichk = Inewvc(Iposabs) - Nocc(Irrp,Ispn)
                Ibgn = Nocc(Irrp,Ispn) + 1
                Iend = Nocc(Irrp,Ispn) + ichk

#ifdef _DEBUG_LVLM
       write(6,"(a,1x,3I4)") "Ichk,Ibgn,Iend", Ichk,Ibgn, Iend
#endif 
                If (Ichk .gt. 1) then
                    If (Ios .eq. 0) then
                        Ioff_bgn  = Nbas*(Ibgn-1) + 1
                        Ioff_end  = Nbas*(Iend-1) + 1
                        Call Dswap(Nbas, EvecA(Ioff_end), 1, 
     &                             EvecA(Ioff_bgn), 1)
                        Call Dswap(Nbas, EvecB(Ioff_end), 1, 
     &                             EvecB(Ioff_bgn), 1)
                   
                       Itmp = Ilocate(Ibgn)
                       Ilocate(Ibgn) = Ilocate(Iend) 
                       Ilocate(Iend) = Itmp
                     Else
                       If (ispn .eq. 1) then
                          Ioff_bgn  = Nbas*(Ibgn-1) + 1
                          Ioff_end  = Nbas*(Iend-1) + 1
                          Call Dswap(Nbas, EvecA(Ioff_end), 1, 
     &                               EvecA(Ioff_bgn), 1)
                   
                          Itmp = Ilocate(Ibgn)
                          Ilocate(Ibgn) = Ilocate(Iend) 
                          Ilocate(Iend) = Itmp
                        Elseif (ispn .eq. 2) then 
                          Ioff_bgn  = Nbas*(Ibgn-1) + 1
                          Ioff_end  = Nbas*(Iend-1) + 1
                          Call Dswap(Nbas, EvecA(Ioff_end), 1,
     &                               EvecB(Ioff_bgn), 1)
                          Itmp = Ilocate(Ibgn)
                          Ilocate(Ibgn) = Ilocate(Iend)
                          Ilocate(Iend) = Itmp
                        Endif 

                     Endif 

                 Endif 
                 Nocc(Irrp, Ispn) = Nocc(Irrp,Ispn) + 1
         Endif 
           
      Enddo 

      Call Putrec(20, "JOBARC", "SCFEVCA0", Nbas*Nbas*Iintfp,
     &            EvecA)
      If (Iuhf .Ne. 0) Call Putrec(20, "JOBARC", "SCFEVCB0", 
     &                      Nbas*Nbas*Iintfp, EvecB)
      Call Icopy(16, Nocc_save, 1, Nocc, 1)
 
      Return
      End

