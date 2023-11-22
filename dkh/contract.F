      Subroutine Contract(Dkh_pints,Ccoef_ishell,Ccoef_jshell,Iprimfns,
     +                    Jprimfns,Icntrfns,Jcntrfns,Iprim_4shell,
     +                    Jprim_4shell,Icntr_4shell,Jcntr_4shell,
     +                    Tot_cnt_fns,Tot_prm_fns,Idegen,Jdegen,
     +                                              Dkh_cints)

      Implicit Double Precision (A-H, O-Z)

      Integer Tot_cnt_fns
      Dimension Dkh_pints(Idegen*Iprim_4shell,Jdegen*Jprim_4shell)
      Dimension Ccoef_ishell(Iprim_4shell,Icntr_4shell)
      Dimension Ccoef_jshell(Jprim_4shell,Jcntr_4shell)
      Dimension Dkh_cints(Idegen*Icntr_4shell,Jdegen*Jcntr_4shell)

#ifdef _DEBUG_LVLM
      Write(6,*)
      write(6,"(1x,a,a)") "Printing from Contract: The primitve DKH",
     +           " ints per shell block."
      call output(Dkh_pints,1,Idegen*Iprim_4shell,1,Jdegen*Jprim_4shell,
     +            Idegen*Iprim_4shell,Jdegen*Jprim_4shell,1)
      write(6,*) "The contraction coefs."
      write(6,"(6(1x,f10.6))")((CCoef_ishell(i,j),i=1,Iprim_4shell),
     +                          j=1,Icntr_4shell)
      write(6,*)
      write(6,"(6(1x,f10.6))")((CCoef_jshell(i,j),i=1,Jprim_4shell),
     +                          j=1,Jcntr_4shell)
      Write(6,*) 
#endif 
      Call Dzero(Dkh_cints,Tot_cnt_fns)

      Int = 0

      Do Jcint = 1, Jcntr_4shell
         Do Icint = 1, Icntr_4shell
            Int = Int + 1

            Do Jprim = 1, Jprim_4shell
               Do Iprim = 1, Iprim_4shell

                  Do J = 1, Jdegen
                     J_cint = (Jcint-1) * Jdegen + J
                     J_prim = (Jprim-1) * Jdegen + J
                     Do I = 1, Idegen
                        I_cint = (Icint-1) * Idegen + I
                        I_prim = (Iprim-1) * Idegen + I

                     Weight = Ccoef_ishell(Iprim,Icint) *
     +                        Ccoef_jshell(Jprim,Jcint)

                     Dkh_cints(I_cint,J_cint) =

     +                        Dkh_cints(I_cint,J_cint) +
     +                        Dkh_pints(I_prim,J_prim) * Weight

#ifdef _DEBUG_LVLM
      write(6,"(1x,5i4,5F10.5)") Iprm,Jprm,Icint,jcint,int,
     +     Ccoef_ishell(Iprim,Icint),Ccoef_jshell(Jprim,Jcint),
     +     Weight,Dkh_pints(I_prim,J_prim),Dkh_cints(I_cint,J_cint)
#endif
                     Enddo
                  Enddo
C
                Enddo
            Enddo
C
         Enddo
      Enddo 
C
#ifdef _DEBUG_LVLM
      Write(6,*)
      Write(6,"(a,2i4)") "The contracted DKH ints per shell block"
      call output(Dkh_cints,1,Idegen*Jdegen,1,Icntr_4shell*Jcntr_4shell,
     +            Idegen*Jdegen,Icntr_4shell*Jcntr_4shell,1)
#endif 
      
      Return
      End 
     
