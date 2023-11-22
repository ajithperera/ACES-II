      Subroutine Pccd_proc_ovvo(Gov,Gvo,Work,Maxcor,Nocc,Nvrt)

      Implicit Double Precision(A-H,O-Z)
    
      Dimension Gvo(Nocc,Nvrt)
      Dimension Gov(Nvrt,Nocc)
      Dimension Work(Maxcor)

      Data Ione /1/

      I000 = Ione
      I010 = I000 + Nocc*Nvrt
      Iend = I010 + Nocc*Nvrt

      Call Transp(Gov,Work(I000),Nvrt,Nocc)
      Call Dcopy(Nocc*Nvrt,Gvo,1,Work(I010),1)
      Call Dcopy(Nocc*Nvrt,Work(I000),1,Gvo,1)

      Call Transp(Work(I010),Work(I000),Nocc,Nvrt)
      Call Dcopy(Nocc*Nvrt,Work(I000),1,Gov,1)

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") "G2(vo)=F(vv)*Dvo contribution"
      call output(Gvo,1,Nvrt,1,Nocc,Nvrt,Nocc,1)
      Write(6,*)
      Write(6,"(a)") "G2(ov)=F(oo)*D(vo)(t) contribution"
      call output(Gov,1,Nocc,1,Nvrt,Nocc,Nvrt,1)
#endif
      Return
      End
