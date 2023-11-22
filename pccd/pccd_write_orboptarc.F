      Subroutine Pccd_write_orboptarc(Hess_vo,Hess_oo,Hess_vv,Grad_vo,
     +                                Grad_ov,Grad_oo,Grad_vv,Nsize_vo,
     +                                Nsize_oo,Nsize_vv,Lenvo,Lenoo,
     +                                Lenvv,Opt_unit)


      Implicit Double Precision(A-H,O-Z)
      Integer Opt_unit 

      Dimension Hess_vo(Nsize_vo),Hess_oo(Nsize_oo) 
      Dimension Hess_vv(Nsize_vv)
      Dimension Grad_ov(Lenvo),Grad_vo(Lenvo)
      Dimension Grad_oo(Lenoo),Grad_vv(Lenvv)

      Write(Opt_unit,*) (Hess_vo(i),i=1,Nsize_vo)
      Write(Opt_unit,*)
      Write(Opt_unit,*) (Hess_oo(i),i=1,Nsize_oo)
      Write(Opt_unit,*)
      Write(Opt_unit,*) (Hess_vv(i),i=1,Nsize_vv)
      Write(Opt_unit,*)
      Write(Opt_unit,*) (Grad_vo(i),i=1,Lenvo)
      Write(Opt_unit,*)
      Write(Opt_unit,*) (Grad_ov(i),i=1,Lenvo)
      Write(Opt_unit,*)
      Write(Opt_unit,*) (Grad_oo(i),i=1,Lenoo)
      Write(Opt_unit,*)
      Write(Opt_unit,*) (Grad_vv(i),i=1,Lenvv)

      Return
      End

