      Subroutine Osc_fin(Osc,Osc_x,Osc_y,Osc_z,Jroot_1,Irrepy_1,Iroot_1,
     +                   Irrepx_1,Length)

      Implicit Double Precision(A-H,O-z)

      Dimension Osc_x(4,8,4,8,4,8,4,8)
      Dimension Osc_y(4,8,4,8,4,8,4,8)
      Dimension Osc_z(4,8,4,8,4,8,4,8)
      Dimension Osc(1024)

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a,1x,F15.8)") "Osc_x: ",
     +                         Osc_x(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +                         Iroot_1,Irrepx_1,Jroot_1,Irrepy_1)
      Write(6,"(a,1x,F15.8)") "Osc_y: ",
     +                         Osc_y(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +                         Iroot_1,Irrepx_1,Jroot_1,Irrepy_1)
      Write(6,"(a,1x,F15.8)") "Osc_z: ",
     +                         Osc_z(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +                         Iroot_1,Irrepx_1,Jroot_1,Irrepy_1)
#endif 

      Osc(Length) = Osc(Length) +
     +              (Osc_x(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +                     Iroot_1,Irrepx_1,Jroot_1,Irrepy_1) + 
     +               Osc_y(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +                     Iroot_1,Irrepx_1,Jroot_1,Irrepy_1) +
     +               Osc_z(Jroot_1,Irrepy_1,Iroot_1,Irrepx_1,
     +                     Iroot_1,Irrepx_1,Jroot_1,Irrepy_1))

      Return
      End


